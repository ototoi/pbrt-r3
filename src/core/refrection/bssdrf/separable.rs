use crate::core::pbrt::*;

pub type BSSRDFMaterialRawPointer = *const dyn Material;

pub fn sw(eta: Float, w: &Vector3f) -> Spectrum {
    let c = 1.0 - 2.0 * fresnel_moment1(1.0 / eta);
    return Spectrum::from((1.0 - fr_dielectric(cos_theta(w), 1.0, eta)) / (c * PI));
}

pub struct SeparableBSSRDFAdapter {
    eta: Float,
    mode: TransportMode,
}

impl SeparableBSSRDFAdapter {
    pub fn new(eta: Float, mode: TransportMode) -> Self {
        SeparableBSSRDFAdapter { eta, mode }
    }
    pub fn sw(&self, w: &Vector3f) -> Spectrum {
        sw(self.eta, w)
    }
}

impl BxDF for SeparableBSSRDFAdapter {
    fn f(&self, _wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let f = self.sw(wi);
        if self.mode == TransportMode::Radiance {
            let eta = self.eta;
            return f * (eta * eta);
        } else {
            return f;
        }
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        sample: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        self.sample_f_default(wo, sample)
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        self.pdf_default(wo, wi)
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_REFLECTION | BSDF_DIFFUSE;
    }

    fn to_string(&self) -> String {
        return format!(
            "SeparableBSSRDFAdapter {:?} {:?}",
            self.get_type(),
            self.mode
        );
    }
}

pub trait SeparableBSSRDF {
    fn projection_axis(&self, u1: Float) -> (Vector3f, Vector3f, Vector3f, Float);
    fn get_p(&self) -> Point3f;
    fn get_time(&self) -> Float;
    fn get_material(&self) -> BSSRDFMaterialRawPointer;
    fn sp(&self, si: &SurfaceInteraction) -> Spectrum;
    fn sample_sr(&self, ch: usize, u: Float) -> Float;
    fn pdf_sp(&self, si: &SurfaceInteraction) -> Float;
    fn sample_sp(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        _arena: &mut MemoryArena,
    ) -> Option<(Spectrum, SurfaceInteraction, Float)> {
        let (vx, vy, vz, u1) = self.projection_axis(u1);
        let n_samples = Spectrum::N_SAMPLES;
        // Choose spectral channel for BSSRDF sampling
        let ch = usize::clamp((u1 * n_samples as Float) as usize, 0, n_samples - 1);
        let u1 = u1 * n_samples as Float - ch as Float;

        // Sample BSSRDF profile in polar coordinates
        let r = self.sample_sr(ch, u2[0]);
        if r < 0.0 {
            return None;
        }
        let phi = 2.0 * PI * u2[1];

        // Compute BSSRDF profile bounds and intersection height
        let r_max = self.sample_sr(ch, 0.999);
        if r_max < 0.0 {
            return None;
        }
        if r >= r_max {
            return None;
        }
        let l = 2.0 * Float::sqrt(r_max * r_max - r * r);

        assert!(l > 1e-6);
        //println!("l:{:?}", l);

        // Compute BSSRDF sampling ray segment
        let p = self.get_p() + r * (vx * Float::cos(phi) + vy * Float::sin(phi)) - l * vz * 0.5;
        let time = self.get_time();
        let mut base = SurfaceInteraction {
            p: p,
            time: time,
            ..Default::default()
        };
        let p_target = p + l * vz;

        // Intersect BSSRDF sampling ray against the scene geometry

        let mut chain = Vec::new();
        loop {
            let distance = Vector3f::distance(&base.p, &p_target);
            if distance <= MACHINE_EPSILON {
                break;
            }
            let r = base.spawn_ray_to_point(&p_target);

            if let Some(si) = scene.intersect(&r) {
                base = si.clone();
                if let Some(primitive) = si.primitive.as_ref() {
                    let primitive = primitive.upgrade().unwrap();
                    let primitive = primitive.as_ref();
                    if let Some(material) = primitive.get_material() {
                        let material = material.as_ref();
                        let ptr = material as BSSRDFMaterialRawPointer;
                        // Append admissible intersection to _IntersectionChain_
                        #[allow(clippy::vtable_address_comparisons)]
                        if std::ptr::eq(ptr, self.get_material()) {
                            chain.push(si);
                        }
                    }
                }
            } else {
                break;
            }
        }

        // Randomly choose one of several intersections during BSSRDF sampling
        if chain.is_empty() {
            return None;
        }
        let n_found = chain.len();
        //println!("n_found:{:?}", n_found);
        let selected = usize::clamp((u1 * n_found as Float) as usize, 0, n_found - 1);
        let pi = chain[selected].clone();

        // Compute sample PDF and return the spatial BSSRDF term $\Sp$
        let pdf = self.pdf_sp(&pi) / (n_found as Float);
        let f = self.sp(&pi);
        return Some((f, pi, pdf));
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BaseSeparableBSSRDF {
    pub base: BaseBSSRDF,
    pub ns: Normal3f,
    pub ss: Vector3f,
    pub ts: Vector3f,
    pub material: BSSRDFMaterialRawPointer,
    pub mode: TransportMode,
}

impl BaseSeparableBSSRDF {
    pub fn new(
        po: &SurfaceInteraction,
        eta: Float,
        material: BSSRDFMaterialRawPointer,
        mode: TransportMode,
    ) -> Self {
        let ns = po.shading.n;
        let ss = po.shading.dpdu.normalize();
        let ts = Vector3f::cross(&ns, &ss);
        BaseSeparableBSSRDF {
            base: BaseBSSRDF::new(po, eta),
            ns,
            ss,
            ts,
            material,
            mode,
        }
    }

    pub fn projection_axis(&self, u1: Float) -> (Vector3f, Vector3f, Vector3f, Float) {
        if u1 < 0.5 {
            let vx = self.ss;
            let vy = self.ts;
            let vz = self.ns;
            let u1 = u1 * 2.0;
            return (vx, vy, vz, u1);
        } else if u1 < 0.75 {
            let vx = self.ts;
            let vy = self.ns;
            let vz = self.ss;
            let u1 = (u1 - 0.5) * 4.0;
            return (vx, vy, vz, u1);
        } else {
            let vx = self.ns;
            let vy = self.ss;
            let vz = self.ts;
            let u1 = (u1 - 0.75) * 4.0;
            return (vx, vy, vz, u1);
        }
    }
}
