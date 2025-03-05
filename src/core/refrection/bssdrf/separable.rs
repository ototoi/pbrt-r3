use crate::core::pbrt::*;
use std::sync::Arc;

pub type BSSRDFMaterialRawPointer = *const dyn Material;

pub trait SeparableBSSRDF: BSSRDF {
    fn sr(&self, d: Float) -> Spectrum;
    fn sample_sr(&self, ch: usize, u: Float) -> Float;
    fn pdf_sr(&self, ch: usize, r: Float) -> Float;
    fn sw(&self, w: &Vector3f) -> Spectrum;
    fn get_mode(&self) -> TransportMode;
    fn get_eta(&self) -> Float;
}

pub struct SeparableBSSRDFAdapter {
    bssrdf: Arc<dyn SeparableBSSRDF>,
}

impl SeparableBSSRDFAdapter {
    pub fn new(bssrdf: Arc<dyn SeparableBSSRDF>) -> Self {
        SeparableBSSRDFAdapter { bssrdf }
    }
}

impl BxDF for SeparableBSSRDFAdapter {
    fn f(&self, _wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let f = self.bssrdf.sw(wi);
        if self.bssrdf.get_mode() == TransportMode::Radiance {
            let eta = self.bssrdf.get_eta();
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
            self.bssrdf.get_mode()
        );
    }
}

/*
pub trait SeparableBSSRDF: BSSRDF {
    fn s_as_separable(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum {
        let eta = self.get_eta();
        let wo = self.get_wo();
        let ft = fr_dielectric(cos_theta(&wo), 1.0, eta);
        return (self.sp(pi) * self.sw(wi)) * (1.0 - ft);
    }

    fn sample_s_as_separable(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        arena: &mut MemoryArena,
    ) -> Option<(Spectrum, SurfaceInteraction, Float)> {
        if let Some((sp, mut si, pdf)) = self.sample_sp(scene, u1, u2, arena) {
            if !sp.is_black() {
                let mut b = arena.alloc_bsdf(&si, 1.0);
                let adapter: Arc<dyn BxDF> =
                    Arc::new(SeparableBSSRDFAdapter::new(self.as_separable()));
                b.add(&adapter);
                let bsdf = Arc::new(b);
                si.bsdf = Some(bsdf);
                si.wo = si.shading.n;
            }
            return Some((sp, si, pdf));
        }
        return None;
    }

    fn sw(&self, w: &Vector3f) -> Spectrum;

    fn sp(&self, pi: &SurfaceInteraction) -> Spectrum {
        let p = self.get_p();
        return self.sr(Vector3f::distance(&p, &pi.p));
    }

    

    

    // SeparableBSSRDF Interface
    fn sr(&self, d: Float) -> Spectrum;
    fn sample_sr(&self, ch: usize, u: Float) -> Float;
    fn pdf_sr(&self, ch: usize, r: Float) -> Float;

    fn get_mode(&self) -> TransportMode;
    fn get_eta(&self) -> Float;
    fn get_p(&self) -> Point3f;
    fn get_wo(&self) -> Vector3f;
}
*/

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
}

/* 
    fn sp(&self, pi: &SurfaceInteraction) -> Spectrum {
        let p = self.base.p;
        return self.sr(Vector3f::distance(&p, &pi.p));
    }

    fn sample_sp(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        _arena: &mut MemoryArena,
    ) -> Option<(Spectrum, SurfaceInteraction, Float)> {
        let (vx, vy, vz, u1) = self.projection_axis(u1);
        let n_samples = 3; //TODO
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
        let p = self.base.p + r * (vx * Float::cos(phi) + vy * Float::sin(phi)) - l * vz * 0.5;
        let time = self.base.time;
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
                        if std::ptr::eq(ptr, self.material) {
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

    pub fn pdf_sp(&self, pi: &SurfaceInteraction) -> Float {
        // Express $\pti-\pto$ and $\bold{n}_i$ with respect to local coordinates at
        // $\pto$
        let po = &self.base;
        let d = po.p - pi.p;

        assert!(pi.n.length() > 0.0);

        let d_local = Vector3f::new(
            Vector3f::dot(&self.ss, &d),
            Vector3f::dot(&self.ts, &d),
            Vector3f::dot(&self.ns, &d),
        );
        let n_local = Vector3f::new(
            Vector3f::dot(&self.ss, &pi.n),
            Vector3f::dot(&self.ts, &pi.n),
            Vector3f::dot(&self.ns, &pi.n),
        );
        // Compute BSSRDF profile radius under projection along each axis
        let r_proj = [
            Float::sqrt(d_local.y * d_local.y + d_local.z * d_local.z),
            Float::sqrt(d_local.z * d_local.z + d_local.x * d_local.x),
            Float::sqrt(d_local.x * d_local.x + d_local.y * d_local.y),
        ];

        // Return combined probability from all BSSRDF sampling strategies
        let mut pdf = 0.0;
        let axis_prob = [0.25, 0.25, 0.5];
        let n_samples = 3; //TODO
        let ch_prob = 1.0 / (n_samples as Float);
        for axis in 0..3 {
            for ch in 0..n_samples {
                pdf += self.pdf_sr(ch, r_proj[axis])
                    * Float::abs(n_local[axis])
                    * ch_prob
                    * axis_prob[axis];
            }
        }
        return pdf;
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

impl BSSRDF for BaseSeparableBSSRDF {
    fn s(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum {
        let eta = self.base.eta;
        let wo = self.base.wo;
        let ft = fr_dielectric(cos_theta(&wo), 1.0, eta);
        return (self.sp(pi) * self.sw(wi)) * (1.0 - ft);
    }

    fn sample_s(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        arena: &mut MemoryArena,
    ) -> Option<(Spectrum, SurfaceInteraction, Float)> {
        if let Some((sp, mut si, pdf)) = self.sample_sp(scene, u1, u2, arena) {
            if !sp.is_black() {
                let mut b = arena.alloc_bsdf(&si, 1.0);
                let adapter: Arc<dyn BxDF> =
                    Arc::new(SeparableBSSRDFAdapter::new(self.as_separable()));
                b.add(&adapter);
                let bsdf = Arc::new(b);
                si.bsdf = Some(bsdf);
                si.wo = si.shading.n;
            }
            return Some((sp, si, pdf));
        }
        return None;
    }
}

impl SeparableBSSRDF for BaseSeparableBSSRDF {
    fn sr(&self, _d: Float) -> Spectrum {
        unimplemented!()
    }
    fn sample_sr(&self, _ch: usize, _u: Float) -> Float {
        unimplemented!()
    }
    fn pdf_sr(&self, _ch: usize, _r: Float) -> Float  {
        unimplemented!()
    }
    fn sw(&self, _w: &Vector3f) -> Spectrum {
        unimplemented!()
    }

    fn get_mode(&self) -> TransportMode {
        self.mode
    }

    fn get_eta(&self) -> Float {
        self.base.eta
    }
}

    */
