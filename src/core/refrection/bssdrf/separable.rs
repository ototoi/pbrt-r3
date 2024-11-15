use crate::core::pbrt::*;
use std::sync::Arc;

pub type BSSRDFMaterialRawPointer = *const dyn Material;

pub struct SeparableBSSRDFAdapter {
    bssrdf: BaseSeparableBSSRDF,
}

impl SeparableBSSRDFAdapter {
    pub fn new(bssrdf: &BaseSeparableBSSRDF) -> Self {
        SeparableBSSRDFAdapter { bssrdf: *bssrdf }
    }
}

impl BxDF for SeparableBSSRDFAdapter {
    fn f(&self, _wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let f = self.bssrdf.sw(wi);
        if self.bssrdf.mode == TransportMode::Radiance {
            let eta = self.bssrdf.base.eta;
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
            self.bssrdf.mode
        );
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

    pub fn sw(&self, w: &Vector3f) -> Spectrum {
        return self.base.sw(w);
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

pub trait SeparableBSSRDF: BSSRDF {
    fn s_as_separable(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum {
        let separable = self.as_separable();
        let eta = separable.base.eta;
        let base = &separable.base;
        let ft = fr_dielectric(cos_theta(&base.wo), 1.0, eta);
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

    fn sw(&self, w: &Vector3f) -> Spectrum {
        let separable = self.as_separable();
        return separable.base.sw(w);
    }

    fn sp(&self, pi: &SurfaceInteraction) -> Spectrum {
        let separable = self.as_separable();
        let base = &separable.base;
        return self.sr(Vector3f::distance(&base.p, &pi.p));
    }

    fn sample_sp(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        _arena: &mut MemoryArena,
    ) -> Option<(Spectrum, SurfaceInteraction, Float)> {
        let separable = self.as_separable();
        let (vx, vy, vz, u1) = separable.projection_axis(u1);
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
        let p = separable.base.p + r * (vx * Float::cos(phi) + vy * Float::sin(phi)) - l * vz * 0.5;
        let time = separable.base.time;
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
                        if std::ptr::eq(ptr, separable.material) {
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

    fn pdf_sp(&self, pi: &SurfaceInteraction) -> Float {
        // Express $\pti-\pto$ and $\bold{n}_i$ with respect to local coordinates at
        // $\pto$
        let separable = self.as_separable();
        let po = &separable.base;
        let d = po.p - pi.p;

        assert!(pi.n.length() > 0.0);

        let d_local = Vector3f::new(
            Vector3f::dot(&separable.ss, &d),
            Vector3f::dot(&separable.ts, &d),
            Vector3f::dot(&separable.ns, &d),
        );
        let n_local = Vector3f::new(
            Vector3f::dot(&separable.ss, &pi.n),
            Vector3f::dot(&separable.ts, &pi.n),
            Vector3f::dot(&separable.ns, &pi.n),
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

    // SeparableBSSRDF Interface
    fn sr(&self, d: Float) -> Spectrum;
    fn sample_sr(&self, ch: usize, u: Float) -> Float;
    fn pdf_sr(&self, ch: usize, r: Float) -> Float;

    //
    fn as_separable(&self) -> &BaseSeparableBSSRDF;
}
