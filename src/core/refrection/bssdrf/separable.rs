use crate::core::pbrt::*;
use std::sync::Arc;

pub type BSSRDFMaterialRawPointer = *const dyn Material;

pub trait SeparableBSSRDF {
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

impl SeparableBSSRDF for BaseSeparableBSSRDF {
    fn sw(&self, w: &Vector3f) -> Spectrum {
        let eta = self.base.eta;
        let c = 1.0 - 2.0 * fresnel_moment1(1.0 / eta);
        return Spectrum::from((1.0 - fr_dielectric(cos_theta(w), 1.0, eta)) / (c * PI));
    }

    fn get_mode(&self) -> TransportMode {
        return self.mode;
    }

    fn get_eta(&self) -> Float {
        return self.base.eta;
    }
}
