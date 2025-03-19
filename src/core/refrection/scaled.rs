use crate::core::pbrt::*;
use crate::core::refrection::*;
use crate::core::spectrum::*;

use std::sync::Arc;

pub struct ScaledBxDF {
    bxdf: Arc<dyn BxDF>,
    scale: Spectrum,
}

impl ScaledBxDF {
    pub fn new(bxdf: &Arc<dyn BxDF>, scale: &Spectrum) -> Self {
        ScaledBxDF {
            bxdf: Arc::clone(bxdf),
            scale: *scale,
        }
    }
}

impl BxDF for ScaledBxDF {
    fn rho(&self, wo: &Vector3f, samples: &[Point2f]) -> Spectrum {
        let bxdf = self.bxdf.as_ref();
        return self.scale * bxdf.rho(wo, samples);
    }

    fn rho2(&self, samples: &[(Point2f, Point2f)]) -> Spectrum {
        let bxdf = self.bxdf.as_ref();
        return self.scale * bxdf.rho2(samples);
    }

    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let bxdf = self.bxdf.as_ref();
        return self.scale * bxdf.f(wo, wi);
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        sample: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        let bxdf = self.bxdf.as_ref();
        if let Some((spec, wi, pdf, t)) = bxdf.sample_f(wo, sample) {
            return Some((self.scale * spec, wi, pdf, t));
        }
        return None;
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        let bxdf = self.bxdf.as_ref();
        return bxdf.pdf(wo, wi);
    }

    fn get_type(&self) -> BxDFType {
        let bxdf = self.bxdf.as_ref();
        return bxdf.get_type();
    }

    fn to_string(&self) -> String {
        return format!("ScaledBxDF {:?}", self.get_type());
    }
}
