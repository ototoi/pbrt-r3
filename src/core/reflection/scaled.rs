use crate::core::base::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;

pub struct ScaledBxDF {
    bxdf: Box<BxDFEnum>,
    scale: Spectrum,
}

impl ScaledBxDF {
    pub fn new(bxdf: BxDFEnum, scale: &Spectrum) -> Self {
        ScaledBxDF {
            bxdf: Box::new(bxdf),
            scale: *scale,
        }
    }
}

impl BxDF for ScaledBxDF {
    fn rho(&self, wo: &Vector3f, samples: &[Point2f]) -> Spectrum {
        return self.scale * self.bxdf.rho(wo, samples);
    }

    fn rho2(&self, samples: &[(Point2f, Point2f)]) -> Spectrum {
        return self.scale * self.bxdf.rho2(samples);
    }

    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        return self.scale * self.bxdf.f(wo, wi);
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        sample: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        if let Some((spec, wi, pdf, t)) = self.bxdf.sample_f(wo, sample) {
            return Some((self.scale * spec, wi, pdf, t));
        }
        return None;
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        return self.bxdf.pdf(wo, wi);
    }

    fn get_type(&self) -> BxDFType {
        return self.bxdf.get_type();
    }

    fn to_string(&self) -> String {
        return format!("ScaledBxDF {:?}", self.get_type());
    }
}
