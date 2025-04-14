use crate::core::base::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;

pub struct DiffuseBxDF {
    pub r: Spectrum,
}

impl DiffuseBxDF {
    pub fn new(r: &Spectrum) -> Self {
        DiffuseBxDF { r: *r }
    }
}

impl BxDF for DiffuseBxDF {
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        return self.r * INV_PI;
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        sample: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        self.sample_f_default(wo, sample)
    }

    fn rho(&self, _wo: &Vector3f, _samples: &[Point2f]) -> Spectrum {
        return self.r;
    }

    fn rho2(&self, _samples: &[(Point2f, Point2f)]) -> Spectrum {
        return self.r;
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        self.pdf_default(wo, wi)
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_REFLECTION | BSDF_DIFFUSE;
    }
    fn to_string(&self) -> String {
        return format!("[ DiffuseBxDF R: {:?} ]", self.r);
    }
}