use crate::core::pbrt::*;

pub struct LambertianReflection {
    pub r: Spectrum,
}

impl LambertianReflection {
    pub fn new(r: &Spectrum) -> Self {
        LambertianReflection { r: *r }
    }
}

impl BxDF for LambertianReflection {
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
        return format!("[ LambertianReflection R: {:?} ]", self.r);
    }
}

pub struct LambertianTransmission {
    pub t: Spectrum,
}

impl LambertianTransmission {
    pub fn new(t: &Spectrum) -> Self {
        LambertianTransmission { t: *t }
    }
}

impl BxDF for LambertianTransmission {
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        return self.t * INV_PI;
    }
    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        let mut wi = cosine_sample_hemisphere(u);
        if wo.z > 0.0 {
            wi.z *= -1.0;
        }
        let pdf = self.pdf(wo, &wi);
        if pdf > 0.0 {
            let f = self.f(wo, &wi);
            return Some((f, wi, pdf, 0));
        } else {
            return None;
        }
    }
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !same_hemisphere(wo, wi) {
            return abs_cos_theta(wi);
        } else {
            return 0.0;
        }
    }
    fn rho(&self, _wo: &Vector3f, _samples: &[Point2f]) -> Spectrum {
        return self.t;
    }
    fn rho2(&self, _samples: &[(Point2f, Point2f)]) -> Spectrum {
        return self.t;
    }
    fn get_type(&self) -> BxDFType {
        return BSDF_TRANSMISSION | BSDF_DIFFUSE;
    }
    fn to_string(&self) -> String {
        return format!("[ LambertianTransmission T: {:?} ]", self.t);
    }
}
