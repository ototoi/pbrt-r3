use crate::core::geometry::*;
use crate::core::material::*;
use crate::core::pbrt::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;

pub struct SpecularReflection {
    r: Spectrum,
    fresnel: Box<dyn Fresnel>,
}

impl SpecularReflection {
    pub fn new(r: &Spectrum, fresnel: Box<dyn Fresnel>) -> Self {
        return SpecularReflection {
            r: *r,
            fresnel: fresnel,
        };
    }
}

impl BxDF for SpecularReflection {
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        return Spectrum::zero();
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        _u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        let wi = Vector3f::new(-wo.x, -wo.y, wo.z);
        let pdf = 1.0;
        let cos_i = cos_theta(&wi);
        let fresnel = self.fresnel.as_ref();
        let spc = fresnel.evaluate(cos_i);
        let spc = (spc * self.r) / abs_cos_theta(&wi);
        return Some((spc, wi, pdf, 0));
    }

    fn pdf(&self, _wo: &Vector3f, _wi: &Vector3f) -> Float {
        return 0.0;
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_REFLECTION | BSDF_SPECULAR;
    }

    fn to_string(&self) -> String {
        return format!("[ SpecularReflection R: {:?} ]", self.r);
    }
}

pub struct SpecularTransmission {
    t: Spectrum,
    eta_a: Float,
    eta_b: Float,
    fresnel: FresnelDielectric,
    mode: TransportMode,
}

impl SpecularTransmission {
    pub fn new(t: &Spectrum, eta_a: Float, eta_b: Float, mode: TransportMode) -> Self {
        SpecularTransmission {
            t: *t,
            eta_a,
            eta_b,
            fresnel: FresnelDielectric::new(eta_a, eta_b),
            mode,
        }
    }
}

impl BxDF for SpecularTransmission {
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        return Spectrum::zero();
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        _u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        // Figure out which $\eta$ is incident and which is transmitted
        let entering = cos_theta(wo) > 0.0;
        let eta_a = self.eta_a;
        let eta_b = self.eta_b;
        let eta_i = if entering { eta_a } else { eta_b };
        let eta_t = if entering { eta_b } else { eta_a };

        // Compute ray direction for specular transmission
        if let Some(wi) = refract(
            wo,
            &face_forward(&Normal3f::new(0.0, 0.0, 1.0), wo),
            eta_i / eta_t,
        ) {
            let pdf = 1.0;
            let mut ft = self.t * (Spectrum::one() - self.fresnel.evaluate(cos_theta(&wi)));
            // Account for non-symmetry with transmission to different medium
            if self.mode == TransportMode::Radiance {
                ft *= (eta_i * eta_i) / (eta_t * eta_t);
            }
            let f = ft / abs_cos_theta(&wi);
            return Some((f, wi, pdf, 0));
        } else {
            return None;
        }
    }
    fn pdf(&self, _wo: &Vector3f, _wi: &Vector3f) -> Float {
        return 0.0;
    }
    fn get_type(&self) -> BxDFType {
        return BSDF_TRANSMISSION | BSDF_SPECULAR;
    }
    fn to_string(&self) -> String {
        return format!(
            "[ SpecularTransmission T: {:?} etaA: {} etaB: {} fresnel: {} mode: {:?}]",
            self.t,
            self.eta_a,
            self.eta_b,
            self.fresnel.to_string(),
            self.mode
        );
    }
}
