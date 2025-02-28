use crate::core::pbrt::*;

pub type BxDFType = u32;

pub const BSDF_REFLECTION: BxDFType = 1 << 0; //1
pub const BSDF_TRANSMISSION: BxDFType = 1 << 1; //2
pub const BSDF_DIFFUSE: BxDFType = 1 << 2; //4
pub const BSDF_GLOSSY: BxDFType = 1 << 3; //8
pub const BSDF_SPECULAR: BxDFType = 1 << 4; //16
pub const BSDF_ALL: BxDFType =
    BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR;

pub trait BxDF {
    fn matches_flags(&self, t: BxDFType) -> bool {
        let tp = self.get_type();
        return (tp & t) == tp;
    }

    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum;

    fn sample_f(
        &self,
        wo: &Vector3f,
        sample: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)>;

    fn rho(&self, wo: &Vector3f, samples: &[Point2f]) -> Spectrum {
        let r = samples
            .iter()
            .map(|u| -> Spectrum {
                if let Some((f, wi, pdf, _)) = self.sample_f(wo, u) {
                    if pdf > 0.0 {
                        return f * (abs_cos_theta(&wi) / pdf);
                    }
                }
                return Spectrum::zero();
            })
            .fold(Spectrum::zero(), |a, b| a + b)
            * (1.0 / (samples.len() as Float));
        return r;
    }

    fn rho2(&self, samples: &[(Point2f, Point2f)]) -> Spectrum {
        let r = samples
            .iter()
            .map(|(u1, u2)| -> Spectrum {
                let wo = uniform_sample_hemisphere(u1);
                let pdfo = uniform_hemisphere_pdf();
                if let Some((f, wi, pdfi, _)) = self.sample_f(&wo, u2) {
                    if pdfi > 0.0 {
                        return f * (abs_cos_theta(&wi) * abs_cos_theta(&wo) / (pdfo * pdfi));
                    }
                }
                return Spectrum::zero();
            })
            .fold(Spectrum::zero(), |a, b| a + b)
            * (1.0 / (PI * samples.len() as Float));
        return r;
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float;

    //fn to_string(&self) -> String;
    fn to_string(&self) -> String {
        return format!("{:?}", self.get_type());
    }
    fn get_type(&self) -> BxDFType;

    //default implementations

    fn sample_f_default(
        &self,
        wo: &Vector3f,
        sample: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        let mut wi = cosine_sample_hemisphere(sample);
        if wo.z < 0.0 {
            wi.z *= -1.0;
        }
        let pdf = self.pdf(wo, &wi);
        let spc = self.f(wo, &wi);
        return Some((spc, wi, pdf, 0));
    }

    fn pdf_default(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        return if same_hemisphere(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0
        };
    }
}
