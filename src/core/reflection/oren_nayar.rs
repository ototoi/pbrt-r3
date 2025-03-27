use crate::core::base::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;

pub struct OrenNayar {
    pub r: Spectrum,
    pub a: Float,
    pub b: Float,
}

#[inline]
fn radians(x: Float) -> Float {
    return x * (PI / 180.0);
}

impl OrenNayar {
    pub fn new(r: &Spectrum, sigma: Float) -> Self {
        let sigma = radians(sigma);
        let sigma2 = sigma * sigma;
        let a = 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33)));
        let b = 0.45 * sigma2 / (sigma2 + 0.09);
        OrenNayar { r: *r, a, b }
    }
}

impl BxDF for OrenNayar {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let sin_theta_i = sin_theta(wi);
        let sin_theta_o = sin_theta(wo);
        // Compute cosine term of Oren-Nayar model
        let mut max_cos = 0.0;
        if sin_theta_i > 1e-4 && sin_theta_o > 1e-4 {
            let sin_phi_i = sin_phi(wi);
            let cos_phi_i = cos_phi(wi);
            let sin_phi_o = sin_phi(wo);
            let cos_phi_o = cos_phi(wo);
            let d_cos = cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;
            max_cos = d_cos.max(0.0);
        }

        // Compute sine and tangent terms of Oren-Nayar model
        let (sin_alpha, tan_beta) = if abs_cos_theta(wi) > abs_cos_theta(wo) {
            (sin_theta_o, sin_theta_i / abs_cos_theta(wi))
        } else {
            (sin_theta_i, sin_theta_o / abs_cos_theta(wo))
        };

        return self.r * INV_PI * (self.a + self.b * max_cos * sin_alpha * tan_beta);
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        self.sample_f_default(wo, u)
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        self.pdf_default(wo, wi)
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_REFLECTION | BSDF_DIFFUSE;
    }

    fn to_string(&self) -> String {
        return format!("OrenNayar {:?}", self.get_type());
    }
}
