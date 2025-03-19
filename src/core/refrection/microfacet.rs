use crate::core::distribution::*;
use crate::core::geometry::*;
use crate::core::material::*;
use crate::core::pbrt::*;
use crate::core::refrection::*;
use crate::core::spectrum::*;

pub struct MicrofacetReflection {
    r: Spectrum,
    distribution: Box<dyn MicrofacetDistribution>,
    fresnel: Box<dyn Fresnel>,
}

impl MicrofacetReflection {
    pub fn new(
        r: &Spectrum,
        distribution: Box<dyn MicrofacetDistribution>,
        fresnel: Box<dyn Fresnel>,
    ) -> Self {
        MicrofacetReflection {
            r: *r,
            distribution: distribution,
            fresnel: fresnel,
        }
    }
}

impl BxDF for MicrofacetReflection {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let cos_theta_o = abs_cos_theta(wo);
        let cos_theta_i = abs_cos_theta(wi);
        let mut wh = *wi + *wo;
        // Handle degenerate cases for microfacet reflection
        if cos_theta_i == 0.0 || cos_theta_o == 0.0 {
            return Spectrum::zero();
        }
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::zero();
        }
        wh = wh.normalize();

        // pbrt-r3:
        debug_assert!(Float::is_finite(wh.x));
        debug_assert!(Float::is_finite(wh.y));
        debug_assert!(Float::is_finite(wh.z));
        // pbrt-r3:

        // For the Fresnel call, make sure that wh is in the same hemisphere
        // as the surface normal, so that TIR is handled correctly.
        let fresnel = self.fresnel.as_ref();
        let distribution = self.distribution.as_ref();
        let f = fresnel.evaluate(Vector3f::dot(
            wi,
            &face_forward(&wh, &Vector3f::new(0.0, 0.0, 1.0)),
        ));
        let r = self.r;
        return r
            * f
            * (distribution.d(&wh) * distribution.g(wo, wi) / (4.0 * cos_theta_i * cos_theta_o));
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        // Sample microfacet orientation $\wh$ and reflected direction $\wi$
        if wo.z == 0.0 {
            return None;
        }
        let distribution = self.distribution.as_ref();
        let wh = distribution.sample_wh(wo, u);
        if Vector3f::dot(wo, &wh) < 0.0 {
            return None; // Should be rare
        }
        let wi = reflect(wo, &wh);
        if !same_hemisphere(wo, &wi) {
            return None;
        }
        // Compute PDF of _wi_ for microfacet reflection
        let pdf = distribution.pdf(wo, &wh) / (4.0 * Vector3f::dot(wo, &wh));
        if pdf == 0.0 {
            return None;
        } else {
            let spec = self.f(wo, &wi);
            return Some((spec, wi, pdf, 0));
        }
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !same_hemisphere(wo, wi) {
            return 0.0;
        }
        let wh = (*wo + *wi).normalize();
        // pbrt-r3:
        if Vector3f::dot(wo, &wh) < 0.0 {
            return 0.0; // Should be rare
        }
        // pbrt-r3:

        let distribution = self.distribution.as_ref();
        return distribution.pdf(wo, &wh) / (4.0 * Vector3f::dot(wo, &wh));
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_REFLECTION | BSDF_GLOSSY;
    }

    fn to_string(&self) -> String {
        return format!("[ MicrofacetReflection R: {:?} ]", self.r);
    }
}

pub struct MicrofacetTransmission {
    t: Spectrum,
    distribution: Box<dyn MicrofacetDistribution>,
    eta_a: Float,
    eta_b: Float,
    fresnel: FresnelDielectric,
    mode: TransportMode,
}

impl MicrofacetTransmission {
    pub fn new(
        t: &Spectrum,
        distribution: Box<dyn MicrofacetDistribution>,
        eta_a: Float,
        eta_b: Float,
        mode: TransportMode,
    ) -> Self {
        MicrofacetTransmission {
            t: *t,
            distribution: distribution,
            eta_a,
            eta_b,
            fresnel: FresnelDielectric::new(eta_a, eta_b),
            mode,
        }
    }
}

impl BxDF for MicrofacetTransmission {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        if same_hemisphere(wo, wi) {
            return Spectrum::zero();
        }
        let cos_theta_o = cos_theta(wo);
        let cos_theta_i = cos_theta(wi);
        if cos_theta_i == 0.0 || cos_theta_o == 0.0 {
            return Spectrum::zero();
        }
        let eta_a = self.eta_a;
        let eta_b = self.eta_b;

        // Compute $\wh$ from $\wo$ and $\wi$ for microfacet transmission
        let eta = if cos_theta_o > 0.0 {
            eta_b / eta_a
        } else {
            eta_a / eta_b
        };
        let mut wh = (*wo + (*wi * eta)).normalize();
        if wh.z < 0.0 {
            wh = -wh;
        }
        let wo_wh = Vector3f::dot(wo, &wh);
        let wi_wh = Vector3f::dot(wi, &wh);

        if wo_wh * wi_wh > 0.0 {
            return Spectrum::zero();
        }

        let f = self.fresnel.evaluate(Vector3f::dot(wo, &wh));
        let sqrt_denom = wo_wh + eta * wi_wh;
        let factor = if self.mode == TransportMode::Radiance {
            1.0 / eta
        } else {
            1.0
        };
        let distribution = self.distribution.as_ref();
        let d = Float::abs(
            distribution.d(&wh)
                * distribution.g(wo, wi)
                * eta
                * eta
                * Vector3f::abs_dot(wi, &wh)
                * Vector3f::abs_dot(wo, &wh)
                * factor
                * factor
                / (cos_theta_i * cos_theta_o * sqrt_denom * sqrt_denom),
        );
        return ((Spectrum::one() - f) * self.t) * d;
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        if wo.z == 0.0 {
            return None;
        }
        let distribution = self.distribution.as_ref();
        let wh = distribution.sample_wh(wo, u);
        if Vector3f::dot(wo, &wh) < 0.0 {
            return None;
        }
        let cos_theta_o = cos_theta(wo);
        let eta_a = self.eta_a;
        let eta_b = self.eta_b;
        let eta = if cos_theta_o > 0.0 {
            eta_a / eta_b
        } else {
            eta_b / eta_a
        };

        if let Some(wi) = refract(wo, &wh, eta) {
            let pdf = self.pdf(wo, &wi);
            if pdf > 0.0 {
                let f = self.f(wo, &wi);
                return Some((f, wi, pdf, 0));
            }
        }
        return None;
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if same_hemisphere(wo, wi) {
            return 0.0;
        }

        let cos_theta_o = cos_theta(wo);
        let eta_a = self.eta_a;
        let eta_b = self.eta_b;
        let eta = if cos_theta_o > 0.0 {
            eta_b / eta_a
        } else {
            eta_a / eta_b
        };

        let wh = (*wo + (*wi * eta)).normalize();

        let wo_wh = Vector3f::dot(wo, &wh);
        let wi_wh = Vector3f::dot(wi, &wh);
        if wo_wh * wi_wh > 0.0 {
            return 0.0;
        }

        // Compute change of variables _dwh\_dwi_ for microfacet transmission
        let sqrt_denom = wo_wh + eta * wi_wh;
        let dwh_dwi = Float::abs((eta * eta * wi_wh) / (sqrt_denom * sqrt_denom));

        let distribution = self.distribution.as_ref();
        return distribution.pdf(wo, &wh) * dwh_dwi;
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_TRANSMISSION | BSDF_GLOSSY;
    }

    fn to_string(&self) -> String {
        return format!("[ MicrofacetTransmission T: {:?} ]", self.t);
    }
}
