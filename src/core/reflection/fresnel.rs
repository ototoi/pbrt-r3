use crate::core::geometry::*;
use crate::core::material::*;
use crate::core::pbrt::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;

use std::mem::swap;

pub trait Fresnel {
    fn evaluate(&self, cos_i: Float) -> Spectrum;
    fn to_string(&self) -> String {
        return "[Fresnel unknown]".to_string();
    }
}

pub fn fr_dielectric(cos_theta_i: Float, eta_i: Float, eta_t: Float) -> Float {
    let mut cos_theta_i = Float::clamp(cos_theta_i, -1.0, 1.0);
    let mut eta_i = eta_i;
    let mut eta_t = eta_t;
    // Potentially swap indices of refraction
    let entering = cos_theta_i > 0.0;
    if !entering {
        swap(&mut eta_i, &mut eta_t);
        cos_theta_i = Float::abs(cos_theta_i);
    }

    // Compute _cosTheta_t_ using Snell's law
    let sin_theta_i = Float::sqrt(Float::max(0.0, 1.0 - cos_theta_i * cos_theta_i));
    let sin_theta_t = eta_i / eta_t * sin_theta_i;

    // Handle total internal reflection
    if sin_theta_t >= 1.0 {
        return 1.0;
    }
    let cos_theta_t = Float::sqrt(Float::max(0.0, 1.0 - sin_theta_t * sin_theta_t));
    let rparl = ((eta_t * cos_theta_i) - (eta_i * cos_theta_t))
        / ((eta_t * cos_theta_i) + (eta_i * cos_theta_t));
    let rperp = ((eta_i * cos_theta_i) - (eta_t * cos_theta_t))
        / ((eta_i * cos_theta_i) + (eta_t * cos_theta_t));

    return (rparl * rparl + rperp * rperp) / 2.0;
}

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
pub fn fr_conductor(
    cos_theta_i: Float,
    eta_i: &Spectrum,
    eta_t: &Spectrum,
    k: &Spectrum,
) -> Spectrum {
    let cos_theta_i = Float::clamp(cos_theta_i, -1.0, 1.0);
    let eta = *eta_t / *eta_i;
    let etak = *k / *eta_i;

    let cos_theta_i2 = cos_theta_i * cos_theta_i;
    let sin_theta_i2 = 1.0 - cos_theta_i2;
    let sin_theta_i2_2 = sin_theta_i2 * sin_theta_i2;
    let eta2 = eta * eta;
    let etak2 = etak * etak;

    let cos_theta_i2 = Spectrum::from(cos_theta_i2);
    let sin_theta_i2 = Spectrum::from(sin_theta_i2);
    let sin_theta_i2_2 = Spectrum::from(sin_theta_i2_2);

    let t0 = eta2 - etak2 - sin_theta_i2;
    let a2plusb2 = Spectrum::sqrt(&(t0 * t0 + eta2 * etak2 * 4.0));
    let t1 = a2plusb2 + cos_theta_i2;
    let a = Spectrum::sqrt(&((a2plusb2 + t0) * 0.5));
    let t2 = Spectrum::from(cos_theta_i) * a * 2.0;
    let rs = (t1 - t2) / (t1 + t2);

    let t3 = cos_theta_i2 * a2plusb2 + sin_theta_i2_2;
    let t4 = t2 * sin_theta_i2;
    let rp = rs * (t3 - t4) / (t3 + t4);

    // pbrt-r3:
    assert!(rp.is_valid());
    assert!(rs.is_valid());
    // pbrt-r3:

    return (rp + rs) * 0.5;
}

pub struct FresnelConductor {
    eta_i: Spectrum,
    eta_t: Spectrum,
    k: Spectrum,
}

impl FresnelConductor {
    pub fn new(eta_i: &Spectrum, eta_t: &Spectrum, k: &Spectrum) -> Self {
        FresnelConductor {
            eta_i: *eta_i,
            eta_t: *eta_t,
            k: *k,
        }
    }
}

impl Fresnel for FresnelConductor {
    fn evaluate(&self, cos_i: Float) -> Spectrum {
        return fr_conductor(Float::abs(cos_i), &self.eta_i, &self.eta_t, &self.k);
    }
}

pub struct FresnelDielectric {
    eta_i: Float,
    eta_t: Float,
}

impl FresnelDielectric {
    pub fn new(eta_i: Float, eta_t: Float) -> Self {
        FresnelDielectric { eta_i, eta_t }
    }
}

impl Fresnel for FresnelDielectric {
    fn evaluate(&self, cos_i: Float) -> Spectrum {
        let r = fr_dielectric(cos_i, self.eta_i, self.eta_t);
        return Spectrum::from(r);
    }
}

pub struct FresnelSpecular {
    r: Spectrum,
    t: Spectrum,
    eta_a: Float,
    eta_b: Float,
    mode: TransportMode,
}

impl FresnelSpecular {
    pub fn new(
        r: &Spectrum,
        t: &Spectrum,
        eta_a: Float,
        eta_b: Float,
        mode: TransportMode,
    ) -> Self {
        return FresnelSpecular {
            r: *r,
            t: *t,
            eta_a,
            eta_b,
            mode,
        };
    }
}

impl BxDF for FresnelSpecular {
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        return Spectrum::zero();
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        let eta_a = self.eta_a;
        let eta_b = self.eta_b;
        let f = fr_dielectric(cos_theta(wo), eta_a, eta_b);
        if u[0] < f {
            // Compute specular reflection for _FresnelSpecular_

            // Compute perfect specular reflection direction
            let wi = Vector3f::new(-wo.x, -wo.y, wo.z);
            let pdf = f;
            let sampled_type = BSDF_SPECULAR | BSDF_REFLECTION;
            let k = f / abs_cos_theta(&wi);
            let spec = self.r * k;
            return Some((spec, wi, pdf, sampled_type));
        } else {
            // Compute specular transmission for _FresnelSpecular_

            // Figure out which $\eta$ is incident and which is transmitted
            let entering = cos_theta(wo) > 0.0;
            let (eta_i, eta_t) = if entering {
                (eta_a, eta_b)
            } else {
                (eta_b, eta_a)
            };

            // Compute ray direction for specular transmission
            if let Some(wi) = refract(
                wo,
                &face_forward(&Normal3f::new(0.0, 0.0, 1.0), wo),
                eta_i / eta_t,
            ) {
                let mut ft = self.t * (1.0 - f);
                // Account for non-symmetry with transmission to different medium
                if self.mode == TransportMode::Radiance {
                    ft *= (eta_i * eta_i) / (eta_t * eta_t);
                }
                let sampled_type = BSDF_SPECULAR | BSDF_TRANSMISSION;
                let pdf = 1.0 - f;
                let spec = ft / abs_cos_theta(&wi);
                return Some((spec, wi, pdf, sampled_type));
            } else {
                return None;
            }
        }
    }

    fn pdf(&self, _wo: &Vector3f, _wi: &Vector3f) -> Float {
        return 0.0;
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR;
    }

    fn to_string(&self) -> String {
        return format!("[ FresnelSpecular T: {:?} ]", self.t);
    }
}

pub struct FresnelNoOp {}

impl FresnelNoOp {
    pub fn new() -> Self {
        FresnelNoOp {}
    }
}

impl Fresnel for FresnelNoOp {
    fn evaluate(&self, _cos_i: Float) -> Spectrum {
        return Spectrum::one();
    }
}

impl Default for FresnelNoOp {
    fn default() -> Self {
        Self::new()
    }
}
