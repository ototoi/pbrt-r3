use super::microfacet::MicrofacetDistribution;
use crate::core::pbrt::*;

fn trowbridge_reitz_sample_11(cos_theta: Float, u1: Float, u2: Float) -> (Float, Float) {
    /* Special case (normal incidence) */
    if cos_theta > 0.9999 {
        let r = Float::sqrt(u1 / (1.0 - u1));
        let phi = 2.0 * PI * u2; //6.28318530718
        let sin_phi = Float::sin(phi);
        let cos_phi = Float::cos(phi);
        let slopex = r * cos_phi;
        let slopey = r * sin_phi;
        return (slopex, slopey);
    }

    let sin_theta = Float::sqrt(Float::max(0.0, 1.0 - cos_theta * cos_theta));
    let tan_theta = sin_theta / cos_theta;
    //let cot_theta_i = 1.0 / tan_theta;

    let a = 1.0 / tan_theta;
    let g1 = 2.0 / (1.0 + Float::sqrt(1.0 + 1.0 / (a * a)));

    // sample slope_x
    let a = 2.0 * u1 / g1 - 1.0;
    let tmp = Float::min(1e10, 1.0 / (a * a - 1.0));
    let b = tan_theta;
    let d = Float::sqrt(Float::max(b * b * tmp * tmp - (a * a - b * b) * tmp, 0.0));
    let slope_x_1 = b * tmp - d;
    let slope_x_2 = b * tmp + d;
    let slope_x = if a < 0.0 || slope_x_2 > 1.0 / tan_theta {
        slope_x_1
    } else {
        slope_x_2
    };

    // sample slope_y
    let (s, u2) = if u2 > 0.5 {
        (1.0, 2.0 * (u2 - 0.5))
    } else {
        (-1.0, 2.0 * (0.5 - u2))
    };
    let z = (u2 * (u2 * (u2 * 0.27385 - 0.73369) + 0.46341))
        / (u2 * (u2 * (u2 * 0.093073 + 0.309420) - 1.000000) + 0.597999);
    let slope_y = s * z * Float::sqrt(1.0 + slope_x * slope_x);

    assert!(!Float::is_infinite(slope_x));
    assert!(!Float::is_nan(slope_x));
    assert!(!Float::is_infinite(slope_y));
    assert!(!Float::is_nan(slope_y));

    return (slope_x, slope_y);
}

pub fn trowbridge_reitz_sample(
    wi: &Vector3f,
    alpha_x: Float,
    alpha_y: Float,
    u1: Float,
    u2: Float,
) -> Vector3f {
    // 1. stretch wi
    let wi_stretched = Vector3f::new(alpha_x * wi.x, alpha_y * wi.y, wi.z).normalize();

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    let (mut slope_x, mut slope_y) = trowbridge_reitz_sample_11(cos_theta(&wi_stretched), u1, u2);

    // 3. rotate
    let tmp = cos_phi(&wi_stretched) * slope_x - sin_phi(&wi_stretched) * slope_y;
    slope_y = sin_phi(&wi_stretched) * slope_x + cos_phi(&wi_stretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x *= alpha_x;
    slope_y *= alpha_y;

    assert!(!Float::is_infinite(slope_x));
    assert!(!Float::is_nan(slope_x));
    assert!(!Float::is_infinite(slope_y));
    assert!(!Float::is_nan(slope_y));

    // 5. compute normal
    return Vector3f::new(-slope_x, -slope_y, 1.0).normalize();
}

pub struct TrowbridgeReitzDistribution {
    alphax: Float,
    alphay: Float,
    samplevis: bool,
}

impl TrowbridgeReitzDistribution {
    pub fn new(alphax: Float, alphay: Float, samplevis: bool) -> Self {
        let alphax = Float::max(0.001, alphax);
        let alphay = Float::max(0.001, alphay);
        TrowbridgeReitzDistribution {
            alphax,
            alphay,
            samplevis,
        }
    }

    pub fn roughness_to_alpha(roughness: Float) -> Float {
        let roughness = Float::max(roughness, 1e-3);
        let x = Float::ln(roughness);
        return 1.62142
            + 0.819955 * x
            + 0.1734 * x * x
            + 0.0171201 * x * x * x
            + 0.000640711 * x * x * x * x;
    }
}

fn sample_wh_helper(alphax: Float, alphay: Float, u1: Float, u2: Float) -> (Float, Float, Float) {
    // pbrt-r3:
    let u1 = u1.clamp(1e-6, 1.0 - 1e-6);
    let u2 = u2.clamp(1e-6, 1.0 - 1e-6);
    //
    if alphax == alphay {
        let tan_theta_2 = alphax * alphax * u1 / (1.0 - u1);
        let phi = 2.0 * PI * u2;
        let cos_theta = 1.0 / Float::sqrt(1.0 + tan_theta_2);
        let sin_theta = Float::sqrt(Float::max(0.0, 1.0 - cos_theta * cos_theta));
        return (sin_theta, cos_theta, phi);
    } else {
        let mut phi = Float::atan(alphay / alphax * Float::tan(2.0 * PI * u2 + 0.5 * PI));
        if u2 > 0.5 {
            phi += PI;
        }
        let sin_phi = Float::sin(phi);
        let cos_phi = Float::cos(phi);
        let alphax2 = alphax * alphax;
        let alphay2 = alphay * alphay;
        let alpha2 = 1.0 / (cos_phi * cos_phi / alphax2 + sin_phi * sin_phi / alphay2);
        let tan_theta_2 = alpha2 * u1 / (1.0 - u1);
        let cos_theta = 1.0 / Float::sqrt(1.0 + tan_theta_2);
        let sin_theta = Float::sqrt(Float::max(0.0, 1.0 - cos_theta * cos_theta));
        return (sin_theta, cos_theta, phi);
    }
}

impl MicrofacetDistribution for TrowbridgeReitzDistribution {
    fn d(&self, wh: &Vector3f) -> Float {
        let alphax = self.alphax;
        let alphay = self.alphay;
        let tan_2_theta = tan_2_theta(wh);
        if Float::is_infinite(tan_2_theta) {
            return 0.0;
        }
        let cos_2_theta = cos_2_theta(wh);
        let cos_4_theta = cos_2_theta * cos_2_theta;
        let e =
            (cos_2_phi(wh) / (alphax * alphax) + sin_2_phi(wh) / (alphay * alphay)) * tan_2_theta;
        let e2 = (1.0 + e) * (1.0 + e);
        return 1.0 / (PI * alphax * alphay * cos_4_theta * e2);
    }

    fn lambda(&self, w: &Vector3f) -> Float {
        let abs_tan_theta = Float::abs(tan_theta(w));
        if Float::is_infinite(abs_tan_theta) {
            return 0.0;
        }
        let alphax = self.alphax;
        let alphay = self.alphay;
        // Compute _alpha_ for direction _w_
        let alpha = Float::sqrt(cos_2_phi(w) * alphax * alphax + sin_2_phi(w) * alphay * alphay);
        let alpha_2_tan_2_theta = (alpha * abs_tan_theta) * (alpha * abs_tan_theta);
        return (-1.0 + Float::sqrt(1.0 + alpha_2_tan_2_theta)) / 2.0;
    }

    fn sample_wh(&self, wo: &Vector3f, u: &Vector2f) -> Vector3f {
        if !self.samplevis {
            let (sin_theta, cos_theta, phi) =
                sample_wh_helper(self.alphax, self.alphay, u[0], u[1]);

            // pbrt-r3:
            assert!(!Float::is_infinite(sin_theta));
            assert!(!Float::is_nan(sin_theta));
            assert!(!Float::is_infinite(cos_theta));
            assert!(!Float::is_nan(cos_theta));
            assert!(!Float::is_infinite(phi));
            assert!(!Float::is_nan(phi));

            let mut wh = spherical_direction(sin_theta, cos_theta, phi).normalize();
            if !same_hemisphere(wo, &wh) {
                wh = -wh;
            }
            return wh;
        } else {
            // Sample visible area of normals for Beckmann distribution
            let flip = wo.z < 0.0;
            let wo = if flip { -*wo } else { *wo };
            let mut wh = trowbridge_reitz_sample(&wo, self.alphax, self.alphay, u[0], u[1]);
            if flip {
                wh = -wh;
            }
            return wh;
        }
    }

    fn pdf(&self, wo: &Vector3f, wh: &Vector3f) -> Float {
        if self.samplevis {
            return self.d(wh) * self.g1(wo) * Vector3f::abs_dot(wo, wh) / abs_cos_theta(wo);
        } else {
            return self.d(wh) * abs_cos_theta(wh);
        }
    }
}
