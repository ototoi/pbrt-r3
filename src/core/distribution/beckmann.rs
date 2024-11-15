use crate::core::pbrt::*;

fn beckmann_sample_11(cos_theta_i: Float, u1: Float, u2: Float) -> (Float, Float) {
    /* Special case (normal incidence) */
    if cos_theta_i > 0.9999 {
        let r = Float::sqrt(-Float::ln(1.0 - u1));
        let sin_phi = Float::sin(2.0 * PI * u2);
        let cos_phi = Float::cos(2.0 * PI * u2);
        let slopex = r * cos_phi;
        let slopey = r * sin_phi;
        return (slopex, slopey);
    }

    /* The original inversion routine from the paper contained
    discontinuities, which causes issues for QMC integration
    and techniques like Kelemen-style MLT. The following code
    performs a numerical inversion with better behavior */
    let sin_theta_i = Float::sqrt(Float::max(0.0, 1.0 - cos_theta_i * cos_theta_i));
    let tan_theta_i = sin_theta_i / cos_theta_i;
    let cot_theta_i = 1.0 / tan_theta_i;

    /* Search interval -- everything is parameterized
    in the Erf() domain */
    let mut a = -1.0;
    let mut c = erf(cot_theta_i);
    let sample_x = Float::max(u1, 1e-6);

    /* Start with a good initial guess */
    // Float b = (1-sample_x) * a + sample_x * c;

    /* We can do better (inverse of an approximation computed in
     * Mathematica) */
    let theta_i = Float::acos(cos_theta_i);
    let fit = 1.0 + theta_i * (-0.876 + theta_i * (0.4265 - 0.0594 * theta_i));
    let mut b = c - (1.0 + c) * Float::powf(1.0 - sample_x, fit);

    /* Normalization factor for the CDF */
    //const SQRT_PI_INV: Float = 1.0 / Float::sqrt(PI);
    const SQRT_PI_INV: Float = INV_SQRT_PI;
    let normalization =
        1.0 / (1.0 + c + SQRT_PI_INV * tan_theta_i * Float::exp(-cot_theta_i * cot_theta_i));

    for _ in 0..10 {
        /* Bisection criterion -- the oddly-looking
        Boolean expression are intentional to check
        for NaNs at little additional cost */
        if !(b >= a && b <= c) {
            b = 0.5 * (a + c);
        }

        /* Evaluate the CDF and its derivative
        (i.e. the density function) */
        let inv_erf = erf_inv(b);
        let value = normalization
            * (1.0 + b + SQRT_PI_INV * tan_theta_i * Float::exp(-inv_erf * inv_erf))
            - sample_x;
        let derivative = normalization * (1.0 - inv_erf * tan_theta_i);

        if Float::abs(value) < 1e-5 {
            break;
        }

        /* Update bisection intervals */
        if value > 0.0 {
            c = b;
        } else {
            a = b;
        }
        b -= value / derivative;
    }

    /* Now convert back into a slope value */
    let slope_x = erf_inv(b);

    /* Simulate Y component */
    let slope_y = erf_inv(2.0 * Float::max(u2, 1e-6) - 1.0);

    assert!(!Float::is_infinite(slope_x));
    assert!(!Float::is_nan(slope_x));
    assert!(!Float::is_infinite(slope_y));
    assert!(!Float::is_nan(slope_y));

    return (slope_x, slope_y);
}

pub fn beckmann_sample(
    wi: &Vector3f,
    alpha_x: Float,
    alpha_y: Float,
    u1: Float,
    u2: Float,
) -> Vector3f {
    // 1. stretch wi
    let wi_stretched = Vector3f::new(alpha_x * wi.x, alpha_y * wi.y, wi.z).normalize();

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    let (mut slope_x, mut slope_y) = beckmann_sample_11(cos_theta(&wi_stretched), u1, u2);

    // 3. rotate
    let tmp = cos_phi(&wi_stretched) * slope_x - sin_phi(&wi_stretched) * slope_y;
    slope_y = sin_phi(&wi_stretched) * slope_x + cos_phi(&wi_stretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x *= alpha_x;
    slope_y *= alpha_y;

    // 5. compute normal
    return Vector3f::new(-slope_x, -slope_y, 1.0).normalize();
}

pub struct BeckmannDistribution {}
