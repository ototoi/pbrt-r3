// Imported from bsdfs.cpp

use pbrt_r3::core::pbrt::*;
use pbrt_r3::shapes::*;

use std::sync::Arc;

static CHI2_THETA_RES: usize = 10;
static CHI2_PHI_RES: usize = 2 * CHI2_THETA_RES;
static CHI2_SAMPLECOUNT: usize = 1000000;
static CHI2_RUNS: usize = 5;

/* Define an recursive lambda function for integration over subintervals */
fn integrate(
    a: Float,
    b: Float,
    c: Float,
    fa: Float,
    fb: Float,
    fc: Float,
    f: &dyn Fn(Float) -> Float,
    i: Float,
    eps: Float,
    depth: i32,
) -> Float {
    /* Evaluate the function at two intermediate points */
    let d = 0.5 * (a + b);
    let e = 0.5 * (b + c);
    let fd = f(d);
    let fe = f(e);
    /* Simpson integration over each subinterval */
    let h = c - a;
    let i0 = (1.0 / 12.0) * h * (fa + 4.0 * fd + fb);
    let i1 = (1.0 / 12.0) * h * (fb + 4.0 * fe + fc);
    let ip = i0 + i1;
    /* Stopping criterion from J.N. Lyness (1969)
    "Notes on the adaptive Simpson quadrature routine" */
    if depth <= 0 || Float::abs(ip - i) < 15.0 * eps {
        // Richardson extrapolation
        return ip + (ip - i) / 15.0;
    }
    return integrate(a, d, b, fa, fd, fb, f, i0, 0.5 * eps, depth - 1)
        + integrate(b, e, c, fb, fe, fc, f, i1, 0.5 * eps, depth - 1);
}

fn adaptive_simpson(
    f: &dyn Fn(Float) -> Float,
    x0: Float,
    x1: Float,
    eps: Float,
    depth: i32,
) -> Float {
    let a = x0;
    let b = 0.5 * (x0 + x1);
    let c = x1;
    let fa = f(a);
    let fb = f(b);
    let fc = f(c);
    let i = (1.0 / 6.0) * (c - a) * (fa + 4.0 * fb + fc);
    return integrate(a, b, c, fa, fb, fc, f, i, eps, depth);
}

/// Adaptive Simpson integration over an 1D interval
fn adaptive_simpson_2d(
    f: &dyn Fn(Float, Float) -> Float,
    x0: Float,
    y0: Float,
    x1: Float,
    y1: Float,
    eps: Float,
    depth: i32,
) -> Float {
    /* Lambda function that integrates over the X axis */
    let integrate_x = |y: Float| -> Float {
        return adaptive_simpson(&|x| f(x, y), x0, x1, eps, depth);
    };
    let value = adaptive_simpson(&integrate_x, y0, y1, eps, depth);
    return value;
}

/// Generate a histogram of the BSDF density function via MC sampling
fn frequency_table(
    bsdf: &BSDF,
    wo: &Vector3f,
    rng: &mut RNG,
    sample_count: usize,
    theta_res: usize,
    phi_res: usize,
) -> Vec<Float> {
    let mut table = vec![0.0; (theta_res * phi_res) as usize];
    let factor_theta = theta_res as Float / PI;
    let factor_phi = phi_res as Float / (2.0 * PI);
    for _ in 0..sample_count {
        let sample = Point2f::new(rng.uniform_float(), rng.uniform_float());
        if let Some((f, wi, pdf, flags)) = bsdf.sample_f(wo, &sample, BSDF_ALL) {
            if f.is_black() || pdf == 0.0 || (flags & BSDF_SPECULAR) != 0 {
                continue;
            }
            let wi_l = bsdf.world_to_local(&wi);
            let x = Float::clamp(wi_l.z, -1.0, 1.0).acos() * factor_theta;
            let mut y = Float::atan2(wi_l.y, wi_l.x) * factor_phi;
            if y < 0.0 {
                y += 2.0 * PI;
            }
            let theta_bin = usize::clamp(Float::floor(x) as usize, 0, theta_res - 1);
            let phi_bin = usize::clamp(Float::floor(y) as usize, 0, phi_res - 1);

            let theta_bin = theta_bin.min(theta_res - 1);
            let phi_bin = phi_bin.min(phi_res - 1);
            let index = theta_bin * phi_res + phi_bin;
            table[index] += 1.0;
        }
    }
    return table;
}

// Numerically integrate the probability density function over rectangles in
// spherical coordinates.
fn integrate_frequency_table(
    bsdf: &BSDF,
    wo: &Vector3f,
    sample_count: usize,
    theta_res: usize,
    phi_res: usize,
) -> Vec<Float> {
    let mut table = vec![0.0; (theta_res * phi_res) as usize];
    let factor_theta = PI / theta_res as Float;
    let factor_phi = (2.0 * PI) / phi_res as Float;
    let mut index = 0;
    for i in 0..theta_res {
        for j in 0..phi_res {
            table[index] = sample_count as Float
                * adaptive_simpson_2d(
                    &|theta, phi| -> Float {
                        let cos_theta = theta.cos();
                        let sin_theta = theta.sin();
                        let cos_phi = phi.cos();
                        let sin_phi = phi.sin();
                        let wi = spherical_direction(
                            sin_theta * cos_phi,
                            sin_theta * sin_phi,
                            cos_theta,
                        );
                        let wi_w = bsdf.local_to_world(&wi);
                        let pdf = bsdf.pdf(wo, &wi_w, BSDF_ALL);
                        return pdf * sin_theta;
                    },
                    i as Float * factor_theta,
                    j as Float * factor_phi,
                    (i + 1) as Float * factor_theta,
                    (j + 1) as Float * factor_phi,
                    1e-6,
                    6,
                );
            index += 1;
        }
    }
    return table;
}

fn chi2_test(
    frequencies: &[Float],
    exp_frequencies: &[Float],
    theta_res: usize,
    phi_res: usize,
    sample_count: usize,
) -> Result<(), String> {
    let mut sum = 0.0;
    for i in 0..frequencies.len() {
        let diff = frequencies[i] - exp_frequencies[i];
        sum += diff * diff / exp_frequencies[i];
    }
    Ok(())
}

fn test_bsdf(create_bsdf: &(dyn Fn(&mut BSDF, &mut MemoryArena)), description: &str) {
    let mut arena = MemoryArena::new();

    let theta_res = CHI2_THETA_RES;
    let phi_res = CHI2_PHI_RES;
    let sample_count = CHI2_SAMPLECOUNT;

    let mut rng = RNG::new();

    let t = Transform::rotate_x(-90.0);
    let t_inv = t.inverse();

    let mut bsdf = None;
    {
        let reverse_orientation = false;
        let disk = Disk::new(&t, &t_inv, reverse_orientation, 0.0, 1.0, 0.0, 360.0);
        let origin = Point3f::new(0.1, 1.0, 0.0);
        let direction = Vector3f::new(0.0, -1.0, 0.0);
        let r = Ray::new(&origin, &direction, Float::INFINITY, 0.0);
        if let Some((_t_hit, isect)) = disk.intersect(&r) {
            let mut b = BSDF::new(&isect, 1.0);
            create_bsdf(&mut b, &mut arena);
            bsdf = Some(b);
        }
    }
    if let Some(bsdf) = bsdf.as_ref() {
        for _ in 0..CHI2_RUNS {
            /* Randomly pick an outgoing direction on the hemisphere */
            let sample = Point2f::new(rng.uniform_float(), rng.uniform_float());
            let wo_l = cosine_sample_hemisphere(&sample);
            let wo = bsdf.world_to_local(&wo_l);
            let frequencies =
                frequency_table(bsdf, &wo, &mut rng, sample_count, theta_res, phi_res);
            //
            let exp_frequencies =
                integrate_frequency_table(bsdf, &wo, sample_count, theta_res, phi_res);
            if let Err(e) = chi2_test(
                &frequencies,
                &exp_frequencies,
                theta_res,
                phi_res,
                sample_count,
            ) {
                panic!("{}: {}", description, e);
            }
        }
    }
}

fn create_lambertian(bsdf: &mut BSDF, _arena: &mut MemoryArena) {
    let kd = Spectrum::from(1.0);
    let r: Arc<dyn BxDF> = Arc::new(LambertianReflection::new(&kd));
    bsdf.add(&r);
}

#[test]
fn bsdf_sampling_lambertian() {
    test_bsdf(&create_lambertian, "Lambertian");
}
