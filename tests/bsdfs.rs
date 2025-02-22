// Imported from bsdfs.cpp

use pbrt_r3::core::pbrt::*;
use pbrt_r3::shapes::*;

use libm::lgamma;
use std::sync::Arc;

/* The null hypothesis will be rejected when the associated
p-value is below the significance level specified here. */
static CHI2_SLEVEL: Float = 0.01;

/* Resolution of the frequency table discretization. The azimuthal
resolution is twice this value. */
static CHI2_THETA_RES: usize = 10;
static CHI2_PHI_RES: usize = 2 * CHI2_THETA_RES;

/* Number of MC samples to compute the observed frequency table */
static CHI2_SAMPLECOUNT: usize = 1000000;

/* Minimum expected bin frequency. The chi^2 test does not
work reliably when the expected frequency in a cell is
low (e.g. less than 5), because normality assumptions
break down in this case. Therefore, the implementation
will merge such low-frequency cells when they fall below
the threshold specified here. */
static CHI2_MINFREQ: Float = 5.0;

/* Each provided BSDF will be tested for a few different
incident directions. The value specified here determines
how many tests will be executed per BSDF */
static CHI2_RUNS: usize = 5;

/// Regularized lower incomplete gamma function (based on code from Cephes)
fn rl_gamma(a: f64, x: f64) -> f64 {
    const EPSILON: f64 = 0.000000000000001;
    const BIG: f64 = 4503599627370496.0;
    const BIG_INV: f64 = 2.22044604925031308085e-16;

    if a < 0.0 || x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return 0.0;
    }

    let ax = (a * f64::ln(x)) - x - lgamma(a);
    if ax < 709.78271289338399 {
        return if a < x { 1.0 } else { 0.0 };
    }

    if x <= 1.0 || x <= a {
        let mut r2 = a;
        let mut c2 = 1.0;
        let mut ans2 = 1.0;

        loop {
            r2 += 1.0;
            c2 *= x / r2;
            ans2 += c2;
            if (c2 / ans2) <= EPSILON {
                break;
            }
        }

        return f64::exp(ax) * ans2 / a;
    }

    let mut c = 0;
    let mut y = 1.0 - a;
    let mut z = x + y + 1.0;
    let mut p3 = 1.0;
    let mut q3 = x;
    let mut p2 = x + 1.0;
    let mut q2 = z * x;
    let mut ans = p2 / q2;

    let mut err;

    loop {
        c += 1;
        y += 1.0;
        z += 2.0;
        let yc = y * c as f64;
        let p = p2 * z - p3 * yc;
        let q = q2 * z - q3 * yc;

        if q != 0.0 {
            let nextans = p / q;
            err = f64::abs((ans - nextans) / nextans);
            ans = nextans;
        } else {
            err = 1.0;
        }

        // shift
        p3 = p2;
        p2 = p;
        q3 = q2;
        q2 = q;

        // normalize fraction when the numerator becomes large
        if f64::abs(p) > BIG {
            p3 *= BIG_INV;
            p2 *= BIG_INV;
            q3 *= BIG_INV;
            q2 *= BIG_INV;
        }
        if err <= EPSILON {
            break;
        }
    }

    return 1.0 - (f64::exp(ax) * ans);
}

/// Chi^2 distribution cumulative distribution function
fn chi2_cdf(x: f64, dof: i32) -> f64 {
    if dof < 1 || x < 0.0 {
        return 0.0;
    } else if dof == 2 {
        return 1.0 - f64::exp(-0.5 * x);
    } else {
        return 1.0 - rl_gamma((dof as f64) / 2.0, x / 2.0);
    }
}

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
    _theta_res: usize,
    _phi_res: usize,
    sample_count: usize,
    min_exp_frequency: Float,
    significance_level: Float,
    num_tests: usize,
) -> Result<(), String> {
    let mut cells: Vec<(Float, usize)> = exp_frequencies
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i))
        .collect();
    cells.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    assert!(cells[0].0 <= cells[1].0);

    /* Compute the Chi^2 statistic and pool cells as necessary */
    let mut pooled_frequencies = 0.0;
    let mut pooled_exp_frequencies = 0.0;
    let mut chsq = 0.0;
    let mut _pooled_cells = 0;
    let mut dof = 0;

    for c in cells.iter() {
        if exp_frequencies[c.1] == 0.0 {
            if frequencies[c.1] > sample_count as Float * 1e-5 {
                let result = format!(
                    "Encountered {} samples in a c with expected frequency 0. Rejecting the null hypothesis!",
                    frequencies[c.1]
                );
                return Err(result);
            }
        } else if exp_frequencies[c.1] < min_exp_frequency {
            /* Pool cells with low expected frequencies */
            pooled_frequencies += frequencies[c.1];
            pooled_exp_frequencies += exp_frequencies[c.1];
            _pooled_cells += 1;
        } else if pooled_exp_frequencies > 0.0 && pooled_exp_frequencies < min_exp_frequency {
            /*
            Keep on pooling cells until a sufficiently high
            expected frequency is achieved.
            */
            pooled_frequencies += frequencies[c.1];
            pooled_exp_frequencies += exp_frequencies[c.1];
            _pooled_cells += 1;
        } else {
            let diff = frequencies[c.1] - exp_frequencies[c.1];
            chsq += (diff * diff) / exp_frequencies[c.1];
            dof += 1;
        }
    }

    if pooled_exp_frequencies > 0.0 || pooled_frequencies > 0.0 {
        let diff = pooled_exp_frequencies - pooled_frequencies;
        chsq += (diff * diff) / pooled_exp_frequencies;
        dof += 1;
    }

    /* All parameters are assumed to be known, so there is no
    additional DF reduction due to model parameters */
    dof -= 1;

    if dof <= 0 {
        let result = format!("The number of degrees of freedom {} is too low!", dof);
        return Err(result);
    }

    /* Probability of obtaining a test statistic at least
    as extreme as the one observed under the assumption
    that the distributions match */
    let pval = 1.0 - chi2_cdf(chsq as f64, dof) as Float;

    /* Apply the Sidak correction term, since we'll be conducting multiple
    independent
    hypothesis tests. This accounts for the fact that the probability of a
    failure
    increases quickly when several hypothesis tests are run in sequence. */
    let alpha = 1.0 - Float::powf(1.0 - significance_level, 1.0 / num_tests as Float);

    if pval < alpha || !pval.is_finite() {
        let result = format!(
            "Rejected the null hypothesis (p-value = {}, significance level = {}",
            pval, alpha
        );
        return Err(result);
    } else {
        return Ok(());
    }
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
        for k in 0..CHI2_RUNS {
            /* Randomly pick an outgoing direction on the hemisphere */
            let sample = Point2f::new(rng.uniform_float(), rng.uniform_float());
            let wo_l = cosine_sample_hemisphere(&sample);
            let wo = bsdf.local_to_world(&wo_l);
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
                CHI2_MINFREQ,
                CHI2_SLEVEL,
                CHI2_RUNS,
            ) {
                assert!(false, "{}: {}, iteration {}", description, e, k);
            }
        }
    }
}

fn create_lambertian(bsdf: &mut BSDF, _arena: &mut MemoryArena) {
    let kd = Spectrum::from(1.0);
    let r: Arc<dyn BxDF> = Arc::new(LambertianReflection::new(&kd));
    bsdf.add(&r);
}

fn create_microfacet(
    bsdf: &mut BSDF,
    _arena: &mut MemoryArena,
    beckmann: bool,
    samplevisible: bool,
    roughx: Float,
    roughy: Float,
) {
    let ks = Spectrum::from(1.0);
    //MicrofacetDistribution* distrib;
    let distrib: Box<dyn MicrofacetDistribution> = if beckmann {
        let alphax = BeckmannDistribution::roughness_to_alpha(roughx);
        let alphay = BeckmannDistribution::roughness_to_alpha(roughy);
        Box::new(BeckmannDistribution::new(alphax, alphay, samplevisible))
    } else {
        let alphax = TrowbridgeReitzDistribution::roughness_to_alpha(roughx);
        let alphay = TrowbridgeReitzDistribution::roughness_to_alpha(roughy);
        Box::new(TrowbridgeReitzDistribution::new(
            alphax,
            alphay,
            samplevisible,
        ))
    };
    let fresnel = Box::new(FresnelNoOp::new());
    let bxdf: Arc<dyn BxDF> = Arc::new(MicrofacetReflection::new(&ks, distrib, fresnel));
    bsdf.add(&bxdf);
}

fn create_fresnel_blend(
    bsdf: &mut BSDF,
    _arena: &mut MemoryArena,
    beckmann: bool,
    samplevisible: bool,
    roughx: Float,
    roughy: Float,
) {
    let d = Spectrum::from(0.5);
    let s = Spectrum::from(0.5);
    let distrib: Box<dyn MicrofacetDistribution> = if beckmann {
        let alphax = BeckmannDistribution::roughness_to_alpha(roughx);
        let alphay = BeckmannDistribution::roughness_to_alpha(roughy);
        Box::new(BeckmannDistribution::new(alphax, alphay, samplevisible))
    } else {
        let alphax = TrowbridgeReitzDistribution::roughness_to_alpha(roughx);
        let alphay = TrowbridgeReitzDistribution::roughness_to_alpha(roughy);
        Box::new(TrowbridgeReitzDistribution::new(
            alphax,
            alphay,
            samplevisible,
        ))
    };
    let bxdf: Arc<dyn BxDF> = Arc::new(FresnelBlend::new(&d, &s, distrib));
    bsdf.add(&bxdf);
}

#[test]
fn bsdf_sampling_lambertian() {
    test_bsdf(&create_lambertian, "Lambertian");
}

#[test]
fn bsdf_sampling_beckmann_va_0p5() {
    test_bsdf(
        &|bsdf, arena| {
            create_microfacet(bsdf, arena, true, true, 0.5, 0.5);
        },
        "Beckmann, visible area sample, alpha = 0.5",
    );
}

#[test]
fn bsdf_sampling_tr_va_0p5() {
    test_bsdf(
        &|bsdf, arena| {
            create_microfacet(bsdf, arena, false, true, 0.5, 0.5);
        },
        "Trowbridge-Reitz, visible area sample, alpha = 0.5",
    );
}

#[test]
fn bsdf_sampling_beckmann_std_0p5() {
    test_bsdf(
        &|bsdf, arena| {
            create_microfacet(bsdf, arena, true, false, 0.5, 0.5);
        },
        "Beckmann, std sample, alpha = 0.5",
    );
}

#[test]
fn bsdf_sampling_tr_std_0p5() {
    test_bsdf(
        &|bsdf, arena| {
            create_microfacet(bsdf, arena, false, false, 0.5, 0.5);
        },
        "Trowbridge-Reitz, std sample, alpha = 0.5",
    );
}

#[test]
fn bsdf_sampling_beckmann_va_0p2_0p1() {
    test_bsdf(
        &|bsdf, arena| {
            create_microfacet(bsdf, arena, true, true, 0.2, 0.1);
        },
        "Beckmann, visible area sample, alpha = 0.2/0.1",
    );
}

#[test]
fn bsdf_sampling_tr_va_0p3_0p15() {
    test_bsdf(
        &|bsdf, arena| {
            create_microfacet(bsdf, arena, false, true, 0.3, 0.15);
        },
        "Trowbridge-Reitz, visible area sample, alpha = 0.3/0.15",
    );
}

#[test]
fn bsdf_sampling_beckmann_std_0p2_0p1() {
    test_bsdf(
        &|bsdf, arena| {
            create_microfacet(bsdf, arena, true, false, 0.2, 0.1);
        },
        "Beckmann, std sample, alpha = 0.2/0.1",
    );
}

#[test]
fn bsdf_sampling_tr_std_0p2_0p1() {
    test_bsdf(
        &|bsdf, arena| {
            create_microfacet(bsdf, arena, false, false, 0.2, 0.1);
        },
        "Trowbridge-Reitz, std sample, alpha = 0.2/0.1",
    );
}

#[test]
fn bsdf_sampling_beckmann_va_0p4_0p3() {
    test_bsdf(
        &|bsdf, arena| {
            create_fresnel_blend(bsdf, arena, true, true, 0.4, 0.3);
        },
        "Fresnel blend Beckmann, visible area sample, alpha = 0.4/0.3",
    );
}

#[test]
fn bsdf_sampling_tr_va_0p3() {
    test_bsdf(
        &|bsdf, arena| {
            create_fresnel_blend(bsdf, arena, false, true, 0.3, 0.3);
        },
        "Fresnel blend Trowbridge-Reitz, visible area sample, alpha = 0.3",
    );
}

#[test]
fn bsdf_sampling_beckmann_std_0p2() {
    test_bsdf(
        &|bsdf, arena| {
            create_fresnel_blend(bsdf, arena, true, false, 0.2, 0.2);
        },
        "Fresnel blend Beckmann, std sample, alpha = 0.2",
    );
}

#[test]
fn bsdf_sampling_tr_std_0p05_0p1() {
    test_bsdf(
        &|bsdf, arena| {
            create_fresnel_blend(bsdf, arena, false, false, 0.05, 0.1);
        },
        "Fresnel blend Trowbridge-Reitz, std sample, alpha = 0.05/0.1",
    );
}
