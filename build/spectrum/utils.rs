use super::config::*;

#[inline]
pub fn lerp(t: f32, v1: f32, v2: f32) -> f32 {
    return (1.0 - t) * v1 + t * v2;
}

pub fn average_spectrum_samples(
    lambda: &[f32],
    vals: &[f32],
    lambda_start: f32,
    lambda_end: f32,
) -> f32 {
    let n = lambda.len();
    for i in 0..n - 1 {
        assert!(lambda[i] <= lambda[i + 1]);
    }
    assert!(lambda_start <= lambda_end);

    // Handle cases with out-of-bounds range or single sample only
    if lambda_end as f32 <= lambda[0] {
        return vals[0];
    }
    if lambda_start as f32 >= lambda[n - 1] {
        return vals[n - 1];
    }
    if n == 1 {
        return vals[0];
    }
    let mut sum: f32 = 0.0;
    // Add contributions of constant segments before/after samples
    if lambda_start < lambda[0] {
        sum += vals[0] * (lambda[0] - lambda_start);
    }
    if lambda_end > lambda[n - 1] {
        sum += vals[n - 1] * (lambda_end - lambda[n - 1]);
    }

    // Advance to first relevant wavelength segment
    let mut i = 0;
    while lambda_start > lambda[i + 1] {
        i += 1;
    }
    assert!((i + 1) <= n);

    // Loop over wavelength sample segments and add contributions
    let interp = |w: f32, i: usize| {
        return lerp(
            (w - lambda[i]) / (lambda[i + 1] - lambda[i]),
            vals[i],
            vals[i + 1],
        );
    };
    while i + 1 < n && lambda_end >= lambda[i] {
        let seg_lambda_start = f32::max(lambda_start, lambda[i]);
        let seg_lambda_end = f32::min(lambda_end, lambda[i + 1]);
        sum += 0.5
            * (interp(seg_lambda_start, i) + interp(seg_lambda_end, i))
            * (seg_lambda_end - seg_lambda_start);
        i += 1;
    }
    return sum / (lambda_end - lambda_start);
}

pub fn sample_spectrum(lambda: &[f32], vals: &[f32]) -> [f32; SPECTRAL_SAMPLES] {
    let mut x: [f32; SPECTRAL_SAMPLES] = [0.0; SPECTRAL_SAMPLES];
    for i in 0..SPECTRAL_SAMPLES {
        let wl0 = lerp(
            (i as f32) / (SPECTRAL_SAMPLES as f32),
            SAMPLED_LAMBDA_START as f32,
            SAMPLED_LAMBDA_END as f32,
        );
        let wl1 = lerp(
            ((i + 1) as f32) / (SPECTRAL_SAMPLES as f32),
            SAMPLED_LAMBDA_START as f32,
            SAMPLED_LAMBDA_END as f32,
        );
        x[i] = average_spectrum_samples(lambda, vals, wl0, wl1);
    }
    return x;
}
