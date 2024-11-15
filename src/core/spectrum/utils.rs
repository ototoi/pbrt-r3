include!(concat!(env!("OUT_DIR"), "/spectrum_utils.rs"));
use crate::core::misc::find_interval;

pub fn spectrum_samples_sorted(lambda: &[f32], _vals: &[f32]) -> bool {
    for i in 0..(lambda.len() - 1) {
        if lambda[i] > lambda[i + 1] {
            return false;
        }
    }
    return true;
}

pub fn sort_spectrum_samples(lambda: &mut [f32], vals: &mut [f32]) {
    let n = lambda.len();
    let mut sort_vec: Vec<(f32, f32)> = Vec::with_capacity(n);
    for i in 0..n {
        sort_vec.push((lambda[i], vals[i]));
    }
    sort_vec.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for i in 0..n {
        lambda[i] = sort_vec[i].0;
        vals[i] = sort_vec[i].1;
    }
}

pub fn interpolate_spectrum_samples(lambda: &[f32], vals: &[f32], l: f32) -> f32 {
    //for (int i = 0; i < n - 1; ++i) CHECK_GT(lambda[i + 1], lambda[i]);
    let n = lambda.len();
    if l <= lambda[0] {
        return vals[0];
    }
    if l >= lambda[n - 1] {
        return vals[n - 1];
    }
    let offset = find_interval(lambda, &|v, index| -> bool {
        return v[index] <= l;
    });
    //CHECK(l >= lambda[offset] && l <= lambda[offset + 1]);
    let t = (l - lambda[offset]) / (lambda[offset + 1] - lambda[offset]);
    return lerp(t, vals[offset], vals[offset + 1]);
}
