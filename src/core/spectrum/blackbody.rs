use crate::core::pbrt::*;

pub fn blackbody(lambda: &[f64], t: f64) -> Vec<f64> {
    let n = lambda.len();
    let mut le = vec![0.0; n];
    if t <= 0.0 {
        return le;
    }
    const C: f64 = 299792458.0;
    const H: f64 = 6.62606957e-34;
    const KB: f64 = 1.3806488e-23;
    for i in 0..n {
        let l = lambda[i] * 1e-9;
        let lambda5 = (l * l) * (l * l) * l;
        let le_i = (2.0 * H * C * C) / (lambda5 * (f64::exp((H * C) / (l * KB * t)) - 1.0));
        le[i] = le_i;
    }
    return le;
}

pub fn blackbody_normalized(lambda: &[Float], t: Float) -> Vec<Float> {
    let n = lambda.len();
    let t = t as f64;
    let lambda: Vec<f64> = lambda.iter().map(|x| *x as f64).collect();
    let mut le = blackbody(&lambda, t);
    let lambda_max = [2.8977721e-3 / t * 1e9];
    let max_l = blackbody(&lambda_max, t);
    for i in 0..n {
        le[i] /= max_l[0];
    }
    let le: Vec<Float> = le.iter().map(|x| *x as Float).collect();
    return le;
}
