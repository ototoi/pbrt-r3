pub use crate::core::pbrt::*;

pub fn shuffle_array<T: Copy>(samp: &mut [T], count: usize, dim: u32, rng: &mut RNG) {
    let dim = dim as usize;
    for i in 0..count {
        let other = i + rng.uniform_uint32_threshold((count - i) as u32) as usize;
        for j in 0..dim {
            let a = dim * i + j;
            let b = dim * other + j;
            samp.swap(a, b);
            //std::mem::swap(&mut samp[a], &mut samp[b]);
        }
    }
}

pub fn stratified_sample_1d(samples: &mut [Float], nsamples: usize, rng: &mut RNG, jitter: bool) {
    let inv_nsamples = 1.0 / (nsamples as Float);
    for i in 0..nsamples {
        let delta = if jitter { rng.uniform_float() } else { 0.5 };
        samples[i] = Float::min(((i as Float) + delta) * inv_nsamples, ONE_MINUS_EPSILON);
    }
}

pub fn stratified_sample_2d(
    samples: &mut [Point2f],
    nx: usize,
    ny: usize,
    rng: &mut RNG,
    jitter: bool,
) {
    let dx = 1.0 / (nx as Float);
    let dy = 1.0 / (ny as Float);
    for y in 0..ny {
        for x in 0..ny {
            let i = y * nx + x;
            let jx = if jitter { rng.uniform_float() } else { 0.5 };
            let jy = if jitter { rng.uniform_float() } else { 0.5 };
            let xx = Float::min(((x as Float) + jx) * dx, ONE_MINUS_EPSILON);
            let yy = Float::min(((y as Float) + jy) * dy, ONE_MINUS_EPSILON);
            samples[i] = Point2f::new(xx, yy);
        }
    }
}

pub fn latin_hypercube(samples: &mut [Float], nsamples: usize, ndim: usize, rng: &mut RNG) {
    // Generate LHS samples along diagonal
    let inv_nsamples = 1.0 / (nsamples as Float);
    for i in 0..nsamples {
        for j in 0..ndim {
            let sj = (i as Float + (rng.uniform_float())) * inv_nsamples;
            samples[ndim * i + j] = Float::min(sj, ONE_MINUS_EPSILON);
        }
    }
    // Permute LHS samples in each dimension
    for i in 0..ndim {
        for j in 0..nsamples {
            let other = j + rng.uniform_uint32_threshold((nsamples - j) as u32) as usize;
            let a = ndim * j + i;
            let b = ndim * other + i;
            samples.swap(a, b);
        }
    }
}

pub fn latin_hypercube_1d(samples: &mut [Float], nsamples: usize, rng: &mut RNG) {
    latin_hypercube(samples, nsamples, 1, rng);
}

pub fn latin_hypercube_2d(samples: &mut [Point2f], nsamples: usize, rng: &mut RNG) {
    let mut v = vec![0.0; samples.len() * 2];
    for i in 0..samples.len() {
        v[2 * i + 0] = samples[i].x;
        v[2 * i + 1] = samples[i].y;
    }
    latin_hypercube(&mut v, nsamples, 2, rng);
    for i in 0..samples.len() {
        samples[i].x = v[2 * i + 0];
        samples[i].y = v[2 * i + 1];
    }
}

//LatinHypercube
//RejectionSampleDisk
pub fn uniform_sample_hemisphere(u: &Point2f) -> Vector3f {
    let z = u[0];
    let r = Float::sqrt(Float::max(0.0, 1.0 - z * z));
    let phi = 2.0 * PI * u[1];
    return Vector3f::new(r * Float::cos(phi), r * Float::sin(phi), z);
}

pub fn uniform_hemisphere_pdf() -> Float {
    return INV_2_PI;
}

#[inline]
pub fn uniform_sample_sphere(u: &Point2f) -> Vector3f {
    let z = 1.0 - 2.0 * u[0];
    let r = Float::sqrt(Float::max(0.0, 1.0 - z * z));
    let phi = 2.0 * PI * u[1];
    return Vector3f::new(r * Float::cos(phi), r * Float::sin(phi), z);
}

#[inline]
pub fn uniform_sphere_pdf() -> Float {
    return INV_2_PI;
}

pub fn uniform_sample_triangle(u: &Point2f) -> Point2f {
    let su0 = Float::sqrt(u[0]);
    return Point2f::new(1.0 - su0, u[1] * su0);
}

pub fn concentric_sample_disk(u: &Point2f) -> Point2f {
    // Map uniform random numbers to $[-1,1]^2$
    let u_offset = *u * 2.0 - Vector2f::new(1.0, 1.0);

    // Handle degeneracy at the origin
    if u_offset.x == 0.0 && u_offset.y == 0.0 {
        return Point2f::zero();
    }

    // Apply concentric mapping to point
    if Float::abs(u_offset.x) > Float::abs(u_offset.y) {
        let r = u_offset.x;
        let theta = PI_OVER_4 * (u_offset.y / u_offset.x);
        return Point2f::new(r * Float::cos(theta), r * Float::sin(theta));
    } else {
        let r = u_offset.y;
        let theta = PI_OVER_2 - PI_OVER_4 * (u_offset.x / u_offset.y);
        return Point2f::new(r * Float::cos(theta), r * Float::sin(theta));
    }
}

#[inline]
pub fn uniform_cone_pdf(cos_theta_max: Float) -> Float {
    return 1.0 / (2.0 * PI * (1.0 - cos_theta_max));
}

pub fn uniform_sample_cone(u: &Point2f, cos_theta_max: Float) -> Vector3f {
    let cos_theta = (1.0 - u[0]) + u[0] * cos_theta_max;
    let sin_theta = Float::sqrt(1.0 - cos_theta * cos_theta);
    let phi = u[1] * 2.0 * PI;
    return Vector3f::new(
        Float::cos(phi) * sin_theta,
        Float::sin(phi) * sin_theta,
        cos_theta,
    );
}

//UniformSampleCone
//UniformSampleCone
//UniformConePdf

pub fn cosine_sample_hemisphere(u: &Point2f) -> Vector3f {
    let d = concentric_sample_disk(u);
    let z = Float::sqrt(Float::max(0.0, 1.0 - d.x * d.x - d.y * d.y));
    return Vector3f::new(d.x, d.y, z);
}

pub fn cosine_hemisphere_pdf(cos_theta: Float) -> Float {
    return cos_theta * INV_PI;
}

//BalanceHeuristic

#[inline]
pub fn power_heuristic(nf: i32, f_pdf: Float, ng: i32, g_pdf: Float) -> Float {
    let f = nf as Float * f_pdf;
    let g = ng as Float * g_pdf;
    return (f * f) / (f * f + g * g);
}
