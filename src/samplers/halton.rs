use crate::core::lowdiscrepancy::primes::{PRIMES, PRIME_SUMS};

use crate::core::pbrt::*;
use std::cell::Cell;
use std::sync::Arc;
use std::sync::OnceLock;
use std::sync::RwLock;

const K_MAX_RESOLUTION: i32 = 128;
static RADICAL_INVERSE_PERMUTATIONS: OnceLock<Vec<u16>> = OnceLock::new();

fn get_radical_inverse_permutations() -> &'static Vec<u16> {
    let v = RADICAL_INVERSE_PERMUTATIONS.get_or_init(|| -> Vec<u16> {
        let mut v = Vec::new();
        let mut rng = RNG::new();
        compute_radical_inverse_permutations(&mut v, &mut rng);
        return v;
    });
    return v;
}

fn math_mod(a: i64, b: i64) -> u64 {
    let result = a - (a / b) * b;
    return if result < 0 {
        (result + b) as u64
    } else {
        result as u64
    };
}

fn extended_gcd(a: u64, b: u64) -> (i64, i64) {
    if b == 0 {
        return (1, 0);
    }
    let d = (a / b) as i64;
    let (xp, yp) = extended_gcd(b, a % b);
    return (yp, xp - (d * yp));
}

fn multiplicative_inverse(a: i64, n: i64) -> u64 {
    let (x, _y) = extended_gcd(a as u64, n as u64);
    return math_mod(x, n);
}

#[derive(Debug, Default, PartialEq, Clone)]
pub struct HaltonSampler {
    base: BaseGlobalSampler,
    base_scales: [i32; 2],
    base_exponents: [i32; 2],
    sample_stride: i32,
    mult_inverse: [i32; 2],
    pixel_for_offset: Cell<Point2i>,
    offset_for_current_pixel: Cell<i64>,
    sample_at_pixel_center: bool,
}

impl HaltonSampler {
    pub fn new(samples_per_pixel: u32, sample_bounds: &Bounds2i, sample_at_center: bool) -> Self {
        let _ = get_radical_inverse_permutations();

        // Find radical inverse base scales and exponents that cover sampling area
        let res = sample_bounds.max - sample_bounds.min;
        const BASES: [i32; 2] = [2, 3];

        let mut base_scales = [1, 1];
        let mut base_exponents = [1, 1];
        for i in 0..2 {
            //let base = if i == 0 { 2 } else { 3 };
            let base = BASES[i];
            let mut scale = 1;
            let mut exp = 0;
            while scale < i32::min(res[i], K_MAX_RESOLUTION) {
                scale *= base;
                exp += 1;
            }
            base_scales[i] = scale;
            base_exponents[i] = exp;
        }

        // Compute stride in samples for visiting each pixel area
        let sample_stride = base_scales[0] * base_scales[1];
        let mult_inverse = [
            multiplicative_inverse(base_scales[1] as i64, base_scales[0] as i64) as i32,
            multiplicative_inverse(base_scales[0] as i64, base_scales[1] as i64) as i32,
        ];

        let pixel_for_offset = Cell::new(Point2i::new(i32::MAX, i32::MAX));
        let offset_for_current_pixel = Cell::new(0);
        let base = BaseGlobalSampler::new(samples_per_pixel);

        HaltonSampler {
            base,
            base_scales,
            base_exponents,
            sample_stride,
            mult_inverse,
            pixel_for_offset,
            offset_for_current_pixel,
            sample_at_pixel_center: sample_at_center,
        }
    }

    fn permutation_for_dimension(dim: u32) -> &'static [u16] {
        let table_size = PRIMES.len();
        if dim as usize >= table_size {
            panic!("HaltonSampler can only sample {} dimensions.", table_size);
        }
        let start = PRIME_SUMS[dim as usize] as usize;
        let rip = get_radical_inverse_permutations();
        return &rip[start..];
    }
}

impl GlobalSampler for HaltonSampler {
    fn get_index_for_sample(&self, sample_num: i64) -> i64 {
        let sample_stride = self.sample_stride as i64;
        if self.base.base.current_pixel != self.pixel_for_offset.get() {
            // Compute Halton sample offset for _currentPixel_
            self.offset_for_current_pixel.set(0);
            if sample_stride > 1 {
                let pm = [
                    math_mod(
                        self.base.base.current_pixel[0] as i64,
                        K_MAX_RESOLUTION as i64,
                    ),
                    math_mod(
                        self.base.base.current_pixel[1] as i64,
                        K_MAX_RESOLUTION as i64,
                    ),
                ];
                const BASES: [u64; 2] = [2, 3];
                for i in 0..2 {
                    let dim_offset =
                        inverse_radical_inverse(BASES[i], pm[i], self.base_exponents[i] as usize);
                    let addition = (dim_offset
                        * ((self.sample_stride / self.base_scales[i]) * self.mult_inverse[i])
                            as u64) as i64;
                    self.offset_for_current_pixel
                        .set(self.offset_for_current_pixel.get() + addition);
                }
                let offset_for_current_pixel = self.offset_for_current_pixel.get() % sample_stride;
                self.offset_for_current_pixel.set(offset_for_current_pixel);
            }
            self.pixel_for_offset.set(self.base.base.current_pixel);
        }
        return self.offset_for_current_pixel.get() + sample_num * sample_stride;
    }

    fn sample_dimension(&self, index: i64, dim: u32) -> Float {
        if self.sample_at_pixel_center && (dim == 0 || dim == 1) {
            return 0.5;
        }
        return match dim {
            0 => radical_inverse(dim, (index >> self.base_exponents[0]) as u64),
            1 => radical_inverse(dim, (index / self.base_scales[1] as i64) as u64),
            _ => scrambled_radical_inverse(
                dim,
                index as u64,
                HaltonSampler::permutation_for_dimension(dim),
            ),
        };
    }

    fn get_base(&mut self) -> &mut BaseGlobalSampler {
        return &mut self.base;
    }
}

impl Sampler for HaltonSampler {
    fn start_pixel(&mut self, p: &Point2i) {
        self.base.base.start_pixel(p);
        self.base.dimension = 0;
        self.base.interval_sample_index = self.get_index_for_sample(0);
        self.base.array_end_dim = BaseGlobalSampler::array_start_dim
            + (self.base.base.sample_array1d.len() + 2 * self.base.base.sample_array2d.len())
                as u32;
        let samples1d_array_sizes = self.base.base.samples1d_array_sizes.len();
        // Compute 1D array samples for _GlobalSampler_
        {
            let l = samples1d_array_sizes;
            for i in 0..l {
                let n_samples =
                    self.base.base.samples1d_array_sizes[i] * self.base.base.samples_per_pixel;
                for j in 0..n_samples {
                    let index = self.get_index_for_sample(j as i64);
                    self.base.base.sample_array1d[i][j as usize] =
                        self.sample_dimension(index, self.base.array_end_dim + i as u32);
                }
            }
        }
        // Compute 2D array samples for _GlobalSampler_
        {
            let dim = BaseGlobalSampler::array_start_dim + samples1d_array_sizes as u32;
            let l = self.base.base.samples2d_array_sizes.len();
            for i in 0..l {
                let n_samples =
                    self.base.base.samples2d_array_sizes[i] * self.base.base.samples_per_pixel;
                for j in 0..n_samples {
                    let index = self.get_index_for_sample(j as i64);
                    self.base.base.sample_array2d[i][j as usize].x =
                        self.sample_dimension(index, dim);
                    self.base.base.sample_array2d[i][j as usize].y =
                        self.sample_dimension(index, dim + 1);
                }
            }
        }
    }

    fn start_next_sample(&mut self) -> bool {
        self.base.dimension = 0;
        self.base.interval_sample_index =
            self.get_index_for_sample((self.base.base.current_pixel_sample_index + 1) as i64);
        return self.base.base.start_next_sample();
    }

    fn set_sample_number(&mut self, sample_num: u32) -> bool {
        self.base.dimension = 0;
        self.base.interval_sample_index = self.get_index_for_sample(sample_num as i64);
        return self.base.base.set_sample_number(sample_num);
    }

    fn get_1d(&mut self) -> Float {
        if self.base.dimension >= BaseGlobalSampler::array_start_dim
            && self.base.dimension < self.base.array_end_dim
        {
            self.base.dimension = self.base.array_end_dim;
        }
        let d = self.base.dimension;
        let x = self.sample_dimension(self.base.interval_sample_index, d);
        self.base.dimension += 1;
        return x;
    }

    fn get_2d(&mut self) -> Vector2f {
        if self.base.dimension + 1 >= BaseGlobalSampler::array_start_dim
            && self.base.dimension < self.base.array_end_dim
        {
            self.base.dimension = self.base.array_end_dim;
        }
        let d = self.base.dimension;
        let x = self.sample_dimension(self.base.interval_sample_index, d);
        let y = self.sample_dimension(self.base.interval_sample_index, d + 1);
        let p = Point2f::new(x, y);
        self.base.dimension += 2;
        return p;
    }

    fn request_1d_array(&mut self, n: u32) {
        self.base.request_1d_array(n);
    }
    fn request_2d_array(&mut self, n: u32) {
        self.base.request_2d_array(n);
    }
    fn get_1d_array(&mut self, n: u32) -> Option<Vec<Float>> {
        return self.base.get_1d_array(n);
    }
    fn get_2d_array(&mut self, n: u32) -> Option<Vec<Vector2f>> {
        return self.base.get_2d_array(n);
    }

    fn clone_with_seed(&self, _seed: u32) -> Arc<RwLock<dyn Sampler>> {
        return Arc::new(RwLock::new(self.clone()));
    }
    fn get_samples_per_pixel(&self) -> u32 {
        return self.base.base.samples_per_pixel;
    }
}

pub fn create_halton_sampler(
    params: &ParamSet,
    sample_bounds: &Bounds2i,
) -> Result<Arc<RwLock<dyn Sampler>>, PbrtError> {
    let ns = params.find_one_int("pixelsamples", 16) as u32;
    let sample_at_center = params.find_one_bool("samplepixelcenter", false);

    // pbrt-r3
    let _ = get_radical_inverse_permutations();
    // pbrt-r3

    return Ok(Arc::new(RwLock::new(HaltonSampler::new(
        ns,
        sample_bounds,
        sample_at_center,
    ))));
}
