use crate::core::pbrt::*;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, PartialEq, Default, Clone)]
pub struct SobolSampler {
    pub base: BaseGlobalSampler,
    //dimension: u32,
    //dimension_1d: u32,
    //dimension_2d: u32,
    sample_bounds: Bounds2i,
    resolution: u32,
    log2_resolution: u32,
}

impl SobolSampler {
    pub fn new(samples_per_pixel: u64, sample_bounds: &Bounds2i) -> Self {
        let diagonal = sample_bounds.diagonal();
        let resolution = (round_up_pow2(i32::max(diagonal.x, diagonal.y))) as u32;
        let log2_resolution = log2int(resolution);
        SobolSampler {
            base: BaseGlobalSampler::new(round_up_pow2(samples_per_pixel as i32) as u32),
            //dimension: 0,
            //dimension_1d: 0,
            //dimension_2d: 0,
            sample_bounds: *sample_bounds,
            resolution,
            log2_resolution,
        }
    }
}

impl Sampler for SobolSampler {
    fn start_pixel(&mut self, p: &Point2i) {
        self.base.base.start_pixel(p);
        self.base.dimension = 0;
        //self.dimension_1d = 0;
        //self.dimension_2d = 0;
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
        //self.dimension_1d = 0;
        //self.dimension_2d = 0;
        self.base.interval_sample_index =
            self.get_index_for_sample((self.base.base.current_pixel_sample_index + 1) as i64);
        return self.base.base.start_next_sample();
    }

    fn set_sample_number(&mut self, sample_num: u32) -> bool {
        self.base.dimension = 0;
        //self.dimension_1d = 0;
        //self.dimension_2d = 0;
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

        //if self.dimension_1d >= BaseGlobalSampler::array_start_dim
        //    && self.dimension_1d < self.base.array_end_dim
        //{
        //    self.dimension_1d = self.base.array_end_dim;
        //}
        //let d = self.dimension_1d;
        //let x = self.sample_dimension(self.base.interval_sample_index, d);
        //self.dimension_1d += 1;
        return x;
    }

    fn get_2d(&mut self) -> Point2f {
        if self.base.dimension + 1 >= BaseGlobalSampler::array_start_dim
            && self.base.dimension < self.base.array_end_dim
        {
            self.base.dimension = self.base.array_end_dim;
        }
        let d = self.base.dimension;
        //println!("d:{}/{}", d, self.base.array_end_dim);
        let x = self.sample_dimension(self.base.interval_sample_index, d);
        let y = self.sample_dimension(self.base.interval_sample_index, d + 1);
        let p = Point2f::new(x, y);
        self.base.dimension += 2;

        //if self.dimension_2d + 1 >= BaseGlobalSampler::array_start_dim
        //    && self.dimension_2d < self.base.array_end_dim
        //{
        //    self.dimension_2d = self.base.array_end_dim;
        //}
        //let d = self.dimension_2d;
        //println!("d:{}/{}", d, self.base.array_end_dim);
        //let x = self.sample_dimension(self.base.interval_sample_index, d);
        //let y = self.sample_dimension(self.base.interval_sample_index, d + 1);
        //let p = Point2f::new(x, y);
        //self.dimension_2d += 2;
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

unsafe impl Sync for SobolSampler {}

impl GlobalSampler for SobolSampler {
    fn get_index_for_sample(&self, sample_num: i64) -> i64 {
        let p = self.base.base.current_pixel - self.sample_bounds.min;
        return sobol_interval_to_index(self.log2_resolution, sample_num as u64, &p) as i64;
    }
    fn sample_dimension(&self, index: i64, dim: u32) -> Float {
        let mut s = sobol_sample(index, dim, 0);
        // Remap Sobol$'$ dimensions used for pixel samples
        if dim == 0 || dim == 1 {
            s = s * (self.resolution as Float) + self.sample_bounds.min[dim as usize] as Float;
            s = Float::clamp(
                s - self.base.base.current_pixel[dim as usize] as Float,
                0.0,
                ONE_MINUS_EPSILON,
            );
            s = Float::fract(s);
        }
        return s;
    }

    fn get_base(&mut self) -> &mut BaseGlobalSampler {
        return &mut self.base;
    }
}

pub fn create_sobol_sampler(
    params: &ParamSet,
    sample_bounds: &Bounds2i,
) -> Result<Arc<RwLock<dyn Sampler>>, PbrtError> {
    let nsamp = params.find_one_int("pixelsamples", 16) as u64;
    return Ok(Arc::new(RwLock::new(SobolSampler::new(
        nsamp,
        sample_bounds,
    ))));
}
