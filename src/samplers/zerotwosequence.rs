use crate::core::pbrt::*;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, PartialEq, Default, Clone)]
pub struct ZeroTwoSequenceSampler {
    base: BasePixelSampler,
}

impl ZeroTwoSequenceSampler {
    pub fn new(samples_per_pixel: u32, n_sampled_dimensions: u32) -> Self {
        let samples_per_pixel_pow2 = round_up_pow2(samples_per_pixel as i32) as u32;
        if !is_power_of_2(samples_per_pixel as i32) {
            log::warn!(
                "Pixel samples being rounded up to power of 2 (from \" {} \" to \" {} \").",
                samples_per_pixel,
                samples_per_pixel_pow2
            );
        }

        ZeroTwoSequenceSampler {
            base: BasePixelSampler::new(samples_per_pixel_pow2, n_sampled_dimensions),
        }
    }
}

impl Sampler for ZeroTwoSequenceSampler {
    fn start_pixel(&mut self, p: &Point2i) {
        // Generate 1D and 2D pixel sample components using $(0,2)$-sequence
        for i in 0..self.base.samples1d.len() {
            van_der_corput(
                1,
                self.base.base.samples_per_pixel,
                &mut self.base.samples1d[i],
                &mut self.base.rng,
            );
        }
        for i in 0..self.base.samples2d.len() {
            sobol_2d(
                1,
                self.base.base.samples_per_pixel,
                &mut self.base.samples2d[i],
                &mut self.base.rng,
            );
        }
        // Generate 1D and 2D array samples using $(0,2)$-sequence
        for i in 0..self.base.base.samples1d_array_sizes.len() {
            van_der_corput(
                self.base.base.samples1d_array_sizes[i],
                self.base.base.samples_per_pixel,
                &mut self.base.base.sample_array1d[i],
                &mut self.base.rng,
            );
        }
        for i in 0..self.base.base.samples2d_array_sizes.len() {
            sobol_2d(
                self.base.base.samples2d_array_sizes[i],
                self.base.base.samples_per_pixel,
                &mut self.base.base.sample_array2d[i],
                &mut self.base.rng,
            );
        }
        self.base.start_pixel(p);
    }

    fn start_next_sample(&mut self) -> bool {
        return self.base.start_next_sample();
    }

    fn set_sample_number(&mut self, sample_num: u32) -> bool {
        return self.base.set_sample_number(sample_num);
    }

    fn get_1d(&mut self) -> Float {
        return self.base.get_1d();
    }

    fn get_2d(&mut self) -> Vector2f {
        return self.base.get_2d();
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

    fn clone_with_seed(&self, seed: u32) -> Arc<RwLock<dyn Sampler>> {
        let mut ss = self.clone();
        ss.base.rng.set_sequence(seed as u64);
        return Arc::new(RwLock::new(ss));
    }

    fn get_samples_per_pixel(&self) -> u32 {
        return self.base.base.samples_per_pixel;
    }
}

impl PixelSampler for ZeroTwoSequenceSampler {}

pub fn create_zerotwosequence_sampler(
    params: &ParamSet,
) -> Result<Arc<RwLock<dyn Sampler>>, PbrtError> {
    let mut nsamp = params.find_one_int("pixelsamples", 16) as u32;
    let sd = params.find_one_int("dimensions", 4) as u32;
    {
        let options = PbrtOptions::get();
        if options.quick_render {
            nsamp = 1;
        }
    }
    return Ok(Arc::new(RwLock::new(ZeroTwoSequenceSampler::new(
        nsamp, sd,
    ))));
}
