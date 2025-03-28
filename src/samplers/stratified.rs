use crate::core::prelude::*;

use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, PartialEq, Default, Clone)]
pub struct StratifiedSampler {
    base: BasePixelSampler,
    x_pixel_samples: u32,
    y_pixel_samples: u32,
    jitter_samples: bool,
}

impl StratifiedSampler {
    pub fn new(
        x_pixel_samples: u32,
        y_pixel_samples: u32,
        jitter_samples: bool,
        n_sampled_dimensions: u32,
    ) -> Self {
        let samples_per_pixel = x_pixel_samples * y_pixel_samples;
        StratifiedSampler {
            base: BasePixelSampler::new(samples_per_pixel, n_sampled_dimensions),
            x_pixel_samples,
            y_pixel_samples,
            jitter_samples,
        }
    }
}

impl Sampler for StratifiedSampler {
    fn start_pixel(&mut self, p: &Point2i) {
        let _p = ProfilePhase::new(Prof::StartPixel);

        // Generate single stratified samples for the pixel
        let samples_per_pixel = (self.base.base.samples_per_pixel) as usize;
        for i in 0..self.base.samples1d.len() {
            stratified_sample_1d(
                &mut self.base.samples1d[i],
                samples_per_pixel,
                &mut self.base.rng,
                self.jitter_samples,
            );
            let count = self.base.samples1d[i].len();
            shuffle_array(&mut self.base.samples1d[i], count, 1, &mut self.base.rng);
        }
        for i in 0..self.base.samples2d.len() {
            stratified_sample_2d(
                &mut self.base.samples2d[i],
                self.x_pixel_samples as usize,
                self.y_pixel_samples as usize,
                &mut self.base.rng,
                self.jitter_samples,
            );
            let count = self.base.samples2d[i].len();
            shuffle_array(&mut self.base.samples2d[i], count, 1, &mut self.base.rng);
        }

        // Generate arrays of stratified samples for the pixel
        for i in 0..self.base.base.samples1d_array_sizes.len() {
            for j in 0..samples_per_pixel {
                let count = self.base.base.samples1d_array_sizes[i] as usize;
                let start = j * count;
                stratified_sample_1d(
                    &mut self.base.base.sample_array1d[i][start..],
                    count,
                    &mut self.base.rng,
                    self.jitter_samples,
                );
                shuffle_array(
                    &mut self.base.base.sample_array1d[i][start..],
                    count,
                    1,
                    &mut self.base.rng,
                );
            }
        }

        for i in 0..self.base.base.samples2d_array_sizes.len() {
            for j in 0..samples_per_pixel {
                let count = self.base.base.samples2d_array_sizes[i] as usize;
                let start = j * count;
                latin_hypercube_2d(
                    &mut self.base.base.sample_array2d[i][start..],
                    count,
                    &mut self.base.rng,
                );
            }
        }
        self.base.base.start_pixel(p);
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

impl PixelSampler for StratifiedSampler {}

pub fn create_stratified_sampler(params: &ParamSet) -> Result<Arc<RwLock<dyn Sampler>>, PbrtError> {
    let jitter = params.find_one_bool("jitter", true);
    let mut xsamp = params.find_one_int("xsamples", 4) as u32;
    let mut ysamp = params.find_one_int("ysamples", 4) as u32;
    let sd = params.find_one_int("dimensions", 4) as u32;
    {
        let options = PbrtOptions::get();
        if options.quick_render {
            xsamp = 1;
            ysamp = 1;
        }
    }
    return Ok(Arc::new(RwLock::new(StratifiedSampler::new(
        xsamp, ysamp, jitter, sd,
    ))));
}
