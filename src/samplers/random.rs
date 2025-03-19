use crate::core::camera::*;
use crate::core::distribution::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::lightdistrib::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::rng::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use crate::core::texture::*;
use crate::core::transform::*;

use std::sync::Arc;
use std::sync::RwLock;

pub struct RandomSampler {
    base: BaseSampler,
    rng: RNG,
}

impl RandomSampler {
    pub fn new(ns: u32, seed: u64) -> Self {
        RandomSampler {
            base: BaseSampler::new(ns),
            rng: RNG::new_sequence(seed),
        }
    }
}

impl Sampler for RandomSampler {
    fn get_1d(&mut self) -> Float {
        let _p = ProfilePhase::new(Prof::GetSample);

        return self.rng.uniform_float();
    }
    fn get_2d(&mut self) -> Vector2f {
        let _p = ProfilePhase::new(Prof::GetSample);

        return Vector2f::new(self.rng.uniform_float(), self.rng.uniform_float());
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

    fn start_pixel(&mut self, p: &Point2i) {
        let _p = ProfilePhase::new(Prof::StartPixel);

        for i in 0..self.base.sample_array1d.len() {
            for j in 0..self.base.sample_array1d[i].len() {
                self.base.sample_array1d[i][j] = self.rng.uniform_float();
            }
        }
        for i in 0..self.base.sample_array2d.len() {
            for j in 0..self.base.sample_array2d[i].len() {
                self.base.sample_array2d[i][j] =
                    Vector2f::new(self.rng.uniform_float(), self.rng.uniform_float());
            }
        }
        self.base.start_pixel(p);
    }

    fn start_next_sample(&mut self) -> bool {
        return self.base.start_next_sample();
    }

    fn clone_with_seed(&self, seed: u32) -> Arc<RwLock<dyn Sampler>> {
        return Arc::new(RwLock::new(RandomSampler {
            base: self.base.clone(),
            rng: RNG::new_sequence(seed as u64),
        }));
    }

    fn get_samples_per_pixel(&self) -> u32 {
        return self.base.samples_per_pixel;
    }
}

pub fn create_random_sampler(params: &ParamSet) -> Result<Arc<RwLock<dyn Sampler>>, PbrtError> {
    let ns = params.find_one_int("pixelsamples", 4) as u32;
    return Ok(Arc::new(RwLock::new(RandomSampler::new(ns, 0))));
}
