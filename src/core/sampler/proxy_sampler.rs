use std::sync::Arc;
use std::sync::RwLock;

use super::sampler::Sampler;
use crate::core::pbrt::*;

#[derive(Clone)]
pub struct ProxySampler {
    pub sampler: Arc<RwLock<dyn Sampler>>,
}

impl ProxySampler {
    pub fn new(s: &Arc<RwLock<dyn Sampler>>) -> Self {
        ProxySampler {
            sampler: Arc::clone(s),
        }
    }
}

impl Sampler for ProxySampler {
    fn start_pixel(&mut self, p: &Point2i) {
        let mut s = self.sampler.as_ref().write().unwrap();
        s.start_pixel(p);
    }
    fn get_1d(&mut self) -> Float {
        let mut s = self.sampler.as_ref().write().unwrap();
        return s.get_1d();
    }
    fn get_2d(&mut self) -> Point2f {
        let mut s = self.sampler.as_ref().write().unwrap();
        return s.get_2d();
    }
    fn request_1d_array(&mut self, n: u32) {
        let mut s = self.sampler.as_ref().write().unwrap();
        s.request_1d_array(n);
    }
    fn request_2d_array(&mut self, n: u32) {
        let mut s = self.sampler.as_ref().write().unwrap();
        s.request_2d_array(n);
    }
    fn get_1d_array(&mut self, n: u32) -> Option<Vec<Float>> {
        let mut s = self.sampler.as_ref().write().unwrap();
        return s.get_1d_array(n);
    }
    fn get_2d_array(&mut self, n: u32) -> Option<Vec<Vector2f>> {
        let mut s = self.sampler.as_ref().write().unwrap();
        return s.get_2d_array(n);
    }
    fn round_count(&self, n: u32) -> u32 {
        let s = self.sampler.as_ref().read().unwrap();
        return s.round_count(n);
    }
    fn get_camera_sample(&mut self, p: &Point2i) -> CameraSample {
        let mut s = self.sampler.as_ref().write().unwrap();
        return s.get_camera_sample(p);
    }
    fn start_next_sample(&mut self) -> bool {
        let mut s = self.sampler.as_ref().write().unwrap();
        return s.start_next_sample();
    }
    fn clone_with_seed(&self, seed: u32) -> Arc<RwLock<dyn Sampler>> {
        return Arc::new(RwLock::new(ProxySampler {
            sampler: self.sampler.as_ref().read().unwrap().clone_with_seed(seed),
        }));
    }
    fn set_sample_number(&mut self, sample_num: u32) -> bool {
        let mut s = self.sampler.as_ref().write().unwrap();
        return s.set_sample_number(sample_num);
    }
    fn get_samples_per_pixel(&self) -> u32 {
        let s = self.sampler.as_ref().read().unwrap();
        return s.get_samples_per_pixel();
    }
}

unsafe impl Sync for ProxySampler {}
unsafe impl Send for ProxySampler {}
