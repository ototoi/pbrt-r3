use crate::core::base::*;
use crate::core::camera::*;

use std::sync::Arc;
use std::sync::RwLock;
/*
    virtual void StartPixel(const Point2i &p);
    virtual Float Get1D() = 0;
    virtual Point2f Get2D() = 0;
    CameraSample GetCameraSample(const Point2i &pRaster);
    void Request1DArray(int n);
    void Request2DArray(int n);
    virtual int RoundCount(int n) const { return n; }
    const Float *Get1DArray(int n);
    const Point2f *Get2DArray(int n);
    virtual bool StartNextSample();
    virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
    virtual bool SetSampleNumber(int64_t sampleNum);

*/

pub trait Sampler {
    fn start_pixel(&mut self, _p: &Point2i) {}
    fn get_1d(&mut self) -> Float;
    fn get_2d(&mut self) -> Point2f;
    fn get_camera_sample(&mut self, p: &Point2i) -> CameraSample {
        let p_film = Point2f::new(p.x as Float, p.y as Float) + self.get_2d();
        let p_lens = self.get_2d();
        let time = self.get_1d();
        //let _ = self.get_1d();
        CameraSample {
            p_film,
            p_lens,
            time,
        }
    }
    fn request_1d_array(&mut self, _n: u32);
    fn request_2d_array(&mut self, _n: u32);
    fn get_1d_array(&mut self, n: u32) -> Option<Vec<Float>>;
    fn get_2d_array(&mut self, n: u32) -> Option<Vec<Vector2f>>;
    fn round_count(&self, n: u32) -> u32 {
        return n;
    }
    fn start_next_sample(&mut self) -> bool {
        return false;
    }
    fn clone_with_seed(&self, seed: u32) -> Arc<RwLock<dyn Sampler>>;
    fn set_sample_number(&mut self, _sample_num: u32) -> bool {
        return false;
    }
    fn get_samples_per_pixel(&self) -> u32;
}
