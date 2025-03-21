use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::spectrum::*;

pub struct SplatImage {
    pub cropped_pixel_bounds: Bounds2i,
    pub width: usize,
    pub height: usize,
    pub splat_buffer: Vec<[Float; 3]>,
}

impl SplatImage {
    pub fn new(cropped_pixel_bounds: &Bounds2i) -> Self {
        let xres = (cropped_pixel_bounds.max.x - cropped_pixel_bounds.min.x) as usize;
        let yres = (cropped_pixel_bounds.max.y - cropped_pixel_bounds.min.y) as usize;
        SplatImage {
            cropped_pixel_bounds: *cropped_pixel_bounds,
            width: xres,
            height: yres,
            splat_buffer: vec![[0.0; 3]; xres * yres],
        }
    }

    pub fn add_splat(&mut self, p: &Vector2f, v: &Spectrum) {
        let _p = ProfilePhase::new(Prof::SplatFilm);

        let pi = Point2i::new(p.x.floor() as i32, p.y.floor() as i32);
        if !self.cropped_pixel_bounds.inside_exclusive(&pi) {
            return;
        }

        let pi = Point2i::new(
            pi.x - self.cropped_pixel_bounds.min.x,
            pi.y - self.cropped_pixel_bounds.min.y,
        );

        let offset = pi.y as usize * self.width + pi.x as usize;
        let xyz = v.to_xyz();
        self.splat_buffer[offset][0] += xyz[0];
        self.splat_buffer[offset][1] += xyz[1];
        self.splat_buffer[offset][2] += xyz[2];
    }

    pub fn clear(&mut self) {
        for i in 0..self.splat_buffer.len() {
            self.splat_buffer[i] = [0.0; 3];
        }
    }

    /*
    pub fn add_splat(&mut self, p: Point2i, v: [Float; 3]) {
        let offset = p.y as usize * p.x as usize;
        self.splat_buffer[offset][0] += v[0];
        self.splat_buffer[offset][1] += v[1];
        self.splat_buffer[offset][2] += v[2];
    }

    pub fn get_splat(&self, p: Point2i) -> [Float; 3] {
        let offset = p.y as usize * p.x as usize;
        self.splat_buffer[offset]
    }
    pub fn add_splat_image(&mut self, other: &SplatImage) {
        for i in 0..self.splat_buffer.len() {
            self.splat_buffer[i][0] += other.splat_buffer[i][0];
            self.splat_buffer[i][1] += other.splat_buffer[i][1];
            self.splat_buffer[i][2] += other.splat_buffer[i][2];
        }
    }
    pub fn scale_splat_image(&mut self, s: Float) {
        for i in 0..self.splat_buffer.len() {
            self.splat_buffer[i][0] *= s;
            self.splat_buffer[i][1] *= s;
            self.splat_buffer[i][2] *= s;
        }
    }
    */
}
