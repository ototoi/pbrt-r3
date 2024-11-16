use crate::core::pbrt::*;

use std::sync::{Arc, RwLock};

pub const FILTER_TABLE_WIDTH: usize = 16;
pub const FT_W: usize = FILTER_TABLE_WIDTH;
pub const FT_SZ: usize = FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct FilmTilePixel {
    pub contrib_sum: Spectrum,
    pub filter_weight_sum: Float,
}

impl FilmTilePixel {
    pub fn zero() -> Self {
        FilmTilePixel {
            contrib_sum: Spectrum::zero(),
            filter_weight_sum: 0.0,
        }
    }
}

pub struct FilmTile {
    pub pixel_bounds: Bounds2i,
    pub filter_radius: Vector2f,
    pub inv_filter_radius: Vector2f,
    pub filter_table: Arc<RwLock<[f32; FT_SZ]>>,
    pub pixels: Vec<FilmTilePixel>,
    pub max_sample_luminance: Float,
}

impl FilmTile {
    pub fn new(
        pixel_bounds: &Bounds2i,
        filter_radius: &Vector2f,
        filter_table: &Arc<RwLock<[f32; FT_SZ]>>,
        max_sample_luminance: Float,
    ) -> Self {
        let inv_filter_radius = Vector2f::new(1.0 / filter_radius.x, 1.0 / filter_radius.y);
        return FilmTile {
            pixel_bounds: *pixel_bounds,
            filter_radius: *filter_radius,
            inv_filter_radius,
            filter_table: Arc::clone(filter_table),
            pixels: vec![FilmTilePixel::zero(); pixel_bounds.area() as usize],
            max_sample_luminance,
        };
    }

    //fn i2f(v: &Vector2i) -> Vector2f {
    //    Vector2f::new(v.x as Float, v.y as Float)
    //}

    fn f2i(v: &Vector2f) -> Vector2i {
        Vector2i::new(v.x as i32, v.y as i32)
    }

    fn ceil(v: &Vector2f) -> Vector2i {
        Vector2i::new(v.x.ceil() as i32, v.y.ceil() as i32)
    }

    fn floor(v: &Vector2f) -> Vector2i {
        Vector2i::new(v.x.floor() as i32, v.y.floor() as i32)
    }

    fn min(a: &Vector2i, b: &Vector2i) -> Vector2i {
        return Vector2i::new(i32::min(a.x, b.x), i32::min(a.y, b.y));
    }

    fn max(a: &Vector2i, b: &Vector2i) -> Vector2i {
        return Vector2i::new(i32::max(a.x, b.x), i32::max(a.y, b.y));
    }

    pub fn add_sample(&mut self, p_film: &Point2f, l: &Spectrum, sample_weight: Float) {
        self.add_sample_filter(p_film, l, sample_weight);
    }

    pub fn add_sample_filter(&mut self, p_film: &Point2f, l: &Spectrum, sample_weight: Float) {
        //println!("x:{}, y:{}", p_film.x, p_film.y);
        let mut l = l.clone();
        if l.y() > self.max_sample_luminance {
            l *= self.max_sample_luminance / l.y();
        }
        let filter_table = self.filter_table.read().unwrap();
        let filter_radius = self.filter_radius;
        let inv_filter_radius = self.inv_filter_radius;
        let p_film_discrete = *p_film;

        //let p0 = Self::floor(&(p_film_discrete - self.filter_radius));
        //let p1 = Self::ceil(&(p_film_discrete + self.filter_radius));
        let p0 = Self::floor(&(p_film_discrete - Point2f::from(0.5) - self.filter_radius)); //0.5-0.5-0.5 -> floor(-0.5) -> -1
        let p1 = Self::ceil(&(p_film_discrete - Point2f::from(0.5) + self.filter_radius)); //0.5-0.5+0.5 -> ceil(0.5) -> +1

        assert!(p1.x - p0.x > 0);
        assert!(p1.y - p0.y > 0);

        let p0 = Self::max(&p0, &self.pixel_bounds.min);
        let p1 = Self::min(&p1, &self.pixel_bounds.max);

        let delta = p1 - p0;
        if delta.x <= 0 || delta.y <= 0 {
            return;
        }

        let filter_table_size = FT_W;
        let mut ifx: Vec<i32> = vec![0; delta.x as usize];
        let lx = inv_filter_radius.x * (filter_table_size - 1) as f32;
        for x in p0.x..p1.x {
            let d = x as f32 + 0.5 - p_film_discrete.x;
            let id = if -filter_radius.x <= d && d < filter_radius.x {
                i32::min(
                    f32::floor(d.abs() * lx) as i32,
                    (filter_table_size - 1) as i32,
                )
            } else {
                -1
            };
            ifx[(x - p0.x) as usize] = id;
        }

        let mut ify: Vec<i32> = vec![0; delta.y as usize];
        let ly = inv_filter_radius.y * (filter_table_size - 1) as f32;
        for y in p0.y..p1.y {
            let d = y as f32 + 0.5 - p_film_discrete.y;
            let id = if -filter_radius.y <= d && d < filter_radius.y {
                i32::min(
                    f32::floor(d.abs() * ly) as i32,
                    (filter_table_size - 1) as i32,
                )
            } else {
                -1
            };
            ify[(y - p0.y) as usize] = id;
        }
        let mut weights = vec![0.0; ifx.len() * ify.len()];
        for y in p0.y..p1.y {
            let yidx = (y - p0.y) as usize;
            let iy = ify[yidx];
            if iy == -1 {
                continue;
            }
            for x in p0.x..p1.x {
                let xidx = (x - p0.x) as usize;
                let ix = ifx[xidx];
                if ix == -1 {
                    continue;
                }
                let offset = (iy * filter_table_size as i32 + ix) as usize;
                let filter_weight = filter_table[offset];
                weights[(yidx * ifx.len() + xidx) as usize] = filter_weight;
            }
        }
        {
            let sum = weights.iter().sum::<f32>();
            if sum <= 0.0 {
                return;
            }
            let isum = 1.0 / sum;
            for w in weights.iter_mut() {
                *w *= isum;
            }
        }

        let width = (self.pixel_bounds.max.x - self.pixel_bounds.min.x) as usize;
        for y in p0.y..p1.y {
            for x in p0.x..p1.x {
                let xidx = (x - p0.x) as usize;
                let yidx = (y - p0.y) as usize;
                let filter_weight = weights[yidx * ifx.len() + xidx];

                let xx = x - self.pixel_bounds.min.x;
                let yy = y - self.pixel_bounds.min.y;

                assert!(xx >= 0);
                assert!(yy >= 0);

                let xx = xx as usize;
                let yy = yy as usize;

                let pindex: usize = yy * width + xx;
                self.pixels[pindex].contrib_sum += l * sample_weight * filter_weight;
                self.pixels[pindex].filter_weight_sum += filter_weight;
            }
        }
    }

    pub fn add_sample_single(&mut self, p_film: &Point2f, l: &Spectrum, _sample_weight: Float) {
        let mut l = l.clone();
        if l.y() > self.max_sample_luminance {
            l *= self.max_sample_luminance / l.y();
        }

        let p = Self::f2i(p_film);
        let x = p.x;
        let y = p.y;

        if x < self.pixel_bounds.min.x || self.pixel_bounds.max.x <= x {
            return;
        }
        if y < self.pixel_bounds.min.y || self.pixel_bounds.max.y <= y {
            return;
        }
        //println!("y:{:?}", y);
        let width = (self.pixel_bounds.max.x - self.pixel_bounds.min.x) as usize;

        let xx = (x - self.pixel_bounds.min.x) as usize;
        let yy = (y - self.pixel_bounds.min.y) as usize;

        let pindex = yy * width + xx;
        self.pixels[pindex].contrib_sum += l;
        self.pixels[pindex].filter_weight_sum += 1.0;
    }

    pub fn get_pixel_index(&self, p: &Vector2i) -> usize {
        let width = self.pixel_bounds.max.x - self.pixel_bounds.min.x;
        let x = p.x - self.pixel_bounds.min.x;
        let y = p.y - self.pixel_bounds.min.y;
        return (y * width + x) as usize;
    }

    pub fn get_pixel_bounds(&self) -> Bounds2i {
        return self.pixel_bounds;
    }
}
