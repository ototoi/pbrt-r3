use super::functions::*;
use super::mipmap_weight_lut::MIPMAP_WEIGHT_LUT;
use super::texture::*;
use crate::core::error::PbrtError;
use crate::core::geometry::*;
use crate::core::imageio::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use std::mem::size_of;

use log::*;
use rayon::prelude::*;
use std::fmt::Debug;
use std::{ops::Deref, vec};

thread_local!(static N_EWA_LOOKUPS: StatCounter = StatCounter::new("Texture/EWA lookups"));
thread_local!(static N_TRILERP_LOOKUPS: StatCounter = StatCounter::new("Texture/Trilinear lookups"));
thread_local!(static MIP_MAP_MEMORY: StatMemoryCounter = StatMemoryCounter::new("Memory/Texture MIP maps"));

pub struct F32MIPMapImage {
    pub resolution: (usize, usize),
    pub data: Vec<f32>,
}

impl F32MIPMapImage {
    pub fn new(data: Vec<f32>, resolution: (usize, usize)) -> Self {
        F32MIPMapImage { resolution, data }
    }
}

impl From<(&[f32], (usize, usize))> for F32MIPMapImage {
    fn from(v: (&[f32], (usize, usize))) -> Self {
        F32MIPMapImage {
            resolution: v.1,
            data: v.0.to_vec(),
        }
    }
}

impl From<(&[f64], (usize, usize))> for F32MIPMapImage {
    fn from(v: (&[f64], (usize, usize))) -> Self {
        let data: Vec<f32> = v.0.iter().map(|x| *x as f32).collect();
        F32MIPMapImage {
            resolution: v.1,
            data: data,
        }
    }
}

impl From<(&[RGBSpectrum], (usize, usize))> for F32MIPMapImage {
    fn from(v: (&[RGBSpectrum], (usize, usize))) -> Self {
        let mut data = vec![0.0; 3 * v.0.len()];
        for i in 0..v.0.len() {
            let c = v.0[i].to_rgb();
            data[3 * i + 0] = c[0] as f32;
            data[3 * i + 1] = c[1] as f32;
            data[3 * i + 2] = c[2] as f32;
        }
        F32MIPMapImage {
            resolution: v.1,
            data,
        }
    }
}

pub trait MIPMapImage<T> {
    fn lookup(&self, i: usize) -> T;
    fn get_width(&self) -> usize;
    fn get_height(&self) -> usize;
    fn as_data(&self) -> &F32MIPMapImage;
}

impl MIPMapImage<Float> for F32MIPMapImage {
    fn lookup(&self, i: usize) -> Float {
        return self.data[i] as Float;
    }
    fn get_width(&self) -> usize {
        self.resolution.0
    }
    fn get_height(&self) -> usize {
        self.resolution.1
    }
    fn as_data(&self) -> &F32MIPMapImage {
        return self;
    }
}

impl MIPMapImage<RGBSpectrum> for F32MIPMapImage {
    fn lookup(&self, i: usize) -> RGBSpectrum {
        let r = self.data[3 * i + 0] as Float;
        let g = self.data[3 * i + 1] as Float;
        let b = self.data[3 * i + 2] as Float;
        return RGBSpectrum::rgb_from_rgb(&[r, g, b]);
    }
    fn get_width(&self) -> usize {
        self.resolution.0
    }
    fn get_height(&self) -> usize {
        self.resolution.1
    }
    fn as_data(&self) -> &F32MIPMapImage {
        return self;
    }
}

struct ResampleWeight {
    first_texel: i32,
    weight: [f32; 4],
}

fn math_mod(a: i32, b: i32) -> i32 {
    let result = a - (a / b) * b;
    return if result < 0 { result + b } else { result };
}

fn downsample_half(img: &[f32], channels: usize, width: usize, height: usize) -> Vec<f32> {
    let nc = channels;
    let nw = width / 2;
    let nh = height;
    /*
    let mut new_img = vec![0.0; channels * nw * nh];
    for y in 0..nh {
        for x in 0..nw {
            for c in 0..nc {
                let src0 = y * width + (2 * x + 0);
                let src1 = y * width + (2 * x + 1);
                let dst = y * nw + x;
                new_img[nc * dst + c] = img[nc * src0 + c] * 0.5 + img[nc * src1 + c] * 0.5;
            }
        }
    }
    */
    let new_img_: Vec<_> = (0..nh)
        .into_par_iter()
        .map(|y| {
            let mut v = vec![0.0; channels * nw];
            for x in 0..nw {
                for c in 0..nc {
                    let i0 = y * width + (2 * x + 0);
                    let i1 = y * width + (2 * x + 1);
                    v[nc * x + c] = img[channels * i0 + c] * 0.5 + img[channels * i1 + c] * 0.5;
                }
            }
            return v;
        })
        .collect();
    let mut new_img = vec![0.0; nc * nw * nh];
    for y in 0..nh {
        let line = &new_img_[y];
        let d0 = y * nc * nw;
        let d1 = d0 + nc * nw;
        let s0 = 0;
        let s1 = s0 + nc * nw;
        new_img[d0..d1].copy_from_slice(&line[s0..s1]);
        /*
        for x in 0..nw {
            for c in 0..nc {
                let i = y * nw + x;
                new_img[nc * i + c] = line[nc * x + c];
            }
        }
        */
    }

    return new_img;
}

fn transpose_image(img: Vec<f32>, channels: usize, width: usize, height: usize) -> Vec<f32> {
    let mut dst = img.clone();
    for y in 0..height {
        for x in 0..width {
            for c in 0..channels {
                let s = channels * (y * width + x) + c;
                let d = channels * (x * height + y) + c;
                dst[d] = img[s];
            }
        }
    }
    return dst;
}

fn resample_weights(old_res: usize, new_res: usize) -> Vec<ResampleWeight> {
    let mut wt = Vec::with_capacity(new_res);
    const FILTER_WIDTH: Float = 2.0;
    for i in 0..new_res {
        let center: Float = (i as Float + 0.5) * (old_res as Float / new_res as Float);
        let first_texel = Float::floor((center - FILTER_WIDTH) + 0.5);
        let mut weight = [0.0; 4];
        for j in 0..4 {
            let pos = first_texel + (j as Float) + 0.5;
            weight[j] = lanczos((pos - center) / FILTER_WIDTH, FILTER_WIDTH) as f32;
        }
        // Normalize filter weights for texel resampling
        let inv_sum_wts = 1.0 / (weight[0] + weight[1] + weight[2] + weight[3]);
        for j in 0..4 {
            weight[j] *= inv_sum_wts;
        }
        let resampled = ResampleWeight {
            first_texel: (first_texel as i32),
            weight,
        };
        wt.push(resampled);
    }
    return wt;
}

fn resample_image(
    img: &[f32],
    channels: usize,
    resolution: (usize, usize),
    swrap_mode: ImageWrap,
    twrap_mode: ImageWrap,
) -> ((usize, usize), Vec<f32>) {
    if is_power_of_2(resolution.0 as u32) && is_power_of_2(resolution.1 as u32) {
        return (resolution, Vec::from(img));
    } else {
        let res_pow2 = (
            round_up_pow2(resolution.0 as u32) as usize,
            round_up_pow2(resolution.1 as u32) as usize,
        );
        let ratio = (res_pow2.0 * res_pow2.1) as Float / (resolution.0 * resolution.1) as Float;
        info!(
            "Resampling MIPMap from {:?} to {:?}. Ratio= {}",
            resolution, res_pow2, ratio
        );
        let mut resampled_image = vec![0.0; channels * res_pow2.0 * res_pow2.1];
        {
            let s_weights = resample_weights(resolution.0, res_pow2.0);
            for t in 0..resolution.1 {
                for s in 0..res_pow2.0 {
                    for j in 0..4 {
                        let mut orig_s = s_weights[s].first_texel + j;
                        if swrap_mode == ImageWrap::Repeat {
                            orig_s = math_mod(orig_s as i32, resolution.0 as i32);
                        } else if swrap_mode == ImageWrap::Clamp {
                            orig_s = i32::clamp(orig_s as i32, 0, resolution.0 as i32 - 1);
                        }
                        if orig_s >= 0 && orig_s < resolution.0 as i32 {
                            let w = s_weights[s].weight[j as usize];
                            let src = t * resolution.0 + orig_s as usize;
                            let dst = (t * res_pow2.0 + s) as usize;
                            for c in 0..channels {
                                resampled_image[channels * dst + c] += img[channels * src + c] * w;
                            }
                        }
                    }
                }
            }
        }

        {
            let t_weights = resample_weights(resolution.1, res_pow2.1);
            for s in 0..res_pow2.0 {
                let mut buffer = vec![0.0; channels * resolution.1];
                for t in 0..resolution.1 {
                    let src = (t * res_pow2.0 + s) as usize;
                    let dst = t;
                    for c in 0..channels {
                        buffer[channels * dst + c] = resampled_image[channels * src + c];
                    }
                }

                for t in 0..res_pow2.1 {
                    let mut l = [0.0; 3];
                    for j in 0..4 {
                        let mut orig_t = t_weights[t].first_texel + j;
                        if twrap_mode == ImageWrap::Repeat {
                            orig_t = math_mod(orig_t as i32, resolution.1 as i32);
                        } else if twrap_mode == ImageWrap::Clamp {
                            orig_t = i32::clamp(orig_t as i32, 0, resolution.1 as i32 - 1);
                        }

                        if orig_t >= 0 && orig_t < resolution.1 as i32 {
                            let w = t_weights[t as usize].weight[j as usize];
                            let src = orig_t as usize;
                            for c in 0..channels {
                                l[c] += buffer[channels * src + c] * w;
                            }
                        }
                    }
                    let dst = (t * res_pow2.0 + s) as usize;
                    for c in 0..channels {
                        resampled_image[channels * dst + c] = l[c];
                    }
                }
            }
        }
        {
            for i in 0..resampled_image.len() {
                resampled_image[i] = resampled_image[i].max(0.0);
            }
        }
        return (res_pow2, resampled_image);
    }
}

fn make_pyramid<T>(
    pyramid: &mut Vec<Box<dyn MIPMapImage<T>>>,
    channels: usize,
    resolution: (usize, usize),
    data: Vec<f32>,
) where
    F32MIPMapImage: MIPMapImage<T>,
{
    let total = resolution.0 * resolution.1;
    if total != 1 {
        let c = channels;
        let w = resolution.0;
        let h = resolution.1;
        let (s_data, w, h) = if w > 1 {
            let s_data = downsample_half(&data, c, w, h);
            (s_data, w / 2, h)
        } else {
            (data.clone(), w, h)
        };
        let (s_data, w, h) = if h > 1 {
            let s_data = transpose_image(s_data, c, w, h);
            let s_data = downsample_half(&s_data, c, h, w);
            let s_data = transpose_image(s_data, c, h / 2, w);
            (s_data, w, h / 2)
        } else {
            (s_data, w, h)
        };

        let image: Box<dyn MIPMapImage<T>> = Box::new(F32MIPMapImage { data, resolution });
        pyramid.push(image);
        make_pyramid(pyramid, c, (w, h), s_data);
    } else {
        let image: Box<dyn MIPMapImage<T>> = Box::new(F32MIPMapImage { data, resolution });
        pyramid.push(image);
    }
}

fn log2int_(x: usize) -> usize {
    return f32::ceil(f32::log(x as f32, 2.0)) as usize;
}

fn make_mipimages<T>(
    data: &[T],
    resolution: (usize, usize),
    swrap_mode: ImageWrap,
    twrap_mode: ImageWrap,
) -> Vec<Box<dyn MIPMapImage<T>>>
where
    T: Clone + Debug,
    F32MIPMapImage: MIPMapImage<T>,
    F32MIPMapImage: for<'a> From<(&'a [T], (usize, usize))>,
{
    let mipdata = F32MIPMapImage::from((data, resolution));
    let channels = mipdata.data.len() / (resolution.0 * resolution.1);

    let (resolution, data) =
        resample_image(&mipdata.data, channels, resolution, swrap_mode, twrap_mode);
    /*
    {
        let filename = "mipmap.exr";
        let _ = write_spectrum_mipmap_image(filename, &mipdata);
        let filename = "mipmap_resampled.exr";
        let _ = write_spectrum_mipmap_image(filename, &F32MIPMapImage::new(data.clone(), resolution));
    }
    */

    // Initialize levels of MIPMap from image
    let n_levels = 1 + log2int_(usize::max(resolution.0, resolution.1));
    let mut pyramid = Vec::with_capacity(n_levels);

    //let t0 = time::Instant::now();
    make_pyramid(&mut pyramid, channels, resolution, data);

    //for i in 0..pyramid.len() {
    //    let filename = format!("mipmap_{}.exr", i);
    //    let _ = write_spectrum_mipmap_image(&filename, &pyramid[i].as_data());
    //}

    return pyramid;
}

pub struct MIPMap<T> {
    pub do_trilinear: bool,
    pub max_anisotropy: Float,
    pub swrap_mode: ImageWrap,
    pub twrap_mode: ImageWrap,
    pub pyramid: Vec<Box<dyn MIPMapImage<T>>>,
}

impl<
        T: Default + Debug + Copy + std::ops::Add<T, Output = T> + std::ops::Mul<Float, Output = T>,
    > MIPMap<T>
where
    F32MIPMapImage: MIPMapImage<T>,
    F32MIPMapImage: for<'a> From<(&'a [T], (usize, usize))>,
{
    pub fn texel_static(
        image: &dyn MIPMapImage<T>,
        s: i32,
        t: i32,
        swrap_mode: ImageWrap,
        twrap_mode: ImageWrap,
    ) -> T {
        let w = image.get_width() as i32;
        let h = image.get_height() as i32;
        let mut s = s;
        let mut t = t;
        match swrap_mode {
            ImageWrap::Repeat => {
                //s = math_mod(s as i32, w);
                s &= w - 1;
            }
            ImageWrap::Clamp => {
                s = i32::clamp(s, 0, w - 1);
            }
            _ => {
                if s < 0 || w <= s {
                    return T::default();
                }
            }
        }
        match twrap_mode {
            ImageWrap::Repeat => {
                //t = math_mod(t as i32, h);
                t &= h - 1;
            }
            ImageWrap::Clamp => {
                t = i32::clamp(t, 0, h - 1);
            }
            _ => {
                if t < 0 || h <= t {
                    return T::default();
                }
            }
        }
        let index = (t * w + s) as usize;
        return image.lookup(index);
    }

    pub fn lerp(t: Float, a: T, b: T) -> T {
        return a * (1.0 - t) + b * t;
    }

    pub fn new(
        resolution: &Point2i,
        data: &[T],
        do_trilinear: bool,
        max_anisotropy: Float,
        swrap_mode: ImageWrap,
        twrap_mode: ImageWrap,
    ) -> Self {
        let _p = ProfilePhase::new(Prof::MIPMapCreation);

        let resolution = (resolution.x as usize, resolution.y as usize);
        let pyramid = make_mipimages(data, resolution, swrap_mode, twrap_mode);
        let mip = Self::make_from_pyramid(
            pyramid,
            do_trilinear,
            max_anisotropy,
            swrap_mode,
            twrap_mode,
        );
        return mip;
    }

    pub fn make_from_pyramid(
        pyramid: Vec<Box<dyn MIPMapImage<T>>>,
        do_trilinear: bool,
        max_anisotropy: Float,
        swrap_mode: ImageWrap,
        twrap_mode: ImageWrap,
    ) -> Self {
        {
            let total_consumption: usize = pyramid
                .iter()
                .map(|p| {
                    let resolution = p.as_data().resolution;
                    resolution.0 * resolution.1 * size_of::<T>()
                })
                .sum();
            MIP_MAP_MEMORY.with(|m| {
                m.add(total_consumption);
            });
        }

        let mip = MIPMap::<T> {
            do_trilinear,
            max_anisotropy,
            swrap_mode,
            twrap_mode,
            pyramid,
        };
        return mip;
    }

    pub fn width(&self) -> usize {
        return self.pyramid[0].get_width();
    }

    pub fn height(&self) -> usize {
        return self.pyramid[0].get_height();
    }

    pub fn levels(&self) -> usize {
        return self.pyramid.len();
    }

    pub fn texel(&self, level: usize, s: i32, t: i32) -> T {
        let image = &(self.pyramid[level]);
        let image = image.deref();
        return Self::texel_static(image, s, t, self.swrap_mode, self.twrap_mode);
    }

    pub fn lookup(&self, st: &Point2f, width: Float) -> T {
        N_TRILERP_LOOKUPS.with(|n| n.inc());
        let _p = ProfilePhase::new(Prof::TexFiltTrilerp);

        let max_level = (self.levels() - 1) as Float;
        let level = max_level + Float::log2(Float::max(width, 1e-8));
        if level < 0.0 {
            return self.triangle(0, st);
        } else if level >= max_level {
            return self.texel(self.levels() - 1, 0, 0);
        } else {
            let i_level = Float::floor(level) as usize;
            let delta = (level - i_level as Float).clamp(0.0, 1.0);
            let a = self.triangle(i_level, st);
            let b = self.triangle(i_level + 1, st);
            return Self::lerp(delta, a, b);
        }
    }

    pub fn lookup_delta(&self, st: &Point2f, dst0: &Vector2f, dst1: &Vector2f) -> T {
        if self.do_trilinear {
            let width = Float::max(
                Float::max(Float::abs(dst0[0]), Float::abs(dst0[1])),
                Float::max(Float::abs(dst1[0]), Float::abs(dst1[1])),
            );
            return self.lookup(st, width);
        } else {
            N_EWA_LOOKUPS.with(|n| n.inc());
            let _p = ProfilePhase::new(Prof::TexFiltEWA);
            let mut dst0 = *dst0;
            let mut dst1 = *dst1;
            // Compute ellipse minor and major axes
            if dst0.length_squared() < dst1.length_squared() {
                std::mem::swap(&mut dst0, &mut dst1);
            }
            let major_length = dst0.length(); //longest axis
            let mut minor_length = dst1.length(); //shortest axis

            // Clamp ellipse eccentricity if too large
            if minor_length * self.max_anisotropy < major_length && minor_length > 0.0 {
                let scale = major_length / (minor_length * self.max_anisotropy);
                dst1 *= scale;
                minor_length *= scale;
            }

            if minor_length <= 0.0 {
                return self.lookup(st, 0.0);
            }

            // Choose level of detail for EWA lookup and perform EWA filtering
            let lod = Float::max(
                0.0,
                self.levels() as Float - 1.0 + Float::log2(minor_length),
            );
            let ilod = lod.floor() as usize;
            return Self::lerp(
                lod - ilod as Float,
                self.ewa(ilod, st, dst0, dst1),
                self.ewa(ilod + 1, st, dst0, dst1),
            );
        }
    }

    pub fn ewa(&self, level: usize, st: &Point2f, dst0: Vector2f, dst1: Vector2f) -> T {
        if level >= self.levels() {
            return self.texel(self.levels() - 1, 0, 0);
        }
        // Convert EWA coordinates to appropriate scale for level
        let st = Point2f::new(
            st[0] * self.pyramid[level].get_width() as Float - 0.5,
            st[1] * self.pyramid[level].get_height() as Float - 0.5,
        );
        let dst0 = Vector2f::new(
            dst0[0] * self.pyramid[level].get_width() as Float,
            dst0[1] * self.pyramid[level].get_height() as Float,
        );
        let dst1 = Vector2f::new(
            dst1[0] * self.pyramid[level].get_width() as Float,
            dst1[1] * self.pyramid[level].get_height() as Float,
        );

        // Compute ellipse coefficients to bound EWA filter region
        let a = dst0[1] * dst0[1] + dst1[1] * dst1[1] + 1.0;
        let b = -2.0 * (dst0[0] * dst0[1] + dst1[0] * dst1[1]);
        let c = dst0[0] * dst0[0] + dst1[0] * dst1[0] + 1.0;
        let inv_f = 1.0 / (a * c - b * b * 0.25);

        let a = a * inv_f;
        let b = b * inv_f;
        let c = c * inv_f;

        // Compute the ellipse's $(s,t)$ bounding box in texture space
        let det = -b * b + 4.0 * a * c;
        let inv_det = 1.0 / det;
        let u_sqrt = (det * c).sqrt();
        let v_sqrt = (det * a).sqrt();
        let s0 = Float::ceil(st[0] - 2.0 * inv_det * u_sqrt) as i32;
        let s1 = Float::floor(st[0] + 2.0 * inv_det * u_sqrt) as i32;
        let t0 = Float::ceil(st[1] - 2.0 * inv_det * v_sqrt) as i32;
        let t1 = Float::floor(st[1] + 2.0 * inv_det * v_sqrt) as i32;

        // Scan over ellipse bound and compute quadratic equation

        let mut sum = T::default();
        let mut sum_wts = 0.0;
        for it in t0..=t1 {
            let tt = it as Float - st[1];
            for is in s0..=s1 {
                let ss = is as Float - st[0];
                // Compute squared radius and filter texel if inside ellipse
                let r2 = a * ss * ss + b * ss * tt + c * tt * tt;
                if r2 < 1.0 {
                    let index = usize::min(
                        (r2 * MIPMAP_WEIGHT_LUT.len() as Float) as usize,
                        MIPMAP_WEIGHT_LUT.len() - 1,
                    );
                    let weight = MIPMAP_WEIGHT_LUT[index];
                    sum = sum + self.texel(level, is, it) * weight;
                    sum_wts += weight;
                }
            }
        }
        return sum * (1.0 / sum_wts);
    }

    pub fn triangle(&self, level: usize, st: &Point2f) -> T {
        let level = usize::clamp(level, 0, self.levels() - 1);
        let w = self.pyramid[level].get_width();
        let h = self.pyramid[level].get_height();
        let s = st[0] * w as Float - 0.5;
        let t = st[1] * h as Float - 0.5;
        let s0 = Float::floor(s) as i32;
        let t0 = Float::floor(t) as i32;
        let ds = s - s0 as Float;
        let dt = t - t0 as Float;
        let a = self.texel(level, s0 + 0, t0 + 0) * ((1.0 - ds) * (1.0 - dt));
        let b = self.texel(level, s0 + 0, t0 + 1) * ((1.0 - ds) * dt);
        let c = self.texel(level, s0 + 1, t0 + 0) * (ds * (1.0 - dt));
        let d = self.texel(level, s0 + 1, t0 + 1) * (ds * dt);
        return a + b + c + d;
    }
}

#[allow(dead_code)]
fn write_spectrum_mipmap_image(path: &str, image: &F32MIPMapImage) -> Result<(), PbrtError> {
    let resolution = image.resolution;
    let w = resolution.0;
    let h = resolution.1;
    // let channels = image.data.len() / (w * h) as usize;
    let mut img: Vec<Float> = vec![0.0; w * h * 3];

    for i in 0..(w * h) {
        let c = image.data[i];
        img[3 * i + 0] = c as Float;
        img[3 * i + 1] = c as Float;
        img[3 * i + 2] = c as Float;
    }

    return write_image(
        &path,
        &img,
        &Bounds2i::from(((0, 0), (w as i32, h as i32))),
        &Point2i::new(w as i32, h as i32),
    );
}

/*
pub fn write_spectrum_mipmap(name: &str, mipmap: &MIPMap<RGBSpectrum>) -> Result<(), PbrtError> {
    for i in 0..mipmap.pyramid.len() {
        let resolution = mipmap.pyramid[i].resolution;
        //println!("res:{:?}", resolution);
        let data = &(mipmap.pyramid[i].data);
        let mut img: Vec<f32> = vec![0.0; data.len() * 3];
        for j in 0..data.len() {
            let c = data[j].to_rgb();
            img[3 * j + 0] = c[0];
            img[3 * j + 1] = c[1];
            img[3 * j + 2] = c[2];
        }
        let path = String::from(name) + &format!("._{}.png", i);
        let w = resolution[0];
        let h = resolution[1];
        write_image(&path, &img, &Bounds2i::from(((0, 0), (w, h))), &resolution)?;
    }
    Ok(())
}
*/

pub fn create_float_mipmap(
    resolution: &Point2i,
    data: &[Float],
) -> Result<MIPMap<Float>, PbrtError> {
    let mip = MIPMap::<Float>::new(
        resolution,
        data,
        false,
        8.0,
        ImageWrap::Repeat,
        ImageWrap::Clamp,
    );
    return Ok(mip);
}

pub fn create_spectrum_mipmap(
    resolution: &Point2i,
    data: &[RGBSpectrum],
) -> Result<MIPMap<RGBSpectrum>, PbrtError> {
    let mip = MIPMap::<RGBSpectrum>::new(
        resolution,
        data,
        false,
        8.0,
        ImageWrap::Repeat,
        ImageWrap::Clamp,
    );
    return Ok(mip);
}
