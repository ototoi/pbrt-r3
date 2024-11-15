use crate::core::pbrt::*;
use image::*;

fn convert_from_luma8(img: &image::GrayImage, gamma: bool) -> (Vec<RGBSpectrum>, Point2i) {
    let (width, height) = img.dimensions();
    let mut spcs = vec![RGBSpectrum::zero(); (width * height) as usize];
    for y in 0..height {
        for x in 0..width {
            let index = (y * width + x) as usize;
            let pixel = img[(x, y)];
            let mut r = pixel[0] as f32 / 255.0;

            if gamma {
                r = inverse_gamma_correct(r);
            }

            spcs[index] = RGBSpectrum::rgb_from_rgb(&[r, r, r]);
        }
    }
    return (spcs, Point2i::from((width as i32, height as i32)));
}

fn convert_from_lumaa8(img: &image::GrayAlphaImage, gamma: bool) -> (Vec<RGBSpectrum>, Point2i) {
    let (width, height) = img.dimensions();
    let mut spcs = vec![RGBSpectrum::zero(); (width * height) as usize];
    for y in 0..height {
        for x in 0..width {
            let index = (y * width + x) as usize;
            let pixel = img[(x, y)];
            let mut r = pixel[0] as f32 / 255.0;
            let _a = pixel[1] as f32 / 255.0; // Ignore alpha channel

            if gamma {
                r = inverse_gamma_correct(r);
            }

            spcs[index] = RGBSpectrum::rgb_from_rgb(&[r, r, r]);
        }
    }
    return (spcs, Point2i::from((width as i32, height as i32)));
}

fn convert_from_rgb8(img: &image::RgbImage, gamma: bool) -> (Vec<RGBSpectrum>, Point2i) {
    let (width, height) = img.dimensions();
    let mut spcs = vec![RGBSpectrum::zero(); (width * height) as usize];
    for y in 0..height {
        for x in 0..width {
            let index = (y * width + x) as usize;
            let pixel = img[(x, y)];
            let mut r = pixel[0] as f32 / 255.0;
            let mut g = pixel[1] as f32 / 255.0;
            let mut b = pixel[2] as f32 / 255.0;

            if gamma {
                r = inverse_gamma_correct(r);
                g = inverse_gamma_correct(g);
                b = inverse_gamma_correct(b);
            }

            spcs[index] = RGBSpectrum::rgb_from_rgb(&[r, g, b]);
        }
    }
    return (spcs, Point2i::from((width as i32, height as i32)));
}

fn convert_from_rgba8(img: &image::RgbaImage, gamma: bool) -> (Vec<RGBSpectrum>, Point2i) {
    let (width, height) = img.dimensions();
    let mut spcs = vec![RGBSpectrum::zero(); (width * height) as usize];
    for y in 0..height {
        for x in 0..width {
            let index = (y * width + x) as usize;
            let pixel = img[(x, y)];
            let mut r = pixel[0] as f32 / 255.0;
            let mut g = pixel[1] as f32 / 255.0;
            let mut b = pixel[2] as f32 / 255.0;

            if gamma {
                r = inverse_gamma_correct(r);
                g = inverse_gamma_correct(g);
                b = inverse_gamma_correct(b);
            }

            spcs[index] = RGBSpectrum::rgb_from_rgb(&[r, g, b]);
        }
    }
    return (spcs, Point2i::from((width as i32, height as i32)));
}

fn convert_from_rgb32f(img: &image::Rgb32FImage) -> (Vec<RGBSpectrum>, Point2i) {
    let (width, height) = img.dimensions();
    let mut spcs = vec![RGBSpectrum::zero(); (width * height) as usize];
    for y in 0..height {
        for x in 0..width {
            let index = (y * width + x) as usize;
            let pixel = img.get_pixel(x, y);
            let r = pixel[0];
            let g = pixel[1];
            let b = pixel[2];
            spcs[index] = RGBSpectrum::new(r, g, b);
        }
    }
    return (spcs, Point2i::from((width as i32, height as i32)));
}

fn convert_from_rgba32f(img: &image::Rgba32FImage) -> (Vec<RGBSpectrum>, Point2i) {
    let (width, height) = img.dimensions();
    let mut spcs = vec![RGBSpectrum::zero(); (width * height) as usize];
    for y in 0..height {
        for x in 0..width {
            let index = (y * width + x) as usize;
            let pixel = img.get_pixel(x, y);
            let r = pixel[0];
            let g = pixel[1];
            let b = pixel[2];
            spcs[index] = RGBSpectrum::new(r, g, b);
        }
    }
    return (spcs, Point2i::from((width as i32, height as i32)));
}

pub fn read_image_gamma_correct(
    name: &str,
    gamma: bool,
) -> Result<(Vec<RGBSpectrum>, Point2i), PbrtError> {
    let r: Result<DynamicImage, ImageError> = image::open(name);
    match r {
        Ok(dimg) => match dimg {
            DynamicImage::ImageLuma8(img) => {
                return Ok(convert_from_luma8(&img, gamma));
            }
            DynamicImage::ImageLumaA8(img) => {
                return Ok(convert_from_lumaa8(&img, gamma));
            }
            DynamicImage::ImageRgb8(img) => {
                return Ok(convert_from_rgb8(&img, gamma));
            }
            DynamicImage::ImageRgba8(img) => {
                return Ok(convert_from_rgba8(&img, gamma));
            }
            DynamicImage::ImageRgb32F(img) => {
                return Ok(convert_from_rgb32f(&img));
            }
            DynamicImage::ImageRgba32F(img) => {
                return Ok(convert_from_rgba32f(&img));
            }
            _ => {
                return Err(PbrtError::from("This file is not supported."));
            }
        },
        Err(e) => {
            return Err(PbrtError::from(e.to_string()));
        }
    };
}

pub fn read_image(name: &str) -> Result<(Vec<RGBSpectrum>, Point2i), PbrtError> {
    return read_image_gamma_correct(name, false);
}

fn convert_f32_from_rgba32f(img: &image::Rgb32FImage) -> (Vec<f32>, Point2i) {
    let (width, height) = img.dimensions();
    let mut values = vec![0.0; (3 * width * height) as usize];
    for y in 0..height {
        for x in 0..width {
            let index = (y * width + x) as usize;
            let pixel = img.get_pixel(x, y);
            let r = pixel[0];
            let g = pixel[1];
            let b = pixel[2];
            values[3 * index + 0] = r;
            values[3 * index + 1] = g;
            values[3 * index + 2] = b;
        }
    }
    return (values, Point2i::from((width as i32, height as i32)));
}

pub fn read_cache_image(name: &str) -> Result<(Vec<f32>, Point2i), PbrtError> {
    let r: Result<DynamicImage, ImageError> = image::open(name);
    match r {
        Ok(dimg) => match dimg {
            DynamicImage::ImageRgb32F(img) => {
                return Ok(convert_f32_from_rgba32f(&img));
            }
            _ => {
                return Err(PbrtError::from("This file is not supported."));
            }
        },
        Err(e) => {
            return Err(PbrtError::from(e.to_string()));
        }
    };
}
