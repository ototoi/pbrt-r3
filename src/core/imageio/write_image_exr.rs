use crate::core::error::*;
use crate::core::geometry::*;
use crate::core::pbrt::*;

use image::*;

pub fn write_image_exr(
    name: &str,
    rgb: &[Float],
    output_bounds: &Bounds2i,
    _total_resolution: &Point2i,
) -> Result<(), PbrtError> {
    let resolution = output_bounds.diagonal();
    let mut float_img: Vec<f32> = vec![0.0; (resolution.x * resolution.y * 3) as usize];
    {
        let x0: u32 = output_bounds.min.x as u32;
        let x1: u32 = output_bounds.max.x as u32;
        let y0: u32 = output_bounds.min.y as u32;
        let y1: u32 = output_bounds.max.y as u32;

        let width = x1 - x0;

        for y in y0..y1 {
            let yy = y - y0;
            for x in x0..x1 {
                let xx = x - x0;
                let index: usize = (yy * width + xx) as usize;
                float_img[3 * index + 0] = rgb[3 * index + 0] as f32;
                float_img[3 * index + 1] = rgb[3 * index + 1] as f32;
                float_img[3 * index + 2] = rgb[3 * index + 2] as f32;
            }
        }
    }
    let img = Rgb32FImage::from_vec(resolution.x as u32, resolution.y as u32, float_img).unwrap();
    match img.save(name) {
        Ok(()) => {
            return Ok(());
        }
        Err(e) => {
            return Err(PbrtError::from(e));
        }
    }
}
