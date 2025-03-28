use super::bilerp::*;
use super::checkerboard::*;
use super::constant::*;
use super::dots::*;
use super::fbm::*;
use super::imagemap::*;
use super::marble::*;
use super::mix::*;
use super::normal::*;
use super::scale::*;
use super::uv::*;
use super::windy::*;
use super::wrinkled::*;
use crate::core::prelude::*;

use std::sync::Arc;

pub fn create_float_texture(
    name: &str,
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    match name {
        "constant" => {
            return create_constant_float_texture(tex2world, tp);
        }
        "scale" => {
            return create_scale_float_texture(tex2world, tp);
        }
        "mix" => {
            return create_mix_float_texture(tex2world, tp);
        }
        "bilerp" => {
            return create_bilerp_float_texture(tex2world, tp);
        }
        "imagemap" => {
            return create_image_float_texture(tex2world, tp);
        }
        "checkerboard" => {
            return create_checkerboard_float_texture(tex2world, tp);
        }
        "dots" => {
            return create_dots_float_texture(tex2world, tp);
        }
        "fbm" => {
            return create_fbm_float_texture(tex2world, tp);
        }
        "wrinkled" => {
            return create_wrinkled_float_texture(tex2world, tp);
        }
        //"marble" => {
        //    return create_marble_float_texture(tex2world, tp);
        //}
        "windy" => {
            return create_windy_float_texture(tex2world, tp);
        }
        _ => {
            let msg = format!("Float texture \"{}\" unknown.", name);
            return Err(PbrtError::warning(&msg));
        }
    }
}

pub fn create_spectrum_texture(
    name: &str,
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    match name {
        "constant" => {
            return create_constant_spectrum_texture(tex2world, tp);
        }
        "scale" => {
            return create_scale_spectrum_texture(tex2world, tp);
        }
        "mix" => {
            return create_mix_spectrum_texture(tex2world, tp);
        }
        "bilerp" => {
            return create_bilerp_spectrum_texture(tex2world, tp);
        }
        "imagemap" => {
            return create_image_spectrum_texture(tex2world, tp);
        }
        "uv" => {
            return create_uv_spectrum_texture(tex2world, tp);
        }
        "checkerboard" => {
            return create_checkerboard_spectrum_texture(tex2world, tp);
        }
        "dots" => {
            return create_dots_spectrum_texture(tex2world, tp);
        }
        "fbm" => {
            return create_fbm_spectrum_texture(tex2world, tp);
        }
        "wrinkled" => {
            return create_wrinkled_spectrum_texture(tex2world, tp);
        }
        "marble" => {
            return create_marble_spectrum_texture(tex2world, tp);
        }
        "windy" => {
            return create_windy_spectrum_texture(tex2world, tp);
        }
        "normal" => {
            return create_normal_texture(tex2world, tp);
        }
        _ => {
            let msg = format!("Spectrum texture \"{}\" unknown.", name);
            return Err(PbrtError::warning(&msg));
        }
    }
}
