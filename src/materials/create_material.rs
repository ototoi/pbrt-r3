use log::warn;

use super::fourier::*;
use super::glass::*;
use super::hair::*;
use super::kdsubsurface::*;
use super::matte::*;
use super::metal::*;
use super::mirror::*;
use super::plastic::*;
use super::substrate::*;
use super::subsurface::*;
use super::translucent::*;
use super::uber::*;

use crate::core::pbrt::*;

use std::sync::Arc;

pub fn create_material(name: &str, mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    match name {
        "matte" => {
            return create_matte_material(mp);
        }
        "plastic" => {
            return create_plastic_material(mp);
        }
        "translucent" => {
            return create_translucent_material(mp);
        }
        "glass" => {
            return create_glass_material(mp);
        }
        "mirror" => {
            return create_mirror_material(mp);
        }
        "hair" => {
            return create_hair_material(mp);
        }
        "mix" => {
            let s = format!("Unable to create mix material \"{}\"", name);
            return Err(PbrtError::error(&s));
        }
        "metal" => {
            return create_metal_material(mp);
        }
        "substrate" => {
            return create_substrate_material(mp);
        }
        "subsurface" => {
            return create_subsurface_material(mp);
        }
        "kdsubsurface" => {
            return create_kdsubsurface_material(mp);
        }
        "uber" => {
            return create_uber_material(mp);
        }
        "fourier" => {
            return create_fourier_material(mp);
        }
        _ => {
            if name != "" {
                warn!("Material \"{}\" unknown. Using \"matte\".", name);
            }
            return create_matte_material(mp);
        }
    }
}
