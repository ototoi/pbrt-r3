use super::disney::*;
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
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use crate::core::texture::*;

use std::sync::Arc;

use log::*;

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
        "disney" => {
            return create_disney_material(mp);
        }
        _ => {
            if name != "" {
                warn!("Material \"{}\" unknown. Using \"matte\".", name);
            }
            return create_matte_material(mp);
        }
    }
}
