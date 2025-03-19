use super::environment::*;
use super::orthographic::*;
use super::perspective::*;
use super::realistic::*;
use crate::core::camera::*;
use crate::core::error::PbrtError;
use crate::core::film::*;
use crate::core::medium::*;
use crate::core::param_set::*;
use crate::core::transform::*;

use std::sync::Arc;
use std::sync::RwLock;

pub fn create_camera(
    name: &str,
    params: &ParamSet,
    cam_to_world: &AnimatedTransform,
    film: &Arc<RwLock<Film>>,
    medium_interface: &MediumInterface,
) -> Result<Arc<dyn Camera>, PbrtError> {
    let medium = medium_interface.get_outside();
    match name {
        "perspective" => create_perspective_camera(params, cam_to_world, film, &medium),
        "realistic" => create_realistic_camera(params, cam_to_world, film, &medium),
        "orthographic" => create_orthographic_camera(params, cam_to_world, film, &medium),
        "environment" => create_environment_camera(params, cam_to_world, film, &medium),
        _ => {
            return Err(PbrtError::warning(&format!("Camera \"{}\" unknown.", name)));
        }
    }
}
