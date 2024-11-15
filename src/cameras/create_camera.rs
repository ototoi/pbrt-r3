use super::environment::*;
use super::orthographic::*;
use super::perspective::*;
use super::realistic::*;
use crate::core::pbrt::*;

use std::sync::Arc;
use std::sync::RwLock;

pub fn create_camera(
    name: &str,
    params: &ParamSet,
    cam_to_world: &AnimatedTransform,
    film: &Arc<RwLock<Film>>,
    medium_interface: &MediumInterface,
) -> Result<Arc<dyn Camera>, PbrtError> {
    match name {
        "perspective" => {
            create_perspective_camera(params, cam_to_world, film, &medium_interface.outside)
        }
        "realistic" => {
            create_realistic_camera(params, cam_to_world, film, &medium_interface.outside)
        }
        "orthographic" => {
            create_orthographic_camera(params, cam_to_world, film, &medium_interface.outside)
        }
        "environment" => {
            create_environment_camera(params, cam_to_world, film, &medium_interface.outside)
        }
        _ => {
            return Err(PbrtError::warning(&format!("Camera \"{}\" unknown.", name)));
        }
    }
}
