use crate::core::prelude::*;

use super::diffuse::*;
use super::distant::*;
use super::goniometric::*;
use super::infinite::*;
use super::point::*;
use super::projection::*;
use super::spot::*;

use std::sync::Arc;

pub fn create_light(
    name: &str,
    light2world: &Transform,
    medium_interface: &MediumInterface,
    params: &ParamSet,
) -> Result<Arc<dyn Light>, PbrtError> {
    let medium = medium_interface.get_outside();
    if name == "point" {
        return create_point_light(light2world, &medium, params);
    } else if name == "spot" {
        return create_spot_light(light2world, &medium, params);
    } else if name == "goniometric" {
        return create_goniometric_light(light2world, &medium, params);
    } else if name == "projection" {
        return create_projection_light(light2world, &medium, params);
    } else if name == "distant" {
        return create_distant_light(light2world, params);
    } else if name == "infinite" || name == "exinfinite" {
        return create_infinite_light(light2world, params);
    } else {
        let msg = format!("Light \"{}\" unknown.", name);
        return Err(PbrtError::warning(&msg));
    }
}

pub fn create_area_light(
    name: &str,
    light2world: &Transform,
    medium_interface: &MediumInterface,
    params: &ParamSet,
    shape: &Arc<dyn Shape>,
) -> Result<Arc<dyn Light>, PbrtError> {
    if name == "area" || name == "diffuse" {
        let medium = medium_interface.get_outside();
        return create_diffuse_area_light(light2world, &medium, params, shape);
    } else {
        return Err(PbrtError::warning(&format!(
            "Area light \"{}\" unknown.",
            name
        )));
    }
}
