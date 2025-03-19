use super::box_filter::*;
use super::gaussian::*;
use super::mitchell::*;
use super::sinc::*;
use super::triangle::*;
use crate::core::error::PbrtError;
use crate::core::filter::*;
use crate::core::param_set::*;

use std::sync::Arc;

pub fn create_filter(name: &str, params: &ParamSet) -> Result<Arc<dyn Filter>, PbrtError> {
    match name {
        "box" => create_box_filter(params),
        "gaussian" => create_gaussian_filter(params),
        "mitchell" => create_mitchell_filter(params),
        "sinc" => create_sinc_filter(params),
        "triangle" => create_triangle_filter(params),
        _ => {
            let msg = format!("Filter \"{}\" unknown.", name);
            return Err(PbrtError::warning(&msg));
        }
    }
}
