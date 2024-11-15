use crate::core::shape::*;

use std::sync::Arc;

pub fn create_nurbs(
    _o2w: &Transform,
    _w2o: &Transform,
    _reverse_orientation: bool,
    _params: &ParamSet,
) -> Result<Vec<Arc<dyn Shape>>, PbrtError> {
    return Err(PbrtError::error("Not implemented yet."));
}
