use super::*;
use crate::core::pbrt::*;

use std::sync::Arc;

pub fn create_accelerator(
    name: &str,
    prims: &[Arc<dyn Primitive>],
    params: &ParamSet,
) -> Result<Arc<dyn Primitive>, PbrtError> {
    if !prims.is_empty() {
        match name {
            "bvh" => create_bvh_accelerator(prims, params),
            "kdtree" => create_kdtree_accelerator(prims, params),
            "exhaustive" => create_exhaustive_accelerator(prims, params),
            _ => {
                return Err(PbrtError::warning(&format!(
                    "Accelerator \"{}\" unknown.",
                    name
                )));
            }
        }
    } else {
        return Err(PbrtError::error(""));
    }
}
