use super::accel::*;
use super::build::*;
use crate::core::pbrt::*;

use std::sync::Arc;

pub fn create_bvh_accelerator(
    prims: &[Arc<dyn Primitive>],
    params: &ParamSet,
) -> Result<Arc<dyn Primitive>, PbrtError> {
    let split_name = params.find_one_string("splitmethod", "middle");
    let split_method = match &split_name as &str {
        "sah" => SplitMethod::SAH,
        "hlbvh" => SplitMethod::HLBVH,
        "middle" => SplitMethod::Middle,
        "equal" => SplitMethod::EqualCounts,
        _ => SplitMethod::SAH,
    };

    let max_prims_in_node: usize = params.find_one_int("maxnodeprims", 4) as usize;
    return Ok(Arc::new(BVHAccel::new(
        prims,
        max_prims_in_node,
        split_method,
    )));
}
