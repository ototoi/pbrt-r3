use super::node::*;
use super::types::*;

use std::sync::Arc;

pub fn split_equal_counts(
    dim: usize,
    primitive_info: &mut [BVHPrimitiveInfo],
    ordered_indices: &mut Vec<usize>,
    max_prims_in_node: usize,
    split_method: SplitMethod,
) -> Arc<BVHBuildNode> {
    let info = primitive_info;
    info.sort_by(|a, b| a.centroid[dim].partial_cmp(&b.centroid[dim]).unwrap());
    let mid = info.len() / 2;

    let c0 = Some(recursive_build(
        &mut info[0..mid],
        ordered_indices,
        max_prims_in_node,
        split_method,
    ));
    let c1 = Some(recursive_build(
        &mut info[mid..],
        ordered_indices,
        max_prims_in_node,
        split_method,
    ));
    return Arc::new(BVHBuildNode::init_interior(dim, c0, c1));
}
