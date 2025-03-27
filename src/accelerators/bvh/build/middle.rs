use super::equal_counts::*;
use super::node::*;
use super::types::*;
use crate::core::base::*;

use std::sync::Arc;

pub fn split_middle(
    dim: usize,
    p_mid: Float,
    primitive_info: &mut [BVHPrimitiveInfo],
    ordered_indices: &mut Vec<usize>,
    max_prims_in_node: usize,
    split_method: SplitMethod,
) -> Arc<BVHBuildNode> {
    let start = 0;
    let end = primitive_info.len();
    let info = primitive_info;
    info.sort_by(|a, b| a.centroid[dim].partial_cmp(&b.centroid[dim]).unwrap());
    let mid = info.iter().position(|x| p_mid <= x.centroid[dim]).unwrap();
    if mid == start || mid == end {
        return split_equal_counts(dim, info, ordered_indices, max_prims_in_node, split_method);
    } else {
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
}
