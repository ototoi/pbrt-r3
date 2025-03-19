use super::node::*;
use super::types::*;
use crate::core::geometry::*;
use crate::core::pbrt::*;

use std::sync::Arc;

#[derive(Debug, Clone, Copy, Default)]
struct BucketInfo {
    count: i32,
    bounds: Bounds3f,
}

fn get_bounds(primitive_info: &[BVHPrimitiveInfo]) -> Bounds3f {
    let mut bounds = primitive_info[0].bounds;
    for pi in primitive_info.iter() {
        bounds = Bounds3f::union(&bounds, &pi.bounds);
    }
    bounds
}

fn get_centroid_bounds(primitive_info: &[BVHPrimitiveInfo]) -> Bounds3f {
    let mut bounds = Bounds3f::new(&primitive_info[0].centroid, &primitive_info[0].centroid);
    for pi in primitive_info.iter() {
        bounds = Bounds3f::union_p(&bounds, &pi.centroid);
    }
    bounds
}

pub fn split_sah(
    dim: usize,
    primitive_info: &mut [BVHPrimitiveInfo],
    ordered_indices: &mut Vec<usize>,
    max_prims_in_node: usize,
    split_method: SplitMethod,
) -> Arc<BVHBuildNode> {
    // Partition primitives using approximate SAH
    let n_primitives = primitive_info.len();
    if n_primitives == 1 {
        // Create leaf _BVHBuildNode_
        let first_prim_offset = ordered_indices.len();
        for pi in primitive_info.iter() {
            ordered_indices.push(pi.primitive_number);
        }
        return Arc::new(BVHBuildNode::init_leaf(
            first_prim_offset,
            n_primitives,
            &primitive_info[0].bounds,
        ));
    } else if n_primitives == 2 {
        // Partition primitives into equally-sized subsets
        let mid = n_primitives / 2;
        primitive_info.sort_by(|a: &BVHPrimitiveInfo, b: &BVHPrimitiveInfo| {
            a.centroid[dim].partial_cmp(&b.centroid[dim]).unwrap()
        });
        let c0 = Some(recursive_build(
            &mut primitive_info[0..mid],
            ordered_indices,
            max_prims_in_node,
            split_method,
        ));
        let c1 = Some(recursive_build(
            &mut primitive_info[mid..],
            ordered_indices,
            max_prims_in_node,
            split_method,
        ));
        return Arc::new(BVHBuildNode::init_interior(dim, c0, c1));
    } else {
        // Compute bounds of all primitives in BVH node
        let bounds = get_bounds(primitive_info);

        // Compute bound of primitive centroids, choose split dimension _dim_
        let centroid_bounds = get_centroid_bounds(primitive_info);

        // Allocate _BucketInfo_ for SAH partition buckets
        const N_BUCKETS: usize = 12;
        let mut buckets = [BucketInfo::default(); N_BUCKETS];

        // Initialize _BucketInfo_ for SAH partition buckets
        for i in 0..primitive_info.len() {
            let b = (Float::floor(
                N_BUCKETS as Float * centroid_bounds.offset(&primitive_info[i].centroid)[dim],
            ) as i32)
                .min(N_BUCKETS as i32 - 1)
                .max(0) as usize;
            buckets[b].count += 1;
            buckets[b].bounds = Bounds3f::union(&buckets[b].bounds, &primitive_info[i].bounds);
        }

        // Compute costs for splitting after each bucket
        let mut cost = [0.0; N_BUCKETS - 1];
        for i in 0..cost.len() {
            let mut count0 = 0;
            let mut count1 = 0;
            let mut b0 = buckets[i + 0].bounds;
            let mut b1 = buckets[i + 1].bounds;

            for j in 0..=i {
                b0 = Bounds3f::union(&b0, &buckets[j].bounds);
                count0 += buckets[j].count;
            }
            for j in i + 1..N_BUCKETS {
                b1 = Bounds3f::union(&b1, &buckets[j].bounds);
                count1 += buckets[j].count;
            }
            cost[i] = 1.0
                + (count0 as Float * b0.surface_area() + count1 as Float * b1.surface_area())
                    / bounds.surface_area();
        }
        // Find bucket to split at that minimizes SAH metric
        let mut min_cost = cost[0];
        let mut min_cost_split_bucket = 0;
        for i in 1..cost.len() {
            if cost[i] < min_cost {
                min_cost = cost[i];
                min_cost_split_bucket = i;
            }
        }

        // Either create leaf or split primitives at selected SAH
        // bucket
        let leaf_cost = n_primitives as Float;
        if n_primitives > max_prims_in_node || min_cost < leaf_cost {
            let (mut left, mut right): (Vec<BVHPrimitiveInfo>, Vec<BVHPrimitiveInfo>) =
                primitive_info.iter().partition(|pi| -> bool {
                    let b = (Float::floor(
                        N_BUCKETS as Float * centroid_bounds.offset(&pi.centroid)[dim],
                    ) as i32)
                        .min(N_BUCKETS as i32 - 1)
                        .max(0) as usize;
                    b <= min_cost_split_bucket
                });
            let c0 = Some(recursive_build(
                &mut left,
                ordered_indices,
                max_prims_in_node,
                split_method,
            ));
            let c1 = Some(recursive_build(
                &mut right,
                ordered_indices,
                max_prims_in_node,
                split_method,
            ));
            return Arc::new(BVHBuildNode::init_interior(dim, c0, c1));
        } else {
            // Create leaf _BVHBuildNode_
            let first_prim_offset = ordered_indices.len();
            for i in 0..n_primitives {
                ordered_indices.push(primitive_info[i].primitive_number);
            }
            return Arc::new(BVHBuildNode::init_leaf(
                first_prim_offset,
                n_primitives,
                &bounds,
            ));
        }
    }
}
