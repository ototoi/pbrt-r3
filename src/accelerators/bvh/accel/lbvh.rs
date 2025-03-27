use super::super::build::*;
use crate::core::base::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::primitive::*;
use crate::core::profile::*;
use crate::core::stats::*;

use log::*;
use std::sync::Arc;

thread_local!(static TREE_BYTES: StatMemoryCounter = StatMemoryCounter::new("Memory/BVH tree"));

pub struct LinearBVHNode {
    pub bounds: [[Float; 3]; 2], //32*3*2
    pub offset: usize,
    pub n_primitives: usize, // 0 -> interior node
    pub axis: u8,            // interior node: xyz
}

fn to_bounds3f(b: &[[Float; 3]; 2]) -> Bounds3f {
    return Bounds3f::from((
        (b[0][0] as Float, b[0][1] as Float, b[0][2] as Float),
        (b[1][0] as Float, b[1][1] as Float, b[1][2] as Float),
    ));
}

fn to_array_bounds3f(b: &Bounds3f) -> [[Float; 3]; 2] {
    return [
        [b.min.x as Float, b.min.y as Float, b.min.z as Float],
        [b.max.x as Float, b.max.y as Float, b.max.z as Float],
    ];
}

impl LinearBVHNode {
    pub fn new(node: &BVHBuildNode) -> Self {
        /*
        let mut has_children = 0;
        if node.children[0].is_some() {
            has_children |= 1;
        }
        if node.children[1].is_some() {
            has_children |= 2;
        }
        */
        LinearBVHNode {
            bounds: to_array_bounds3f(&node.bounds),
            offset: 0,
            n_primitives: node.n_primitives, // 0 -> interior node
            axis: node.split_axis as u8,     // interior node: xyz
        }
    }
}

fn flatten_bvh_tree(nodes: &mut Vec<LinearBVHNode>, node: &Arc<BVHBuildNode>) -> usize {
    let mut linear_node = LinearBVHNode::new(node);
    let offset = nodes.len();
    if node.n_primitives > 0 {
        linear_node.offset = node.first_prim_offset;
        linear_node.n_primitives = node.n_primitives;
        nodes.push(linear_node);
    } else {
        linear_node.n_primitives = 0;
        nodes.push(linear_node);
        if let Some(c) = node.children[0].as_ref() {
            let _ = flatten_bvh_tree(nodes, c);
        }
        if let Some(c) = node.children[1].as_ref() {
            nodes[offset].offset = flatten_bvh_tree(nodes, c);
        }
    }
    return offset;
}

fn sign(x: Float) -> usize {
    if x.is_sign_negative() {
        1
    } else {
        0
    }
}

pub struct LBVHAccel {
    primitives: Vec<Arc<dyn Primitive>>,
    nodes: Vec<LinearBVHNode>,
}

impl LBVHAccel {
    pub fn new(
        prims: &[Arc<dyn Primitive>],
        max_prims_in_node: usize,
        split_method: SplitMethod,
    ) -> Self {
        let max_prims_in_node = usize::min(max_prims_in_node, 255);
        let mut orderd_prims = Vec::new();
        let root = create_bvh_node(&mut orderd_prims, prims, max_prims_in_node, split_method);

        let total_nodes = root.node_count();
        let allocated_memory = total_nodes * std::mem::size_of::<LinearBVHNode>();
        info!(
            "BVH created with {} nodes for {} primitives ({:.2} MB)",
            total_nodes,
            prims.len(),
            allocated_memory as Float / (1024.0 * 1024.0)
        );

        // Compute representation of depth-first traversal of BVH tree
        let tree_bytes = total_nodes * std::mem::size_of::<LinearBVHNode>()
            + std::mem::size_of::<Self>()
            + prims.len() * std::mem::size_of::<Arc<dyn Primitive>>();
        TREE_BYTES.with(|c| c.add(tree_bytes));

        let mut nodes = Vec::new();
        nodes.reserve(root.node_count());
        flatten_bvh_tree(&mut nodes, &root);
        LBVHAccel {
            primitives: orderd_prims,
            nodes,
        }
    }
}

impl Primitive for LBVHAccel {
    fn world_bound(&self) -> Bounds3f {
        return to_bounds3f(&(self.nodes[0].bounds));
    }

    fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction> {
        let _p = ProfilePhase::new(Prof::AccelIntersect);

        let mut isect = None;
        let mut nodes_to_visit: Vec<(usize, Float, Float)> = Vec::with_capacity(16);
        let mut t_max = r.t_max.get();
        let t0: Float = 0.0;
        let t1: Float = t_max;
        let org = [r.o.x, r.o.y, r.o.z];
        let idir = [
            Float::recip(r.d.x),
            Float::recip(r.d.y),
            Float::recip(r.d.z),
        ];
        let dir_is_neg = [sign(idir[0]), sign(idir[1]), sign(idir[2])];
        nodes_to_visit.push((0, t0, t1));
        while let Some((current_node_index, t0, mut t1)) = nodes_to_visit.pop() {
            t1 = Float::min(t_max, t1);
            if t1 < t0 {
                continue;
            }
            assert!(t0 <= t1);
            let node = &self.nodes[current_node_index];
            let min = &node.bounds[0];
            let max = &node.bounds[1];
            if let Some((t0, t1)) =
                intersect_box_array_i(min, max, &org, &idir, &dir_is_neg, t0, t1)
            {
                if node.n_primitives > 0 {
                    let start = node.offset as usize;
                    let end = start + (node.n_primitives as usize);
                    for i in start..end {
                        let prim = self.primitives[i].as_ref();
                        if let Some(mut isect_n) = prim.intersect(r) {
                            if prim.is_geometric() {
                                isect_n.primitive = Some(Arc::downgrade(&self.primitives[i]));
                            }
                            t_max = r.t_max.get();
                            isect = Some(isect_n);
                        }
                    }
                } else {
                    //let has_children = [(node.has_children & 1) != 0, (node.has_children & 2) != 0];
                    let index_children = [current_node_index + 1, node.offset as usize];
                    let indices = [
                        dir_is_neg[node.axis as usize],
                        1 - dir_is_neg[node.axis as usize],
                    ];
                    //if has_children[indices[1]] {
                    nodes_to_visit.push((index_children[indices[1]], t0, t1));
                    //}
                    //if has_children[indices[0]] {
                    nodes_to_visit.push((index_children[indices[0]], t0, t1));
                    //}
                }
            }
        }
        return isect;
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        let _p = ProfilePhase::new(Prof::AccelIntersectP);

        let mut nodes_to_visit: Vec<(usize, Float, Float)> = Vec::with_capacity(16);
        let t_max = r.t_max.get();
        let t0: Float = 0.0;
        let t1: Float = t_max;
        let org = [r.o.x, r.o.y, r.o.z];
        let idir = [
            Float::recip(r.d.x),
            Float::recip(r.d.y),
            Float::recip(r.d.z),
        ];
        let dir_is_neg = [sign(idir[0]), sign(idir[1]), sign(idir[2])];
        nodes_to_visit.push((0, t0, t1));
        while let Some((current_node_index, t0, t1)) = nodes_to_visit.pop() {
            //t1 = Float::min(t_max, t1);
            //if t1 < t0 {
            //    continue;
            //}
            //assert!(t0 <= t1);
            let node = &self.nodes[current_node_index];
            let min = &node.bounds[0];
            let max = &node.bounds[1];
            if let Some((t0, t1)) =
                intersect_box_array_i(min, max, &org, &idir, &dir_is_neg, t0, t1)
            {
                if node.n_primitives > 0 {
                    let start = node.offset as usize;
                    let end = start + (node.n_primitives as usize);
                    for i in start..end {
                        let p = self.primitives[i].as_ref();
                        if p.intersect_p(r) {
                            return true;
                        }
                    }
                } else {
                    //let has_children = [(node.has_children & 1) != 0, (node.has_children & 2) != 0];
                    let index_children = [current_node_index + 1, node.offset as usize];
                    let indices = [
                        dir_is_neg[node.axis as usize],
                        1 - dir_is_neg[node.axis as usize],
                    ];
                    //if has_children[indices[1]] {
                    nodes_to_visit.push((index_children[indices[1]], t0, t1));
                    //}
                    //if has_children[indices[0]] {
                    nodes_to_visit.push((index_children[indices[0]], t0, t1));
                    //}
                }
            }
        }
        return false;
    }
}

impl Aggregate for LBVHAccel {}
