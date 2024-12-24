use super::super::super::build::*;
use crate::core::pbrt::*;

use std::sync::Arc;

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[derive(Debug, Clone, Copy)]
#[repr(C, align(16))]
pub struct SIMDBVHNode {
    pub bboxes: [[__m128; 3]; 2], //768
    pub children: [usize; 4],
    pub axis_top: u8,
    pub axis_left: u8,
    pub axis_right: u8,
    pub is_leaf: u8,
}

unsafe fn test_aabb(
    bboxes: &[[__m128; 3]; 2], //4boxes : min-max[2] of xyz[3] of boxes[4]
    org: &[__m128; 3],         //ray origin
    idir: &[__m128; 3],        //ray inveresed direction
    sign: &[usize; 3],         //ray xyz direction -> +:0,-:1
    tmin: __m128,              //ray range tmin
    tmax: __m128,              //ray range tmax
) -> u32 {
    let mut tmin = tmin;
    let mut tmax = tmax;
    // x coordinate
    tmin = _mm_max_ps(
        tmin,
        _mm_mul_ps(_mm_sub_ps(bboxes[sign[0]][0], org[0]), idir[0]),
    );
    tmax = _mm_min_ps(
        tmax,
        _mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[0]][0], org[0]), idir[0]),
    );

    // y coordinate
    tmin = _mm_max_ps(
        tmin,
        _mm_mul_ps(_mm_sub_ps(bboxes[sign[1]][1], org[1]), idir[1]),
    );
    tmax = _mm_min_ps(
        tmax,
        _mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[1]][1], org[1]), idir[1]),
    );

    // z coordinate
    tmin = _mm_max_ps(
        tmin,
        _mm_mul_ps(_mm_sub_ps(bboxes[sign[2]][2], org[2]), idir[2]),
    );
    tmax = _mm_min_ps(
        tmax,
        _mm_mul_ps(_mm_sub_ps(bboxes[1 - sign[2]][2], org[2]), idir[2]),
    );

    return _mm_movemask_ps(_mm_cmpge_ps(tmax, tmin)) as u32;
}

const EMPTY_MASK: usize = !0;

#[inline]
fn is_empty(i: usize) -> bool {
    return i == EMPTY_MASK;
}

impl Default for SIMDBVHNode {
    fn default() -> Self {
        unsafe {
            let bboxes = [[_mm_setzero_ps(); 3]; 2];
            let children = [0, 0, 0, 0];
            SIMDBVHNode {
                bboxes,
                children,
                axis_top: 0,
                axis_left: 0,
                axis_right: 0,
                is_leaf: 0,
            }
        }
    }
}

fn flatten_qbvh_tree(nodes: &mut Vec<SIMDBVHNode>, node: &Arc<BVHBuildNode>) -> usize {
    let offset = nodes.len();
    nodes.push(SIMDBVHNode::default());
    if node.n_primitives > 0 {
        nodes[offset].is_leaf = 1;
        nodes[offset].children[0] = node.first_prim_offset;
        nodes[offset].children[1] = node.n_primitives;
    } else {
        let mut indices: [usize; 4] = [0; 4];
        let mut boxes: [[Vector3f; 2]; 4] = [[Vector3f::zero(); 2]; 4];

        let c0 = node.children[0].as_ref().unwrap();
        let c1 = node.children[1].as_ref().unwrap();
        if c0.n_primitives > 0 {
            indices[0] = flatten_qbvh_tree(nodes, c0);
            indices[1] = EMPTY_MASK;
            boxes[0][0] = c0.bounds.min;
            boxes[0][1] = c0.bounds.max;
            //boxes[1][0] = c0.bounds.min;
            //boxes[1][1] = c1.bounds.max;
        } else {
            let c00 = c0.children[0].as_ref().unwrap();
            let c01 = c0.children[1].as_ref().unwrap();
            indices[0] = flatten_qbvh_tree(nodes, c00);
            indices[1] = flatten_qbvh_tree(nodes, c01);

            boxes[0][0] = c00.bounds.min;
            boxes[0][1] = c00.bounds.max;
            boxes[1][0] = c01.bounds.min;
            boxes[1][1] = c01.bounds.max;
        }

        if c1.n_primitives > 0 {
            indices[2] = flatten_qbvh_tree(nodes, c1);
            indices[3] = EMPTY_MASK;

            boxes[2][0] = c1.bounds.min;
            boxes[2][1] = c1.bounds.max;
            //boxes[3][0] = c0.bounds.min;
            //boxes[3][1] = c1.bounds.max;
        } else {
            let c10 = c1.children[0].as_ref().unwrap();
            let c11 = c1.children[1].as_ref().unwrap();
            indices[2] = flatten_qbvh_tree(nodes, c10);
            indices[3] = flatten_qbvh_tree(nodes, c11);

            boxes[2][0] = c10.bounds.min;
            boxes[2][1] = c10.bounds.max;
            boxes[3][0] = c11.bounds.min;
            boxes[3][1] = c11.bounds.max;
        }

        //convert & swizzle
        let mut bboxes: [[[f32; 4]; 3]; 2] = [[[0.0; 4]; 3]; 2];
        for j in 0..3 {
            //xyz
            for k in 0..4 {
                bboxes[0][j][k] = boxes[k][0][j] as f32;
                bboxes[1][j][k] = boxes[k][1][j] as f32;
            }
        }

        unsafe {
            for m in 0..2 {
                for j in 0..3 {
                    let a = bboxes[m][j][0];
                    let b = bboxes[m][j][1];
                    let c = bboxes[m][j][2];
                    let d = bboxes[m][j][3];
                    nodes[offset].bboxes[m][j] = _mm_set_ps(d, c, b, a);
                }
            }
        }
        //for i in 0..4 {
        //    nodes[offset].children[i] = indices[i];
        //}
        nodes[offset].children = indices.clone();

        nodes[offset].axis_top = node.split_axis as u8;
        nodes[offset].axis_left = c0.split_axis as u8;
        nodes[offset].axis_right = c1.split_axis as u8;
    }
    return offset;
}

fn get_sign(x: Float) -> usize {
    if x.is_sign_negative() {
        1
    } else {
        0
    }
}

#[rustfmt::skip]
const ORDER_TABLE: [u32; 128] = [
    0x44444, 0x44444, 0x44444, 0x44444, 0x44444, 0x44444, 0x44444, 0x44444,
    0x44440, 0x44440, 0x44440, 0x44440, 0x44440, 0x44440, 0x44440, 0x44440,
    0x44441, 0x44441, 0x44441, 0x44441, 0x44441, 0x44441, 0x44441, 0x44441,
    0x44401, 0x44401, 0x44410, 0x44410, 0x44401, 0x44401, 0x44410, 0x44410,
    0x44442, 0x44442, 0x44442, 0x44442, 0x44442, 0x44442, 0x44442, 0x44442,
    0x44402, 0x44402, 0x44402, 0x44402, 0x44420, 0x44420, 0x44420, 0x44420,
    0x44412, 0x44412, 0x44412, 0x44412, 0x44421, 0x44421, 0x44421, 0x44421,
    0x44012, 0x44012, 0x44102, 0x44102, 0x44201, 0x44201, 0x44210, 0x44210,
    0x44443, 0x44443, 0x44443, 0x44443, 0x44443, 0x44443, 0x44443, 0x44443,
    0x44403, 0x44403, 0x44403, 0x44403, 0x44430, 0x44430, 0x44430, 0x44430,
    0x44413, 0x44413, 0x44413, 0x44413, 0x44431, 0x44431, 0x44431, 0x44431,
    0x44013, 0x44013, 0x44103, 0x44103, 0x44301, 0x44301, 0x44310, 0x44310,
    0x44423, 0x44432, 0x44423, 0x44432, 0x44423, 0x44432, 0x44423, 0x44432,
    0x44023, 0x44032, 0x44023, 0x44032, 0x44230, 0x44320, 0x44230, 0x44320,
    0x44123, 0x44132, 0x44123, 0x44132, 0x44231, 0x44321, 0x44231, 0x44321,
    0x40123, 0x40132, 0x41023, 0x41032, 0x42301, 0x43201, 0x42310, 0x43210,
];

#[inline]
fn intersect_primitives(primitives: &[Arc<dyn Primitive>], r: &Ray) -> Option<SurfaceInteraction> {
    let mut isect = None;
    for p in primitives.iter() {
        let prim = p.as_ref();
        if let Some(mut isect_n) = prim.intersect(r) {
            if prim.is_geometric() {
                isect_n.primitive = Some(Arc::downgrade(p));
            }
            isect = Some(isect_n);
        }
    }
    return isect;
}

#[inline]
fn intersect_primitives_p(primitives: &[Arc<dyn Primitive>], r: &Ray) -> bool {
    for p in primitives.iter() {
        let prim = p.as_ref();
        if prim.intersect_p(r) {
            return true;
        }
    }
    return false;
}

#[inline]
unsafe fn intersect_simd(
    primitives: &[Arc<dyn Primitive>],
    nodes: &[SIMDBVHNode],
    r: &Ray,
) -> Option<SurfaceInteraction> {
    let mut isect = None;
    let mut nodes_to_visit: Vec<usize> = Vec::with_capacity(16);

    let org: [__m128; 3] = [
        _mm_set1_ps(r.o.x as f32),
        _mm_set1_ps(r.o.y as f32),
        _mm_set1_ps(r.o.z as f32),
    ];
    let idir: [__m128; 3] = [
        _mm_set1_ps(r.d.x.recip() as f32),
        _mm_set1_ps(r.d.y.recip() as f32),
        _mm_set1_ps(r.d.z.recip() as f32),
    ];

    let sign = [get_sign(r.d.x), get_sign(r.d.y), get_sign(r.d.z)];

    let t = r.t_max.get();
    let tmin = _mm_set1_ps(0.0);
    let mut tmax = _mm_set1_ps(t as f32);

    nodes_to_visit.push(0);
    while let Some(current_node_index) = nodes_to_visit.pop() {
        assert!(!is_empty(current_node_index));
        let node = &nodes[current_node_index];
        if node.is_leaf == 0 {
            let hit_mask = test_aabb(&node.bboxes, &org, &idir, &sign, tmin, tmax) as usize;
            if hit_mask != 0 {
                let node_idx = (sign[node.axis_top as usize] << 2)
                    | (sign[node.axis_left as usize] << 1)
                    | (sign[node.axis_right as usize]);
                let mut order = ORDER_TABLE[hit_mask * 8 + node_idx];
                while (order & 0x4) == 0 {
                    let cidx = node.children[(order & 0x3) as usize];
                    if !is_empty(cidx) {
                        nodes_to_visit.push(cidx);
                    }
                    order >>= 4;
                }
            }
        } else {
            let start = node.children[0];
            let end = start + node.children[1];
            if let Some(isect_n) = intersect_primitives(&primitives[start..end], r) {
                let t = r.t_max.get();
                tmax = _mm_set1_ps(t as f32);
                isect = Some(isect_n);
            }
        }
    }
    return isect;
}

#[inline]
unsafe fn intersect_simd_p(
    primitives: &[Arc<dyn Primitive>],
    nodes: &[SIMDBVHNode],
    r: &Ray,
) -> bool {
    let mut nodes_to_visit: Vec<usize> = Vec::with_capacity(16);

    let org: [__m128; 3] = [
        _mm_set1_ps(r.o.x as f32),
        _mm_set1_ps(r.o.y as f32),
        _mm_set1_ps(r.o.z as f32),
    ];
    let idir: [__m128; 3] = [
        _mm_set1_ps(r.d.x.recip() as f32),
        _mm_set1_ps(r.d.y.recip() as f32),
        _mm_set1_ps(r.d.z.recip() as f32),
    ];

    let sign = [get_sign(r.d.x), get_sign(r.d.y), get_sign(r.d.z)];

    let t = r.t_max.get();
    let tmin = _mm_set1_ps(0.0);
    let tmax = _mm_set1_ps(t as f32);

    nodes_to_visit.push(0);
    while let Some(current_node_index) = nodes_to_visit.pop() {
        assert!(!is_empty(current_node_index));
        let node = &nodes[current_node_index];
        if node.is_leaf == 0 {
            let hit_mask = test_aabb(&node.bboxes, &org, &idir, &sign, tmin, tmax) as usize;
            if hit_mask != 0 {
                let node_idx = (sign[node.axis_top as usize] << 2)
                    | (sign[node.axis_left as usize] << 1)
                    | (sign[node.axis_right as usize]);
                let mut order = ORDER_TABLE[hit_mask * 8 + node_idx];
                while (order & 0x4) == 0 {
                    let cidx = node.children[(order & 0x3) as usize];
                    if !is_empty(cidx) {
                        nodes_to_visit.push(cidx);
                    }
                    order >>= 4;
                }
            }
        } else {
            let start = node.children[0];
            let end = start + node.children[1];
            if intersect_primitives_p(&primitives[start..end], r) {
                return true;
            }
        }
    }
    return false;
}

pub struct QBVHAccel {
    primitives: Vec<Arc<dyn Primitive>>,
    nodes: Vec<SIMDBVHNode>,
    bounds: Bounds3f,
}

impl QBVHAccel {
    pub fn new(
        prims: &[Arc<dyn Primitive>],
        max_prims_in_node: usize,
        split_method: SplitMethod,
    ) -> Self {
        let _p = ProfilePhase::new(Prof::AccelConstruction);

        let max_prims_in_node = usize::min(max_prims_in_node, 255);
        let mut orderd_prims = Vec::new();
        let root = create_bvh_node(&mut orderd_prims, prims, max_prims_in_node, split_method);
        let mut nodes = Vec::new();
        nodes.reserve(root.node_count());
        flatten_qbvh_tree(&mut nodes, &root);
        QBVHAccel {
            primitives: orderd_prims,
            nodes,
            bounds: root.bounds,
        }
    }
}

impl Primitive for QBVHAccel {
    fn world_bound(&self) -> Bounds3f {
        return self.bounds;
    }

    fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction> {
        let _p = ProfilePhase::new(Prof::AccelIntersect);

        unsafe {
            return intersect_simd(&self.primitives, &self.nodes, r);
        }
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        let _p = ProfilePhase::new(Prof::AccelIntersectP);

        unsafe {
            return intersect_simd_p(&self.primitives, &self.nodes, r);
        }
    }
}

impl Aggregate for QBVHAccel {}
