use super::super::super::build::*;
use crate::core::pbrt::*;

use std::sync::Arc;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

#[derive(Debug, Clone, Copy)]
#[repr(C, align(16))]
pub struct SIMDBVHNode {
    pub bboxes: [[float32x4_t; 3]; 2], //768
    pub children: [usize; 4],
    pub axis_top: u8,
    pub axis_left: u8,
    pub axis_right: u8,
    pub is_leaf: u8,
}

/*
uint32_t _mm_movemask_ps(float32x4_t x) {
  uint32x4_t mmA = vandq_u32(
    vreinterpretq_u32_f32(x), (uint32x4_t) {0x1, 0x2, 0x4, 0x8}); // [0 1 2 3]
  uint32x4_t mmB = vextq_u32(mmA, mmA, 2);                        // [2 3 0 1]
  uint32x4_t mmC = vorrq_u32(mmA, mmB);                           // [0+2 1+3 0+2 1+3]
  uint32x4_t mmD = vextq_u32(mmC, mmC, 3);                        // [1+3 0+2 1+3 0+2]
  uint32x4_t mmE = vorrq_u32(mmC, mmD);                           // [0+1+2+3 ...]
  return vgetq_lane_u32(mmE, 0);
}
*/
#[inline]
unsafe fn get_mask(x: uint32x4_t) -> u32 {
    let mask_list: [u32; 4] = [0x1, 0x2, 0x4, 0x8];
    let mm_a = vandq_u32(x, vld1q_u32(&mask_list as *const u32)); // [0 1 2 3]
    let mm_b = vextq_u32(mm_a, mm_a, 2); // [2 3 0 1]
    let mm_c = vorrq_u32(mm_a, mm_b); // [0+2 1+3 0+2 1+3]
    let mm_d = vextq_u32(mm_c, mm_c, 3); // [1+3 0+2 1+3 0+2]
    let mm_e = vorrq_u32(mm_c, mm_d); // [0+1+2+3 ...]
    return vgetq_lane_u32(mm_e, 0);
}

unsafe fn test_aabb(
    bboxes: &[[float32x4_t; 3]; 2], //4boxes : min-max[2] of xyz[3] of boxes[4]
    org: &[float32x4_t; 3],         //ray origin
    idir: &[float32x4_t; 3],        //ray inveresed direction
    sign: &[usize; 3],              //ray xyz direction -> +:0,-:1
    tmin: float32x4_t,              //ray range tmin
    tmax: float32x4_t,              //ray range tmax
) -> u32 {
    let mut tmin = tmin;
    let mut tmax = tmax;
    // x coordinate
    tmin = vmaxq_f32(
        tmin,
        vmulq_f32(vsubq_f32(bboxes[sign[0]][0], org[0]), idir[0]),
    );
    tmax = vminq_f32(
        tmax,
        vmulq_f32(vsubq_f32(bboxes[1 - sign[0]][0], org[0]), idir[0]),
    );
    // y coordinate
    tmin = vmaxq_f32(
        tmin,
        vmulq_f32(vsubq_f32(bboxes[sign[1]][1], org[1]), idir[1]),
    );
    tmax = vminq_f32(
        tmax,
        vmulq_f32(vsubq_f32(bboxes[1 - sign[1]][1], org[1]), idir[1]),
    );
    // z coordinate
    tmin = vmaxq_f32(
        tmin,
        vmulq_f32(vsubq_f32(bboxes[sign[2]][2], org[2]), idir[2]),
    );
    tmax = vminq_f32(
        tmax,
        vmulq_f32(vsubq_f32(bboxes[1 - sign[2]][2], org[2]), idir[2]),
    );
    return get_mask(vcgeq_f32(tmax, tmin));
}

const EMPTY_MASK: usize = !0;

#[inline]
fn is_empty(i: usize) -> bool {
    return i == EMPTY_MASK;
}

impl Default for SIMDBVHNode {
    fn default() -> Self {
        unsafe {
            let bboxes = [[vdupq_n_f32(0.0); 3]; 2];
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
                    let l = [a, b, c, d];
                    nodes[offset].bboxes[m][j] = vld1q_f32(&l as *const f32);
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

fn get_sign(x: Float) -> usize {
    if x.is_sign_negative() {
        1
    } else {
        0
    }
}

fn is_valid_hitmask(
    order: u32,
    sign: &[usize; 3],
    axis_top: u8,
    axis_left: u8,
    axis_right: u8,
) -> bool {
    let org_order = order;
    let mut order = order;
    let mut indices = Vec::new();
    while (order & 0x4) == 0 {
        let cidx = order & 0x3;
        indices.push(cidx);
        order >>= 4;
    }
    indices.reverse();

    let sign_top = sign[axis_top as usize];
    let sign_left = sign[axis_left as usize];
    let sign_right = sign[axis_right as usize];
    let mut predicts = vec![0, 1, 2, 3];
    if sign_left != 0 {
        predicts.swap(0, 1);
    }
    if sign_right != 0 {
        predicts.swap(2, 3);
    }
    if sign_top != 0 {
        predicts.swap(0, 2);
        predicts.swap(1, 3);
    }
    let org_predicts = predicts.clone();
    for i in 0..indices.len() {
        let idx = indices[i];
        if predicts.len() == 0 {
            println!(
                "order: {:X}, indices: {:?}, predicts: {:?}, sign: {:?}, axis: [{}, {}, {}]",
                org_order, indices, org_predicts, sign, axis_top, axis_left, axis_right
            );
            return false;
        }
        while predicts.len() > 0 {
            let first = *predicts.first().unwrap();
            if first == idx {
                predicts.remove(0);
                break;
            } else {
                predicts.remove(0);
            }
        }
    }
    return true;
}

#[inline]
unsafe fn intersect_simd(
    primitives: &[Arc<dyn Primitive>],
    nodes: &[SIMDBVHNode],
    r: &Ray,
) -> Option<SurfaceInteraction> {
    let mut isect = None;
    let mut nodes_to_visit: Vec<usize> = Vec::with_capacity(16);

    let org: [float32x4_t; 3] = [
        vdupq_n_f32(r.o.x as f32),
        vdupq_n_f32(r.o.y as f32),
        vdupq_n_f32(r.o.z as f32),
    ];

    let idir: [float32x4_t; 3] = [
        vdupq_n_f32(r.d.x.recip() as f32),
        vdupq_n_f32(r.d.y.recip() as f32),
        vdupq_n_f32(r.d.z.recip() as f32),
    ];

    let sign = [get_sign(r.d.x), get_sign(r.d.y), get_sign(r.d.z)];

    let t = r.t_max.get();
    let tmin = vdupq_n_f32(0.0 as f32);
    let mut tmax = vdupq_n_f32(t as f32);

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

                if cfg!(debug_assertions) {
                    assert!(is_valid_hitmask(
                        order,
                        &sign,
                        node.axis_top,
                        node.axis_left,
                        node.axis_right
                    ));
                }

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
                tmax = vdupq_n_f32(t as f32);
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

    let org: [float32x4_t; 3] = [
        vdupq_n_f32(r.o.x as f32),
        vdupq_n_f32(r.o.y as f32),
        vdupq_n_f32(r.o.z as f32),
    ];

    let idir: [float32x4_t; 3] = [
        vdupq_n_f32(r.d.x.recip() as f32),
        vdupq_n_f32(r.d.y.recip() as f32),
        vdupq_n_f32(r.d.z.recip() as f32),
    ];

    let sign = [get_sign(r.d.x), get_sign(r.d.y), get_sign(r.d.z)];

    let t = r.t_max.get();
    let tmin = vdupq_n_f32(0.0 as f32);
    let tmax = vdupq_n_f32(t as f32);

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
        unsafe {
            return intersect_simd(&self.primitives, &self.nodes, r);
        }
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        unsafe {
            return intersect_simd_p(&self.primitives, &self.nodes, r);
        }
    }
}

impl Aggregate for QBVHAccel {}

#[cfg(all(test, target_arch = "aarch64"))]
mod tests {
    use super::*;
    use std::arch::aarch64::*;

    #[test]
    fn test_001() {
        unsafe {
            let tmin = vld1q_f32(&[0.0f32, 1.0, 0.0, 1.0] as *const f32);
            let tmax = vld1q_f32(&[1.0f32, 0.0, 1.0, 1.0] as *const f32);
            // tmin <= tmax
            let cmp = vcleq_f32(tmin, tmax);
            let mut cmp_u: [u32; 4] = [0; 4]; //std::mem::transmute(cmp);
            vst1q_u32(&mut cmp_u as *mut u32, vandq_u32(cmp, vdupq_n_u32(0x1)));

            print!("cmp: {:?}\n", cmp);
            print!("cmp_u: {:?}\n", cmp_u);
            assert_eq!(cmp_u, [1, 0, 1, 1]);

            let mask = get_mask(cmp);
            print!("mask: {:?}\n", mask);
            assert_eq!(mask, 0x1 | 0x4 | 0x8);
        }
    }
}
