use crate::core::shape::*;
use std::sync::Arc;

use super::create_triangle_mesh;

// doesn't handle flat out discontinuities in the curve...

#[derive(Debug, Default, Copy, Clone)]
struct Homogeneous3 {
    x: Float,
    y: Float,
    z: Float,
    w: Float,
}
impl Homogeneous3 {
    pub fn new(x: Float, y: Float, z: Float, w: Float) -> Self {
        Self { x, y, z, w }
    }
}

fn knot_offset(knot: &[Float], order: usize, t: Float) -> usize {
    let first_knot = order - 1;
    assert!(first_knot < knot.len());

    let mut knot_offset = first_knot as usize;
    while t >= knot[knot_offset + 1] {
        knot_offset += 1;
    }
    assert!(t >= knot[knot_offset] && t < knot[knot_offset + 1]);
    return knot_offset;
}

fn nurbs_evaluate(
    order: usize,
    knot: &[Float],
    cp: &[Homogeneous3],
    cp_stride: usize,
    t: Float,
) -> (Homogeneous3, Vector3f) {
    let knot_offset = knot_offset(knot, order, t);
    let knot = &knot[knot_offset..];

    assert!(knot_offset >= order - 1);
    let cp_offset = knot_offset - order + 1;
    assert!(cp_offset < cp.len());

    let mut cp_work = vec![Homogeneous3::default(); order as usize];
    for i in 0..order {
        cp_work[i as usize] = cp[((cp_offset + i) * cp_stride) as usize];
    }

    let order = order as usize;
    for i in 0..order - 2 {
        for j in 0..order - 1 - i {
            let i = i as usize;
            let j = j as usize;

            let alpha = (knot[1 + j] - t) / (knot[1 + j] - knot[j + 2 - order + i]);
            assert!(alpha >= 0.0 && alpha <= 1.0);
            cp_work[j].x = alpha * cp_work[j].x + (1.0 - alpha) * cp_work[j + 1].x;
            cp_work[j].y = alpha * cp_work[j].y + (1.0 - alpha) * cp_work[j + 1].y;
            cp_work[j].z = alpha * cp_work[j].z + (1.0 - alpha) * cp_work[j + 1].z;
            cp_work[j].w = alpha * cp_work[j].w + (1.0 - alpha) * cp_work[j + 1].w;
        }
    }

    let alpha = (knot[1] - t) / (knot[1] - knot[0]);
    assert!(alpha >= 0.0 && alpha <= 1.0);
    let x = alpha * cp_work[0].x + (1.0 - alpha) * cp_work[1].x;
    let y = alpha * cp_work[0].y + (1.0 - alpha) * cp_work[1].y;
    let z = alpha * cp_work[0].z + (1.0 - alpha) * cp_work[1].z;
    let w = alpha * cp_work[0].w + (1.0 - alpha) * cp_work[1].w;
    let val = Homogeneous3::new(x, y, z, w);

    let factor = (order - 1) as Float / (knot[1] - knot[0]);
    let dx = factor * (cp_work[1].x - cp_work[0].x);
    let dy = factor * (cp_work[1].y - cp_work[0].y);
    let dz = factor * (cp_work[1].z - cp_work[0].z);
    let dw = factor * (cp_work[1].w - cp_work[0].w);

    let dx = (dx / val.w) - (val.x * dw / (val.w * val.w));
    let dy = (dy / val.w) - (val.y * dw / (val.w * val.w));
    let dz = (dz / val.w) - (val.z * dw / (val.w * val.w));
    let deriv = Vector3f::new(dx, dy, dz);

    return (val, deriv);
}

fn nurbs_evaluate_surface(
    u_order: usize,
    u_knot: &[Float],
    u_cp: usize,
    u: Float,
    v_order: usize,
    v_knot: &[Float],
    v_cp: usize,
    v: Float,
    cp: &[Homogeneous3],
) -> (Point3f, Vector3f, Vector3f) {
    let mut iso = vec![Homogeneous3::default(); usize::max(u_order as usize, v_order as usize)];

    let u_offset = knot_offset(u_knot, u_order, u);
    assert!(u_offset >= u_order - 1);
    let u_first_cp = u_offset - u_order + 1;
    let u_first_cp = u_first_cp as usize;
    let u_order = u_order as usize;
    for i in 0..u_order {
        iso[i] = nurbs_evaluate(v_order, v_knot, &cp[(u_first_cp + i) as usize..], u_cp, v).0;
    }

    let v_offset = knot_offset(v_knot, v_order, v);
    assert!(v_offset >= v_order - 1);
    let v_first_cp = v_offset - v_order + 1;
    assert!(v_first_cp < v_cp);
    let v_first_cp = v_first_cp as usize;
    let v_order = v_order as usize;

    let (p, dpdu) = nurbs_evaluate(u_order, u_knot, &iso[v_first_cp..], 1, u);
    for i in 0..v_order {
        iso[i] = nurbs_evaluate(u_order, u_knot, &cp[((v_first_cp + i) * u_cp)..], 1, u).0;
    }
    let (_p, dpdv) = nurbs_evaluate(v_order, v_knot, &iso, 1, v);

    return (Point3f::new(p.x / p.w, p.y / p.w, p.z / p.w), dpdu, dpdv);
}

fn get_points(params: &ParamSet) -> Vec<Float> {
    let p = params.get_floats("P");
    if !p.is_empty() {
        return p;
    }
    let p = params.get_floats("Pw");
    if !p.is_empty() {
        return p;
    }
    return Vec::new();
}

pub fn create_nurbs(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
) -> Result<Vec<Arc<dyn Shape>>, PbrtError> {
    let nu = params.find_one_int("nu", -1);
    if nu == -1 {
        return Err(PbrtError::error(
            "Must provide number of control points \"nu\" with NURBS shape.",
        ));
    }
    let uorder = params.find_one_int("uorder", -1);
    if uorder == -1 {
        return Err(PbrtError::error(
            "Must provide u order \"uorder\" with NURBS shape.",
        ));
    }
    let uknots = params.get_floats("uknots");
    if uknots.is_empty() {
        return Err(PbrtError::error(
            "Must provide u knot vector \"uknots\" with NURBS shape.",
        ));
    }
    let nu = nu as usize;
    let uorder = uorder as usize;
    let nuknots = uknots.len();
    if nuknots != nu + uorder {
        let msg = format!(
            "Number of knots in u knot vector {} doesn't match sum of number of u control points {} and u order {}.",
            nuknots, nu, uorder
        );
        return Err(PbrtError::error(&msg));
    }

    let nv = params.find_one_int("nv", -1);
    if nv == -1 {
        return Err(PbrtError::error(
            "Must provide number of control points \"nv\" with NURBS shape.",
        ));
    }
    let vorder = params.find_one_int("vorder", -1);
    if vorder == -1 {
        return Err(PbrtError::error(
            "Must provide v order \"vorder\" with NURBS shape.",
        ));
    }
    let vknots = params.get_floats("vknots");
    if vknots.is_empty() {
        return Err(PbrtError::error(
            "Must provide v knot vector \"vknots\" with NURBS shape.",
        ));
    }
    let nv = nv as usize;
    let vorder = vorder as usize;
    let nvknots = vknots.len();
    if nvknots != nv + vorder {
        let msg = format!(
            "Number of knots in v knot vector {} doesn't match sum of number of v control points {} and v order {}.",
            nvknots, nv, vorder
        );
        return Err(PbrtError::error(&msg));
    }

    let p = get_points(params);
    if p.is_empty() {
        return Err(PbrtError::error(
            "Must provide control points via \"P\" or \"Pw\" parameter to NURBS shape.",
        ));
    }
    let is_homogeneous; //
    let mut npts = p.len();
    if npts % 3 == 0 {
        is_homogeneous = false;
        npts /= 3;
    } else if npts % 4 == 0 {
        is_homogeneous = true;
        npts /= 4;
    } else {
        return Err(PbrtError::error(
            "Number of control points must be multiple of 3 or 4.",
        ));
    }

    if npts != nu * nv {
        let msg = format!(
            "Number of control points {} doesn't match nu * nv = {} * {} = {}.",
            npts,
            nu,
            nv,
            nu * nv
        );
        return Err(PbrtError::error(&msg));
    }

    let mut pw = vec![Homogeneous3::default(); npts];
    if is_homogeneous {
        for i in 0..npts {
            pw[i] = Homogeneous3::new(p[4 * i], p[4 * i + 1], p[4 * i + 2], p[4 * i + 3]);
        }
    } else {
        for i in 0..npts {
            pw[i] = Homogeneous3::new(p[3 * i], p[3 * i + 1], p[3 * i + 2], 1.0);
        }
    }

    let u0: f32 = uknots[uorder as usize - 1];
    let u1: f32 = uknots[nu];
    let v0: f32 = vknots[vorder as usize - 1];
    let v1: f32 = vknots[nv];

    let u0 = params.find_one_float("u0", u0);
    let u1 = params.find_one_float("u1", u1);
    let v0 = params.find_one_float("v0", v0);
    let v1 = params.find_one_float("v1", v1);

    // Compute NURBS dicing rates
    let diceu = params.find_one_int("diceu", 30) as usize;
    let dicev = params.find_one_int("dicev", 30) as usize;

    let mut ueval = vec![0.0; diceu];
    let mut veval = vec![0.0; dicev];
    let mut eval_ps = vec![Point3f::default(); diceu * dicev];
    let mut eval_ns = vec![Normal3f::default(); diceu * dicev];
    let mut uvs = vec![Point2f::default(); diceu * dicev];
    for i in 0..diceu {
        ueval[i] = lerp(i as Float / (diceu - 1) as Float, u0, u1);
    }
    for i in 0..dicev {
        veval[i] = lerp(i as Float / (dicev - 1) as Float, v0, v1);
    }
    for v in 0..dicev {
        for u in 0..diceu {
            let uu = ueval[u];
            let vv = veval[v];
            uvs[v * diceu + u] = Point2f::new(uu, vv);
            let (p, dpdu, dpdv) =
                nurbs_evaluate_surface(uorder, &uknots, nu, uu, vorder, &vknots, nv, vv, &pw);

            eval_ps[v * diceu + u] = p;
            eval_ns[v * diceu + u] = Vector3f::cross(&dpdu, &dpdv).normalize();
        }
    }

    // Generate points-polygons mesh
    let ntris = 2 * (diceu - 1) * (dicev - 1);
    let mut vertex_indices = vec![0; 3 * ntris];
    let mut index = 0;
    let vn = |u: u32, v: u32| -> u32 { v * diceu as u32 + u };
    for v in 0..dicev - 1 {
        for u in 0..diceu - 1 {
            let u = u as u32;
            let v = v as u32;
            vertex_indices[index] = vn(u, v);
            vertex_indices[index + 1] = vn(u + 1, v);
            vertex_indices[index + 2] = vn(u + 1, v + 1);
            index += 3;

            vertex_indices[index] = vn(u, v);
            vertex_indices[index + 1] = vn(u + 1, v + 1);
            vertex_indices[index + 2] = vn(u, v + 1);
            index += 3;
        }
    }
    let params = ParamSet::new();
    let mesh = create_triangle_mesh(
        o2w,
        w2o,
        reverse_orientation,
        vertex_indices,
        eval_ps,
        Vec::new(),
        eval_ns,
        uvs,
        &params,
    );
    return Ok(mesh);
}
