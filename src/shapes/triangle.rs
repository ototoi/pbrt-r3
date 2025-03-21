use super::alphamask::*;
use crate::core::error::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::sampling::*;
use crate::core::shape::*;
use crate::core::stats::*;
use crate::core::texture::*;

use std::collections::HashMap;
use std::sync::Arc;

thread_local!(static TESTS: StatPercent = StatPercent::new("Intersections/Ray-triangle intersection tests"));
thread_local!(static TRI_MESH_BYTES: StatMemoryCounter = StatMemoryCounter::new("Memory/Triangle meshes"));

pub struct TriangleMesh {
    pub object_to_world: Transform,
    pub world_to_object: Transform,
    pub reverse_orientation: bool,
    pub swaps_handedness: bool,
    pub two_sided: bool,
    //pub vertex_indices: Vec<usize>,
    pub p: Vec<Point3f>,
    pub s: Vec<Vector3f>,
    pub n: Vec<Vector3f>,
    pub uv: Vec<Point2f>,
    //pub face_index: Vec<usize>,
}

/*
const Transform &ObjectToWorld, int nTriangles,
const int *vertexIndices, int nVertices, const Point3f *P,
const Vector3f *S, const Normal3f *N, const Point2f *uv

*/
const MACHINE_EPSILON: Float = Float::EPSILON * 0.5;
const GAMMA2: Float = (2.0 * MACHINE_EPSILON) / (1.0 - (2.0 * MACHINE_EPSILON));
const GAMMA3: Float = (3.0 * MACHINE_EPSILON) / (1.0 - (3.0 * MACHINE_EPSILON));
const GAMMA5: Float = (5.0 * MACHINE_EPSILON) / (1.0 - (5.0 * MACHINE_EPSILON));
const GAMMA6: Float = (6.0 * MACHINE_EPSILON) / (1.0 - (5.0 * MACHINE_EPSILON));
const GAMMA7: Float = (7.0 * MACHINE_EPSILON) / (1.0 - (7.0 * MACHINE_EPSILON));
const TRI: [usize; 4] = [0, 1, 2, 0];

impl TriangleMesh {
    pub fn new(
        object_to_world: &Transform,
        reverse_orientation: bool,
        two_sided: bool,
        //vertex_indices: Vec<usize>,
        p: Vec<Point3f>,
        s: Vec<Vector3f>,
        n: Vec<Normal3f>,
        uv: Vec<Point2f>,
    ) -> Self {
        let p: Vec<Point3f> = p
            .iter()
            .map(|p| -> Point3f {
                return object_to_world.transform_point(p);
            })
            .collect();
        let s: Vec<Vector3f> = s
            .iter()
            .map(|s| -> Vector3f {
                return object_to_world.transform_vector(s);
            })
            .collect();
        let n: Vec<Normal3f> = n
            .iter()
            .map(|n| -> Normal3f {
                return object_to_world.transform_normal(n);
            })
            .collect();
        let swaps_handedness = object_to_world.swaps_handedness();
        TriangleMesh {
            object_to_world: *object_to_world,
            world_to_object: object_to_world.inverse(),
            reverse_orientation,
            swaps_handedness,
            two_sided,

            //vertex_indices: vertex_indices,
            p,
            s,
            n,
            uv,
        }
    }

    pub fn calc_normal(&self, dpdu: &Vector3f, dpdv: &Vector3f) -> Normal3f {
        let mut n = Vector3f::cross(dpdu, dpdv).normalize();
        if self.reverse_orientation ^ self.swaps_handedness {
            n *= -1.0;
        }
        return n;
    }
}

pub struct Triangle {
    pub mesh: Arc<TriangleMesh>,
    pub v: [u32; 3],
    pub face_index: usize,
}

impl Triangle {
    pub fn new(mesh: &Arc<TriangleMesh>, v: &[u32; 3], face_index: usize) -> Self {
        TRI_MESH_BYTES.with(|s| {
            s.add(std::mem::size_of::<Triangle>());
        });

        Triangle {
            mesh: Arc::clone(mesh),
            v: *v,
            face_index,
        }
    }

    pub fn calc_normal(&self, dpdu: &Vector3f, dpdv: &Vector3f) -> Normal3f {
        return self.mesh.as_ref().calc_normal(dpdu, dpdv);
    }

    pub fn get_uvs(&self) -> [Vector2f; 3] {
        let mesh = self.mesh.as_ref();
        if !mesh.uv.is_empty() {
            return [
                mesh.uv[self.v[0] as usize],
                mesh.uv[self.v[1] as usize],
                mesh.uv[self.v[2] as usize],
            ];
        } else {
            return [
                Point2f::new(0.0, 0.0),
                Point2f::new(1.0, 0.0),
                Point2f::new(1.0, 1.0),
            ];
        }
    }

    fn get_dpdu_dpdv_from_uv(
        &self,
        p0: &Vector3f,
        p1: &Vector3f,
        p2: &Vector3f,
    ) -> Option<([Vector2f; 3], Vector3f, Vector3f)> {
        // Handle the case where there are no UVs
        // pbrt-r3
        //if self.mesh.uv.is_empty() {
        //    return None;
        //}
        // pbrt-r3
        let uv = self.get_uvs();
        // Compute deltas for triangle partial derivatives
        let duv02 = uv[0] - uv[2];
        let duv12 = uv[1] - uv[2];
        let dp02 = *p0 - *p2;
        let dp12 = *p1 - *p2;
        let determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
        let degenerate_uv = Float::abs(determinant) < 1e-8;
        if !degenerate_uv {
            let invdet = 1.0 / determinant;
            let dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
            let dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
            if Vector3f::cross(&dpdu, &dpdv).length_squared() <= 0.0 {
                return None;
            } else {
                return Some((uv, dpdu, dpdv));
            }
        }
        return None;
    }

    pub fn get_dpdu_dpdv(
        &self,
        p0: &Vector3f,
        p1: &Vector3f,
        p2: &Vector3f,
    ) -> Option<([Vector2f; 3], Vector3f, Vector3f)> {
        if let Some((uv, dpdu, dpdv)) = self.get_dpdu_dpdv_from_uv(p0, p1, p2) {
            return Some((uv, dpdu, dpdv));
        } else {
            // Handle zero determinant for triangle partial derivative matrix
            let ng = Vector3f::cross(&(*p2 - *p0), &(*p1 - *p0));
            if ng.length_squared() <= 0.0 {
                // The triangle is actually degenerate; the intersection is
                // bogus.
                return None;
            } else {
                let uv = self.get_uvs();
                let (dpdu, dpdv) = coordinate_system(&ng.normalize());
                return Some((uv, dpdu, dpdv));
            }
        }
    }
}

fn union3(p0: Point3f, p1: Point3f, p2: Point3f) -> Bounds3f {
    let a = [[p1.x, p1.y, p1.z], [p2.x, p2.y, p2.z]];
    let mut min: [Float; 3] = [p0.x, p0.y, p0.z];
    let mut max: [Float; 3] = [p0.x, p0.y, p0.z];
    for j in 0..2 {
        for i in 0..3 {
            min[i] = Float::min(min[i], a[j][i]);
            max[i] = Float::max(max[i], a[j][i]);
        }
    }
    return Bounds3f::from(((min[0], min[1], min[2]), (max[0], max[1], max[2])));
}

impl Shape for Triangle {
    fn object_bound(&self) -> Bounds3f {
        let mesh = self.mesh.as_ref();
        let world_to_object = mesh.world_to_object;
        let i0 = self.v[0] as usize;
        let i1 = self.v[1] as usize;
        let i2 = self.v[2] as usize;
        let p0 = world_to_object.transform_point(&mesh.p[i0]);
        let p1 = world_to_object.transform_point(&mesh.p[i1]);
        let p2 = world_to_object.transform_point(&mesh.p[i2]);
        return union3(p0, p1, p2);
    }

    fn world_bound(&self) -> Bounds3f {
        let mesh = self.mesh.as_ref();
        let i0 = self.v[0] as usize;
        let i1 = self.v[1] as usize;
        let i2 = self.v[2] as usize;
        let p0 = mesh.p[i0];
        let p1 = mesh.p[i1];
        let p2 = mesh.p[i2];
        return union3(p0, p1, p2);
    }

    fn intersect(&self, r: &Ray) -> Option<(Float, SurfaceInteraction)> {
        let _p = ProfilePhase::new(Prof::TriIntersect);
        TESTS.with(|stat| {
            stat.add_denom(1);
        });

        let mesh = self.mesh.as_ref();
        let i0 = self.v[0] as usize;
        let i1 = self.v[1] as usize;
        let i2 = self.v[2] as usize;
        let p0 = mesh.p[i0];
        let p1 = mesh.p[i1];
        let p2 = mesh.p[i2];

        let dp02 = p0 - p2;
        let dp12 = p1 - p2;
        let mut n = Vector3f::cross(&dp02, &dp12);
        if mesh.reverse_orientation ^ mesh.swaps_handedness {
            n *= -1.0;
        }

        if !mesh.two_sided {
            if Vector3f::dot(&n, &r.d) >= 0.0 {
                return None;
            }
        }

        n = n.normalize();

        let mut p0t = p0 - r.o;
        let mut p1t = p1 - r.o;
        let mut p2t = p2 - r.o;

        let kz = max_dimension(&(r.d.abs())) as usize;
        let kx = TRI[kz + 1];
        let ky = TRI[kx + 1];

        let d = permute(&r.d, kx, ky, kz);
        p0t = permute(&p0t, kx, ky, kz);
        p1t = permute(&p1t, kx, ky, kz);
        p2t = permute(&p2t, kx, ky, kz);

        let sx = -d.x / d.z;
        let sy = -d.y / d.z;
        let sz = 1.0 / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;

        // Compute edge function coefficients _e0_, _e1_, and _e2_
        let mut e0 = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1 = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2 = p0t.x * p1t.y - p0t.y * p1t.x;

        if e0 == 0.0 || e1 == 0.0 || e2 == 0.0 {
            let p2txp1ty = p2t.x as f64 * p1t.y as f64;
            let p2typ1tx = p2t.y as f64 * p1t.x as f64;
            e0 = (p2typ1tx - p2txp1ty) as Float;
            let p0txp2ty = p0t.x as f64 * p2t.y as f64;
            let p0typ2tx = p0t.y as f64 * p2t.x as f64;
            e1 = (p0typ2tx - p0txp2ty) as Float;
            let p1txp0ty = p1t.x as f64 * p0t.y as f64;
            let p1typ0tx = p1t.y as f64 * p0t.x as f64;
            e2 = (p1typ0tx - p1txp0ty) as Float;
        }

        // Perform triangle edge and determinant tests
        if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0) && (e0 > 0.0 || e1 > 0.0 || e2 > 0.0) {
            return None;
        }
        let det = e0 + e1 + e2;
        if det == 0.0 {
            return None;
        }

        // Compute scaled hit distance to triangle and test against ray $t$ range
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;

        let t_max = r.t_max.get();

        let t_scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        #[allow(clippy::if_same_then_else)]
        if det < 0.0 && (t_scaled >= 0.0 || t_scaled < t_max * det) {
            return None;
        } else if det > 0.0 && (t_scaled <= 0.0 || t_scaled > t_max * det) {
            return None;
        }

        // Compute barycentric coordinates and $t$ value for triangle intersection
        let inv_det = 1.0 / det;
        let b0 = e0 * inv_det;
        let b1 = e1 * inv_det;
        let b2 = e2 * inv_det;
        let t = t_scaled * inv_det;

        // Ensure that computed triangle $t$ is conservatively greater than zero

        // Compute $\delta_z$ term for triangle $t$ error bounds
        let max_zt = max_component(&Vector3f::new(p0t.z, p1t.z, p2t.z).abs());
        let delta_z = GAMMA3 * max_zt;

        // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        let max_xt = max_component(&Vector3f::new(p0t.x, p1t.x, p2t.x).abs());
        let max_yt = max_component(&Vector3f::new(p0t.y, p1t.y, p2t.y).abs());
        let delta_x = GAMMA5 * (max_xt + max_zt);
        let delta_y = GAMMA5 * (max_yt + max_zt);

        // Compute $\delta_e$ term for triangle $t$ error bounds
        let delta_e = 2.0 * (GAMMA2 * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);
        // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        let max_e = max_component(&Vector3f::new(e0, e1, e2).abs());
        let delta_t = 3.0
            * (GAMMA3 * max_e * max_zt + delta_e * max_zt + delta_z * max_e)
            * Float::abs(inv_det);
        if t <= delta_t {
            return None;
        }

        let (uv, dpdu, dpdv) = self.get_dpdu_dpdv(&p0, &p1, &p2)?;
        // Compute error bounds for triangle intersection
        let x_abs_sum = Float::abs(b0 * p0.x) + Float::abs(b1 * p1.x) + Float::abs(b2 * p2.x);
        let y_abs_sum = Float::abs(b0 * p0.y) + Float::abs(b1 * p1.y) + Float::abs(b2 * p2.y);
        let z_abs_sum = Float::abs(b0 * p0.z) + Float::abs(b1 * p1.z) + Float::abs(b2 * p2.z);
        let p_error = GAMMA7 * Vector3f::new(x_abs_sum, y_abs_sum, z_abs_sum);

        // Interpolate $(u,v)$ parametric coordinates and hit point
        let p_hit = b0 * p0 + b1 * p1 + b2 * p2;
        let uv_hit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

        // Fill in _SurfaceInteraction_ from triangle hit
        let mut isect = SurfaceInteraction::new(
            &p_hit,
            &p_error,
            &uv_hit,
            &(-r.d),
            &n,
            &dpdu,
            &dpdv,
            &Normal3f::zero(),
            &Normal3f::zero(),
            r.time,
            self.face_index as u32,
        );
        isect.shading.n = n;

        if !mesh.n.is_empty() || !mesh.s.is_empty() {
            let mut ns = isect.n;
            if !mesh.n.is_empty() {
                let nns = b0 * mesh.n[i0] + b1 * mesh.n[i1] + b2 * mesh.n[i2];
                if nns.length_squared() > 0.0 {
                    ns = nns.normalize();
                }
            }

            let mut ss = isect.dpdu.normalize();
            if !mesh.s.is_empty() {
                let nns = b0 * mesh.s[i0] + b1 * mesh.s[i1] + b2 * mesh.s[i2];
                if nns.length_squared() > 0.0 {
                    ss = nns.normalize();
                }
            }

            //ns x ss -> ts
            //ss x ts -> ns
            //ts x ns -> ss
            //let mut ts = Vector3f::cross(&ss, &ns); //pbrt-v3
            let mut ts = Vector3f::cross(&ns, &ss); //zx->y //pbrt-r3
            if ts.length_squared() > 0.0 {
                ts = ts.normalize();
                ss = Vector3f::cross(&ts, &ns).normalize(); //yz->x
            } else {
                let (ss1, ts1) = coordinate_system(&ns);
                ss = ss1;
                ts = ts1;
            }

            let mut dndu = Vector3f::zero();
            let mut dndv = Vector3f::zero();
            if !mesh.n.is_empty() {
                let duv02 = uv[0] - uv[2];
                let duv12 = uv[1] - uv[2];
                let dn1 = mesh.n[i0] - mesh.n[i2];
                let dn2 = mesh.n[i1] - mesh.n[i2];
                let determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
                let degenerate_uv = Float::abs(determinant) < 1e-8;
                if degenerate_uv {
                    // We can still compute dndu and dndv, with respect to the
                    // same arbitrary coordinate system we use to compute dpdu
                    // and dpdv when this happens. It's important to do this
                    // (rather than giving up) so that ray differentials for
                    // rays reflected from triangles with degenerate
                    // parameterizations are still reasonable.
                    let dn =
                        Vector3f::cross(&(mesh.n[i2] - mesh.n[i0]), &(mesh.n[i1] - mesh.n[i0]));
                    if dn.length_squared() == 0.0 {
                        //
                    } else {
                        let (dnu, dnv) = coordinate_system(&dn);
                        dndu = dnu;
                        dndv = dnv;
                    }
                } else {
                    let inv_det = 1.0 / determinant;
                    dndu = (duv12[1] * dn1 - duv02[1] * dn2) * inv_det;
                    dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * inv_det;
                }
            }

            if mesh.reverse_orientation {
                ts *= -1.0;
            }

            isect.set_shading_geometry(&ss, &ts, &dndu, &dndv, true);
        }

        TESTS.with(|stat| {
            stat.add_num(1);
        });
        return Some((t, isect));
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        let _p = ProfilePhase::new(Prof::TriIntersectP);
        TESTS.with(|stat| {
            stat.add_denom(1);
        });

        let mesh = self.mesh.as_ref();
        let i0 = self.v[0] as usize;
        let i1 = self.v[1] as usize;
        let i2 = self.v[2] as usize;
        let p0 = mesh.p[i0];
        let p1 = mesh.p[i1];
        let p2 = mesh.p[i2];

        let dp02 = p0 - p2;
        let dp12 = p1 - p2;
        let mut n = Vector3f::cross(&dp02, &dp12);
        if mesh.reverse_orientation ^ mesh.swaps_handedness {
            n *= -1.0;
        }

        if !mesh.two_sided {
            if Vector3f::dot(&n, &r.d) >= 0.0 {
                return false;
            }
        }

        let mut p0t = p0 - r.o;
        let mut p1t = p1 - r.o;
        let mut p2t = p2 - r.o;

        let kz = max_dimension(&(r.d.abs())) as usize;
        let kx = TRI[kz + 1];
        let ky = TRI[kx + 1];

        let d = permute(&r.d, kx, ky, kz);
        p0t = permute(&p0t, kx, ky, kz);
        p1t = permute(&p1t, kx, ky, kz);
        p2t = permute(&p2t, kx, ky, kz);

        let sx = -d.x / d.z;
        let sy = -d.y / d.z;
        let sz = 1.0 / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;

        // Compute edge function coefficients _e0_, _e1_, and _e2_
        let mut e0 = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1 = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2 = p0t.x * p1t.y - p0t.y * p1t.x;

        if e0 == 0.0 || e1 == 0.0 || e2 == 0.0 {
            let p2txp1ty = p2t.x as f64 * p1t.y as f64;
            let p2typ1tx = p2t.y as f64 * p1t.x as f64;
            e0 = (p2typ1tx - p2txp1ty) as Float;
            let p0txp2ty = p0t.x as f64 * p2t.y as f64;
            let p0typ2tx = p0t.y as f64 * p2t.x as f64;
            e1 = (p0typ2tx - p0txp2ty) as Float;
            let p1txp0ty = p1t.x as f64 * p0t.y as f64;
            let p1typ0tx = p1t.y as f64 * p0t.x as f64;
            e2 = (p1typ0tx - p1txp0ty) as Float;
        }

        // Perform triangle edge and determinant tests
        if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0) && (e0 > 0.0 || e1 > 0.0 || e2 > 0.0) {
            return false;
        }
        let det = e0 + e1 + e2;
        if det == 0.0 {
            return false;
        }

        // Compute scaled hit distance to triangle and test against ray $t$ range
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;

        let t_max = r.t_max.get();

        let t_scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        #[allow(clippy::if_same_then_else)]
        if det < 0.0 && (t_scaled >= 0.0 || t_scaled < t_max * det) {
            return false;
        } else if det > 0.0 && (t_scaled <= 0.0 || t_scaled > t_max * det) {
            return false;
        }

        // Compute barycentric coordinates and $t$ value for triangle intersection
        let inv_det = 1.0 / det;
        //let b0 = e0 * inv_det;
        //let b1 = e1 * inv_det;
        //let b2 = e2 * inv_det;
        let t = t_scaled * inv_det;

        // Ensure that computed triangle $t$ is conservatively greater than zero

        // Compute $\delta_z$ term for triangle $t$ error bounds
        let max_zt = max_component(&Vector3f::new(p0t.z, p1t.z, p2t.z).abs());
        let delta_z = GAMMA3 * max_zt;

        // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        let max_xt = max_component(&Vector3f::new(p0t.x, p1t.x, p2t.x).abs());
        let max_yt = max_component(&Vector3f::new(p0t.y, p1t.y, p2t.y).abs());
        let delta_x = GAMMA5 * (max_xt + max_zt);
        let delta_y = GAMMA5 * (max_yt + max_zt);

        // Compute $\delta_e$ term for triangle $t$ error bounds
        let delta_e = 2.0 * (GAMMA2 * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);
        // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        let max_e = max_component(&Vector3f::new(e0, e1, e2).abs());
        let delta_t = 3.0
            * (GAMMA3 * max_e * max_zt + delta_e * max_zt + delta_z * max_e)
            * Float::abs(inv_det);
        if t <= delta_t {
            return false;
        }

        TESTS.with(|stat| {
            stat.add_num(1);
        });
        return true;
    }

    fn area(&self) -> Float {
        let mesh = self.mesh.as_ref();
        let i0 = self.v[0] as usize;
        let i1 = self.v[1] as usize;
        let i2 = self.v[2] as usize;
        let p0 = mesh.p[i0];
        let p1 = mesh.p[i1];
        let p2 = mesh.p[i2];
        return 0.5 * Vector3f::cross(&(p1 - p0), &(p2 - p0)).length();
    }

    fn sample(&self, u: &Point2f) -> Option<(Interaction, Float)> {
        let b = uniform_sample_triangle(u);
        // Get triangle vertices in _p0_, _p1_, and _p2_
        let mesh = self.mesh.as_ref();
        let i0 = self.v[0] as usize;
        let i1 = self.v[1] as usize;
        let i2 = self.v[2] as usize;
        let p0 = mesh.p[i0];
        let p1 = mesh.p[i1];
        let p2 = mesh.p[i2];
        let p = b[0] * p0 + b[1] * p1 + (1.0 - b[0] - b[1]) * p2;
        // Compute surface normal for sampled point on triangle
        let mut n = Vector3f::cross(&(p1 - p0), &(p2 - p0)).normalize();
        // Ensure correct orientation of the geometric normal; follow the same
        // approach as was used in Triangle::Intersect().
        if !mesh.n.is_empty() {
            let ns = b[0] * mesh.n[i0] + b[1] * mesh.n[i1] + (1.0 - b[0] - b[1]) * mesh.n[i2];
            n = face_forward(&n, &ns);
        } else if mesh.reverse_orientation ^ mesh.swaps_handedness {
            n *= -1.0;
        }
        // Compute error bounds for sampled point on triangle
        let p_abs_sum = Vector3f::abs(&(b[0] * p0))
            + Vector3f::abs(&(b[1] * p1))
            + Vector3f::abs(&((1.0 - b[0] - b[1]) * p2));
        let p_error = GAMMA6 * Vector3f::new(p_abs_sum.x, p_abs_sum.y, p_abs_sum.z);
        let pdf = 1.0 / self.area();
        let it = Interaction::from_surface_sample(&p, &p_error, &n);
        return Some((it, pdf));
    }

    fn sample_from(&self, ref_: &Interaction, u: &Point2f) -> Option<(Interaction, Float)> {
        let (intr, pdf) = self.sample(u)?;
        assert!(intr.is_surface_interaction());
        let wi = intr.get_p() - ref_.get_p();
        if wi.length_squared() <= 0.0 {
            return None;
        } else {
            assert!(intr.get_n().length() > 0.0);
            let wi = wi.normalize();

            //pbrt-r3
            {
                let mesh = self.mesh.as_ref();
                if !mesh.two_sided {
                    if Vector3::dot(&intr.get_n(), &-wi) <= 0.0 {
                        return None;
                    }
                }
            }
            //pbrt-r3

            // Convert from area measure, as returned by the Sample() call
            // above, to solid angle measure.
            let pdf = pdf * Vector3f::distance_squared(&ref_.get_p(), &intr.get_p())
                / Vector3f::abs_dot(&intr.get_n(), &-wi);
            if pdf <= 0.0 || pdf.is_infinite() {
                return None;
            }
            return Some((intr, pdf));
        }
    }
}

pub fn get_alpha_texture(
    params: &ParamSet,
    float_textures: &FloatTextureMap,
) -> Option<AlphaMaskInfo> {
    if let Some(textures) = params.get_textures_ref("alpha") {
        if textures.len() >= 1 {
            let alpha_tex_name = textures[0].clone();
            if let Some(tex) = float_textures.get(&alpha_tex_name) {
                return Some(AlphaMaskInfo::Texture {
                    texture: Arc::clone(tex),
                });
            }
        }
    } else if let Some(alpha) = params.get_floats_ref("alpha") {
        if alpha.len() > 0 {
            return Some(AlphaMaskInfo::Value { value: alpha[0] });
        }
    }
    return None;
}

pub fn get_shadow_alpha_texture(
    params: &ParamSet,
    float_textures: &FloatTextureMap,
) -> Option<AlphaMaskInfo> {
    if let Some(textures) = params.get_textures_ref("shadowalpha") {
        if textures.len() >= 1 {
            let alpha_tex_name = textures[0].clone();
            if let Some(tex) = float_textures.get(&alpha_tex_name) {
                return Some(AlphaMaskInfo::Texture {
                    texture: Arc::clone(tex),
                });
            }
        }
    } else if let Some(alpha) = params.get_floats_ref("shadowalpha") {
        if alpha.len() > 0 {
            return Some(AlphaMaskInfo::Value { value: alpha[0] });
        }
    }
    return None;
}

pub fn create_triangle_mesh(
    o2w: &Transform,
    _: &Transform,
    reverse_orientation: bool,
    vertex_indices: Vec<u32>,
    p: Vec<Point3f>,
    s: Vec<Vector3f>,
    n: Vec<Vector3f>,
    uv: Vec<Point2f>,
    params: &ParamSet,
) -> Vec<Arc<dyn Shape>> {
    let two_sided = params.find_one_bool("twosided", true);
    let mesh = Arc::new(TriangleMesh::new(
        o2w,
        reverse_orientation,
        two_sided,
        p,
        s,
        n,
        uv,
    ));
    let n_triangles = vertex_indices.len() / 3;
    let mut tris: Vec<Arc<dyn Shape>> = Vec::with_capacity(n_triangles);
    for i in 0..n_triangles {
        let v: [u32; 3] = [
            vertex_indices[3 * i + 0],
            vertex_indices[3 * i + 1],
            vertex_indices[3 * i + 2],
        ];
        let tri = Arc::new(Triangle::new(&mesh, &v, i));
        if tri.as_ref().area() > 1e-16 {
            tris.push(tri);
        }
    }
    return tris;
}

fn is_fillable_uv(vertex_indices: &[u32], vertices_length: usize) -> bool {
    let mut index_check = vec![0; vertices_length];
    let fsz = vertex_indices.len() / 3;
    for i in 0..fsz {
        let f = 3 * i;
        let vi = [
            vertex_indices[f + 0] as usize,
            vertex_indices[f + 1] as usize,
            vertex_indices[f + 2] as usize,
        ];
        for j in 0..3 {
            if index_check[vi[j]] == 0 || index_check[vi[j]] == j {
                index_check[vi[j]] = j;
            } else {
                return false;
            }
        }
    }
    return true;
}

type FloatTextureMap = HashMap<String, Arc<dyn Texture<Float>>>;

pub fn create_triangle_mesh_shape(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
    float_textures: &FloatTextureMap,
) -> Result<Vec<Arc<dyn Shape>>, PbrtError> {
    let mut vertex_indices = Vec::new();
    let mut p: Vec<Point3f> = Vec::new();
    let mut s: Vec<Vector3f> = Vec::new();
    let mut n: Vec<Normal3f> = Vec::new();
    let mut uv: Vec<Vector2f> = Vec::new();

    if let Some(vi) = params.get_ints_ref("indices") {
        vertex_indices.resize(vi.len(), 0);
        for i in 0..vi.len() {
            vertex_indices[i] = vi[i] as u32;
        }
    }

    if let Some(ps) = params.get_points_ref("P") {
        let sz = ps.len() / 3;
        p.resize(sz, Point3f::zero());
        for i in 0..sz {
            p[i] = Point3f::new(ps[3 * i + 0], ps[3 * i + 1], ps[3 * i + 2]);
        }
    }

    if let Some(ps) = params.get_floats_ref("uv") {
        let sz = ps.len() / 2;
        uv.resize(sz, Vector2::zero());
        for i in 0..sz {
            uv[i] = Vector2::new(ps[2 * i + 0], ps[2 * i + 1]);
        }
    } else if let Some(ps) = params.get_floats_ref("st") {
        let sz = ps.len() / 2;
        uv.resize(sz, Vector2::zero());
        for i in 0..sz {
            uv[i] = Vector2::new(ps[2 * i + 0], ps[2 * i + 1]);
        }
    } else if !vertex_indices.is_empty() {
        // pbrt-r3:
        if is_fillable_uv(&vertex_indices, p.len()) {
            let tri_uv = [
                Point2f::new(0.0, 0.0),
                Point2f::new(1.0, 0.0),
                Point2f::new(1.0, 1.0),
            ];
            uv.resize(p.len(), Vector2::zero());
            let fsz = vertex_indices.len() / 3;
            for i in 0..fsz {
                let f = 3 * i;
                let vi = [
                    vertex_indices[f + 0] as usize,
                    vertex_indices[f + 1] as usize,
                    vertex_indices[f + 2] as usize,
                ];
                for j in 0..3 {
                    if uv[vi[j]].x == 0.0 && uv[vi[j]].y == 0.0 {
                        uv[vi[j]] = tri_uv[j];
                    }
                }
            }
        }
        // pbrt-r3:
    }

    if let Some(ps) = params.get_points_ref("S") {
        let sz = ps.len() / 3;
        s.resize(sz, Vector3::zero());
        for i in 0..sz {
            s[i] = Vector3f::new(ps[3 * i + 0], ps[3 * i + 1], ps[3 * i + 2]);
        }
    }

    if let Some(ps) = params.get_points_ref("N") {
        let sz = ps.len() / 3;
        n.resize(sz, Normal3f::zero());
        for i in 0..sz {
            n[i] = Normal3f::new(ps[3 * i + 0], ps[3 * i + 1], ps[3 * i + 2]);
        }
    }

    if !vertex_indices.is_empty() && !p.is_empty() {
        let mut mesh = create_triangle_mesh(
            o2w,
            w2o,
            reverse_orientation,
            vertex_indices,
            p,
            s,
            n,
            uv,
            params,
        );

        let alpha_mask_info = get_alpha_texture(params, float_textures);
        let shadow_alpha_mask_info = get_shadow_alpha_texture(params, float_textures);
        if alpha_mask_info.is_some() || shadow_alpha_mask_info.is_some() {
            for i in 0..mesh.len() {
                mesh[i] = Arc::new(AlphaMaskShape::new(
                    &mesh[i],
                    &alpha_mask_info,
                    &shadow_alpha_mask_info,
                ));
            }
        }

        return Ok(mesh);
    } else {
        return Err(PbrtError::from("Invalid mesh"));
    }
}
