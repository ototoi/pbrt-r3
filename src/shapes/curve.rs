use crate::core::pbrt::constants::*;
use crate::core::pbrt::*;

use std::sync::Arc;

thread_local!(static TESTS: StatPercent = StatPercent::new("Intersections/Ray-curve intersection tests"));
thread_local!(static REFINEMENT_LEVEL: StatIntDistribution = StatIntDistribution::new("Intersections/Curve refinement level"));
thread_local!(static N_CURVES: StatCounter = StatCounter::new("Scene/Curves"));
thread_local!(static N_SPLIT_CURVES: StatCounter = StatCounter::new("Scene/Split curves"));

// CurveType Declarations
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CurveType {
    Flat,
    Cylinder,
    Ribbon,
}

#[derive(Debug, Clone, Copy)]
pub struct CurveNormal {
    pub n: [Normal3f; 2],
    pub normal_angle: Float,
    pub inv_sin_normal_angle: Float,
}

// CurveCommon Declarations
#[derive(Debug, Clone)]
pub struct CurveCommon {
    pub base: BaseShape,
    pub t: CurveType,
    pub cp_obj: [Point3f; 4],
    pub width: [Float; 2],
    pub normals: Option<CurveNormal>,
}

// Curve Declarations
#[derive(Debug, Clone)]
pub struct Curve {
    common: Arc<CurveCommon>,
    u_min: Float,
    u_max: Float,
}

// Curve Utility Functions
#[inline]
fn lerp3(t: Float, p0: &Point3f, p1: &Point3f) -> Point3f {
    return Point3f::new(
        lerp(t, p0.x, p1.x),
        lerp(t, p0.y, p1.y),
        lerp(t, p0.z, p1.z),
    );
}

#[inline]
fn blossom_bezier(cp: &[Point3f], u0: Float, u1: Float, u2: Float) -> Point3f {
    let a = [
        lerp3(u0, &cp[0], &cp[1]),
        lerp3(u0, &cp[1], &cp[2]),
        lerp3(u0, &cp[2], &cp[3]),
    ];
    let b = [lerp3(u1, &a[0], &a[1]), lerp3(u1, &a[1], &a[2])];
    return lerp3(u2, &b[0], &b[1]);
}

#[inline]
fn subdivide_bezier(cp: &[Point3f]) -> [Point3f; 7] {
    let cp_split = [
        cp[0],
        (cp[0] + cp[1]) * 0.5,
        (cp[0] + 2.0 * cp[1] + cp[2]) * 0.25,
        (cp[0] + 3.0 * cp[1] + 3.0 * cp[2] + cp[3]) * 0.125,
        (cp[1] + 2.0 * cp[2] + cp[3]) * 0.25,
        (cp[2] + cp[3]) * 0.5,
        cp[3],
    ];
    return cp_split;
}

#[inline]
fn eval_bezier(cp: &[Point3f], u: Float) -> (Point3f, Vector3f) {
    let cp1 = [
        lerp3(u, &cp[0], &cp[1]),
        lerp3(u, &cp[1], &cp[2]),
        lerp3(u, &cp[2], &cp[3]),
    ];
    let cp2 = [lerp3(u, &cp1[0], &cp1[1]), lerp3(u, &cp1[1], &cp1[2])];
    let cp3 = lerp3(u, &cp2[0], &cp2[1]);
    let delta = cp2[1] - cp2[0];
    let deriv = if delta.length_squared() > 0.0 {
        3.0 * delta
    } else {
        // For a cubic Bezier, if the first three control points (say) are
        // coincident, then the derivative of the curve is legitimately (0,0,0)
        // at u=0.  This is problematic for us, though, since we'd like to be
        // able to compute a surface normal there.  In that case, just punt and
        // take the difference between the first and last control points, which
        // ain't great, but will hopefully do.
        cp[3] - cp[0]
    };
    return (cp3, deriv);
}

#[inline]
fn transform_lookat(e: &Point3f, l: &Point3f, u: &Vector3f) -> Transform {
    return Transform::look_at(e.x, e.y, e.z, l.x, l.y, l.z, u.x, u.y, u.z);
}

impl CurveCommon {
    pub fn new(
        o2w: &Transform,
        w2o: &Transform,
        reverse_orientation: bool,
        c: &[Point3f],
        w0: Float,
        w1: Float,
        t: CurveType,
        norm: &Option<Vec<Normal3f>>,
    ) -> Self {
        let cp_obj: Vec<Point3f> = c.to_vec();
        let cp_obj: [Point3f; 4] = cp_obj.try_into().unwrap();
        let width = [w0, w1];
        let normals: Option<CurveNormal> = match norm {
            Some(n) => {
                let n0 = n[0].normalize();
                let n1 = n[1].normalize();
                let normal_angle = Float::acos(Float::clamp(Normal3f::dot(&n0, &n1), 0.0, 1.0));
                let inv_sin_normal_angle = 1.0 / Float::sin(normal_angle);
                Some(CurveNormal {
                    n: [n0, n1],
                    normal_angle,
                    inv_sin_normal_angle,
                })
            }
            None => None,
        };
        N_CURVES.with(|n| n.inc());
        CurveCommon {
            base: BaseShape::new(o2w, w2o, reverse_orientation),
            t,
            cp_obj,
            width,
            normals,
        }
    }
}

fn get_bounds(cp: &[Point3f]) -> Bounds3f {
    let min = cp.iter().fold(cp[0], |a, b| -> Vector3f {
        return Vector3f::new(
            Float::min(a[0], b[0]),
            Float::min(a[1], b[1]),
            Float::min(a[2], b[2]),
        );
    });
    let max = cp.iter().fold(cp[0], |a, b| -> Vector3f {
        return Vector3f::new(
            Float::max(a[0], b[0]),
            Float::max(a[1], b[1]),
            Float::max(a[2], b[2]),
        );
    });
    return Bounds3f::from(((min[0], min[1], min[2]), (max[0], max[1], max[2])));
}

fn recursive_intersect(
    common: &CurveCommon,
    is_shadow: bool,
    ray: &Ray,
    cp: &[Point3f],
    ray_to_object: &Transform,
    u0: Float,
    u1: Float,
    t_max: Float,
    depth: i32,
) -> Option<(Float, Option<SurfaceInteraction>)> {
    let mut t_max = t_max;
    let width = &common.width;
    let t = common.t;
    let normals = &common.normals;
    let ray_length = ray.d.length();

    if depth > 0 {
        // Split curve segment into sub-segments and test for intersection
        let cp_split = subdivide_bezier(cp);
        // For each of the two segments, see if the ray's bounding box
        // overlaps the segment before recursively checking for
        // intersection with it.
        let u = [u0, (u0 + u1) / 2.0, u1];
        let mut t_hit = None;
        // Pointer to the 4 control poitns for the current segment.
        for seg in 0..2 {
            let z_max = ray_length * t_max;
            //0:0,1,2,3
            //1:3,4,5,6
            let cps = &cp_split[(3 * seg)..];
            let max_width = Float::max(
                lerp(u[seg], width[0], width[1]),
                lerp(u[seg + 1], width[0], width[1]),
            );
            let max_radius = 0.5 * max_width;

            // As above, check y first, since it most commonly lets us exit
            // out early.
            let max_y = Float::max(
                Float::max(cps[0].y, cps[1].y),
                Float::max(cps[2].y, cps[3].y),
            );
            let min_y = Float::min(
                Float::min(cps[0].y, cps[1].y),
                Float::min(cps[2].y, cps[3].y),
            );
            if ((max_y + max_radius) < 0.0) || ((min_y - max_radius) > 0.0) {
                continue;
            }

            let max_x = Float::max(
                Float::max(cps[0].x, cps[1].x),
                Float::max(cps[2].x, cps[3].x),
            );
            let min_x = Float::min(
                Float::min(cps[0].x, cps[1].x),
                Float::min(cps[2].x, cps[3].x),
            );
            if ((max_x + max_radius) < 0.0) || ((min_x - max_radius) > 0.0) {
                continue;
            }

            let max_z = Float::max(
                Float::max(cps[0].z, cps[1].z),
                Float::max(cps[2].z, cps[3].z),
            );
            let min_z = Float::min(
                Float::min(cps[0].z, cps[1].z),
                Float::min(cps[2].z, cps[3].z),
            );
            if ((max_z + max_radius) < 0.0) || ((min_z - max_radius) > z_max) {
                continue;
            }

            {
                let t_t_hit = recursive_intersect(
                    common,
                    is_shadow,
                    ray,
                    cps,
                    ray_to_object,
                    u[seg],
                    u[seg + 1],
                    t_max,
                    depth - 1,
                );
                if let Some(t_t_hit) = t_t_hit {
                    // If we found an intersection and this is a shadow ray,
                    // we can exit out immediately.
                    if is_shadow {
                        return Some(t_t_hit);
                    }
                    //if t_t_hit.0 < t_max {
                    t_max = Float::min(t_max, t_t_hit.0);
                    t_hit = Some(t_t_hit);
                    //}
                }
            }
        }
        return t_hit;
    } else {
        // Intersect ray with curve segment

        // Test ray against segment endpoint boundaries

        // Test sample point against tangent perpendicular at curve start
        let edge = (cp[1].y - cp[0].y) * -cp[0].y + cp[0].x * (cp[0].x - cp[1].x);
        if edge < 0.0 {
            return None;
        }

        // Test sample point against tangent perpendicular at curve end
        let edge = (cp[2].y - cp[3].y) * -cp[3].y + cp[3].x * (cp[3].x - cp[2].x);
        if edge < 0.0 {
            return None;
        }

        // Compute line $w$ that gives minimum distance to sample point
        let segment_direction = Point2f::new(cp[3][0], cp[3][1]) - Point2f::new(cp[0][0], cp[0][1]);
        let denom = segment_direction.length_squared();
        if denom == 0.0 {
            return None;
        }
        let w = Vector2f::dot(&-Vector2f::new(cp[0][0], cp[0][1]), &segment_direction) / denom;
        // Compute $u$ coordinate of curve intersection point and _hitWidth_
        let u = Float::clamp(lerp(w, u0, u1), u0, u1);

        let mut hit_width = lerp(u, width[0], width[1]);
        let mut n_hit = Normal3f::new(1.0, 0.0, 0.0);
        if t == CurveType::Ribbon {
            if let Some(n) = normals {
                // Scale _hitWidth_ based on ribbon orientation
                let sin0 = Float::sin((1.0 - u) * n.normal_angle) * n.inv_sin_normal_angle;
                let sin1 = Float::sin(u * n.normal_angle) * n.inv_sin_normal_angle;
                n_hit = (sin0 * n.n[0] + sin1 * n.n[1]).normalize();
                hit_width *= Vector3f::abs_dot(&n_hit, &ray.d) / ray_length;
            }
        }

        // Test intersection point against curve width
        let z_max = ray_length * t_max;
        let (pc, dpcdw) = eval_bezier(cp, Float::clamp(w, 0.0, 1.0));
        let pt_curve_dist2 = pc.x * pc.x + pc.y * pc.y;
        if pt_curve_dist2 > (hit_width * hit_width * 0.25) {
            return None;
        }
        if pc.z < 0.0 || pc.z > z_max {
            return None;
        }

        // Compute $v$ coordinate of curve intersection point
        let pt_curve_dist = Float::sqrt(pt_curve_dist2);
        let edge_func = dpcdw.x * -pc.y + pc.x * dpcdw.y;
        let v = if edge_func > 0.0 {
            0.5 + pt_curve_dist / hit_width
        } else {
            0.5 - pt_curve_dist / hit_width
        };

        // Compute hit _t_ and partial derivatives for curve intersection

        // FIXME: this tHit isn't quite right for ribbons...
        let t_hit = pc.z / ray_length;
        //println!("{}", t_hit);

        TESTS.with(|stat| stat.add_num(1));

        if !is_shadow {
            // Compute error bounds for curve intersection
            let p_error_coef = 2.0;
            let p_error = Vector3f::new(
                p_error_coef * hit_width,
                p_error_coef * hit_width,
                p_error_coef * hit_width,
            );
            //let p_error = ray_to_object.transform_vector(&p_error);//object space;

            // Compute $\dpdu$ and $\dpdv$ for curve intersection
            let (_, dpdu) = eval_bezier(&common.cp_obj, u); //object space

            let dpdv; // = Vector3f::zero();
            if t == CurveType::Ribbon {
                dpdv = Vector3f::normalize(&Vector3f::cross(&n_hit, &dpdu)) * hit_width;
            } else {
                let object_to_ray = ray_to_object.inverse();
                // Compute curve $\dpdv$ for flat and cylinder curves
                let dpdu_plane = object_to_ray.transform_vector(&dpdu).normalize(); //ray space
                let p_plane = Vector3f::new(0.0, 0.0, 1.0);
                let mut dpdv_plane = Vector3f::cross(&p_plane, &dpdu_plane);
                //let mut dpdv_plane =
                //    Vector3f::normalize(&Vector3f::new(-dpdu_plane.y, dpdu_plane.x, 0.0))
                //        * hit_width;
                if t == CurveType::Cylinder {
                    // Rotate _dpdvPlane_ to give cylindrical appearance
                    let theta = lerp(v, -90.0, 90.0);
                    let rot = Transform::rotate(-theta, dpdu_plane.x, dpdu_plane.y, dpdu_plane.z);
                    dpdv_plane = rot.transform_vector(&dpdv_plane);
                }
                dpdv = ray_to_object.transform_vector(&dpdv_plane).normalize(); //object space
            }
            let p = ray.o + t_hit * ray.d;
            let uv = Point2f::new(u, v); //u, v
            let wo = -ray.d;
            let mut n = common.base.calc_normal(&dpdu, &dpdv);
            //if t == CurveType::Cylinder {
            //    n = (p - p_center).normalize();
            //} else {
            if Vector3f::dot(&ray.d, &n) > 0.0 {
                n *= -1.0;
            }
            //}
            let si = SurfaceInteraction::new(
                &p,
                &p_error,
                &uv,
                &wo,
                &n,
                &dpdu,
                &dpdv,
                &Vector3f::zero(),
                &Vector3f::zero(),
                ray.time,
                0,
            );
            let si = common
                .base
                .object_to_world
                .transform_surface_interaction(&si);
            return Some((t_hit, Some(si)));
        } else {
            return Some((t_hit, None));
        }
    }
}

fn log2(x: Float) -> i32 {
    if x < 1.0 {
        return 0;
    } else {
        return ((Float::log2(x) + 0.5) * 10.0) as i32;
    }
}

impl Curve {
    pub fn new(common: &Arc<CurveCommon>, u_min: Float, u_max: Float) -> Self {
        Curve {
            common: common.clone(),
            u_min,
            u_max,
        }
    }

    fn intersect_common(
        &self,
        r: &Ray,
        is_shadow: bool,
    ) -> Option<(Float, Option<SurfaceInteraction>)> {
        TESTS.with(|stat| stat.add_denom(1));

        let u_min = self.u_min;
        let u_max = self.u_max;
        let common = self.common.as_ref();
        // Transform _Ray_ to object space
        let (ray, _o_err, _d_err) = common.base.world_to_object.transform_ray(r);
        //let ray_length = ray.d.length();
        let t_max = ray.t_max.get();
        //println!("rl:{}, tmax:{}", ray_length, t_max);

        // Compute object-space control points for curve segment, _cpObj_
        let cp_obj = [
            blossom_bezier(&common.cp_obj, u_min, u_min, u_min),
            blossom_bezier(&common.cp_obj, u_min, u_min, u_max),
            blossom_bezier(&common.cp_obj, u_min, u_max, u_max),
            blossom_bezier(&common.cp_obj, u_max, u_max, u_max),
        ];

        // Project curve control points to plane perpendicular to ray

        // Be careful to set the "up" direction passed to LookAt() to equal the
        // vector from the first to the last control points.  In turn, this
        // helps orient the curve to be roughly parallel to the x axis in the
        // ray coordinate system.
        //
        // In turn (especially for curves that are approaching stright lines),
        // we get curve bounds with minimal extent in y, which in turn lets us
        // early out more quickly in recursiveIntersect().
        let delta = cp_obj[3] - cp_obj[0];
        let mut dx = Vector3f::cross(&ray.d, &delta);
        if dx.length_squared() == 0.0 {
            // If the ray and the vector between the first and last control
            // points are parallel, dx will be zero.  Generate an arbitrary xy
            // orientation for the ray coordinate system so that intersection
            // tests can proceeed in this unusual case.
            let (_dx, _) = coordinate_system(&ray.d);
            dx = _dx;
        }
        let object_to_ray = transform_lookat(&ray.o, &(ray.o + ray.d), &dx);
        let cp = [
            object_to_ray.transform_point(&cp_obj[0]),
            object_to_ray.transform_point(&cp_obj[1]),
            object_to_ray.transform_point(&cp_obj[2]),
            object_to_ray.transform_point(&cp_obj[3]),
        ];

        // Before going any further, see if the ray's bounding box intersects
        // the curve's bounding box. We start with the y dimension, since the y
        // extent is generally the smallest (and is often tiny) due to our
        // careful orientation of the ray coordinate ysstem above.
        let max_width = Float::max(
            lerp(u_min, common.width[0], common.width[1]),
            lerp(u_max, common.width[0], common.width[1]),
        );
        let max_radius = 0.5 * max_width;

        // Check for non-overlap in y.
        {
            let max_y = Float::max(Float::max(cp[0].y, cp[1].y), Float::max(cp[2].y, cp[3].y));
            let min_y = Float::min(Float::min(cp[0].y, cp[1].y), Float::min(cp[2].y, cp[3].y));
            if ((max_y + max_radius) < 0.0) || ((min_y - max_radius) > 0.0) {
                return None;
            }
        }

        // Check for non-overlap in x.
        {
            let max_x = Float::max(Float::max(cp[0].x, cp[1].x), Float::max(cp[2].x, cp[3].x));
            let min_x = Float::min(Float::min(cp[0].x, cp[1].x), Float::min(cp[2].x, cp[3].x));
            if ((max_x + max_radius) < 0.0) || ((min_x - max_radius) > 0.0) {
                return None;
            }
        }

        // Check for non-overlap in z.
        {
            let ray_length = ray.d.length();
            let z_max = ray_length * t_max;

            let max_z = Float::max(Float::max(cp[0].z, cp[1].z), Float::max(cp[2].z, cp[3].z));
            let min_z = Float::min(Float::min(cp[0].z, cp[1].z), Float::min(cp[2].z, cp[3].z));
            if ((max_z + max_radius) < 0.0) || ((min_z - max_radius) > z_max) {
                return None;
            }
        }

        // Compute refinement depth for curve, _maxDepth_
        let mut l0 = 0.0;
        for i in 0..2 {
            l0 = Float::max(
                l0,
                Float::max(
                    Float::max(
                        Float::abs(cp[i].x - 2.0 * cp[i + 1].x + cp[i + 2].x),
                        Float::abs(cp[i].y - 2.0 * cp[i + 1].y + cp[i + 2].y),
                    ),
                    Float::abs(cp[i].z - 2.0 * cp[i + 1].z + cp[i + 2].z),
                ),
            );
        }
        let eps = Float::max(common.width[0], common.width[1]) * 0.05; // width / 20

        // Compute log base 4 by dividing log2 in half.
        //let r0 = log2(1.41421356237 * 6.0 * l0 / (8.0 * eps)) / 2;
        let r0 = log2(SQRT_2 * 6.0 * l0 / (8.0 * eps)) / 2;
        //let r0 = log2(1.41421356237 * 6.0 * l0 / (8.0 * eps) * 16.0) / 2;//x16

        let max_depth = i32::clamp(r0, 0, 10); //

        REFINEMENT_LEVEL.with(|stat| stat.add(max_depth as u64));

        return recursive_intersect(
            &common,
            is_shadow,
            &ray,
            &cp,
            &Transform::inverse(&object_to_ray),
            u_min,
            u_max,
            t_max,
            max_depth,
        );
    }
}

impl Shape for Curve {
    fn object_bound(&self) -> Bounds3f {
        let u_min = self.u_min;
        let u_max = self.u_max;
        let common = self.common.as_ref();
        let cp_obj = [
            blossom_bezier(&common.cp_obj, u_min, u_min, u_min),
            blossom_bezier(&common.cp_obj, u_min, u_min, u_max),
            blossom_bezier(&common.cp_obj, u_min, u_max, u_max),
            blossom_bezier(&common.cp_obj, u_max, u_max, u_max),
        ];
        let b = get_bounds(&cp_obj);
        let width = [
            lerp(u_min, common.width[0], common.width[1]),
            lerp(u_max, common.width[0], common.width[1]),
        ];
        let radius = Float::max(width[0], width[1]) * 0.5;
        let b = b.expand(radius);
        return b;
    }

    fn world_bound(&self) -> Bounds3f {
        let common = self.common.as_ref();
        let b = self.object_bound();
        let b = common.base.object_to_world.transform_bounds(&b);
        return b;
    }

    fn intersect(&self, r: &Ray) -> Option<(Float, SurfaceInteraction)> {
        if let Some((t_hit, isect)) = self.intersect_common(r, false) {
            return Some((t_hit, isect.unwrap()));
        } else {
            return None;
        }
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        return self.intersect_common(r, true).is_some();
    }

    fn area(&self) -> Float {
        // Compute object-space control points for curve segment, _cpObj_
        let u_min = self.u_min;
        let u_max = self.u_max;
        let common = self.common.as_ref();
        let cp_obj = [
            blossom_bezier(&common.cp_obj, u_min, u_min, u_min),
            blossom_bezier(&common.cp_obj, u_min, u_min, u_max),
            blossom_bezier(&common.cp_obj, u_min, u_max, u_max),
            blossom_bezier(&common.cp_obj, u_max, u_max, u_max),
        ];
        let width0 = lerp(u_min, common.width[0], common.width[1]);
        let width1 = lerp(u_max, common.width[0], common.width[1]);
        let avg_width = (width0 + width1) * 0.5;
        let mut approx_length = 0.0;
        for i in 0..3 {
            approx_length += Vector3f::distance(&cp_obj[i], &cp_obj[i + 1]);
        }
        return approx_length * avg_width;
    }

    fn sample(&self, _u: &Point2f) -> Option<(Interaction, Float)> {
        //LOG(FATAL) << "Curve::Sample not implemented.";
        return None;
    }
}

fn create_curve(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    c: &[Point3f],
    w0: Float,
    w1: Float,
    t: CurveType,
    norm: &Option<Vec<Normal3f>>,
    split_depth: i32,
) -> Result<Vec<Arc<dyn Shape>>, PbrtError> {
    let mut segments: Vec<Arc<dyn Shape>> = Vec::new();
    let common = Arc::new(CurveCommon::new(
        o2w,
        w2o,
        reverse_orientation,
        c,
        w0,
        w1,
        t,
        norm,
    ));
    let n_segments = (1 << split_depth) as usize;
    segments.reserve(n_segments);
    for i in 0..n_segments {
        let u_min = (i as Float) / (n_segments as Float);
        let u_max = ((i + 1) as Float) / (n_segments as Float);
        let curve = Arc::new(Curve::new(&common, u_min, u_max));
        segments.push(curve);
    }
    N_SPLIT_CURVES.with(|c| c.add(n_segments as u64));
    return Ok(segments);
}

fn get_curve_type(s: &str) -> Result<CurveType, PbrtError> {
    return match s {
        "flat" => Ok(CurveType::Flat),
        "ribbon" => Ok(CurveType::Ribbon),
        "cylinder" => Ok(CurveType::Cylinder),
        _ => {
            let msg = format!("Unknown curve type \"{}\".  Using \"cylinder\".", s);
            return Err(PbrtError::error(&msg));
        }
    };
}

fn convert_to_point(f: &[Float]) -> Vec<Point3f> {
    let l = f.len() / 3;
    let mut v = Vec::new();
    for i in 0..l {
        v.push(Vector3f::new(f[3 * i + 0], f[3 * i + 1], f[3 * i + 2]));
    }
    return v;
}

pub fn create_curve_shape(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
) -> Result<Vec<Arc<dyn Shape>>, PbrtError> {
    let width = params.find_one_float("width", 1.0);
    let width0 = params.find_one_float("width0", width);
    let width1 = params.find_one_float("width1", width);

    let degree = params.find_one_int("degree", 3);
    if degree != 2 && degree != 3 {
        let msg = format!(
            "Invalid degree {}: only degree 2 and 3 curves are supported.",
            degree
        );
        return Err(PbrtError::error(&msg));
    }
    let degree = degree as usize;
    let basis = params.find_one_string("basis", "bezier");
    if basis != "bezier" && basis != "bspline" {
        let msg = format!(
            "Invalid basis \"{}\": only \"bezier\" and \"bspline\" are supported.",
            degree
        );
        return Err(PbrtError::error(&msg));
    }
    let cp = params.get_points_ref("P").ok_or("\"P\"")?;
    let cp = convert_to_point(&cp);
    let ncp = cp.len();
    let n_segments: usize = match &basis as &str {
        "bezier" => {
            // After the first segment, which uses degree+1 control points,
            // subsequent segments reuse the last control point of the previous
            // one and then use degree more control points.
            if ((ncp - 1 - degree) % degree) != 0 {
                let msg = format!(
                    "Invalid number of control points {}: for the degree {} Bezier basis {} + n * {} are required, for n >= 0.",
                    ncp, degree, degree +1, degree
                );
                return Err(PbrtError::error(&msg));
            }
            (ncp - 1) / degree
        }
        "bspline" => ncp - degree,
        _ => {
            panic!();
        }
    };

    let curve_type = params.find_one_string("type", "flat");
    let curve_type = get_curve_type(&curve_type)?;

    let n = params.get_points_ref("N");
    /*
    let mut n = if let Some(nn) = n.as_ref() {
        Some(convert_to_point(nn))
    } else {
        None
    };
    */
    let mut n = n.as_ref().map(|nn| convert_to_point(nn));
    if let Some(nn) = n.as_ref() {
        let nnorm = nn.len();
        if curve_type != CurveType::Ribbon {
            //Warning("Curve normals are only used with \"ribbon\" type curves.");
            n = None;
        } else if nnorm != (n_segments + 1) {
            let msg = format!("Invalid number of normals {}: must provide {} normals for ribbon curves with {} segments.", nnorm, n_segments + 1, n_segments);
            return Err(PbrtError::error(&msg));
        }
    } else {
        if curve_type == CurveType::Ribbon {
            return Err(PbrtError::error(
                "Must provide normals \"N\" at curve endpoints with ribbon curves.",
            ));
        }
    }
    let sd = params.find_one_float("splitdepth", 3.0) as i32;
    let sd = params.find_one_int("splitdepth", sd);

    let mut curves: Vec<Arc<dyn Shape>> = Vec::new();
    for seg in 0..n_segments {
        // First, compute the cubic Bezier control points for the current
        // segment and store them in segCpBezier. (It is admittedly
        // wasteful storage-wise to turn b-splines into Bezier segments and
        // wasteful computationally to turn quadratic curves into cubics,
        // but yolo.)
        let cp_base = &cp[(seg * degree)..];
        let mut seg_cp_bezier = vec![
            Point3f::zero(),
            Point3f::zero(),
            Point3f::zero(),
            Point3f::zero(),
        ];
        if basis == "bezier" {
            if degree == 2 {
                // Elevate to degree 3.
                seg_cp_bezier[0] = cp_base[0];
                seg_cp_bezier[1] = lerp3(2.0 / 3.0, &cp_base[0], &cp_base[1]);
                seg_cp_bezier[2] = lerp3(1.0 / 3.0, &cp_base[1], &cp_base[2]);
                seg_cp_bezier[3] = cp_base[2];
            } else {
                //for i in 0..4 {
                //    seg_cp_bezier[i] = cp_base[i];
                //}
                seg_cp_bezier.copy_from_slice(&cp_base[..4]);
            }
        } else {
            // Uniform b-spline.
            if degree == 2 {
                // First compute equivalent Bezier control points via some
                // blossiming.  We have three control points and a uniform
                // knot vector; we'll label the points p01, p12, and p23.
                // We want the Bezier control points of the equivalent
                // curve, which are p11, p12, and p22.
                let p01 = cp_base[0];
                let p12 = cp_base[1];
                let p23 = cp_base[2];

                // We already have p12.
                let p11 = lerp3(0.5, &p01, &p12);
                let p22 = lerp3(0.5, &p12, &p23);

                // Now elevate to degree 3.
                seg_cp_bezier[0] = p11;
                seg_cp_bezier[1] = lerp3(2.0 / 3.0, &p11, &p12);
                seg_cp_bezier[2] = lerp3(1.0 / 3.0, &p12, &p22);
                seg_cp_bezier[3] = p22;
            } else {
                // Otherwise we will blossom from p012, p123, p234, and p345
                // to the Bezier control points p222, p223, p233, and p333.
                // https://people.eecs.berkeley.edu/~sequin/CS284/IMGS/cubicbsplinepoints.gif
                let p012 = cp_base[0];
                let p123 = cp_base[1];
                let p234 = cp_base[2];
                let p345 = cp_base[3];

                let p122 = lerp3(2.0 / 3.0, &p012, &p123);
                let p223 = lerp3(1.0 / 3.0, &p123, &p234);
                let p233 = lerp3(2.0 / 3.0, &p123, &p234);
                let p334 = lerp3(1.0 / 3.0, &p234, &p345);

                let p222 = lerp3(0.5, &p122, &p223);
                let p333 = lerp3(0.5, &p233, &p334);

                seg_cp_bezier[0] = p222;
                seg_cp_bezier[1] = p223;
                seg_cp_bezier[2] = p233;
                seg_cp_bezier[3] = p333;
            }
        }
        let w0 = lerp((seg as Float) / (n_segments as Float), width0, width1);
        let w1 = lerp(((seg + 1) as Float) / (n_segments as Float), width0, width1);
        let mut c = create_curve(
            o2w,
            w2o,
            reverse_orientation,
            &seg_cp_bezier,
            w0,
            w1,
            curve_type,
            &n,
            sd,
        )?;
        curves.append(&mut c);
    }

    return Ok(curves);
}
