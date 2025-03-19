use crate::core::efloat::EFloat;
use crate::core::error::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::shape::*;

#[inline]
fn radians(x: Float) -> Float {
    return x * (PI / 180.0);
}

pub struct Hyperboloid {
    pub base: BaseShape,
    pub p1: Point3f,
    pub p2: Point3f,
    pub z_min: Float,
    pub z_max: Float,
    pub r_max: Float,
    pub phi_max: Float,
    pub ah: Float,
    pub ch: Float,
}

fn compute_implicit_function_coefficients(point1: &Point3f, point2: &Point3f) -> (Float, Float) {
    let mut p1 = *point1;
    let mut p2 = *point2;
    if p2.z == 0.0 {
        std::mem::swap(&mut p1, &mut p2);
    }
    let mut pp = p1;
    loop {
        pp += 2.0 * (p2 - p1);
        let xy1 = pp.x * pp.x + pp.y * pp.y;
        let xy2 = p2.x * p2.x + p2.y * p2.y;
        let ah = (1.0 / xy1 - (pp.z * pp.z) / (xy1 * p2.z * p2.z))
            / (1.0 - (xy2 * pp.z * pp.z) / (xy1 * p2.z * p2.z));
        let ch = (ah * xy2 - 1.0) / (p2.z * p2.z);
        if Float::is_finite(ah) {
            return (ah, ch);
        }
    }
}

#[inline]
fn sqr(x: Float) -> Float {
    x * x
}

#[inline]
fn quad(x: Float) -> Float {
    x * x * x * x
}

impl Hyperboloid {
    pub fn new(
        o2w: &Transform,
        w2o: &Transform,
        reverse_orientation: bool,
        point1: &Point3f,
        point2: &Point3f,
        phi_max: Float,
    ) -> Self {
        let p1 = *point1;
        let p2 = *point2;
        let phi_max = radians(Float::clamp(phi_max, 0.0, 360.0));
        let radius1 = Float::sqrt(p1.x * p1.x + p1.y * p1.y);
        let radius2 = Float::sqrt(p2.x * p2.x + p2.y * p2.y);
        let r_max = Float::max(radius1, radius2);
        let z_min = Float::min(p1.z, p2.z);
        let z_max = Float::max(p1.z, p2.z);

        // Compute implicit function coefficients for hyperboloid
        let (ah, ch) = compute_implicit_function_coefficients(&p1, &p2);

        Hyperboloid {
            base: BaseShape::new(o2w, w2o, reverse_orientation),
            p1,
            p2,
            z_min,
            z_max,
            r_max,
            phi_max,
            ah,
            ch,
        }
    }
}

impl Shape for Hyperboloid {
    fn object_bound(&self) -> Bounds3f {
        let r_max = self.r_max;
        return Bounds3f::new(
            &Point3f::new(-r_max, -r_max, -r_max),
            &Point3f::new(r_max, r_max, r_max),
        );
    }

    fn world_bound(&self) -> Bounds3f {
        return self
            .base
            .object_to_world
            .transform_bounds(&self.object_bound());
    }
    fn intersect(&self, r: &Ray) -> Option<(Float, SurfaceInteraction)> {
        // Transform _Ray_ to object space
        let (ray, o_err, d_err) = self.base.world_to_object.transform_ray(r);

        // Compute quadratic hyperboloid coefficients

        // Initialize _EFloat_ ray coordinate values
        let p1 = self.p1;
        let p2 = self.p2;
        let ah = EFloat::from(self.ah);
        let ch = EFloat::from(self.ch);
        let z_min = self.z_min;
        let z_max = self.z_max;
        let phi_max = self.phi_max;

        let ox = EFloat::from((ray.o.x, o_err.x));
        let oy = EFloat::from((ray.o.y, o_err.y));
        let oz = EFloat::from((ray.o.z, o_err.z));
        let dx = EFloat::from((ray.d.x, d_err.x));
        let dy = EFloat::from((ray.d.y, d_err.y));
        let dz = EFloat::from((ray.d.z, d_err.z));
        let a = ah * (dx * dx) + ah * (dy * dy) - ch * (dz * dz);
        let b = 2.0 * (ah * (dx * ox) + ah * (dy * oy) - ch * (dz * oz));
        let c = ah * (ox * ox) + ah * (oy * oy) - ch * (oz * oz) - EFloat::from(1.0);

        // Solve quadratic equation for _t_ values
        let t_max = ray.t_max.get();

        let (t0, t1) = EFloat::quadratic(a, b, c)?;
        // pbrt-r3:
        if t0.v.is_infinite() || t1.v.is_infinite() {
            return None;
        }

        assert!(t0.v.is_finite());
        assert!(t1.v.is_finite());
        assert!(t0.v <= t1.v);
        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > t_max || t1.lower_bound() <= 0.0 {
            return None;
        }

        // Compute hyperboloid inverse mapping
        let mut t_shape_hit = t0;
        if t_shape_hit.lower_bound() <= 0.0 {
            t_shape_hit = t1;
            if t_max < t_shape_hit.upper_bound() {
                return None;
            }
        }

        // Compute hyperboloid inverse mapping
        let mut p_hit = ray.o + ray.d * Float::from(t_shape_hit);
        let mut v = (p_hit.z - p1.z) / (p2.z - p1.z);
        let pr = (1.0 - v) * p1 + v * p2;
        let mut phi = Float::atan2(
            pr.x * p_hit.y - p_hit.x * pr.y,
            p_hit.x * pr.x + p_hit.y * pr.y,
        );
        if phi < 0.0 {
            phi += 2.0 * PI;
        }

        // Test hyperboloid intersection against clipping parameters
        if p_hit.z < z_min || p_hit.z > z_max || phi > phi_max {
            if t_shape_hit == t1 {
                return None;
            }
            if t1.upper_bound() > t_max {
                return None;
            }
            t_shape_hit = t1;
            // Compute hyperboloid inverse mapping
            p_hit = ray.o + ray.d * Float::from(t_shape_hit);
            v = (p_hit.z - p1.z) / (p2.z - p1.z);
            let pr = (1.0 - v) * p1 + v * p2;
            phi = Float::atan2(
                pr.x * p_hit.y - p_hit.x * pr.y,
                p_hit.x * pr.x + p_hit.y * pr.y,
            );
            if phi < 0.0 {
                phi += 2.0 * PI;
            }
            if p_hit.z < z_min || p_hit.z > z_max || phi > phi_max {
                return None;
            }
        }

        // Find parametric representation of hyperboloid hit
        let u = phi / phi_max;

        // Compute hyperboloid $\dpdu$ and $\dpdv$
        let cos_phi = Float::cos(phi);
        let sin_phi = Float::sin(phi);
        let dpdu = Vector3f::new(-phi_max * p_hit.y, phi_max * p_hit.x, 0.0);
        let dpdv = Vector3f::new(
            (p2.x - p1.x) * cos_phi - (p2.y - p1.y) * sin_phi,
            (p2.x - p1.x) * sin_phi + (p2.y - p1.y) * cos_phi,
            p2.z - p1.z,
        );

        // Compute hyperboloid $\dndu$ and $\dndv$
        let d2pduu = -phi_max * phi_max * Vector3f::new(p_hit.x, p_hit.y, 0.0);
        let d2pduv = phi_max * Vector3f::new(-dpdv.y, dpdv.x, 0.0);
        let d2pdvv = Vector3f::new(0.0, 0.0, 0.0);

        // Compute coefficients for fundamental forms
        #[allow(non_snake_case)]
        let E = Vector3f::dot(&dpdu, &dpdu);
        #[allow(non_snake_case)]
        let F = Vector3f::dot(&dpdu, &dpdv);
        #[allow(non_snake_case)]
        let G = Vector3f::dot(&dpdv, &dpdv);

        let n = self.base.calc_normal(&dpdu, &dpdv);

        let e = Vector3f::dot(&n, &d2pduu);
        let f = Vector3f::dot(&n, &d2pduv);
        let g = Vector3f::dot(&n, &d2pdvv);

        let inv_egf2 = Float::recip(E * G - F * F);
        let dndu = dpdu * ((f * F - e * G) * inv_egf2) + dpdv * ((e * F - f * E) * inv_egf2);
        let dndv = dpdu * ((g * F - f * G) * inv_egf2) + dpdv * ((f * F - g * E) * inv_egf2);

        // Compute error bounds for intersection computed with ray equation
        let px = ox + t_shape_hit * dx;
        let py = oy + t_shape_hit * dy;
        let pz = oz + t_shape_hit * dz;
        let p_error = Vector3f::new(
            px.get_absolute_error() as Float,
            py.get_absolute_error() as Float,
            pz.get_absolute_error() as Float,
        );
        let mut isect = SurfaceInteraction::new(
            &p_hit,
            &p_error,
            &Point2f::new(u, v),
            &(-ray.d),
            &n,
            &dpdu,
            &dpdv,
            &dndu,
            &dndv,
            ray.time,
            0,
        );
        isect = self
            .base
            .object_to_world
            .transform_surface_interaction(&isect);
        return Some((t_shape_hit.into(), isect));
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        // Transform _Ray_ to object space
        let (ray, o_err, d_err) = self.base.world_to_object.transform_ray(r);

        // Compute quadratic hyperboloid coefficients

        // Initialize _EFloat_ ray coordinate values
        let p1 = self.p1;
        let p2 = self.p2;
        let ah = EFloat::from(self.ah);
        let ch = EFloat::from(self.ch);
        let z_min = self.z_min;
        let z_max = self.z_max;
        let phi_max = self.phi_max;

        let ox = EFloat::from((ray.o.x, o_err.x));
        let oy = EFloat::from((ray.o.y, o_err.y));
        let oz = EFloat::from((ray.o.z, o_err.z));
        let dx = EFloat::from((ray.d.x, d_err.x));
        let dy = EFloat::from((ray.d.y, d_err.y));
        let dz = EFloat::from((ray.d.z, d_err.z));
        let a = ah * (dx * dx) + ah * (dy * dy) - ch * (dz * dz);
        let b = 2.0 * (ah * (dx * ox) + ah * (dy * oy) - ch * (dz * oz));
        let c = ah * (ox * ox) + ah * (oy * oy) - ch * (oz * oz) - EFloat::from(1.0);

        // Solve quadratic equation for _t_ values
        let t_max = ray.t_max.get();

        if let Some((t0, t1)) = EFloat::quadratic(a, b, c) {
            // pbrt-r3:
            if t0.v.is_infinite() || t1.v.is_infinite() {
                return false;
            }

            assert!(t0.v.is_finite());
            assert!(t1.v.is_finite());
            assert!(t0.v <= t1.v);
            // Check quadric shape _t0_ and _t1_ for nearest intersection
            if t0.upper_bound() > t_max || t1.lower_bound() <= 0.0 {
                return false;
            }

            // Compute hyperboloid inverse mapping
            let mut t_shape_hit = t0;
            if t_shape_hit.lower_bound() <= 0.0 {
                t_shape_hit = t1;
                if t_max < t_shape_hit.upper_bound() {
                    return false;
                }
            }

            // Compute hyperboloid inverse mapping
            let mut p_hit = ray.o + ray.d * Float::from(t_shape_hit);
            let mut v = (p_hit.z - p1.z) / (p2.z - p1.z);
            let pr = (1.0 - v) * p1 + v * p2;
            let mut phi = Float::atan2(
                pr.x * p_hit.y - p_hit.x * pr.y,
                p_hit.x * pr.x + p_hit.y * pr.y,
            );
            if phi < 0.0 {
                phi += 2.0 * PI;
            }

            // Test hyperboloid intersection against clipping parameters
            if p_hit.z < z_min || p_hit.z > z_max || phi > phi_max {
                if t_shape_hit == t1 {
                    return false;
                }
                if t1.upper_bound() > t_max {
                    return false;
                }
                t_shape_hit = t1;
                // Compute hyperboloid inverse mapping
                p_hit = ray.o + ray.d * Float::from(t_shape_hit);
                v = (p_hit.z - p1.z) / (p2.z - p1.z);
                let pr = (1.0 - v) * p1 + v * p2;
                phi = Float::atan2(
                    pr.x * p_hit.y - p_hit.x * pr.y,
                    p_hit.x * pr.x + p_hit.y * pr.y,
                );
                if phi < 0.0 {
                    phi += 2.0 * PI;
                }
                if p_hit.z < z_min || p_hit.z > z_max || phi > phi_max {
                    return false;
                }
            }
            return true;
        } else {
            return false;
        }
    }

    fn area(&self) -> Float {
        let p1 = self.p1;
        let p2 = self.p2;
        let phi_max = self.phi_max;
        return phi_max / 6.0
            * (2.0 * quad(p1.x) - 2.0 * p1.x * p1.x * p1.x * p2.x
                + 2.0 * quad(p2.x)
                + 2.0
                    * (p1.y * p1.y + p1.y * p2.y + p2.y * p2.y)
                    * (sqr(p1.y - p2.y) + sqr(p1.z - p2.z))
                + p2.x
                    * p2.x
                    * (5.0 * p1.y * p1.y + 2.0 * p1.y * p2.y - 4.0 * p2.y * p2.y
                        + 2.0 * sqr(p1.z - p2.z))
                + p1.x
                    * p1.x
                    * (-4.0 * p1.y * p1.y
                        + 2.0 * p1.y * p2.y
                        + 5.0 * p2.y * p2.y
                        + 2.0 * sqr(p1.z - p2.z))
                - 2.0
                    * p1.x
                    * p2.x
                    * (p2.x * p2.x - p1.y * p1.y + 5.0 * p1.y * p2.y - p2.y * p2.y - p1.z * p1.z
                        + 2.0 * p1.z * p2.z
                        - p2.z * p2.z));
    }

    fn sample(&self, _u: &Point2f) -> Option<(Interaction, Float)> {
        //LOG(FATAL) << "Hyperboloid::Sample not implemented.";
        return None;
    }
}

pub fn create_hyperboloid_shape(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
) -> Result<Hyperboloid, PbrtError> {
    let p1 = params.find_one_point3f("p1", &Point3f::new(0.0, 0.0, 0.0));
    let p2 = params.find_one_point3f("p2", &Point3f::new(1.0, 1.0, 1.0));
    let phimax = params.find_one_float("phimax", 360.0);

    return Ok(Hyperboloid::new(
        o2w,
        w2o,
        reverse_orientation,
        &p1,
        &p2,
        phimax,
    ));
}
