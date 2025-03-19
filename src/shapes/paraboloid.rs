use crate::core::camera::*;
use crate::core::distribution::*;
use crate::core::efloat::EFloat;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::lightdistrib::*;
use crate::core::lowdiscrepancy::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::rng::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use crate::core::texture::*;
use crate::core::transform::*;

#[inline]
fn radians(x: Float) -> Float {
    return x * (PI / 180.0);
}

pub struct Paraboloid {
    pub base: BaseShape,
    pub radius: Float,
    pub z_min: Float,
    pub z_max: Float,
    pub phi_max: Float,
}

impl Paraboloid {
    pub fn new(
        o2w: &Transform,
        w2o: &Transform,
        reverse_orientation: bool,
        radius: Float,
        z_min: Float,
        z_max: Float,
        phi_max: Float,
    ) -> Self {
        let phi_max = radians(Float::clamp(phi_max, 0.0, 360.0));
        Paraboloid {
            base: BaseShape::new(o2w, w2o, reverse_orientation),
            radius,
            z_min,
            z_max,
            phi_max,
        }
    }
}

impl Shape for Paraboloid {
    fn object_bound(&self) -> Bounds3f {
        let radius = self.radius;
        let z_min = self.z_min;
        let z_max = self.z_max;
        return Bounds3f::new(
            &Point3f::new(-radius, -radius, z_min),
            &Point3f::new(radius, radius, z_max),
        );
    }

    fn world_bound(&self) -> Bounds3f {
        return self
            .base
            .object_to_world
            .transform_bounds(&self.object_bound());
    }
    fn intersect(&self, r: &Ray) -> Option<(Float, SurfaceInteraction)> {
        let (ray, o_err, d_err) = self.base.world_to_object.transform_ray(r);

        // Compute quadratic paraboloid coefficients

        // Initialize _EFloat_ ray coordinate values
        let radius = self.radius;
        let z_min = self.z_min;
        let z_max = self.z_max;
        let phi_max = self.phi_max;

        let ox = EFloat::from((ray.o.x, o_err.x));
        let oy = EFloat::from((ray.o.y, o_err.y));
        let oz = EFloat::from((ray.o.z, o_err.z));
        let dx = EFloat::from((ray.d.x, d_err.x));
        let dy = EFloat::from((ray.d.y, d_err.y));
        let dz = EFloat::from((ray.d.z, d_err.z));
        let k = EFloat::from(z_max) / EFloat::from(radius) * EFloat::from(radius);
        let a = k * (dx * dx + dy * dy);
        let b = ((dx * ox + dy * oy) * k * 2.0) - dz;
        let c = k * (ox * ox + oy * oy) - oz;

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

        // Compute paraboloid inverse mapping
        let mut t_shape_hit = t0;
        if t_shape_hit.lower_bound() <= 0.0 {
            t_shape_hit = t1;
            if t_max < t_shape_hit.upper_bound() {
                return None;
            }
        }

        // Compute paraboloid inverse mapping
        let mut p_hit = ray.o + ray.d * Float::from(t_shape_hit);
        let mut phi = Float::atan2(p_hit.y, p_hit.x);
        if phi < 0.0 {
            phi += 2.0 * PI;
        }

        // Test paraboloid intersection against clipping parameters
        if p_hit.z < z_min || p_hit.z > z_max || phi > phi_max {
            if t_shape_hit == t1 {
                return None;
            }
            if t1.upper_bound() > t_max {
                return None;
            }
            t_shape_hit = t1;
            // Compute paraboloid inverse mapping
            p_hit = ray.o + ray.d * Float::from(t_shape_hit);
            phi = Float::atan2(p_hit.y, p_hit.x);
            if phi < 0.0 {
                phi += 2.0 * PI;
            }
            if p_hit.z < z_min || p_hit.z > z_max || phi > phi_max {
                return None;
            }
        }

        // Find parametric representation of paraboloid hit
        let u = phi / phi_max;
        let v = (p_hit.z - z_min) / (z_max - z_min);

        // Compute paraboloid $\dpdu$ and $\dpdv$
        let dpdu = Vector3f::new(-phi_max * p_hit.y, phi_max * p_hit.x, 0.0);
        let dpdv = (z_max - z_min)
            * Vector3f::new(p_hit.x / (2.0 * p_hit.z), p_hit.y / (2.0 * p_hit.z), 1.0);

        // Compute paraboloid $\dndu$ and $\dndv$
        let d2pduu = -phi_max * phi_max * Vector3f::new(p_hit.x, p_hit.y, 0.0);
        let d2pduv = (z_max - z_min)
            * phi_max
            * Vector3f::new(-p_hit.y / (2.0 * p_hit.z), p_hit.x / (2.0 * p_hit.z), 0.0);
        let d2pdvv = -(z_max - z_min)
            * (z_max - z_min)
            * Vector3f::new(
                p_hit.x / (4.0 * p_hit.z * p_hit.z),
                p_hit.y / (4.0 * p_hit.z * p_hit.z),
                0.0,
            );

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
        let (ray, o_err, d_err) = self.base.world_to_object.transform_ray(r);

        // Compute quadratic paraboloid coefficients

        // Initialize _EFloat_ ray coordinate values
        let radius = self.radius;
        let z_min = self.z_min;
        let z_max = self.z_max;
        let phi_max = self.phi_max;

        let ox = EFloat::from((ray.o.x, o_err.x));
        let oy = EFloat::from((ray.o.y, o_err.y));
        let oz = EFloat::from((ray.o.z, o_err.z));
        let dx = EFloat::from((ray.d.x, d_err.x));
        let dy = EFloat::from((ray.d.y, d_err.y));
        let dz = EFloat::from((ray.d.z, d_err.z));
        let k = EFloat::from(z_max) / EFloat::from(radius) * EFloat::from(radius);
        let a = k * (dx * dx + dy * dy);
        let b = ((dx * ox + dy * oy) * k * 2.0) - dz;
        let c = k * (ox * ox + oy * oy) - oz;

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

            // Compute paraboloid inverse mapping
            let mut t_shape_hit = t0;
            if t_shape_hit.lower_bound() <= 0.0 {
                t_shape_hit = t1;
                if t_max < t_shape_hit.upper_bound() {
                    return false;
                }
            }

            // Compute paraboloid inverse mapping
            let mut p_hit = ray.o + ray.d * Float::from(t_shape_hit);
            let mut phi = Float::atan2(p_hit.y, p_hit.x);
            if phi < 0.0 {
                phi += 2.0 * PI;
            }

            // Test paraboloid intersection against clipping parameters
            if p_hit.z < z_min || p_hit.z > z_max || phi > phi_max {
                if t_shape_hit == t1 {
                    return false;
                }
                if t1.upper_bound() > t_max {
                    return false;
                }
                t_shape_hit = t1;
                // Compute paraboloid inverse mapping
                p_hit = ray.o + ray.d * Float::from(t_shape_hit);
                phi = Float::atan2(p_hit.y, p_hit.x);
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
        let radius = self.radius;
        let z_min = self.z_min;
        let z_max = self.z_max;
        let phi_max = self.phi_max;
        let radius2 = radius * radius;
        let k = 4.0 * z_max / radius2;

        return (radius2 * radius2 * phi_max / (12.0 * z_max * z_max))
            * (Float::powf(k * z_max + 1.0, 1.5) - Float::powf(k * z_min + 1.0, 1.5));
    }

    fn sample(&self, _u: &Point2f) -> Option<(Interaction, Float)> {
        //log
        return None;
    }
}

pub fn create_paraboloid_shape(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
) -> Result<Paraboloid, PbrtError> {
    let radius = params.find_one_float("radius", 1.0);
    let mut zmin = params.find_one_float("zmin", 0.0);
    let mut zmax = params.find_one_float("zmax", 1.0);
    let phimax = params.find_one_float("phimax", 360.0);

    if zmin > zmax {
        std::mem::swap(&mut zmin, &mut zmax);
    }

    // pbrt-r3
    if radius == 0.0 {
        let msg = format!(
            "Unable to create paraboloid shape: radius={}, phimax={}",
            radius, phimax
        );
        return Err(PbrtError::error(&msg));
    // pbrt-r3
    } else {
        return Ok(Paraboloid::new(
            o2w,
            w2o,
            reverse_orientation,
            radius,
            zmin,
            zmax,
            phimax,
        ));
    }
}
