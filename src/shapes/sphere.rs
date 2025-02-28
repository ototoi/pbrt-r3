use crate::core::efloat::EFloat;
use crate::core::pbrt::types::*;
use crate::core::shape::*;

const MACHINE_EPSILON: Float = Float::EPSILON * 0.5;
const GAMMA5: Float = (5.0 * MACHINE_EPSILON) / (1.0 - (5.0 * MACHINE_EPSILON));

#[derive(Debug, PartialEq, Clone)]
pub struct Sphere {
    pub base: BaseShape,
    pub radius: Float,
    pub z_min: Float,
    pub z_max: Float,
    pub theta_min: Float,
    pub theta_max: Float,
    pub phi_max: Float,
}

impl Sphere {
    pub fn new(
        o2w: &Transform,
        w2o: &Transform,
        reverse_orientation: bool,
        radius: Float,
        z_min: Float,
        z_max: Float,
        phi_max: Float,
    ) -> Self {
        let z_min = Float::clamp(Float::min(z_min, z_max), -radius, radius);
        let z_max = Float::clamp(Float::max(z_min, z_max), -radius, radius);
        let theta_min = Float::acos(Float::clamp(z_min / radius, -1.0, 1.0));
        let theta_max = Float::acos(Float::clamp(z_max / radius, -1.0, 1.0));
        let phi_max = radians(Float::clamp(phi_max, 0.0, 360.0));
        Sphere {
            base: BaseShape::new(o2w, w2o, reverse_orientation),
            radius,
            z_min,
            z_max,
            theta_min,
            theta_max,
            phi_max,
        }
    }
}

impl Shape for Sphere {
    fn object_bound(&self) -> Bounds3f {
        let radius = self.radius * 1.001;
        let diff = radius - self.radius;
        return Bounds3f::new(
            &Vector3f::new(-radius, -radius, self.z_min - diff),
            &Vector3f::new(radius, radius, self.z_max + diff),
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
        let ox = EFloat::from((ray.o.x, o_err.x));
        let oy = EFloat::from((ray.o.y, o_err.y));
        let oz = EFloat::from((ray.o.z, o_err.z));
        let dx = EFloat::from((ray.d.x, d_err.x));
        let dy = EFloat::from((ray.d.y, d_err.y));
        let dz = EFloat::from((ray.d.z, d_err.z));
        let rad = EFloat::from((self.radius, 0.0));
        let a = dx * dx + dy * dy + dz * dz;
        let b = (dx * ox + dy * oy + dz * oz) * 2.0;
        let c = ox * ox + oy * oy + oz * oz - rad * rad;

        let t_max = ray.t_max.get();

        let (t0, t1) = EFloat::quadratic(a, b, c)?;
        // pbrt-r3:
        if t0.v.is_infinite() || t1.v.is_infinite() {
            return None;
        }

        assert!(t0.v <= t1.v);
        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > t_max || t1.lower_bound() <= 0.0 {
            return None;
        }

        let mut t_shape_hit = t0;
        if t_shape_hit.lower_bound() <= 0.0 {
            t_shape_hit = t1;
            if t_max < t_shape_hit.upper_bound() {
                return None;
            }
        }

        let mut p_hit = ray.o + ray.d * Float::from(t_shape_hit);
        p_hit *= self.radius / Vector3f::distance(&p_hit, &Vector3f::zero());
        if p_hit.x == 0.0 && p_hit.y == 0.0 {
            p_hit.x = 1e-5 * self.radius;
        }
        let mut phi = Float::atan2(p_hit.y, p_hit.x);
        if phi < 0.0 {
            phi += 2.0 * PI;
        }
        if (self.z_min > -self.radius && p_hit.z < self.z_min)
            || (self.z_max < self.radius && p_hit.z > self.z_max)
            || (phi > self.phi_max)
        {
            if t_shape_hit == t1 {
                return None;
            }
            if t1.upper_bound() > t_max {
                return None;
            }
            t_shape_hit = t1;
            p_hit = ray.o + ray.d * Float::from(t_shape_hit);
            p_hit *= self.radius / Vector3f::distance(&p_hit, &Vector3f::zero());
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5 * self.radius;
            }
            phi = Float::atan2(p_hit.y, p_hit.x);
            if phi < 0.0 {
                phi += PI;
            }
            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || (phi > self.phi_max)
            {
                return None;
            }
        }

        let radius = self.radius;
        let theta_min = self.theta_min;
        let theta_max = self.theta_max;
        let dtheta = theta_max - theta_min;
        let phi_max = self.phi_max;

        let u = phi / phi_max;
        let theta = Float::acos(Float::clamp(p_hit.z / radius, -1.0, 1.0));
        let v = (theta - theta_min) / dtheta;
        assert!(u >= 0.0);
        assert!(v >= 0.0);

        let z_radius = Float::sqrt(p_hit.x * p_hit.x + p_hit.y * p_hit.y);
        let inv_z_radius = Float::recip(z_radius);
        let cos_phi = p_hit.x * inv_z_radius;
        let sin_phi = p_hit.y * inv_z_radius;
        let dpdu = Vector3f::new(-phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
        let dpdv = Vector3f::new(
            p_hit.z * cos_phi,
            p_hit.z * sin_phi,
            -radius * Float::sin(theta),
        ) * dtheta;

        let d2pduu = Vector3f::new(p_hit.x, p_hit.y, 0.0) * (-phi_max * phi_max);
        let d2pduv = Vector3f::new(-sin_phi, cos_phi, 0.0) * (p_hit.z * dtheta * phi_max);
        let d2pdvv = Vector3f::new(p_hit.x, p_hit.y, p_hit.z) * (-dtheta * dtheta);

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

        let p_error = GAMMA5 * p_hit.abs();
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
        let ox = EFloat::from((ray.o.x, o_err.x));
        let oy = EFloat::from((ray.o.y, o_err.y));
        let oz = EFloat::from((ray.o.z, o_err.z));
        let dx = EFloat::from((ray.d.x, d_err.x));
        let dy = EFloat::from((ray.d.y, d_err.y));
        let dz = EFloat::from((ray.d.z, d_err.z));
        let rad = EFloat::from((self.radius, 0.0));
        let a = dx * dx + dy * dy + dz * dz;
        let b = (dx * ox + dy * oy + dz * oz) * 2.0;
        let c = ox * ox + oy * oy + oz * oz - rad * rad;

        let t_max = ray.t_max.get();

        if let Some((t0, t1)) = EFloat::quadratic(a, b, c) {
            // pbrt-r3:
            if t0.v.is_infinite() || t1.v.is_infinite() {
                return false;
            }

            if t0.upper_bound() > t_max || t1.lower_bound() <= 0.0 {
                return false;
            }
            let mut t_shape_hit = t0;
            if t_shape_hit.lower_bound() <= 0.0 {
                t_shape_hit = t1;
                if t_max < t_shape_hit.upper_bound() {
                    return false;
                }
            }
            let mut p_hit = ray.o + ray.d * Float::from(t_shape_hit);
            p_hit *= self.radius / Vector3f::distance(&p_hit, &Vector3f::zero());
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5 * self.radius;
            }
            let mut phi = Float::atan2(p_hit.y, p_hit.x);
            if phi < 0.0 {
                phi += 2.0 * PI; //Float
            }
            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || phi > self.phi_max
            {
                if t_shape_hit == t1 {
                    return false;
                }
                if t1.upper_bound() > t_max {
                    return false;
                }
                t_shape_hit = t1;
                p_hit = ray.o + ray.d * Float::from(t_shape_hit);
                p_hit *= self.radius / Vector3f::distance(&p_hit, &Vector3f::zero());
                if p_hit.x == 0.0 && p_hit.y == 0.0 {
                    p_hit.x = 1e-5 * self.radius;
                }
                phi = Float::atan2(p_hit.y, p_hit.x);
                if phi < 0.0 {
                    phi += 2.0 * PI; //Float
                }
                if (self.z_min > -self.radius && p_hit.z < self.z_min)
                    || (self.z_max < self.radius && p_hit.z > self.z_max)
                    || phi > self.phi_max
                {
                    return false;
                }
            }
            return true;
        } else {
            return false;
        }
    }

    fn area(&self) -> Float {
        return self.phi_max * self.radius * (self.z_max - self.z_min);
    }

    fn sample(&self, u: &Point2f) -> Option<(Interaction, Float)> {
        let mut p_obj = Point3f::zero() + self.radius * uniform_sample_sphere(u);
        let mut n = self
            .base
            .object_to_world
            .transform_normal(&p_obj)
            .normalize();
        if self.base.reverse_orientation {
            n *= -1.0;
        }
        p_obj *= self.radius / Point3f::distance(&p_obj, &Point3f::zero());
        let p_obj_error = GAMMA5 * p_obj.abs();
        let (p, p_error) = self
            .base
            .object_to_world
            .transform_point_with_abs_error(&p_obj, &p_obj_error);
        let it = Interaction::from_surface_sample(&p, &p_error, &n);
        let pdf = Float::recip(self.area());
        return Some((it, pdf));
    }

    fn sample_from(&self, ref_: &Interaction, u: &Point2f) -> Option<(Interaction, Float)> {
        let p_center = self.base.object_to_world.transform_point(&Point3f::zero());
        // Sample uniformly on sphere if $\pt{}$ is inside it
        let p_origin = offset_ray_origin(
            &ref_.get_p(),
            &ref_.get_p_error(),
            &ref_.get_n(),
            &(p_center - ref_.get_p()),
        );
        let radius = self.radius;
        if Vector3f::distance_squared(&p_origin, &p_center) <= radius * radius {
            let (intr, pdf) = self.sample(u)?;
            let wi = intr.get_p() - ref_.get_p();
            if wi.length_squared() == 0.0 {
                return None;
            } else {
                // Convert from area measure returned by Sample() call above to
                // solid angle measure.
                let wi = wi.normalize();
                let pdf = pdf * Vector3f::distance_squared(&intr.get_p(), &ref_.get_p())
                    / Vector3f::abs_dot(&intr.get_n(), &-wi);
                if pdf <= 0.0 || pdf.is_infinite() {
                    return None;
                }
                return Some((intr, pdf));
            }
        }
        // Sample sphere uniformly inside subtended cone

        // Compute coordinate system for sphere sampling
        //
        let dc = Vector3f::distance(&ref_.get_p(), &p_center);
        let inv_dc = 1.0 / dc;
        let wc = (p_center - ref_.get_p()) * inv_dc;
        let (wc_x, wc_y) = coordinate_system(&wc);

        // Compute $\theta$ and $\phi$ values for sample in cone
        let sin_theta_max = radius * inv_dc;
        let sin_theta_max2 = sin_theta_max * sin_theta_max;
        let inv_sin_theta_max = 1.0 / sin_theta_max;
        let cos_theta_max = Float::sqrt(Float::max(0.0, 1.0 - sin_theta_max2));

        let pdf = 1.0 / (2.0 * PI * (1.0 - cos_theta_max));
        if pdf <= 0.0 || pdf.is_infinite() {
            return None;
        }

        let mut cos_theta = (cos_theta_max - 1.0) * u[0] + 1.0;
        let mut sin_theta2 = 1.0 - cos_theta * cos_theta;

        if sin_theta_max2 < 0.00068523
        //sin^2(1.5 deg)
        {
            /* Fall back to a Taylor series expansion for small angles, where
            the standard approach suffers from severe cancellation errors */
            sin_theta2 = Float::max(0.0, sin_theta_max2 * u[0]);
            cos_theta = Float::sqrt(1.0 - sin_theta2);
        }

        // Compute angle $\alpha$ from center of sphere to sampled point on surface
        let cos_alpha = sin_theta2 * inv_sin_theta_max
            + cos_theta
                * Float::sqrt(Float::max(
                    0.0,
                    1.0 - sin_theta2 * inv_sin_theta_max * inv_sin_theta_max,
                ));
        let sin_alpha = Float::sqrt(Float::max(0.0, 1.0 - cos_alpha * cos_alpha));
        let phi = u[1] * 2.0 * PI;

        // Compute surface normal and sampled point on sphere
        let n_world = spherical_direction_axes(sin_alpha, cos_alpha, phi, &-wc_x, &-wc_y, &-wc);
        let p_world = p_center + radius * Point3f::from(n_world);

        // Return _Interaction_ for sampled point on sphere
        let p = p_world;
        let p_error = GAMMA5 * Vector3::abs(&p_world);
        let mut n = Normal3f::from(n_world);
        if self.base.reverse_orientation {
            n *= -1.0;
        }
        let it = Interaction::from_surface_sample(&p, &p_error, &n);
        return Some((it, pdf));
    }
}

#[inline]
fn radians(x: Float) -> Float {
    return x * (PI / 180.0);
}
pub fn create_sphere_shape(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
) -> Result<Sphere, PbrtError> {
    let radius = params.find_one_float("radius", 1.0);
    let zmin = params.find_one_float("zmin", -radius);
    let zmax = params.find_one_float("zmax", radius);
    let phimax = params.find_one_float("phimax", 360.0);

    return Ok(Sphere::new(
        o2w,
        w2o,
        reverse_orientation,
        radius,
        zmin,
        zmax,
        phimax,
    ));
}

//------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let mut params = ParamSet::new();
        params.add_float("radius", 1.0);
        let t = Transform::identity();
        let s1 = create_sphere_shape(&t, &t, false, &params).unwrap();
        assert_eq!(s1.radius, 1.0);
    }

    #[test]
    fn test_002() {
        let mut params = ParamSet::new();
        params.add_float("radius", 1.0);
        let t = Transform::identity();
        let s1 = create_sphere_shape(&t, &t, false, &params).unwrap();
        let d = Vector3f::new(0.0, 0.0, 1.0);
        let p1 = s1.intersect_p(&Ray::new(&Point3f::new(0.0, 0.0, -5.0), &d, 1000.0, 0.0));
        let p2 = s1.intersect_p(&Ray::new(&Point3f::new(0.0, 1.0, -5.0), &d, 1000.0, 0.0));
        let p3 = s1.intersect_p(&Ray::new(&Point3f::new(0.0, 1.2, -5.0), &d, 1000.0, 0.0));
        let p4 = s1.intersect_p(&Ray::new(&Point3f::new(1.0, 0.0, -5.0), &d, 1000.0, 0.0));
        let p5 = s1.intersect_p(&Ray::new(&Point3f::new(1.2, 0.0, -5.0), &d, 1000.0, 0.0));

        assert_eq!(p1, true);
        assert_eq!(p2, true);
        assert_eq!(p3, false);
        assert_eq!(p4, true);
        assert_eq!(p5, false);
    }
}
