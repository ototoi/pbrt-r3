use crate::core::prelude::*;

#[inline]
fn radians(x: Float) -> Float {
    return x * (PI / 180.0);
}

pub struct Cone {
    pub base: BaseShape,
    pub radius: Float,
    pub height: Float,
    pub phi_max: Float,
}

impl Cone {
    pub fn new(
        o2w: &Transform,
        w2o: &Transform,
        reverse_orientation: bool,
        height: Float,
        radius: Float,
        phi_max: Float,
    ) -> Self {
        let phi_max = radians(Float::clamp(phi_max, 0.0, 360.0));
        Cone {
            base: BaseShape::new(o2w, w2o, reverse_orientation),
            radius,
            height,
            phi_max,
        }
    }
}

impl Shape for Cone {
    fn object_bound(&self) -> Bounds3f {
        let radius = self.radius;
        let height = Float::max(self.height, MACHINE_EPSILON);
        return Bounds3f::new(
            &Point3f::new(-radius, -radius, 0.0),
            &Point3f::new(radius, radius, height),
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

        // Compute quadratic cone coefficients

        // Initialize _EFloat_ ray coordinate values
        let radius = self.radius;
        let height = self.height;
        let phi_max = self.phi_max;

        let ox = EFloat::from((ray.o.x, o_err.x));
        let oy = EFloat::from((ray.o.y, o_err.y));
        let oz = EFloat::from((ray.o.z, o_err.z));
        let dx = EFloat::from((ray.d.x, d_err.x));
        let dy = EFloat::from((ray.d.y, d_err.y));
        let dz = EFloat::from((ray.d.z, d_err.z));
        let k = EFloat::from(radius) / EFloat::from(height);
        let k = k * k;
        let height_ = EFloat::from(height);
        let a = dx * dx + dy * dy - k * dz * dz;
        let b = (dx * ox + dy * oy - k * dz * (oz - height_)) * 2.0;
        let c = ox * ox + oy * oy - k * (oz - height_) * (oz - height_);

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

        // Compute cone inverse mapping
        let mut t_shape_hit = t0;
        if t_shape_hit.lower_bound() <= 0.0 {
            t_shape_hit = t1;
            if t_max < t_shape_hit.upper_bound() {
                return None;
            }
        }

        // Compute cone inverse mapping
        let mut p_hit = ray.o + ray.d * Float::from(t_shape_hit);
        let mut phi = Float::atan2(p_hit.y, p_hit.x);
        if phi < 0.0 {
            phi += 2.0 * PI;
        }

        // Test cone intersection against clipping parameters
        if p_hit.z < 0.0 || p_hit.z > height || phi > phi_max {
            if t_shape_hit == t1 {
                return None;
            }
            if t1.upper_bound() > t_max {
                return None;
            }
            t_shape_hit = t1;
            // Compute cone inverse mapping
            p_hit = ray.o + ray.d * Float::from(t_shape_hit);
            phi = Float::atan2(p_hit.y, p_hit.x);
            if phi < 0.0 {
                phi += 2.0 * PI;
            }
            if p_hit.z < 0.0 || p_hit.z > height || phi > phi_max {
                return None;
            }
        }

        // Find parametric representation of cone hit
        let u = phi / phi_max;
        let v = p_hit.z / height;

        // Compute cone $\dpdu$ and $\dpdv$
        let dpdu = Vector3f::new(-phi_max * p_hit.y, phi_max * p_hit.x, 0.0);
        let dpdv = Vector3f::new(-p_hit.x / (1.0 - v), -p_hit.y / (1.0 - v), height);

        // Compute cone $\dndu$ and $\dndv$
        let d2pduu = -phi_max * phi_max * Vector3f::new(p_hit.x, p_hit.y, 0.0);
        let d2pduv = phi_max / (1.0 - v) * Vector3f::new(p_hit.y, -p_hit.x, 0.0);
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
        let (ray, o_err, d_err) = self.base.world_to_object.transform_ray(r);

        // Compute quadratic cone coefficients

        // Initialize _EFloat_ ray coordinate values
        let radius = self.radius;
        let height = self.height;
        let phi_max = self.phi_max;

        let ox = EFloat::from((ray.o.x, o_err.x));
        let oy = EFloat::from((ray.o.y, o_err.y));
        let oz = EFloat::from((ray.o.z, o_err.z));
        let dx = EFloat::from((ray.d.x, d_err.x));
        let dy = EFloat::from((ray.d.y, d_err.y));
        let dz = EFloat::from((ray.d.z, d_err.z));
        let k = EFloat::from(radius) / EFloat::from(height);
        let k = k * k;
        let height_ = EFloat::from(height);
        let a = dx * dx + dy * dy - k * dz * dz;
        let b = (dx * ox + dy * oy - k * dz * (oz - height_)) * 2.0;
        let c = ox * ox + oy * oy - k * (oz - height_) * (oz - height_);

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

            // Compute cone inverse mapping
            let mut t_shape_hit = t0;
            if t_shape_hit.lower_bound() <= 0.0 {
                t_shape_hit = t1;
                if t_max < t_shape_hit.upper_bound() {
                    return false;
                }
            }

            // Compute cone inverse mapping
            let mut p_hit = ray.o + ray.d * Float::from(t_shape_hit);
            let mut phi = Float::atan2(p_hit.y, p_hit.x);
            if phi < 0.0 {
                phi += 2.0 * PI;
            }

            // Test cone intersection against clipping parameters
            if p_hit.z < 0.0 || p_hit.z > height || phi > phi_max {
                if t_shape_hit == t1 {
                    return false;
                }
                if t1.upper_bound() > t_max {
                    return false;
                }
                t_shape_hit = t1;
                // Compute cone inverse mapping
                p_hit = ray.o + ray.d * Float::from(t_shape_hit);
                phi = Float::atan2(p_hit.y, p_hit.x);
                if phi < 0.0 {
                    phi += 2.0 * PI;
                }
                if p_hit.z < 0.0 || p_hit.z > height || phi > phi_max {
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
        let height = self.height;
        let phi_max = self.phi_max;
        return radius * Float::sqrt((height * height) + (radius * radius)) * phi_max / 2.0;
    }

    fn sample(&self, _u: &Point2f) -> Option<(Interaction, Float)> {
        //log
        return None;
    }
}

pub fn create_cone_shape(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
) -> Result<Cone, PbrtError> {
    let height = params.find_one_float("height", 1.0);
    let radius = params.find_one_float("radius", 1.0);
    let phimax = params.find_one_float("phimax", 360.0);

    // pbrt-r3
    if height == 0.0 || radius == 0.0 {
        let msg = format!(
            "Unable to create cone shape: height={}, radius={}, phimax={}",
            height, radius, phimax
        );
        return Err(PbrtError::error(&msg));
    // pbrt-r3
    } else {
        return Ok(Cone::new(
            o2w,
            w2o,
            reverse_orientation,
            height,
            radius,
            phimax,
        ));
    }
}
