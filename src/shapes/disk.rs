use crate::core::pbrt::types::*;
use crate::core::shape::*;

#[inline]
fn radians(x: Float) -> Float {
    return x * (PI / 180.0);
}

pub struct Disk {
    pub base: BaseShape,
    pub height: Float,
    pub radius: Float,
    pub inner_radius: Float,
    pub phi_max: Float,
}

impl Disk {
    pub fn new(
        o2w: &Transform,
        w2o: &Transform,
        reverse_orientation: bool,
        height: Float,
        radius: Float,
        inner_radius: Float,
        phi_max: Float,
    ) -> Self {
        let phi_max = radians(Float::clamp(phi_max, 0.0, 360.0));
        Disk {
            base: BaseShape::new(o2w, w2o, reverse_orientation),
            height,
            radius,
            inner_radius,
            phi_max,
        }
    }
}

impl Shape for Disk {
    fn object_bound(&self) -> Bounds3f {
        let radius = self.radius;
        let height = self.height;
        return Bounds3f::new(
            &Point3f::new(-radius, -radius, height - 0.001), //originally, Point3f(-radius, -radius, height)
            &Point3f::new(radius, radius, height + 0.001),
        );
    }
    fn world_bound(&self) -> Bounds3f {
        return self
            .base
            .object_to_world
            .transform_bounds(&self.object_bound());
    }
    fn intersect(&self, r: &Ray) -> Option<(Float, SurfaceInteraction)> {
        let (ray, _o_err, _d_err) = self.base.world_to_object.transform_ray(r);

        // Compute plane intersection for disk

        // Reject disk intersections for rays parallel to the disk's plane
        let t_max = ray.t_max.get();
        let height = self.height;
        if ray.d.z == 0.0 {
            return None;
        }
        let t_shape_hit = (height - ray.o.z) / ray.d.z;
        if t_shape_hit <= 0.0 || t_shape_hit >= t_max {
            return None;
        }

        // See if hit point is inside disk radii and $\phimax$
        let radius = self.radius;
        let inner_radius = self.inner_radius;
        let mut p_hit = ray.o + ray.d * t_shape_hit; //TODO
        let dist2 = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > radius * radius || dist2 < inner_radius * inner_radius {
            return None;
        }

        // Test disk $\phi$ value against $\phimax$
        let phi_max = self.phi_max;
        let mut phi = Float::atan2(p_hit.y, p_hit.x);
        if phi < 0.0 {
            phi += 2.0 * PI;
        }
        if phi > phi_max {
            return None;
        }

        let u = phi / phi_max;
        let r_hit = Float::sqrt(dist2);
        let v = (radius - r_hit) / (radius - inner_radius);
        let dpdu = Vector3f::new(-phi_max * p_hit.y, phi_max * p_hit.x, 0.0);
        let dpdv = Vector3f::new(p_hit.x, p_hit.y, 0.0) * ((inner_radius - radius) / r_hit);
        let dndu = Vector3f::new(0.0, 0.0, 0.0);
        let dndv = Vector3f::new(0.0, 0.0, 0.0);

        let mut n = self.base.calc_normal(&dpdu, &dpdv);
        if Vector3f::dot(&ray.d, &n) > 0.0 {
            n *= -1.0;
        }

        // Refine disk intersection point
        p_hit.z = height;

        // Compute error bounds for disk intersection
        let p_error = Vector3f::new(0.0, 0.0, 0.0);

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
        return Some((t_shape_hit, isect));
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        let (ray, _o_err, _d_err) = self.base.world_to_object.transform_ray(r);

        // Compute plane intersection for disk

        // Reject disk intersections for rays parallel to the disk's plane
        let t_max = ray.t_max.get();
        let height = self.height;
        if ray.d.z == 0.0 {
            return false;
        }
        let t_shape_hit = (height - ray.o.z) / ray.d.z;
        if t_shape_hit <= 0.0 || t_shape_hit >= t_max {
            return false;
        }

        // See if hit point is inside disk radii and $\phimax$
        let radius = self.radius;
        let inner_radius = self.inner_radius;
        let p_hit = ray.o + ray.d * t_shape_hit; //TODO
        let dist2 = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > radius * radius || dist2 < inner_radius * inner_radius {
            return false;
        }

        // Test disk $\phi$ value against $\phimax$
        let phi_max = self.phi_max;
        let mut phi = Float::atan2(p_hit.y, p_hit.x);
        if phi < 0.0 {
            phi += 2.0 * PI;
        }
        if phi > phi_max {
            return false;
        }

        return true;
    }

    fn area(&self) -> Float {
        let radius = self.radius;
        let inner_radius = self.inner_radius;
        let phi_max = self.phi_max;
        return phi_max * 0.5 * (radius * radius - inner_radius * inner_radius);
    }

    fn sample(&self, u: &Point2f) -> Option<(Interaction, Float)> {
        let radius = self.radius;
        let height = self.height;
        let pd = concentric_sample_disk(u);
        let p_obj = Point3f::new(pd.x * radius, pd.y * radius, height);
        let mut n = self
            .base
            .object_to_world
            .transform_normal(&Normal3f::new(0.0, 0.0, 1.0));
        if self.base.reverse_orientation {
            n *= -1.0;
        }
        let (p, p_error) = self
            .base
            .object_to_world
            .transform_point_with_abs_error(&p_obj, &Point3f::new(0.0, 0.0, 0.0));
        let pdf = 1.0 / self.area();
        let it = Interaction::from((p, p_error, n, 0.0));
        return Some((it, pdf));
    }
}

pub fn create_disk_shape(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
) -> Result<Disk, PbrtError> {
    let height = params.find_one_float("height", 0.0);
    let radius = params.find_one_float("radius", 1.0);
    let inner_radius = params.find_one_float("innerradius", 0.0);
    let phimax = params.find_one_float("phimax", 360.0);

    return Ok(Disk::new(
        o2w,
        w2o,
        reverse_orientation,
        height,
        radius,
        inner_radius,
        phimax,
    ));
}
