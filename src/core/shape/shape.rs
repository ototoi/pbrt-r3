pub use crate::core::pbrt::*;

pub trait Shape {
    fn object_bound(&self) -> Bounds3f;
    fn world_bound(&self) -> Bounds3f;
    /* pbrt-r3:
    The below functions - intersect/intersect_p had an argument variable "bool testAlphaTexture" in the original code,
    but since alpha masking has been changed to the responsibility of the class AlphaMaskShape, this variable has been removed. */
    fn intersect(&self, r: &Ray) -> Option<(Float, SurfaceInteraction)>;
    fn intersect_p(&self, r: &Ray) -> bool;
    fn area(&self) -> Float;
    fn sample(&self, u: &Point2f) -> Option<(Interaction, Float)>;
    fn pdf(&self, _inter: &Interaction) -> Float {
        Float::recip(self.area())
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

    fn pdf_from(&self, inter: &Interaction, wi: &Vector3f) -> Float {
        let ray = inter.spawn_ray(wi);
        if let Some((_, isect_light)) = self.intersect(&ray) {
            assert!(isect_light.n.length() > 0.0);

            let pdf = Vector3f::distance_squared(&inter.get_p(), &isect_light.p)
                / (Vector3f::abs_dot(&isect_light.n, &(-*wi)) * self.area());
            if pdf.is_infinite() {
                return 0.0;
            }
            return pdf;
        } else {
            return 0.0;
        }
    }

    fn solid_angle(&self, p: &Point3f, n_samples: i32) -> Float {
        let mut it = BaseInteraction::default();
        it.p = *p;
        it.wo = Vector3f::new(0.0, 0.0, 1.0);
        let inter = Interaction::from(it);
        let mut solid_angle = 0.0;
        for i in 0..n_samples {
            let u = Point2f::new(radical_inverse(0, i as u64), radical_inverse(1, i as u64));
            if let Some((p_shape, pdf)) = self.sample_from(&inter, &u) {
                let r = Ray::new(p, &(p_shape.get_p() - *p), 0.999, 0.0);
                if !self.intersect_p(&r) {
                    solid_angle += 1.0 / pdf
                }
            }
        }
        return solid_angle / n_samples as Float;
    }
}
