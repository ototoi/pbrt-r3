use crate::core::pbrt::*;

use std::sync::Arc;

pub struct AlphaMaskShape {
    shape: Arc<dyn Shape>,
    mask: Arc<dyn Texture<Float>>,
}

impl AlphaMaskShape {
    pub fn new(shape: &Arc<dyn Shape>, mask: &Arc<dyn Texture<Float>>) -> Self {
        AlphaMaskShape {
            shape: Arc::clone(shape),
            mask: Arc::clone(mask),
        }
    }
}

impl Shape for AlphaMaskShape {
    fn object_bound(&self) -> Bounds3f {
        let shape = self.shape.as_ref();
        return shape.object_bound();
    }
    fn world_bound(&self) -> Bounds3f {
        let shape = self.shape.as_ref();
        return shape.world_bound();
    }
    fn intersect(&self, r: &Ray) -> Option<(Float, SurfaceInteraction)> {
        let shape = self.shape.as_ref();
        let t_max = r.t_max.get();
        if let Some((t, si)) = shape.intersect(r) {
            let mask = self.mask.as_ref();
            let a = mask.evaluate(&si);
            if a > 0.0 {
                return Some((t, si));
            } else {
                r.t_max.set(t_max);
                return None;
            }
        }
        return None;
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        let shape = self.shape.as_ref();
        let t_max = r.t_max.get();
        if let Some((_, si)) = shape.intersect(r) {
            let mask = self.mask.as_ref();
            let a = mask.evaluate(&si);
            if a > 0.0 {
                return true;
            } else {
                r.t_max.set(t_max);
                return false;
            }
        }
        return false;
    }
    fn area(&self) -> Float {
        let shape = self.shape.as_ref();
        return shape.area();
    }
    fn pdf(&self, inter: &Interaction) -> Float {
        // Ignore any alpha textures used for trimming the shape when performing
        // this intersection. Hack for the "San Miguel" scene, where this is used
        // to make an invisible area light.
        let shape = self.shape.as_ref();
        return shape.pdf(inter);
    }
    fn pdf_from(&self, inter: &Interaction, wi: &Vector3f) -> Float {
        let shape = self.shape.as_ref();
        return shape.pdf_from(inter, wi);
    }
    fn sample(&self, u: &Point2f) -> Option<(Interaction, Float)> {
        let shape = self.shape.as_ref();
        return shape.sample(u);
    }
    fn sample_from(&self, inter: &Interaction, u: &Point2f) -> Option<(Interaction, Float)> {
        let shape = self.shape.as_ref();
        return shape.sample_from(inter, u);
    }
    fn solid_angle(&self, p: &Point3f, n_samples: i32) -> Float {
        let shape = self.shape.as_ref();
        return shape.solid_angle(p, n_samples);
    }
}
