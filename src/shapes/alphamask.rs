use crate::core::camera::*;
use crate::core::distribution::*;
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

use std::sync::Arc;

#[derive(Clone)]
pub enum AlphaMaskInfo {
    Texture { texture: Arc<dyn Texture<Float>> },
    Value { value: Float },
}

pub struct AlphaMaskShape {
    shape: Arc<dyn Shape>,
    test_intersection: bool,
    test_intersection_p: bool,
    alpha_mask_texture: Option<Arc<dyn Texture<Float>>>,
    shadow_alpha_mask_texture: Option<Arc<dyn Texture<Float>>>,
}

impl AlphaMaskShape {
    pub fn new(
        shape: &Arc<dyn Shape>,
        alpha_mask_info: &Option<AlphaMaskInfo>,
        shadow_alpha_mask_info: &Option<AlphaMaskInfo>,
    ) -> Self {
        let mut test_intersection = true;
        let mut test_intersection_p = true;
        let mut alpha_mask_texture = None;
        let mut shadow_alpha_mask_texture = None;
        if let Some(info) = alpha_mask_info.as_ref() {
            match info {
                AlphaMaskInfo::Texture { texture } => {
                    alpha_mask_texture = Some(Arc::clone(texture));
                }
                AlphaMaskInfo::Value { value: alpha } => {
                    if *alpha <= 0.0 {
                        test_intersection = false;
                        test_intersection_p = false;
                    }
                }
            }
        }
        if let Some(info) = shadow_alpha_mask_info.as_ref() {
            match info {
                AlphaMaskInfo::Texture { texture } => {
                    shadow_alpha_mask_texture = Some(Arc::clone(texture));
                }
                AlphaMaskInfo::Value { value: alpha } => {
                    if *alpha <= 0.0 {
                        test_intersection_p = false;
                    }
                }
            }
        }
        AlphaMaskShape {
            shape: Arc::clone(shape),
            test_intersection,
            test_intersection_p,
            alpha_mask_texture,
            shadow_alpha_mask_texture,
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
        if self.test_intersection {
            let shape = self.shape.as_ref();
            let t_max = r.t_max.get();
            if let Some((t, si)) = shape.intersect(r) {
                if let Some(mask) = self.alpha_mask_texture.as_ref() {
                    let a = mask.evaluate(&si);
                    if a <= 0.0 {
                        r.t_max.set(t_max);
                        return None;
                    }
                }
                return Some((t, si));
            }
        }
        return None;
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        if self.test_intersection_p {
            let shape = self.shape.as_ref();
            let t_max = r.t_max.get();
            if let Some((_, si)) = shape.intersect(r) {
                r.t_max.set(t_max);
                if let Some(mask) = self.alpha_mask_texture.as_ref() {
                    let a = mask.evaluate(&si);
                    if a <= 0.0 {
                        return false;
                    }
                }
                if let Some(mask) = self.shadow_alpha_mask_texture.as_ref() {
                    let a = mask.evaluate(&si);
                    if a <= 0.0 {
                        return false;
                    }
                }
                return true;
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
