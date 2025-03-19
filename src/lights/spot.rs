use crate::core::error::PbrtError;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::medium::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::transform::*;

use std::sync::Arc;

pub struct SpotLight {
    base: BaseLight,
    p_light: Point3f,
    intensity: Spectrum,
    cos_total_width: Float,
    cos_falloff_start: Float,
}

#[inline]
fn radians(x: Float) -> Float {
    return x * (PI / 180.0);
}

impl SpotLight {
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        intensity: &Spectrum,
        total_width: Float,
        falloff_start: Float,
    ) -> Self {
        let base = BaseLight::new(
            LightFlags::DeltaPosition as u32,
            light_to_world,
            medium_interface,
            1,
        );
        let p_light = light_to_world.transform_point(&Point3f::zero());
        let intensity = *intensity;
        let cos_total_width = Float::cos(radians(total_width));
        let cos_falloff_start = Float::cos(radians(falloff_start));
        SpotLight {
            base,
            p_light,
            intensity,
            cos_total_width,
            cos_falloff_start,
        }
    }

    pub fn falloff(&self, w: &Vector3f) -> Float {
        let wl = self.base.world_to_light.transform_vector(w).normalize();
        let cos_theta = wl.z;
        if cos_theta < self.cos_total_width {
            return 0.0;
        }
        if cos_theta >= self.cos_falloff_start {
            return 1.0;
        }
        // Compute falloff inside spotlight cone
        let delta =
            (cos_theta - self.cos_total_width) / (self.cos_falloff_start - self.cos_total_width);
        return (delta * delta) * (delta * delta);
    }
}

impl Light for SpotLight {
    fn sample_li(
        &self,
        inter: &Interaction,
        _u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, VisibilityTester)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let inter_p = inter.get_p();
        let wi = (self.p_light - inter_p).normalize();
        let pdf: Float = 1.0;
        let intensity = self.intensity
            * (self.falloff(&-wi) / Vector3f::distance_squared(&self.p_light, &inter_p));
        let p = self.p_light;
        let inter_light =
            Interaction::from_light_sample(&p, inter.get_time(), &self.base.medium_interface);
        let vis = VisibilityTester::from((inter, &inter_light));
        return Some((intensity, wi, pdf, vis));
    }

    fn power(&self) -> Spectrum {
        return self.intensity
            * (2.0 * PI * (1.0 - 0.5 * (self.cos_falloff_start - self.cos_total_width)));
    }

    fn pdf_li(&self, _inter: &Interaction, _wi: &Vector3f) -> Float {
        return 0.0;
    }

    fn sample_le(
        &self,
        u1: &Point2f,
        _u2: &Point2f,
        time: Float,
    ) -> Option<(Spectrum, Ray, Normal3f, Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let w = uniform_sample_cone(u1, self.cos_total_width);
        let medium = self.base.medium_interface.get_inside();
        let ray = Ray::from((
            &self.p_light,
            &self.base.light_to_world.transform_vector(&w),
            Float::INFINITY,
            time,
            &medium,
        ));
        let n_light = ray.d;
        let pdf_pos = 1.0;
        let pdf_dir = uniform_cone_pdf(self.cos_total_width);
        let spec = self.intensity * self.falloff(&ray.d);
        return Some((spec, ray, n_light, pdf_pos, pdf_dir));
    }

    fn pdf_le(&self, ray: &Ray, _n_light: &Normal3f) -> Option<(Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightPdf);

        if cos_theta(&self.base.world_to_light.transform_vector(&ray.d)) >= self.cos_total_width {
            let pdf_pos = 0.0;
            let pdf_dir = uniform_cone_pdf(self.cos_total_width);
            assert!(pdf_dir > 0.0);
            return Some((pdf_pos, pdf_dir));
        } else {
            return None;
        }
    }

    fn get_light_flags(&self) -> u32 {
        return self.base.flags;
    }

    fn get_sample_count(&self) -> u32 {
        return self.base.n_samples;
    }
}

pub fn create_spot_light(
    l2w: &Transform,
    medium: &Option<Arc<dyn Medium>>,
    params: &ParamSet,
) -> Result<Arc<dyn Light>, PbrtError> {
    let intensity = params.find_one_spectrum("I", &Spectrum::one());
    let sc = params.find_one_spectrum("scale", &Spectrum::one());
    let coneangle = params.find_one_float("coneangle", 30.0);
    let conedelta = params.find_one_float("conedelta", 5.0); //5.0
    let conedelta = params.find_one_float("conedeltaangle", conedelta);
    let from = params.find_one_point3f("from", &Point3f::new(0.0, 0.0, 0.0));
    let to = params.find_one_point3f("to", &Point3f::new(0.0, 0.0, 1.0));
    let dir = (to - from).normalize();
    let (dir, du, dv) = Vector3f::coordinate_system(&dir);
    let dir_to_z = Transform::from(Matrix4x4::new(
        du.x, du.y, du.z, 0.0, dv.x, dv.y, dv.z, 0., dir.x, dir.y, dir.z, 0.0, 0.0, 0.0, 0.0, 1.0,
    ));
    let light2world =
        *l2w * Transform::translate(from.x, from.y, from.z) * Transform::inverse(&dir_to_z);
    return Ok(Arc::new(SpotLight::new(
        &light2world,
        &MediumInterface::from(medium),
        &(intensity * sc),
        coneangle,
        coneangle - conedelta,
    )));
}
