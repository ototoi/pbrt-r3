use crate::core::pbrt::*;

use std::sync::Arc;

pub struct PointLight {
    base: BaseLight,
    p_light: Point3f,
    intensity: Spectrum,
}

impl PointLight {
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        intensity: &Spectrum,
    ) -> Self {
        let base = BaseLight::new(
            LightFlags::DeltaPosition as u32,
            light_to_world,
            medium_interface,
            1,
        );
        let p_light = light_to_world.transform_point(&Vector3f::zero());
        let intensity = *intensity;
        PointLight {
            base,
            p_light,
            intensity,
        }
    }
}

impl Light for PointLight {
    fn sample_li(
        &self,
        inter: &Interaction,
        _u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, VisibilityTester)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let f = self.intensity / Vector3f::distance_squared(&self.p_light, &inter.get_p());
        let wi = (self.p_light - inter.get_p()).normalize();
        let pdf = 1.0;
        let p = self.p_light;
        let inter_light =
            Interaction::from_light_sample(&p, inter.get_time(), &self.base.medium_interface);
        let vis = VisibilityTester::from((inter, &inter_light));
        return Some((f, wi, pdf, vis));
    }

    fn power(&self) -> Spectrum {
        return self.intensity * (4.0 * PI);
    }

    fn pdf_li(&self, _inter: &Interaction, _wi: &Vector3f) -> Float {
        return 0.0;
    }

    fn sample_le(
        &self,
        u1: &Point2f,
        _: &Point2f,
        time: Float,
    ) -> Option<(Spectrum, Ray, Normal3f, Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let medium = self.base.medium_interface.get_inside();
        let f = self.intensity;
        let ray = Ray::from((
            &self.p_light,
            &uniform_sample_sphere(u1),
            Float::INFINITY,
            time,
            &medium,
        ));
        let n = ray.d;
        let pdf = 1.0;
        let pdf_dir = uniform_sphere_pdf();
        return Some((f, ray, n, pdf, pdf_dir));
    }

    fn pdf_le(&self, _ray: &Ray, _n: &Normal3f) -> Option<(Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightPdf);

        let pdf_pos = 0.0;
        let pdf_dir = uniform_sphere_pdf();
        if pdf_dir > 0.0 {
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

pub fn create_point_light(
    light2world: &Transform,
    medium: &Option<Arc<dyn Medium>>,
    params: &ParamSet,
) -> Result<Arc<dyn Light>, PbrtError> {
    let i = params.find_one_spectrum("I", &Spectrum::one());
    let sc = params.find_one_spectrum("scale", &Spectrum::one());
    let p = params.find_one_point3f("from", &Point3f::zero());
    let l2w = Transform::translate(p.x, p.y, p.z) * (*light2world);
    let mi = MediumInterface::from(medium);
    return Ok(Arc::new(PointLight::new(&l2w, &mi, &(i * sc))));
}
