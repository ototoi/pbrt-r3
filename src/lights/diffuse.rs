use crate::core::error::PbrtError;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::medium::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::sampling::*;
use crate::core::shape::*;
use crate::core::spectrum::*;

use std::sync::Arc;

pub struct DiffuseAreaLight {
    base: BaseLight,
    lemit: Spectrum,
    shape: Arc<dyn Shape>,
    two_sided: bool,
    area: Float,
}

impl DiffuseAreaLight {
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        le: &Spectrum,
        n_samples: u32,
        shape: &Arc<dyn Shape>,
        two_sided: bool,
    ) -> Self {
        let base = BaseLight::new(
            LightFlags::Area as u32,
            light_to_world,
            medium_interface,
            n_samples,
        );
        let lemit = *le;
        let shape = shape.clone();
        let area = shape.area();
        DiffuseAreaLight {
            base,
            lemit,
            shape,
            two_sided,
            area,
        }
    }

    fn sample_wo(&self, _u1: &Vector2f, u2: &Vector2f) -> (Vector3f, Float) {
        if self.two_sided {
            let mut u = *u2;
            // Choose a side to sample and then remap u[0] to [0,1] before
            // applying cosine-weighted hemisphere sampling for the chosen side.
            let mut w;
            if u[0] < 0.5 {
                u.x = Float::min(u[0] * 2.0, ONE_MINUS_EPSILON);
                w = cosine_sample_hemisphere(&u);
            } else {
                u.x = Float::min((u[0] - 0.5) * 2.0, ONE_MINUS_EPSILON);
                w = cosine_sample_hemisphere(&u);
                w.z *= -1.0;
            };
            let pdf = 0.5 * cosine_hemisphere_pdf(Float::abs(w.z));
            return (w, pdf);
        } else {
            let w = cosine_sample_hemisphere(u2);
            let pdf = cosine_hemisphere_pdf(w.z);
            return (w, pdf);
        }
    }
}

impl Light for DiffuseAreaLight {
    fn power(&self) -> Spectrum {
        let n = if self.two_sided { 2.0 } else { 1.0 };
        return self.lemit * (n * self.area * PI);
    }

    fn sample_li(
        &self,
        inter: &Interaction,
        u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, VisibilityTester)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let shape = self.shape.as_ref();
        let (mut p_shape, pdf) = shape.sample_from(inter, u)?;
        p_shape.set_medium_interface(&self.base.medium_interface);
        if pdf <= 0.0 || (p_shape.get_p() - inter.get_p()).length_squared() <= 0.0 {
            return None;
        }
        let wi = (p_shape.get_p() - inter.get_p()).normalize();
        let vis = VisibilityTester::from((inter, &p_shape));
        let spec = self.l(&p_shape, &-wi);
        return Some((spec, wi, pdf, vis));
    }

    fn pdf_li(&self, inter: &Interaction, wi: &Vector3f) -> Float {
        let _p = ProfilePhase::new(Prof::LightPdf);

        let shape = self.shape.as_ref();
        return shape.pdf_from(inter, wi);
    }

    fn sample_le(
        &self,
        u1: &Point2f,
        u2: &Point2f,
        _time: Float,
    ) -> Option<(Spectrum, Ray, Normal3f, Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let shape = self.shape.as_ref();
        let (mut p_shape, pdf_pos) = shape.sample(u1)?;
        p_shape.set_medium_interface(&self.base.medium_interface);
        let n = p_shape.get_n();

        // Sample a cosine-weighted outgoing direction _w_ for area light
        let (w, pdf_dir) = self.sample_wo(u1, u2);
        let (v1, v2) = coordinate_system(&n);
        let w = w.x * v1 + w.y * v2 + w.z * n;
        let ray = p_shape.spawn_ray(&w);
        let spec = self.l(&p_shape, &w);
        return Some((spec, ray, n, pdf_pos, pdf_dir));
    }

    fn pdf_le(&self, ray: &Ray, n: &Normal3f) -> Option<(Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightPdf);

        let shape = self.shape.as_ref();
        let n = *n;
        let wo = Vector3f::from(n);
        let medium_interface = self.base.medium_interface.clone();
        let it = Interaction::from((ray.o, n, Vector3f::zero(), wo, ray.time, medium_interface));
        let pdf_pos = shape.pdf(&it);
        let pdf_dir = if self.two_sided {
            0.5 * cosine_hemisphere_pdf(Vector3::abs_dot(&n, &ray.d))
        } else {
            // pbrt-r3: Clamp to 0.0 in case of numerical error.
            //cosine_hemisphere_pdf(Vector3f::dot(n, &ray.d))
            cosine_hemisphere_pdf(Vector3f::dot(&n, &ray.d)).max(0.0)
            // pbrt-r3:
        };
        if pdf_pos > 0.0 || pdf_dir > 0.0 {
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

    fn as_area_light(&self) -> Option<&dyn AreaLight> {
        return Some(self);
    }
}

impl AreaLight for DiffuseAreaLight {
    fn l(&self, inter: &Interaction, w: &Vector3f) -> Spectrum {
        if self.two_sided || Vector3f::dot(&inter.get_n(), w) > 0.0 {
            return self.lemit;
        } else {
            return Spectrum::zero();
        }
    }
}

unsafe impl Sync for DiffuseAreaLight {}
unsafe impl Send for DiffuseAreaLight {}

pub fn create_diffuse_area_light(
    light2world: &Transform,
    medium: &Option<Arc<dyn Medium>>,
    params: &ParamSet,
    shape: &Arc<dyn Shape>,
) -> Result<Arc<dyn Light>, PbrtError> {
    let l = params.find_one_spectrum("L", &Spectrum::one());
    let sc = params.find_one_spectrum("scale", &Spectrum::one());
    let mi = MediumInterface::from(medium);
    let mut n_samples = params.find_one_int("nsamples", 1) as u32;
    let two_sided = params.find_one_bool("twosided", false);
    {
        let options = PbrtOptions::get();
        if options.quick_render {
            n_samples = (n_samples / 4).max(1);
        }
    }
    return Ok(Arc::new(DiffuseAreaLight::new(
        light2world,
        &mi,
        &(l * sc),
        n_samples,
        shape,
        two_sided,
    )));
}
