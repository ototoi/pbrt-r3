use crate::core::error::PbrtError;
use crate::core::geometry::*;
use crate::core::imageio::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::medium::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::sampling::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use crate::core::transform::*;

use std::sync::Arc;

pub struct ProjectionLight {
    base: BaseLight,
    projection_map: MIPMap<RGBSpectrum>,
    p_light: Point3f,
    intensity: Spectrum,
    light_projection: Transform,
    hitcher: Float,
    #[allow(dead_code)]
    yon: Float,
    screen_bounds: Bounds2f,
    cos_total_width: Float,
}

impl ProjectionLight {
    pub fn new(
        light_to_world: &Transform,
        medium: &MediumInterface,
        intensity: &Spectrum,
        projection_map: MIPMap<RGBSpectrum>,
        resolution: &Point2i,
        fov: Float,
    ) -> Self {
        let base = BaseLight::new(LightFlags::DeltaPosition as u32, light_to_world, medium, 1);
        let p_light = light_to_world.transform_point(&Point3f::new(0.0, 0.0, 0.0));
        let intensity = *intensity;
        let aspect = resolution.x as Float / resolution.y as Float;
        let screen_bounds = if aspect > 1.0 {
            Bounds2f::new(&Point2f::new(-aspect, -1.0), &Point2f::new(aspect, 1.0))
        } else {
            Bounds2f::new(
                &Point2f::new(-1.0, -1.0 / aspect),
                &Point2f::new(1.0, 1.0 / aspect),
            )
        };
        let hitcher = 1e-3;
        let yon = 1e30;
        let light_projection = Transform::perspective(fov, hitcher, yon);
        let screen_to_light = light_projection.inverse();
        let p_corner = Point3f::new(screen_bounds.max.x, screen_bounds.max.y, 0.0);
        let w_corner = screen_to_light.transform_point(&p_corner).normalize();
        let cos_total_width = w_corner.z;
        return ProjectionLight {
            base,
            projection_map,
            p_light,
            intensity,
            light_projection,
            hitcher,
            yon,
            screen_bounds,
            cos_total_width,
        };
    }

    fn projection(&self, w: &Vector3f) -> Spectrum {
        let wl = self.base.world_to_light.transform_vector(w);
        if wl.z < self.hitcher {
            return Spectrum::from(0.0);
        }

        // Project point onto projection plane and compute light
        let p = self
            .light_projection
            .transform_point(&Point3f::new(wl.x, wl.y, wl.z));
        if !self.screen_bounds.inside(&Point2f::new(p.x, p.y)) {
            return Spectrum::from(0.0);
        }
        let st = self.screen_bounds.offset(&Point2f::new(p.x, p.y));
        let rgb = self.projection_map.lookup(&st, 0.0);
        let spec = Spectrum::from_rgb(&rgb.to_rgb(), SpectrumType::Illuminant);
        return spec;
    }
}

impl Light for ProjectionLight {
    fn sample_li(
        &self,
        inter: &Interaction,
        _u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, VisibilityTester)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let wi = (self.p_light - inter.get_p()).normalize();
        let pdf = 1.0;

        let p = self.p_light;
        let inter_light =
            Interaction::from_light_sample(&p, inter.get_time(), &self.base.medium_interface);
        let vis = VisibilityTester::from((inter.clone(), inter_light));
        let spec = self.intensity * self.projection(&-wi)
            / Point3f::distance_squared(&self.p_light, &inter.get_p());
        return Some((spec, wi, pdf, vis));
    }

    fn power(&self) -> Spectrum {
        let rgb = self.projection_map.lookup(&Point2f::new(0.5, 0.5), 0.5);
        let spec = Spectrum::from_rgb(&rgb.to_rgb(), SpectrumType::Illuminant);
        return spec * self.intensity * (2.0 * PI * (1.0 - self.cos_total_width));
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

        let v = uniform_sample_cone(u1, self.cos_total_width);
        let ray = Ray::new(
            &self.p_light,
            &self.base.light_to_world.transform_vector(&v),
            Float::INFINITY,
            time,
        );
        let n_light = Normal3f::from(ray.d);
        let pdf_pos = 1.0;
        let pdf_dir = uniform_cone_pdf(self.cos_total_width);
        let spec = self.intensity * self.projection(&ray.d);
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

unsafe impl Sync for ProjectionLight {}
unsafe impl Send for ProjectionLight {}

fn make_mipmap(path: &str) -> Result<(MIPMap<RGBSpectrum>, Point2i), PbrtError> {
    if !path.is_empty() {
        let (mut texels, resolution) = read_image(path)?;
        let total = (resolution.x * resolution.y) as usize;
        for i in 0..total {
            let cc = texels[i].to_rgb();
            let cc = cc.iter().map(|x| x.max(0.0)).collect::<Vec<Float>>(); // pbrt-r3: Clamp negative values to zero
            texels[i] = RGBSpectrum::from(cc);
        }
        let mipmap = create_spectrum_mipmap(&resolution, texels.as_ref())?;
        return Ok((mipmap, resolution));
    } else {
        let c = RGBSpectrum::one();
        let texels: Vec<RGBSpectrum> = vec![c];
        let resolution = Point2i::new(1, 1);
        let mipmap = create_spectrum_mipmap(&resolution, texels.as_ref())?;
        return Ok((mipmap, resolution));
    }
}

pub fn create_projection_light(
    light2world: &Transform,
    medium: &Option<Arc<dyn Medium>>,
    params: &ParamSet,
) -> Result<Arc<dyn Light>, PbrtError> {
    let intensity = params.find_one_spectrum("I", &Spectrum::from(1.0));
    let sc = params.find_one_spectrum("scale", &Spectrum::from(1.0));
    let fov = params.find_one_float("fov", 45.0);
    let texname = params.find_one_filename("mapname", "");
    let (map, resolution) = make_mipmap(&texname)?;
    let mi = MediumInterface::from(medium);
    return Ok(Arc::new(ProjectionLight::new(
        light2world,
        &mi,
        &(intensity * sc),
        map,
        &resolution,
        fov,
    )));
}
