use crate::core::imageio::*;
use crate::core::pbrt::*;

use std::sync::Arc;

pub struct GonioPhotometricLight {
    base: BaseLight,
    p_light: Point3f,
    intensity: Spectrum,
    mipmap: MIPMap<RGBSpectrum>,
}

impl GonioPhotometricLight {
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        intensity: &Spectrum,
        mipmap: MIPMap<RGBSpectrum>,
    ) -> Self {
        let base = BaseLight::new(
            LightFlags::DeltaPosition as u32,
            light_to_world,
            medium_interface,
            1,
        );
        let p_light = light_to_world.transform_point(&Vector3f::zero());
        let intensity = *intensity;
        GonioPhotometricLight {
            base,
            p_light,
            intensity,
            mipmap,
        }
    }

    fn scale(&self, w: &Vector3f) -> Spectrum {
        let wp = self.base.world_to_light.transform_point(w);
        //swap y,z
        let theta = spherical_theta(&wp);
        let phi = spherical_phi(&wp);
        let st = Point2f::new(phi * INV_2_PI, theta * INV_PI);
        let rgb = self.mipmap.lookup(&st, 0.0);
        let spec = Spectrum::from_rgb(&rgb.to_rgb(), SpectrumType::Illuminant);
        return spec;
    }
}

impl Light for GonioPhotometricLight {
    fn sample_li(
        &self,
        inter: &Interaction,
        _u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, VisibilityTester)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let wi = (self.p_light - inter.get_p()).normalize();
        let f = self.intensity
            * self.scale(&-wi)
            * (1.0 / Vector3f::distance_squared(&self.p_light, &inter.get_p()));
        let pdf = 1.0;
        let inter2 = Interaction::from((
            self.p_light,
            inter.get_time(),
            self.base.medium_interface.clone(),
        ));
        let vis = VisibilityTester::from((inter, &inter2));
        return Some((f, wi, pdf, vis));
    }

    fn power(&self) -> Spectrum {
        let rgb = self.mipmap.lookup(&Point2f::new(0.5, 0.5), 0.5);
        let c = Spectrum::from_rgb(&rgb.to_rgb(), SpectrumType::Illuminant);
        return 4.0 * PI * self.intensity * c;
    }

    fn get_light_flags(&self) -> u32 {
        return self.base.flags;
    }

    fn get_sample_count(&self) -> u32 {
        return self.base.n_samples;
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

        let ray = Ray::from((
            &self.p_light,
            &uniform_sample_sphere(u1),
            Float::INFINITY,
            time,
            &self.base.medium_interface.inside,
        ));
        let n_light = ray.d;
        let f = self.intensity * self.scale(&ray.d);
        let pdf_pos = 1.0;
        let pdf_dir = uniform_sphere_pdf();
        return Some((f, ray, n_light, pdf_pos, pdf_dir));
    }

    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f) -> Option<(Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightPdf);

        let pdf_pos = 0.0;
        let pdf_dir = uniform_sphere_pdf();
        if pdf_dir > 0.0 {
            return Some((pdf_pos, pdf_dir));
        } else {
            return None;
        }
    }
}

unsafe impl Sync for GonioPhotometricLight {}
unsafe impl Send for GonioPhotometricLight {}

fn make_mipmap(path: &str) -> Result<(MIPMap<RGBSpectrum>, Point2i), PbrtError> {
    if let Ok((mut texels, resolution)) = read_image(path) {
        //println!("make_mipmap {}", path);
        let total = (resolution.x * resolution.y) as usize;
        for i in 0..total {
            let cc = texels[i].to_rgb();
            // pbrt-r3: Clamp negative values to zero
            let cc = cc.iter().map(|x| x.max(0.0)).collect::<Vec<Float>>();
            // pbrt-rs: Clamp negative values to zero
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

pub fn create_goniometric_light(
    light2world: &Transform,
    medium: &Option<Arc<dyn Medium>>,
    params: &ParamSet,
) -> Result<Arc<dyn Light>, PbrtError> {
    let l = params.find_one_spectrum("L", &Spectrum::one());
    let sc = params.find_one_spectrum("scale", &Spectrum::one());
    let texmap = params.find_one_filename("mapname", "");
    let mi = MediumInterface::from(medium);
    let (mipmap, _resolution) = make_mipmap(&texmap)?;
    return Ok(Arc::new(GonioPhotometricLight::new(
        light2world,
        &mi,
        &(l * sc),
        mipmap,
    )));
}
