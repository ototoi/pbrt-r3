use vector3::Vector3f;

use crate::core::imageio::*;
use crate::core::pbrt::*;

use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, Default, Copy, Clone)]
struct WorldBound {
    pub center: Point3f,
    pub radius: Float,
}

pub struct InfiniteAreaLight {
    base: BaseLight,
    lmap: MIPMap<RGBSpectrum>,
    bound: RwLock<WorldBound>,
    distribution: Distribution2D,
}

impl InfiniteAreaLight {
    pub fn new(
        light_to_world: &Transform,
        n_samples: u32,
        lmap: MIPMap<RGBSpectrum>,
        distribution: Distribution2D,
    ) -> Self {
        let base = BaseLight::new(
            LightFlags::Infinite as u32,
            light_to_world,
            &MediumInterface::new(),
            n_samples,
        );
        InfiniteAreaLight {
            base,
            lmap,
            bound: RwLock::new(WorldBound {
                center: Point3f::new(0.0, 0.0, 0.0),
                radius: Float::INFINITY,
            }),
            distribution,
        }
    }

    fn get_bound(&self) -> WorldBound {
        let bound = self.bound.read().unwrap();
        return *bound;
    }
}

fn make_mipmap(path: &str, l: &Spectrum) -> Result<MIPMap<RGBSpectrum>, PbrtError> {
    if let Ok((mut texels, resolution)) = read_image(path) {
        //println!("make_mipmap {}", path);
        let total = (resolution.x * resolution.y) as usize;
        let c = RGBSpectrum::from(l.to_rgb());
        for i in 0..total {
            let cc = texels[i].to_rgb();
            let cc = cc.iter().map(|x| x.max(0.0)).collect::<Vec<Float>>();// pbrt-r3: Clamp negative values to zero
            texels[i] = RGBSpectrum::from(cc) * c;

            assert!(texels[i].y() >= 0.0);
        }
        return create_spectrum_mipmap(&resolution, texels.as_ref());
    } else {
        let c = RGBSpectrum::from(l.to_rgb());
        let texels: Vec<RGBSpectrum> = vec![c];
        let resolution = Point2i::new(1, 1);
        return create_spectrum_mipmap(&resolution, texels.as_ref());
    }
}

fn make_distribution(lmap: &MIPMap<RGBSpectrum>) -> Result<Distribution2D, PbrtError> {
    // Initialize sampling PDFs for infinite area light

    // Compute scalar-valued image _img_ from environment map
    let width = lmap.width() * 2;
    let height = lmap.height() * 2;
    let mut img = vec![0.0; width * height];
    let fwidth = 0.5 / usize::min(width, height) as Float;
    for v in 0..height {
        let vp = (v as Float + 0.5) / height as Float;
        assert!(vp >= 0.0 && vp <= 1.0);
        let sin_theta = Float::sin(PI * vp);
        assert!(sin_theta >= 0.0);
        for u in 0..width {
            let up = (u as Float + 0.5) / width as Float;
            let y = lmap.lookup(&Point2f::new(up, vp), fwidth).y();
            let y = y.max(0.0); // pbrt-r3: Clamp negative values to zero
            let c = y * sin_theta;
            img[v * width + u] = c;
        }
    }
    let d = Distribution2D::new(&img, width, height);
    return Ok(d);
}

impl Light for InfiniteAreaLight {
    fn preprocess(&self, scene: &Scene) {
        //scene.WorldBound().BoundingSphere(&worldCenter, &worldRadius);
        let (center, radius) = scene.world_bound.bounding_sphere();
        let mut bound = self.bound.write().unwrap();
        bound.center = center;
        bound.radius = radius;
    }

    fn power(&self) -> Spectrum {
        let bound = self.bound.read().unwrap();
        let world_radius = bound.radius;
        let rgb = self.lmap.lookup(&Point2f::new(0.5, 0.5), 0.5);
        return Spectrum::from_rgb(&rgb.to_rgb(), SpectrumType::Illuminant)
            * (PI * world_radius * world_radius);
    }

    fn le(&self, ray: &RayDifferential) -> Spectrum {
        let w = self
            .base
            .world_to_light
            .transform_vector(&ray.ray.d)
            .normalize();
        let s = spherical_phi(&w) * INV_2_PI;
        let t = spherical_theta(&w) * INV_PI;
        let rgb = self.lmap.lookup(&Point2f::new(s, t), 0.0);
        let spc = Spectrum::from_rgb(&rgb.to_rgb(), SpectrumType::Illuminant);
        return spc;
    }

    fn sample_li(
        &self,
        inter: &Interaction,
        u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, VisibilityTester)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        // Find $(u,v)$ sample coordinates in infinite light texture
        let (uv, map_pdf) = self.distribution.sample_continuous(u);
        if map_pdf <= 0.0 {
            return None;
        }
        let bound = self.get_bound();

        // Convert infinite light sample point to direction
        let theta = uv[1] * PI;
        let phi = uv[0] * 2.0 * PI;
        let cos_theta = Float::cos(theta);
        let sin_theta = Float::sin(theta);
        let sin_theta = Float::clamp(sin_theta, 0.0, 1.0); // pbrt-r3: Clamp negative values to zero
        let cos_phi = Float::cos(phi);
        let sin_phi = Float::sin(phi);
        let wi = Vector3f::new(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
        let wi = self.base.light_to_world.transform_vector(&wi);
        let pdf = if sin_theta == 0.0 {
            0.0
        } else {
            map_pdf / (2.0 * PI * PI * sin_theta)
        };
        let world_radius = bound.radius;
        let p = inter.get_p() + wi * (2.0 * world_radius);
        let inter2 = Interaction::from((p, inter.get_time(), self.base.medium_interface.clone()));
        let vis = VisibilityTester::from((inter, &inter2));
        let rgb = self.lmap.lookup(&uv, 0.0);
        let spc = Spectrum::from_rgb(&rgb.to_rgb(), SpectrumType::Illuminant);
        return Some((spc, wi, pdf, vis));
    }

    fn pdf_li(&self, _inter: &Interaction, w: &Vector3f) -> Float {
        let _p = ProfilePhase::new(Prof::LightPdf);

        let wi = self.base.world_to_light.transform_vector(w);
        let theta = spherical_theta(&wi);
        let phi = spherical_phi(&wi);
        let sin_theta = Float::sin(theta);
        let sin_theta = Float::clamp(sin_theta, 0.0, 1.0); // pbrt-r3: Clamp negative values to zero
        assert!(sin_theta >= 0.0);
        if sin_theta == 0.0 {
            return 0.0;
        } else {
            return self
                .distribution
                .pdf(&Point2f::new(phi * INV_2_PI, theta * INV_PI))
                / (2.0 * PI * PI * sin_theta);
        }
    }

    fn sample_le(
        &self,
        u1: &Point2f,
        u2: &Point2f,
        time: Float,
    ) -> Option<(Spectrum, Ray, Normal3f, Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        // Compute direction for infinite light sample ray
        let u = *u1;

        // Find $(u,v)$ sample coordinates in infinite light texture
        let (uv, map_pdf) = self.distribution.sample_continuous(&u);
        if map_pdf <= 0.0 {
            return None;
        }
        let theta = uv[1] * PI;
        let phi = uv[0] * 2.0 * PI;
        let cos_theta = Float::cos(theta);
        let sin_theta = Float::sin(theta);
        let sin_theta = Float::clamp(sin_theta, 0.0, 1.0); // pbrt-r3: Clamp negative values to zero
        let sin_phi = Float::sin(phi);
        let cos_phi = Float::cos(phi);
        let d = -self.base.light_to_world.transform_vector(&Vector3f::new(
            sin_theta * cos_phi,
            sin_theta * sin_phi,
            cos_theta,
        ));
        let n_light = d;
        let (v1, v2) = coordinate_system(&-d);
        let cd = concentric_sample_disk(u2);

        let bound = self.get_bound();
        let world_center: Vector3<f32> = bound.center;
        let world_radius = bound.radius;

        let p_disk = world_center + world_radius * (cd.x * v1 + cd.y * v2);
        let o = p_disk + world_radius * -d;
        let ray = Ray::new(&o, &d, Float::INFINITY, time);

        // Compute _InfiniteAreaLight_ ray PDFs
        let pdf_dir = if sin_theta == 0.0 {
            0.0
        } else {
            map_pdf / (2.0 * PI * PI * sin_theta)
        };
        let pdf_pos = 1.0 / (PI * world_radius * world_radius);
        let rgb = self.lmap.lookup(&uv, 0.0);
        let spc = Spectrum::from_rgb(&rgb.to_rgb(), SpectrumType::Illuminant);
        return Some((spc, ray, n_light, pdf_pos, pdf_dir));
    }

    fn pdf_le(&self, ray: &Ray, _n_light: &Normal3f) -> Option<(Float, Float)> {
        let _ = ProfilePhase::new(Prof::LightPdf);

        let d = -self.base.world_to_light.transform_vector(&ray.d);
        let theta = spherical_theta(&d);
        let phi = spherical_phi(&d);
        let sin_theta = Float::sin(theta);
        let sin_theta = Float::clamp(sin_theta, 0.0, 1.0); // pbrt-r3: Clamp negative values to zero

        let bound = self.get_bound();
        //let world_center = bound.center;
        let world_radius = bound.radius;

        let uv = Point2f::new(phi * INV_2_PI, theta * INV_PI);
        let map_pdf = self.distribution.pdf(&uv);
        let map_pdf = Float::max(0.0, map_pdf);
        let pdf_dir = if sin_theta == 0.0 {
            0.0
        } else {
            map_pdf / (2.0 * PI * PI * sin_theta)
        };
        let pdf_pos = 1.0 / (PI * world_radius * world_radius);
        return Some((pdf_pos, pdf_dir));
    }

    fn is_infinite(&self) -> bool {
        return true;
    }

    fn get_light_flags(&self) -> u32 {
        return self.base.flags;
    }

    fn get_sample_count(&self) -> u32 {
        return self.base.n_samples;
    }
}

unsafe impl Send for InfiniteAreaLight {}
unsafe impl Sync for InfiniteAreaLight {}

pub fn create_infinite_light(
    light2world: &Transform,
    params: &ParamSet,
) -> Result<Arc<dyn Light>, PbrtError> {
    let l = params.find_one_spectrum("L", &Spectrum::one());
    let sc = params.find_one_spectrum("scale", &Spectrum::one());
    let texmap = params.find_one_filename("mapname", "");
    let n_samples = params.find_one_int("samples", params.find_one_int("nsamples", 1)) as u32;
    let lmap = make_mipmap(&texmap, &(l * sc))?;
    let distribution = make_distribution(&lmap)?;
    return Ok(Arc::new(InfiniteAreaLight::new(
        light2world,
        n_samples,
        lmap,
        distribution,
    )));
}
