use crate::core::prelude::*;

use std::sync::Arc;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

#[derive(Debug, Default, Copy, Clone)]
struct WorldBound {
    pub center: Point3f,
    pub radius: Float,
}

pub struct InfiniteAreaLight {
    base: BaseLight,
    lmap: MIPMap<RGBSpectrum>,
    bound: OnceLock<WorldBound>,
    distribution: Distribution2D,
}

static INFINITE_LIGHT_TIMING_ENABLED: OnceLock<bool> = OnceLock::new();
static INFINITE_SAMPLE_LI_CALLS: AtomicU64 = AtomicU64::new(0);
static INFINITE_SAMPLE_LI_TOTAL_NS: AtomicU64 = AtomicU64::new(0);
static INFINITE_SAMPLE_LI_DISTR_NS: AtomicU64 = AtomicU64::new(0);
static INFINITE_SAMPLE_LI_DIR_PDF_NS: AtomicU64 = AtomicU64::new(0);
static INFINITE_SAMPLE_LI_LMAP_NS: AtomicU64 = AtomicU64::new(0);
static INFINITE_SAMPLE_LI_VIS_NS: AtomicU64 = AtomicU64::new(0);

#[inline]
fn infinite_light_timing_enabled() -> bool {
    *INFINITE_LIGHT_TIMING_ENABLED
        .get_or_init(|| std::env::var_os("PBRT_R3_BDPT_TIMING").is_some())
}

pub fn reset_infinite_light_sample_timing_counters() {
    INFINITE_SAMPLE_LI_CALLS.store(0, Ordering::Relaxed);
    INFINITE_SAMPLE_LI_TOTAL_NS.store(0, Ordering::Relaxed);
    INFINITE_SAMPLE_LI_DISTR_NS.store(0, Ordering::Relaxed);
    INFINITE_SAMPLE_LI_DIR_PDF_NS.store(0, Ordering::Relaxed);
    INFINITE_SAMPLE_LI_LMAP_NS.store(0, Ordering::Relaxed);
    INFINITE_SAMPLE_LI_VIS_NS.store(0, Ordering::Relaxed);
}

pub fn report_infinite_light_sample_timing_counters(label: &str) {
    let calls = INFINITE_SAMPLE_LI_CALLS.load(Ordering::Relaxed);
    if calls == 0 {
        return;
    }
    let total = INFINITE_SAMPLE_LI_TOTAL_NS.load(Ordering::Relaxed) as f64 * 1e-9;
    let t_distr = INFINITE_SAMPLE_LI_DISTR_NS.load(Ordering::Relaxed) as f64 * 1e-9;
    let t_dir = INFINITE_SAMPLE_LI_DIR_PDF_NS.load(Ordering::Relaxed) as f64 * 1e-9;
    let t_lmap = INFINITE_SAMPLE_LI_LMAP_NS.load(Ordering::Relaxed) as f64 * 1e-9;
    let t_vis = INFINITE_SAMPLE_LI_VIS_NS.load(Ordering::Relaxed) as f64 * 1e-9;
    let denom = if total > 0.0 { total } else { 1.0 };
    eprintln!(
        "[infinite sample_li timing:{}] calls={} total={:.3}s distr={:.3}s ({:.1}%) dir_pdf={:.3}s ({:.1}%) lmap={:.3}s ({:.1}%) vis={:.3}s ({:.1}%)",
        label,
        calls,
        total,
        t_distr,
        100.0 * t_distr / denom,
        t_dir,
        100.0 * t_dir / denom,
        t_lmap,
        100.0 * t_lmap / denom,
        t_vis,
        100.0 * t_vis / denom
    );
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
            bound: OnceLock::new(),
            distribution,
        }
    }

    fn get_bound(&self) -> WorldBound {
        self.bound.get().copied().unwrap_or(WorldBound {
            center: Point3f::new(0.0, 0.0, 0.0),
            radius: Float::INFINITY,
        })
    }
}

fn make_mipmap(path: &str, l: &Spectrum) -> Result<MIPMap<RGBSpectrum>, PbrtError> {
    if !path.is_empty() {
        let (mut texels, resolution) = read_image(path)?;
        let total = (resolution.x * resolution.y) as usize;
        let c = RGBSpectrum::from(l.to_rgb());
        for i in 0..total {
            let cc = texels[i].to_rgb();
            let cc = cc.iter().map(|x| x.max(0.0)).collect::<Vec<Float>>(); // pbrt-r3: Clamp negative values to zero
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
        let _ = self.bound.set(WorldBound { center, radius });
    }

    fn power(&self) -> Spectrum {
        let world_radius = self.get_bound().radius;
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
        let timing = infinite_light_timing_enabled();
        let t_total = if timing { Some(Instant::now()) } else { None };
        if timing {
            INFINITE_SAMPLE_LI_CALLS.fetch_add(1, Ordering::Relaxed);
        }

        // Find $(u,v)$ sample coordinates in infinite light texture
        let t_distr = if timing { Some(Instant::now()) } else { None };
        let (uv, map_pdf) = self.distribution.sample_continuous(u);
        if let Some(t0) = t_distr {
            INFINITE_SAMPLE_LI_DISTR_NS.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        if map_pdf <= 0.0 {
            return None;
        }
        let bound = self.get_bound();

        // Convert infinite light sample point to direction
        let t_dir = if timing { Some(Instant::now()) } else { None };
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
        if let Some(t0) = t_dir {
            INFINITE_SAMPLE_LI_DIR_PDF_NS.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        let world_radius = bound.radius;
        let p = inter.get_p() + wi * (2.0 * world_radius);
        let t_vis = if timing { Some(Instant::now()) } else { None };
        let inter_light =
            Interaction::from_light_sample(&p, inter.get_time(), &self.base.medium_interface);
        let vis = VisibilityTester::from((inter, &inter_light));
        if let Some(t0) = t_vis {
            INFINITE_SAMPLE_LI_VIS_NS.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        let t_lmap = if timing { Some(Instant::now()) } else { None };
        let rgb = self.lmap.lookup(&uv, 0.0);
        if let Some(t0) = t_lmap {
            INFINITE_SAMPLE_LI_LMAP_NS.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
        let spc = Spectrum::from_rgb(&rgb.to_rgb(), SpectrumType::Illuminant);
        if let Some(t0) = t_total {
            INFINITE_SAMPLE_LI_TOTAL_NS.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
        }
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
        let world_center = bound.center;
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
    let mut n_samples = params.find_one_int("samples", params.find_one_int("nsamples", 1)) as u32;
    let lmap = make_mipmap(&texmap, &(l * sc))?;
    let distribution = make_distribution(&lmap)?;
    {
        let options = PbrtOptions::get();
        if options.quick_render {
            n_samples = (n_samples / 4).max(1);
        }
    }
    return Ok(Arc::new(InfiniteAreaLight::new(
        light2world,
        n_samples,
        lmap,
        distribution,
    )));
}
