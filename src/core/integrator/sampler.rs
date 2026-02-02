use super::integrator::*;
use crate::core::base::*;
use crate::core::camera::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::reflection::*;
use crate::core::sampler::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use crate::core::stats::*;

use std::ops::DerefMut;
use std::path::Path;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::RwLock;

use log::*;
use rayon::prelude::*;

thread_local!(static N_CAMERA_RAYS: StatCounter = StatCounter::new("Integrator/Camera rays traced"));

pub trait SamplerIntegrator: Integrator + Sync {
    fn preprocess(&mut self, _scene: &Scene, _sampler: &mut dyn Sampler) {}
    fn li(
        &self,
        ray: &RayDifferential,
        scene: &Scene,
        sampler: &mut dyn Sampler,
        arena: &mut MemoryArena,
        depth: i32,
    ) -> Spectrum;

    fn specular_reflect(
        &self,
        ray: &RayDifferential,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut dyn Sampler,
        arena: &mut MemoryArena,
        depth: i32,
    ) -> Spectrum {
        if let Some(bsdf) = isect.bsdf.as_ref() {
            let bsdf = bsdf.as_ref();
            let wo = isect.wo;
            let u = sampler.get_2d();
            let t = BSDF_REFLECTION | BSDF_SPECULAR;
            if let Some((f, wi, pdf, _tt)) = bsdf.sample_f(&wo, &u, t) {
                assert!(isect.n.length() > 0.0);
                assert!(isect.shading.n.length() > 0.0);

                let ns = isect.shading.n;
                let wi_ns = Vector3f::abs_dot(&wi, &ns);
                if pdf > 0.0 && !f.is_black() && wi_ns != 0.0 {
                    let mut rd = RayDifferential::from(isect.spawn_ray(&wi));
                    if ray.has_differentials {
                        rd.has_differentials = true;
                        rd.rx_origin = isect.p + isect.dpdx;
                        rd.ry_origin = isect.p + isect.dpdy;

                        let dndx =
                            isect.shading.dndu * isect.dudx + isect.shading.dndv * isect.dvdx;
                        let dndy =
                            isect.shading.dndu * isect.dudy + isect.shading.dndv * isect.dvdy;
                        let dwodx = -ray.rx_direction - wo;
                        let dwody = -ray.ry_direction - wo;

                        let d_dndx = Vector3f::dot(&dwodx, &ns) + Vector3f::dot(&wo, &dndx);
                        let d_dndy = Vector3f::dot(&dwody, &ns) + Vector3f::dot(&wo, &dndy);
                        let wo_ns = Vector3f::dot(&wo, &ns);
                        rd.rx_direction = wi - dwodx + 2.0 * (wo_ns * dndx + d_dndx * ns);
                        rd.ry_direction = wi - dwody + 2.0 * (wo_ns * dndy + d_dndy * ns);
                    }
                    return f * self.li(&rd, scene, sampler, arena, depth + 1) * (wi_ns / pdf);
                }
            }
        }
        return Spectrum::zero();
    }

    fn specular_transmit(
        &self,
        ray: &RayDifferential,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut dyn Sampler,
        arena: &mut MemoryArena,
        depth: i32,
    ) -> Spectrum {
        if let Some(bsdf) = isect.bsdf.as_ref() {
            let bsdf = bsdf.as_ref();
            let wo = isect.wo;
            let u = sampler.get_2d();
            let t = BSDF_TRANSMISSION | BSDF_SPECULAR;
            if let Some((f, wi, pdf, _tt)) = bsdf.sample_f(&wo, &u, t) {
                let mut ns = isect.shading.n;
                let mut wi_ns = Vector3f::abs_dot(&wi, &ns);
                let mut wo_ns = Vector3f::dot(&wo, &ns);
                if pdf > 0.0 && !f.is_black() && wi_ns != 0.0 {
                    let mut rd = RayDifferential::from(isect.spawn_ray(&wi));
                    if ray.has_differentials {
                        rd.has_differentials = true;
                        rd.rx_origin = isect.p + isect.dpdx;
                        rd.ry_origin = isect.p + isect.dpdy;

                        let mut dndx =
                            isect.shading.dndu * isect.dudx + isect.shading.dndv * isect.dvdx;
                        let mut dndy =
                            isect.shading.dndu * isect.dudy + isect.shading.dndv * isect.dvdy;
                        let mut eta = 1.0 / bsdf.eta;
                        if Vector3f::dot(&wo, &ns) < 0.0 {
                            eta = 1.0 / eta;
                            ns = -ns;
                            dndx = -dndx;
                            dndy = -dndy;

                            wi_ns = Vector3f::abs_dot(&wi, &ns);
                            wo_ns = Vector3f::dot(&wo, &ns);
                        }

                        let dwodx = -ray.rx_direction - wo;
                        let dwody = -ray.ry_direction - wo;

                        let d_dndx = Vector3f::dot(&dwodx, &ns) + Vector3f::dot(&wo, &dndx);
                        let d_dndy = Vector3f::dot(&dwody, &ns) + Vector3f::dot(&wo, &dndy);

                        let mu = eta * wo_ns - wi_ns;
                        let dmudx = (eta - (eta * eta * wo_ns) / wi_ns) * d_dndx;
                        let dmudy = (eta - (eta * eta * wo_ns) / wi_ns) * d_dndy;

                        rd.rx_direction = wi - eta * dwodx + (mu * dndx + dmudx * ns);
                        rd.ry_direction = wi - eta * dwody + (mu * dndy + dmudy * ns);
                    }
                    return f * self.li(&rd, scene, sampler, arena, depth + 1) * (wi_ns / pdf);
                }
            }
        }
        return Spectrum::zero();
    }

    //fn get_camera(&self) -> Arc<dyn Camera>;

    fn get_sampler(&self) -> Arc<RwLock<dyn Sampler>>;

    fn get_pixel_bounds(&self) -> Bounds2i;
}

fn validate_radiance_result(l: Spectrum, pixel: &Point2i) -> Spectrum {
    if !l.is_valid() {
        error!(
            "Not-a-number radiance value returned for pixel ({}, {}). Setting to black.",
            pixel.x, pixel.y
        );
        return Spectrum::zero();
    }
    if l.y() < -1e-5 {
        error!(
            "Negative luminance value, {}, returned for pixel ({}, {}). Setting to black.",
            l.y(),
            pixel.x,
            pixel.y
        );
        return Spectrum::zero();
    }
    if l.y().is_infinite() {
        error!(
            "Infinite luminance value returned for pixel ({}, {}). Setting to black.",
            pixel.x, pixel.y
        );
        return Spectrum::zero();
    }
    return l;
}

struct SampleIntegratorCore {}

impl SampleIntegratorCore {
    pub fn get_film_tile(film: &Arc<Mutex<ProxyFilm>>, tile_bounds: &Bounds2i) -> FilmTile {
        let film = film.lock().unwrap();
        return film.get_film_tile(tile_bounds);
    }

    pub fn merge_film_tile(film: &Arc<Mutex<ProxyFilm>>, tile: &FilmTile) {
        let mut film = film.lock().unwrap();
        film.merge_film_tile(tile);
        film.update_display(&tile.get_pixel_bounds());
    }

    pub fn get_filename(film: &Arc<Mutex<ProxyFilm>>) -> String {
        let film = film.lock().unwrap();
        let path = film.get_filename();
        let path = Path::new(&path);
        let filename = path.file_name().unwrap();
        let filename = filename.to_str().unwrap();
        return String::from(filename);
    }

    pub fn render_tile(
        integrator: &dyn SamplerIntegrator,
        scene: &Scene,
        camera: &dyn Camera,
        film: &Arc<Mutex<ProxyFilm>>,
        tile_bounds: &Bounds2i,
        sampler: &Arc<RwLock<dyn Sampler>>,
        reporter: &Arc<Mutex<ProgressReporter>>,
    ) {
        let mut arena = MemoryArena::new();
        let x0 = tile_bounds.min.x;
        let x1 = tile_bounds.max.x;
        let y0 = tile_bounds.min.y;
        let y1 = tile_bounds.max.y;

        let mut sampler = sampler.write().unwrap();

        let ray_scale = (1.0 / sampler.get_samples_per_pixel() as Float).sqrt();

        let mut film_tile = Self::get_film_tile(film, tile_bounds);
        for yy in y0..y1 {
            for xx in x0..x1 {
                let pixel = Point2i::new(xx, yy);
                sampler.start_pixel(&pixel);
                loop {
                    // Initialize _CameraSample_ for current sample
                    let camera_sample = sampler.get_camera_sample(&pixel);

                    // Generate camera ray for current sample
                    if let Some((ray_weight, mut ray)) =
                        camera.generate_ray_differential(&camera_sample)
                    {
                        ray.scale_differentials(ray_scale);

                        N_CAMERA_RAYS.with(|c| c.inc());

                        // Evaluate radiance along camera ray
                        let l = integrator.li(&ray, scene, sampler.deref_mut(), &mut arena, 0);
                        let l = validate_radiance_result(l, &pixel);
                        film_tile.add_sample(&camera_sample.p_film, &l, ray_weight);
                    } else {
                        film_tile.add_sample(&camera_sample.p_film, &Spectrum::zero(), 0.0);
                    }
                    arena.reset();

                    if !sampler.start_next_sample() {
                        break;
                    }
                }
            }
        }
        Self::merge_film_tile(film, &film_tile);
        {
            let mut reporter = reporter.lock().unwrap();
            reporter.update(1);
        }
    }

    pub fn render(
        integrator: &dyn SamplerIntegrator,
        scene: &Scene,
        camera: &dyn Camera,
        film: &Arc<RwLock<Film>>,
        sampler: &Arc<RwLock<dyn Sampler>>,
    ) {
        let mut tile_indices = Vec::new();
        {
            let mut film = film.as_ref().write().unwrap();
            let sample_bounds = film.get_sample_bounds();
            let sample_extent = sample_bounds.diagonal();
            const TILE_SIZE: i32 = 16;
            let n_tiles = Point2i::from((
                (sample_extent.x + TILE_SIZE - 1) / TILE_SIZE,
                (sample_extent.y + TILE_SIZE - 1) / TILE_SIZE,
            ));
            tile_indices.reserve((n_tiles.x * n_tiles.y) as usize);
            for y in 0..n_tiles.y {
                for x in 0..n_tiles.x {
                    let x0 = sample_bounds.min.x + x * TILE_SIZE;
                    let x1 = i32::min(x0 + TILE_SIZE, sample_bounds.max.x);
                    let y0 = sample_bounds.min.y + y * TILE_SIZE;
                    let y1 = i32::min(y0 + TILE_SIZE, sample_bounds.max.y);
                    let tile_bounds = Bounds2i::from(((x0, y0), (x1, y1)));

                    let seed = (y * n_tiles.x + x) as u32;
                    let s = sampler.read().unwrap().clone_with_seed(seed);
                    tile_indices.push((tile_bounds, s));
                }
            }

            film.render_start();
        }

        {
            let proxy_film = Arc::new(Mutex::new(ProxyFilm::new(film)));
            let filename = Self::get_filename(&proxy_film);

            let total = tile_indices.len();
            let reporter = Arc::new(Mutex::new(ProgressReporter::new(total, &filename)));

            {
                tile_indices.par_iter().for_each(|(tile_bounds, sampler)| {
                    Self::render_tile(
                        integrator,
                        scene,
                        camera,
                        &proxy_film,
                        tile_bounds,
                        sampler,
                        &reporter,
                    );
                });
            }

            {
                let mut reporter = reporter.lock().unwrap();
                reporter.done();
            }
        }

        {
            let mut film = film.as_ref().write().unwrap();
            film.render_end();
            film.write_image();
        }
    }
}

pub struct BaseSamplerIntegrator {
    pub camera: Arc<dyn Camera>,
    pub sampler: Arc<RwLock<dyn Sampler>>,
    pub pixel_bounds: Bounds2i,
}

impl BaseSamplerIntegrator {
    pub fn new(
        camera: &Arc<dyn Camera>,
        sampler: &Arc<RwLock<dyn Sampler>>,
        pixel_bounds: &Bounds2i,
    ) -> Self {
        BaseSamplerIntegrator {
            camera: Arc::clone(camera),
            sampler: Arc::clone(sampler),
            pixel_bounds: pixel_bounds.clone(),
        }
    }

    pub fn render(integrator: &mut dyn SamplerIntegrator, scene: &Scene) {
        let acamera = integrator.get_camera();
        let camera = acamera.as_ref();
        let sampler = integrator.get_sampler();
        {
            let mut sampler = sampler.write().unwrap();
            integrator.preprocess(scene, sampler.deref_mut());
        }

        let film = camera.get_film();
        SampleIntegratorCore::render(integrator, scene, camera, &film, &sampler);
    }
}
