use crate::core::prelude::*;

use std::ops::Deref;
use std::sync::Arc;
use std::sync::RwLock;

thread_local!(static PATHS: StatPercent = StatPercent::new("Integrator/Zero-radiance paths"));
thread_local!(static PATH_LENGTH: StatIntDistribution = StatIntDistribution::new("Integrator/Path length"));

pub struct PathIntegrator {
    base: BaseSamplerIntegrator,
    light_distribution: Option<Arc<dyn LightDistribution>>,
    max_depth: i32,
    rr_threshold: Float,
    light_sample_sterategy: String,
}

impl PathIntegrator {
    pub fn new(
        max_depth: i32,
        camera: &Arc<dyn Camera>,
        sampler: &Arc<RwLock<dyn Sampler>>,
        pixel_bounds: &Bounds2i,
        rr_threshold: Float,
        light_sample_sterategy: &str,
    ) -> Self {
        PathIntegrator {
            base: BaseSamplerIntegrator::new(camera, sampler, pixel_bounds),
            light_distribution: None,
            max_depth,
            rr_threshold,
            light_sample_sterategy: String::from(light_sample_sterategy),
        }
    }
}

impl Integrator for PathIntegrator {
    fn render(&mut self, scene: &Scene) {
        BaseSamplerIntegrator::render(self, scene);
    }

    fn get_camera(&self) -> Arc<dyn Camera> {
        return self.base.camera.clone();
    }
}

unsafe impl Sync for PathIntegrator {}

impl SamplerIntegrator for PathIntegrator {
    fn preprocess(&mut self, scene: &Scene, _sampler: &mut dyn Sampler) {
        match create_light_sample_distribution(&self.light_sample_sterategy, scene) {
            Ok(distrib) => {
                self.light_distribution = Some(distrib);
            }
            Err(e) => {
                println!("{:?}", e);
            }
        }
    }

    fn li(
        &self,
        r: &RayDifferential,
        scene: &Scene,
        sampler: &mut dyn Sampler,
        arena: &mut MemoryArena,
        _depth: i32,
    ) -> Spectrum {
        let _p = ProfilePhase::new(Prof::SamplerIntegratorLi);

        let light_distribution = self.light_distribution.as_ref();
        if light_distribution.is_none() {
            return Spectrum::zero();
        }
        let light_distribution = light_distribution.unwrap();

        let max_depth = self.max_depth;
        let mut ray = r.clone();
        let mut eta_scale = 1.0;

        let mut l = Spectrum::zero();
        let mut beta = Spectrum::one();
        let mut specular_bounce = false;
        let mut bounces = 0;
        loop {
            let found_intersection = scene.intersect(&ray.ray);
            if bounces == 0 || specular_bounce {
                if let Some(isect) = found_intersection.as_ref() {
                    let w = -(ray.ray.d);
                    l += beta * isect.le(&w);
                } else {
                    for light in scene.infinite_lights.iter() {
                        l += beta * light.le(&ray);
                    }
                }
                assert!(beta.y() >= 0.0);
                assert!(l.y() >= -1e-5);
            }

            // Terminate path if ray escaped or _maxDepth_ was reached
            if found_intersection.is_none() || bounces >= max_depth {
                break;
            }

            let mut isect = found_intersection.unwrap();
            isect.compute_scattering_functions(&ray, arena, TransportMode::Radiance, true);

            let tisect = Interaction::from(&isect);
            //let isect = tisect.as_surface_interaction().unwrap();

            if isect.bsdf.is_none() {
                ray = isect.spawn_ray(&ray.ray.d).into();
                continue;
            }

            let bsdf = isect.bsdf.as_ref().unwrap();
            let distrib = light_distribution.lookup(&isect.p);

            // assert!(beta.max_component_value() <= 1.0);
            // Sample illumination from lights to find path contribution.
            // (But skip this for perfectly specular BSDFs.)
            if bsdf.num_components(BSDF_ALL & !(BSDF_SPECULAR as u32)) > 0 {
                PATHS.with(|stat| stat.add_denom(1)); //totalPaths

                let ld = beta
                    * uniform_sample_one_light(
                        &tisect,
                        scene,
                        arena,
                        sampler,
                        false,
                        Some(distrib.deref()),
                    );
                if ld.is_black() {
                    PATHS.with(|stat| stat.add_num(1)); //zeroRadiancePaths
                }
                assert!(beta.y() >= 0.0);
                //assert!(ld.y() >= -1e-5); //Should be true
                l += ld;
            }

            // Sample BSDF to get new path direction
            let wo = -ray.ray.d;
            if let Some((f, wi, pdf, flags)) = bsdf.sample_f(&wo, &sampler.get_2d(), BSDF_ALL) {
                if f.is_black() || pdf == 0.0 {
                    break;
                }
                assert!(isect.n.length() > 0.0);
                assert!(isect.shading.n.length() > 0.0);
                assert!(Float::is_finite(beta.y()));

                beta *= f * (Vector3f::abs_dot(&wi, &isect.shading.n) / pdf);
                assert!(Float::is_finite(beta.y()));

                specular_bounce = (flags & BSDF_SPECULAR) != 0;
                let flag_s = (flags & BSDF_SPECULAR) != 0;
                let flag_t = (flags & BSDF_TRANSMISSION) != 0;

                if flag_s && flag_t {
                    let eta = bsdf.eta;
                    // Update the term that tracks radiance scaling for refraction
                    // depending on whether the ray is entering or leaving the
                    // medium.
                    eta_scale *= if Vector3f::dot(&wo, &isect.n) > 0.0 {
                        eta * eta
                    } else {
                        1.0 / (eta * eta)
                    }
                }
                ray = isect.spawn_ray(&wi).into();

                // Account for subsurface scattering, if applicable
                if let Some(bssrdf) = isect.bssrdf.as_ref() {
                    if flag_t {
                        let light_distribution = light_distribution.as_ref();
                        if let Some((s, pi, pdf)) =
                            bssrdf.sample_s(scene, sampler.get_1d(), &sampler.get_2d(), arena)
                        {
                            if s.is_black() || pdf == 0.0 {
                                break;
                            }
                            beta *= s / pdf;

                            let tpi = Interaction::from(&pi);
                            //let pi = tisect.as_surface_interaction().unwrap();

                            // Account for the direct subsurface scattering component
                            let distrib = light_distribution.lookup(&pi.p);
                            l += beta
                                * uniform_sample_one_light(
                                    &tpi,
                                    scene,
                                    arena,
                                    sampler,
                                    false,
                                    Some(distrib.as_ref()),
                                );

                            // Account for the indirect subsurface scattering component
                            if let Some(pi_bsdf) = pi.bsdf.as_ref() {
                                if let Some((f, wi, pdf, flags)) =
                                    pi_bsdf.sample_f(&pi.wo, &sampler.get_2d(), BSDF_ALL)
                                {
                                    assert!(isect.n.length() > 0.0);
                                    assert!(isect.shading.n.length() > 0.0);

                                    if f.is_black() || pdf == 0.0 {
                                        break;
                                    }

                                    beta *= f * (Vector3f::abs_dot(&wi, &pi.shading.n) / pdf);
                                    assert!(Float::is_finite(beta.y()));
                                    specular_bounce = (flags & BSDF_SPECULAR) != 0;
                                    ray = pi.spawn_ray(&wi).into();
                                } else {
                                    break;
                                }
                            }
                        } else {
                            break;
                        }
                    }
                }

                // Possibly terminate the path with Russian roulette.
                // Factor out radiance scaling due to refraction in rrBeta.
                let rr_threshold = self.rr_threshold;
                let rr_beta = beta * eta_scale;
                if rr_beta.max_component_value() < rr_threshold && bounces > 3 {
                    let q = Float::max(0.05, 1.0 - rr_beta.max_component_value());
                    if sampler.get_1d() < q {
                        break;
                    }
                    beta /= 1.0 - q;
                }
            } else {
                break;
            }
            bounces += 1;
        }

        {
            PATH_LENGTH.with(|stat| stat.add(bounces as u64));
        }

        return l;
    }

    fn get_sampler(&self) -> Arc<RwLock<dyn Sampler>> {
        return Arc::clone(&self.base.sampler);
    }

    fn get_pixel_bounds(&self) -> Bounds2i {
        return self.base.pixel_bounds;
    }
}

pub fn create_path_integrator(
    params: &ParamSet,
    sampler: &Arc<RwLock<dyn Sampler>>,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    let max_depth = params.find_one_int("maxdepth", 5);
    let film = camera.as_ref().get_film();
    let pixel_bounds = film.read().unwrap().get_sample_bounds();

    let rr_threshold = params.find_one_float("rrthreshold", 1.0);
    let light_strategy = params.find_one_string("lightsamplestrategy", "spatial");
    return Ok(Arc::new(RwLock::new(PathIntegrator::new(
        max_depth,
        camera,
        sampler,
        &pixel_bounds,
        rr_threshold,
        &light_strategy,
    ))));
}
