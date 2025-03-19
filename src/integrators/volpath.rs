use crate::core::camera::*;
use crate::core::error::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::lightdistrib::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::memory::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::sampler::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use crate::core::stats::*;

use std::sync::Arc;
use std::sync::RwLock;

thread_local!(static PATH_LENGTH: StatIntDistribution = StatIntDistribution::new("Integrator/Path length"));
thread_local!(static VOLUME_INTERACTIONS: StatCounter = StatCounter::new("Integrator/Volume interactions"));
thread_local!(static SURFACE_INTERACTIONS: StatCounter = StatCounter::new("Integrator/Surface interactions"));

pub struct VolPathIntegrator {
    pub base: BaseSamplerIntegrator,
    pub max_depth: u32,
    pub rr_threshold: Float,
    pub light_sample_strategy: String,
    pub light_distribution: Option<Arc<dyn LightDistribution>>,
}

impl VolPathIntegrator {
    pub fn new(
        max_depth: u32,
        camera: &Arc<dyn Camera>,
        sampler: &Arc<RwLock<dyn Sampler>>,
        pixel_bounds: &Bounds2i,
        rr_threshold: Float,
        light_sample_strategy: String,
    ) -> Self {
        let base = BaseSamplerIntegrator::new(camera, sampler, pixel_bounds);
        VolPathIntegrator {
            base,
            max_depth,
            rr_threshold,
            light_sample_strategy,
            light_distribution: None,
        }
    }
}

impl Integrator for VolPathIntegrator {
    fn render(&mut self, scene: &Scene) {
        BaseSamplerIntegrator::render(self, scene);
    }
    fn get_camera(&self) -> Arc<dyn Camera> {
        self.base.camera.clone()
    }
}

unsafe impl Sync for VolPathIntegrator {}

impl SamplerIntegrator for VolPathIntegrator {
    fn preprocess(&mut self, scene: &Scene, _sampler: &mut dyn Sampler) {
        match create_light_sample_distribution(&self.light_sample_strategy, scene) {
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

        let rr_threshold = self.rr_threshold;
        let mut l = Spectrum::zero();
        let mut beta = Spectrum::one();
        let mut ray = r.clone();
        let mut specular_bounce = false;

        let light_distr = self.light_distribution.clone();
        let light_distr = light_distr.unwrap();

        // Added after book publication: etaScale tracks the accumulated effect
        // of radiance scaling due to rays passing through refractive
        // boundaries (see the derivation on p. 527 of the third edition). We
        // track this value in order to remove it from beta when we apply
        // Russian roulette; this is worthwhile, since it lets us sometimes
        // avoid terminating refracted rays that are about to be refracted back
        // out of a medium and thus have their beta value increased.
        let mut eta_scale = 1.0;

        let mut bounces = 0;
        loop {
            // Find next path vertex and accumulate contribution

            // Intersect _ray_ with scene and store intersection in _isect_
            let mut found_intersection = scene.intersect(&ray.ray);

            let mut mi = None;
            // Sample the participating medium, if present
            if let Some(medium) = ray.ray.medium.as_ref() {
                //println!("MediumInteraction: {:?}", medium);
                let (spec, m) = medium.sample(&ray.ray, sampler, arena);
                if let Some(mut m) = m {
                    m.medium_interface = MediumInterface::from(medium);
                    mi = Some(m);
                }
                beta *= spec;
            }
            if beta.is_black() {
                break;
            }

            if let Some(mi) = mi {
                // Terminate path if ray escaped or _maxDepth_ was reached
                if bounces >= self.max_depth {
                    break;
                }

                VOLUME_INTERACTIONS.with(|c| c.inc());

                // Handle scattering at point in medium for volumetric path tracer
                let light_distrib = light_distr.lookup(&mi.p);
                let mit = Interaction::from(&mi);
                let ld = beta
                    * uniform_sample_one_light(
                        &mit,
                        scene,
                        arena,
                        sampler,
                        true,
                        Some(&light_distrib),
                    );
                l += ld;
                let wo = -ray.ray.d;
                let (_pdf, wi) = mi.phase.sample_p(&wo, &sampler.get_2d());
                ray = mi.spawn_ray(&wi).into();
                specular_bounce = false;
            } else {
                SURFACE_INTERACTIONS.with(|c| c.inc());
                // Handle scattering at point on surface for volumetric path tracer

                // Possibly add emitted light at intersection
                if bounces == 0 || specular_bounce {
                    // Add emitted light at path vertex or from the environment
                    if let Some(isect) = found_intersection.as_ref() {
                        l += beta * isect.le(&-ray.ray.d);
                    } else {
                        for light in scene.infinite_lights.iter() {
                            l += beta * light.le(&ray);
                        }
                    }
                }

                // Terminate path if ray escaped or _maxDepth_ was reached
                if found_intersection.is_none() || bounces >= self.max_depth {
                    break;
                }

                // Compute scattering functions and skip over medium boundaries
                let isect = found_intersection.as_mut().unwrap();
                //if let Some(isect) = found_intersection.as_mut() {
                isect.compute_scattering_functions(&ray, arena, TransportMode::Radiance, true);
                if isect.bsdf.is_none() {
                    ray = isect.spawn_ray(&ray.ray.d).into();
                    continue;
                }
                let distrib = light_distr.lookup(&isect.p);

                // Sample illumination from lights to find path contribution.
                {
                    let tisect = Interaction::from(isect as &SurfaceInteraction);
                    let ld = beta
                        * uniform_sample_one_light(
                            &tisect,
                            scene,
                            arena,
                            sampler,
                            true,
                            Some(&distrib),
                        );
                    l += ld;
                }

                // Sample BSDF to get new path direction
                let bsdf = isect.bsdf.as_ref().unwrap();
                let wo = -ray.ray.d;
                if let Some((f, wi, pdf, flags)) = bsdf.sample_f(&wo, &sampler.get_2d(), BSDF_ALL) {
                    if f.is_black() || pdf == 0.0 {
                        break;
                    }
                    beta *= f * (wi.abs_dot(&isect.shading.n) / pdf);
                    //beta *= 0.1;

                    specular_bounce = (flags & BSDF_SPECULAR) != 0;
                    let flag_s = (flags & BSDF_SPECULAR) != 0;
                    let flag_t = (flags & BSDF_TRANSMISSION) != 0;
                    if flag_s && flag_t {
                        let eta = bsdf.eta;
                        // Update the term that tracks radiance scaling for refraction
                        // depending on whether the ray is entering or leaving the
                        // medium.
                        if Vector3f::dot(&wo, &isect.n) > 0.0 {
                            eta_scale *= eta * eta;
                        } else {
                            eta_scale *= 1.0 / (eta * eta);
                        }
                    }

                    ray = isect.spawn_ray(&wi).into();

                    // Account for subsurface scattering, if applicable
                    if let Some(bssrdf) = isect.bssrdf.as_ref() {
                        if flag_t {
                            // Importance sample the BSSRDF
                            if let Some((s, pi, pdf)) = bssrdf.as_ref().sample_s(
                                scene,
                                sampler.get_1d(),
                                &sampler.get_2d(),
                                arena,
                            ) {
                                if s.is_black() || pdf == 0.0 {
                                    break;
                                }
                                beta *= s / pdf;
                                //beta *= 0.1;

                                // Account for the direct subsurface scattering component
                                let distrib = light_distr.lookup(&pi.p);
                                let pit = Interaction::from(&pi);
                                l += beta
                                    * uniform_sample_one_light(
                                        &pit,
                                        scene,
                                        arena,
                                        sampler,
                                        true,
                                        Some(&distrib),
                                    );

                                // Account for the indirect subsurface scattering component
                                if let Some(pi_bsdf) = pi.bsdf.as_ref() {
                                    if let Some((f, wi, pdf, flags)) =
                                        pi_bsdf.sample_f(&pi.wo, &sampler.get_2d(), BSDF_ALL)
                                    {
                                        if f.is_black() || pdf == 0.0 {
                                            break;
                                        }
                                        beta *= f * (wi.abs_dot(&pi.shading.n) / pdf);
                                        //beta *= 0.1;
                                        specular_bounce = (flags & BSDF_SPECULAR) != 0;
                                        ray = pi.spawn_ray(&wi).into();
                                    } else {
                                        break;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                } else {
                    break;
                }
            }

            // Possibly terminate the path with Russian roulette.
            // Factor out radiance scaling due to refraction in rrBeta.
            let rr_beta = beta * eta_scale;
            if rr_beta.max_component_value() < rr_threshold && bounces > 3 {
                let q = Float::max(0.05, 1.0 - rr_beta.max_component_value());
                if sampler.get_1d() < q {
                    break;
                }
                beta /= 1.0 - q;
                assert!(beta.max_component_value().is_finite());
            }

            bounces += 1;
        } // end loop

        PATH_LENGTH.with(|c| c.add(bounces as u64));

        return l;
    }

    fn get_sampler(&self) -> Arc<RwLock<dyn Sampler>> {
        return self.base.sampler.clone();
    }

    fn get_pixel_bounds(&self) -> Bounds2i {
        return self.base.pixel_bounds;
    }
}

pub fn create_volpath_integrator(
    params: &ParamSet,
    sampler: &Arc<RwLock<dyn Sampler>>,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    let max_depth = params.find_one_int("maxdepth", 5) as u32;

    let pixel_bounds = camera.get_film().read().unwrap().get_sample_bounds();

    let rr_threshold = params.find_one_float("rrthreshold", 1.0);
    let light_sample_strategy = params.find_one_string("lightsamplestrategy", "spatial");
    return Ok(Arc::new(RwLock::new(VolPathIntegrator::new(
        max_depth,
        camera,
        sampler,
        &pixel_bounds,
        rr_threshold,
        light_sample_strategy,
    ))));
}
