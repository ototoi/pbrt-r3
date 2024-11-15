use crate::core::pbrt::*;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum DirectLightingLightStrategy {
    UniformSampleAll,
    UniformSampleOne,
}

type LightStrategy = DirectLightingLightStrategy;

pub struct DirectLightingIntegrator {
    base: BaseSamplerIntegrator,
    strategy: LightStrategy,
    max_depth: i32,
    n_light_samples: Vec<u32>,
}

impl DirectLightingIntegrator {
    pub fn new(
        strategy: LightStrategy,
        max_depth: i32,
        camera: &Arc<dyn Camera>,
        sampler: &Arc<RwLock<dyn Sampler>>,
        pixel_bounds: &Bounds2i,
    ) -> Self {
        DirectLightingIntegrator {
            base: BaseSamplerIntegrator::new(camera, sampler, pixel_bounds),
            strategy,
            max_depth,
            n_light_samples: Vec::new(),
        }
    }
}

impl Integrator for DirectLightingIntegrator {
    fn render(&mut self, scene: &Scene) {
        BaseSamplerIntegrator::render(self, scene);
    }

    fn get_camera(&self) -> Arc<dyn Camera> {
        return Arc::clone(&self.base.camera);
    }
}

unsafe impl Sync for DirectLightingIntegrator {}

impl SamplerIntegrator for DirectLightingIntegrator {
    fn preprocess(&mut self, scene: &Scene, sampler: &mut dyn Sampler) {
        if self.strategy == LightStrategy::UniformSampleAll {
            for light in scene.lights.iter() {
                let count = sampler.round_count(light.as_ref().get_sample_count());
                self.n_light_samples.push(count);
            }
            for _ in 0..self.max_depth {
                for j in 0..scene.lights.len() {
                    let n = self.n_light_samples[j];
                    sampler.request_2d_array(n);
                    sampler.request_2d_array(n);
                }
            }
        }
    }
    fn li(
        &self,
        ray: &RayDifferential,
        scene: &Scene,
        sampler: &mut dyn Sampler,
        arena: &mut MemoryArena,
        depth: i32,
    ) -> Spectrum {
        if let Some(mut isect) = scene.intersect(&ray.ray) {
            let wo = isect.wo;

            isect.compute_scattering_functions(ray, arena, TransportMode::Radiance, false);

            let tisect = Interaction::from(isect);
            let isect = tisect.as_surface_interaction().unwrap();

            if isect.bsdf.is_some() {
                let mut l = isect.le(&wo);
                if !scene.lights.is_empty() {
                    match self.strategy {
                        LightStrategy::UniformSampleAll => {
                            l += uniform_sample_all_lights(
                                &tisect,
                                scene,
                                arena,
                                sampler,
                                &self.n_light_samples,
                                false,
                            );
                        }
                        LightStrategy::UniformSampleOne => {
                            l += uniform_sample_one_light(
                                &tisect, scene, arena, sampler, false, None,
                            );
                        }
                    }
                }
                if depth + 1 < self.max_depth {
                    l += self.specular_reflect(ray, isect, scene, sampler, arena, depth);
                    l += self.specular_transmit(ray, isect, scene, sampler, arena, depth);
                }
                return l;
            } else {
                let rd = RayDifferential::from(isect.spawn_ray(&ray.ray.d));
                return self.li(&rd, scene, sampler, arena, depth);
            }
        } else {
            let l = scene
                .lights
                .iter()
                .map(|light| -> Spectrum {
                    return light.as_ref().le(ray);
                })
                .fold(Spectrum::zero(), |a, b| {
                    return a + b;
                });
            return l;
        }
    }

    fn get_sampler(&self) -> Arc<RwLock<dyn Sampler>> {
        return Arc::clone(&self.base.sampler);
    }

    fn get_pixel_bounds(&self) -> Bounds2i {
        return self.base.pixel_bounds;
    }
}

pub fn create_direct_lighting_integrator(
    params: &ParamSet,
    sampler: &Arc<RwLock<dyn Sampler>>,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    let max_depth = params.find_one_int("maxdepth", 5);
    let st = params.find_one_string("strategy", "all");
    let sterategy = if st == "one" {
        LightStrategy::UniformSampleOne
    } else {
        LightStrategy::UniformSampleAll
    };
    let film = camera.as_ref().get_film();
    let pixel_bounds = film.read().unwrap().get_sample_bounds();

    return Ok(Arc::new(RwLock::new(DirectLightingIntegrator::new(
        sterategy,
        max_depth,
        camera,
        sampler,
        &pixel_bounds,
    ))));
}
