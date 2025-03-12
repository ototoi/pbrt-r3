use crate::core::pbrt::*;
use std::sync::Arc;
use std::sync::RwLock;

pub struct WhittedIntegrator {
    pub base: BaseSamplerIntegrator,
    pub max_depth: i32,
}

impl WhittedIntegrator {
    pub fn new(
        max_depth: i32,
        camera: &Arc<dyn Camera>,
        sampler: &Arc<RwLock<dyn Sampler>>,
        pixel_bounds: &Bounds2i,
    ) -> Self {
        WhittedIntegrator {
            base: BaseSamplerIntegrator::new(camera, sampler, pixel_bounds),
            max_depth,
        }
    }
}

impl Integrator for WhittedIntegrator {
    fn render(&mut self, scene: &Scene) {
        BaseSamplerIntegrator::render(self, scene);
    }

    fn get_camera(&self) -> Arc<dyn Camera> {
        return Arc::clone(&self.base.camera);
    }
}

unsafe impl Sync for WhittedIntegrator {}

impl SamplerIntegrator for WhittedIntegrator {
    fn li(
        &self,
        ray: &RayDifferential,
        scene: &Scene,
        sampler: &mut dyn Sampler,
        arena: &mut MemoryArena,
        depth: i32,
    ) -> Spectrum {
        let _p = ProfilePhase::new(Prof::SamplerIntegratorLi);

        if let Some(mut isect) = scene.intersect(&ray.ray) {
            assert!(isect.n.length() > 0.0);
            assert!(isect.shading.n.length() > 0.0);

            let n = isect.shading.n;
            let wo = isect.wo;

            isect.compute_scattering_functions(ray, arena, TransportMode::Radiance, false);

            let tisect = Interaction::from(isect);
            let isect = tisect.as_surface_interaction().unwrap();

            if let Some(bsdf) = isect.bsdf.as_ref() {
                //return Spectrum::one();
                let mut l = scene
                    .lights
                    .iter()
                    .map(|light| -> Spectrum {
                        let lt = light.as_ref();
                        let u = sampler.get_2d();
                        if let Some((li, wi, pdf, visibility)) = lt.sample_li(&tisect, &u) {
                            if !(pdf <= 0.0 || li.is_black()) {
                                let f = bsdf.f(&wo, &wi, BSDF_ALL);
                                if !f.is_black() && visibility.unoccluded(scene) {
                                    //let dot = Float::max(Vector3f::dot(&wi, &n), 0.0);
                                    return f * li * (Vector3f::abs_dot(&wi, &n) / pdf);
                                }
                            }
                        }
                        return Spectrum::zero();
                    })
                    .fold(Spectrum::zero(), |a, b| {
                        return a + b;
                    });
                if depth + 1 < self.max_depth {
                    l += self.specular_reflect(ray, isect, scene, sampler, arena, depth);
                    l += self.specular_transmit(ray, isect, scene, sampler, arena, depth);
                }
                return l;
            }
            return Spectrum::zero();
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

pub fn create_whitted_integrator(
    params: &ParamSet,
    sampler: &Arc<RwLock<dyn Sampler>>,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    let max_depth = params.find_one_int("maxdepth", 5);
    let film = camera.as_ref().get_film();
    let mut pixel_bounds = film.read().unwrap().get_sample_bounds();
    if let Some(pb) = params.get_ints_ref("pixelbounds") {
        if pb.len() >= 4 {
            pixel_bounds = Bounds2i::from(((pb[0], pb[2]), (pb[1], pb[3])));
        }
    }
    return Ok(Arc::new(RwLock::new(WhittedIntegrator::new(
        max_depth,
        camera,
        sampler,
        &pixel_bounds,
    ))));
}
