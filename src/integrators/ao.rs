use crate::core::camera::*;
use crate::core::error::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::material::*;
use crate::core::memory::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::spectrum::*;

use std::sync::Arc;
use std::sync::RwLock;

use log::warn;

pub struct AOIntegrator {
    base: BaseSamplerIntegrator,
    cos_sample: bool,
    n_samples: u32,
}

impl AOIntegrator {
    pub fn new(
        cos_sample: bool,
        ns: u32,
        camera: &Arc<dyn Camera>,
        sampler: &Arc<RwLock<dyn Sampler>>,
        pixel_bounds: &Bounds2i,
    ) -> Self {
        let n_samples = sampler.read().unwrap().round_count(ns);
        if ns != n_samples {
            warn!(
                "Rounding AO samples to {} (from {} specified).",
                n_samples, ns
            );
        }
        sampler.write().unwrap().request_2d_array(n_samples);
        AOIntegrator {
            base: BaseSamplerIntegrator::new(camera, sampler, pixel_bounds),
            cos_sample,
            n_samples,
        }
    }
}

impl Integrator for AOIntegrator {
    fn render(&mut self, scene: &Scene) {
        BaseSamplerIntegrator::render(self, scene);
    }

    fn get_camera(&self) -> Arc<dyn Camera> {
        return self.base.camera.clone();
    }
}

impl SamplerIntegrator for AOIntegrator {
    fn li(
        &self,
        r: &RayDifferential,
        scene: &Scene,
        sampler: &mut dyn Sampler,
        arena: &mut MemoryArena,
        _depth: i32,
    ) -> Spectrum {
        let _p = ProfilePhase::new(Prof::SamplerIntegratorLi);

        let mut l = Spectrum::zero();
        let mut ray = r.clone();

        let cos_sample = self.cos_sample;
        let n_samples = self.n_samples as u32;
        // Intersect _ray_ with scene and store intersection in _isect_
        loop {
            if let Some(mut isect) = scene.intersect(&ray.ray) {
                isect.compute_scattering_functions(&ray, arena, TransportMode::Radiance, true);
                if isect.bsdf.is_none() {
                    ray = isect.spawn_ray(&ray.ray.d).into();
                    continue;
                }
                // Compute coordinate frame based on true geometry, not shading
                // geometry.
                let n = face_forward(&isect.n, &-ray.ray.d);
                let s = isect.dpdu.normalize();
                let t = Vector3f::cross(&isect.n, &s);

                let u = sampler.get_2d_array(n_samples).unwrap();
                for i in 0..u.len() {
                    let (wi, pdf) = if cos_sample {
                        let wi = cosine_sample_hemisphere(&u[i]);
                        let pdf = cosine_hemisphere_pdf(wi.z.abs());
                        (wi, pdf)
                    } else {
                        let wi = uniform_sample_hemisphere(&u[i]);
                        let pdf = uniform_hemisphere_pdf();
                        (wi, pdf)
                    };

                    // Transform wi from local frame to world space.
                    let wi = Vector3f::new(
                        wi.x * s.x + wi.y * t.x + wi.z * n.x,
                        wi.x * s.y + wi.y * t.y + wi.z * n.y,
                        wi.x * s.z + wi.y * t.z + wi.z * n.z,
                    );
                    if !scene.intersect_p(&isect.spawn_ray(&wi)) {
                        l += Spectrum::from(wi.dot(&n) / (pdf * n_samples as Float));
                    }
                }
            }
            break;
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

unsafe impl Sync for AOIntegrator {}

pub fn create_ao_integrator(
    params: &ParamSet,
    sampler: &Arc<RwLock<dyn Sampler>>,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    let pixel_bounds = camera.get_film().read().unwrap().get_sample_bounds();
    let cos_sample = params.find_one_bool("cossample", true);
    let mut n_samples = params.find_one_int("nsamples", 64);
    {
        let options = PbrtOptions::get();
        if options.quick_render {
            n_samples = 1;
        }
    }
    Ok(Arc::new(RwLock::new(AOIntegrator::new(
        cos_sample,
        n_samples as u32,
        camera,
        sampler,
        &pixel_bounds,
    ))))
}
