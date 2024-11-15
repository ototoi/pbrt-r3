use crate::core::pbrt::*;
use std::sync::Arc;

thread_local!(pub static N_INTERSECTION_TESTS: StatCounter = StatCounter::new("Intersections/Regular ray intersection tests"));
thread_local!(pub static N_SHADOW_TESTS: StatCounter = StatCounter::new("Intersections/Shadow ray intersection tests"));

pub struct Scene {
    pub lights: Vec<Arc<dyn Light>>,
    pub infinite_lights: Vec<Arc<dyn Light>>,
    pub aggregate: Arc<dyn Primitive>,
    pub world_bound: Bounds3f,
}

impl Scene {
    pub fn new(aggregate: &Arc<dyn Primitive>, lights: &[Arc<dyn Light>]) -> Self {
        let world_bound = aggregate.world_bound();
        let infinite_lights: Vec<Arc<dyn Light>> =
            lights.iter().filter(|l| l.is_infinite()).cloned().collect();
        let scene = Scene {
            lights: lights.to_vec(),
            infinite_lights,
            aggregate: aggregate.clone(),
            world_bound,
        };
        for light in lights.iter() {
            let l = light.as_ref();
            l.preprocess(&scene);
        }

        return scene;
    }

    pub fn world_bound(&self) -> Bounds3f {
        return self.world_bound;
    }

    pub fn intersect(&self, ray: &Ray) -> Option<SurfaceInteraction> {
        N_INTERSECTION_TESTS.with(|c| c.inc());
        let aggregate = self.aggregate.as_ref();
        return aggregate.intersect(ray);
    }

    pub fn intersect_p(&self, ray: &Ray) -> bool {
        N_SHADOW_TESTS.with(|c| c.inc());
        let aggregate = self.aggregate.as_ref();
        return aggregate.intersect_p(ray);
    }

    pub fn intersect_tr(
        &self,
        r: &Ray,
        sampler: &mut dyn Sampler,
    ) -> (Option<SurfaceInteraction>, Spectrum) {
        let mut ray = r.clone();
        let mut tr = Spectrum::one();
        loop {
            if let Some(it) = self.intersect(&ray) {
                if let Some(medium) = ray.medium.as_ref() {
                    tr *= medium.as_ref().tr(&ray, sampler);
                }

                if let Some(prim) = it.get_primitive() {
                    if prim.get_material().is_some() {
                        return (Some(it), tr);
                    }
                }

                ray = it.spawn_ray(&ray.d);
            } else {
                if let Some(medium) = ray.medium.as_ref() {
                    tr *= medium.as_ref().tr(&ray, sampler);
                }
                return (None, tr);
            }
        }
    }
}

unsafe impl Sync for Scene {}
