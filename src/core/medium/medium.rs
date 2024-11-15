use std::fmt::Debug;

use crate::core::pbrt::*;

pub trait Medium: Debug + Send + Sync {
    fn tr(&self, _ray: &Ray, _sampler: &mut dyn Sampler) -> Spectrum {
        Spectrum::zero()
    }
    fn sample(
        &self,
        _ray: &Ray,
        _sampler: &mut dyn Sampler,
        _arena: &mut MemoryArena,
    ) -> (Spectrum, Option<MediumInteraction>) {
        (Spectrum::one(), None)
    }
}
