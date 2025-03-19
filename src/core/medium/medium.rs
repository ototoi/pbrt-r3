use std::fmt::Debug;

use crate::core::geometry::ray::Ray;
use crate::core::interaction::MediumInteraction;
use crate::core::memory::MemoryArena;
use crate::core::pbrt::*;
use crate::core::sampler::Sampler;
use crate::core::spectrum::*;

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
