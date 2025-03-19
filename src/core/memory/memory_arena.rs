use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::*;
use crate::core::reflection::BSDF;

pub struct MemoryArena {}

impl MemoryArena {
    pub fn new() -> Self {
        MemoryArena {}
    }

    pub fn reset(&mut self) {}

    pub fn alloc_bsdf(&mut self, si: &SurfaceInteraction, eta: Float) -> BSDF {
        BSDF::new(si, eta)
    }
}

unsafe impl Send for MemoryArena {}
unsafe impl Sync for MemoryArena {}
