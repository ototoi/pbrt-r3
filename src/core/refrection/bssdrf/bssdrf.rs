use crate::core::interaction::*;
use crate::core::memory::*;
use crate::core::pbrt::*;
use crate::core::scene::*;
use crate::core::spectrum::*;

use std::fmt::Debug;

pub trait BSSRDF: Debug {
    fn s(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum;
    fn sample_s(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        arena: &mut MemoryArena,
    ) -> Option<(Spectrum, SurfaceInteraction, Float)>;
}

#[derive(Clone, Copy, Default, Debug)]
pub struct BaseBSSRDF {
    //pub po: SurfaceInteraction,
    pub p: Point3f,
    pub wo: Vector3f,
    pub time: Float,
    pub eta: Float,
}

impl BaseBSSRDF {
    pub fn new(po: &SurfaceInteraction, eta: Float) -> Self {
        BaseBSSRDF {
            p: po.p,
            wo: po.wo,
            time: po.time,
            eta,
        }
    }
}
