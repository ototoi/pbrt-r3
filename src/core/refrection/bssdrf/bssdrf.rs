use std::fmt::Debug;

use super::functions::*;
use crate::core::pbrt::*;

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

    pub fn sw(&self, w: &Vector3f) -> Spectrum {
        let eta = self.eta;
        let c = 1.0 - 2.0 * fresnel_moment1(1.0 / eta);
        return Spectrum::from((1.0 - fr_dielectric(cos_theta(w), 1.0, eta)) / (c * PI));
    }
}
