use crate::core::camera::*;
use crate::core::distribution::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::lightdistrib::*;
use crate::core::lowdiscrepancy::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::rng::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use crate::core::texture::*;
use crate::core::transform::*;

use std::sync::Arc;

pub struct NormalTexture {}

impl Texture<Spectrum> for NormalTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        let n = si.shading.n;
        let n = n * 2.0 - Normal3f::new(1.0, 1.0, 1.0);
        return Spectrum::from([n[0], n[1], n[2]]);
    }
}

pub fn create_normal_texture(
    _tex2world: &Transform,
    _tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    return Ok(Arc::new(NormalTexture {}));
}
