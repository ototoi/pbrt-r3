use crate::core::error::*;
use crate::core::interaction::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

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
