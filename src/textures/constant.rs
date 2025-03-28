use crate::core::prelude::*;

use std::sync::Arc;

pub fn create_constant_float_texture(
    _tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let value = tp.find_float("value", 1.0);
    return Ok(Arc::new(ConstantTexture::<Float>::new(&value)));
}

pub fn create_constant_spectrum_texture(
    _tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let value = tp.find_spectrum("value", &Spectrum::one());
    return Ok(Arc::new(ConstantTexture::<Spectrum>::new(&value)));
}
