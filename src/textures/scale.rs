use crate::core::error::*;
use crate::core::interaction::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

use std::sync::Arc;

pub struct ScaleTexture<T1, T2> {
    tex1: Arc<dyn Texture<T1>>,
    tex2: Arc<dyn Texture<T2>>,
}

impl<T1: Copy, T2: Copy> ScaleTexture<T1, T2> {
    pub fn new(tex1: &Arc<dyn Texture<T1>>, tex2: &Arc<dyn Texture<T2>>) -> Self {
        return ScaleTexture::<T1, T2> {
            tex1: Arc::clone(tex1),
            tex2: Arc::clone(tex2),
        };
    }
}

impl<T1: Copy + std::ops::Mul<T2, Output = T2>, T2: Copy> Texture<T2> for ScaleTexture<T1, T2> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T2 {
        let value1 = self.tex1.as_ref().evaluate(si);
        let value2 = self.tex2.as_ref().evaluate(si);
        return value1 * value2;
    }
}

pub fn create_scale_float_texture(
    _tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let tex1 = tp.get_float_texture("tex1", 1.0);
    let tex2 = tp.get_float_texture("tex2", 1.0);
    return Ok(Arc::new(ScaleTexture::<Float, Float>::new(&tex1, &tex2)));
}

pub fn create_scale_spectrum_texture(
    _tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let tex1 = tp.get_spectrum_texture("tex1", &Spectrum::one());
    let tex2 = tp.get_spectrum_texture("tex2", &Spectrum::one());
    return Ok(Arc::new(ScaleTexture::<Spectrum, Spectrum>::new(
        &tex1, &tex2,
    )));
}
