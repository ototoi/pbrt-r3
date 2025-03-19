use crate::core::error::*;
use crate::core::interaction::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

use std::ops::*;
use std::sync::Arc;

pub struct MixTexture<T> {
    tex1: Arc<dyn Texture<T>>,
    tex2: Arc<dyn Texture<T>>,
    amount: Arc<dyn Texture<Float>>,
}

impl<T: Copy> MixTexture<T> {
    pub fn new(
        tex1: &Arc<dyn Texture<T>>,
        tex2: &Arc<dyn Texture<T>>,
        amount: &Arc<dyn Texture<Float>>,
    ) -> Self {
        return MixTexture::<T> {
            tex1: Arc::clone(tex1),
            tex2: Arc::clone(tex2),
            amount: Arc::clone(amount),
        };
    }
}

impl<T: Copy + Add<T, Output = T> + Mul<Float, Output = T>> Texture<T> for MixTexture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let t1 = self.tex1.as_ref().evaluate(si);
        let t2 = self.tex2.as_ref().evaluate(si);
        let amt = self.amount.as_ref().evaluate(si);
        return t1 * (1.0 - amt) + t2 * amt;
    }
}

pub fn create_mix_float_texture(
    _tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let tex1 = tp.get_float_texture("tex1", 1.0);
    let tex2 = tp.get_float_texture("tex2", 1.0);
    let amount = tp.get_float_texture("amount", 0.5);
    return Ok(Arc::new(MixTexture::<Float>::new(&tex1, &tex2, &amount)));
}

pub fn create_mix_spectrum_texture(
    _tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let tex1 = tp.get_spectrum_texture("tex1", &Spectrum::one());
    let tex2 = tp.get_spectrum_texture("tex2", &Spectrum::one());
    let amount = tp.get_float_texture("amount", 0.5);
    return Ok(Arc::new(MixTexture::<Spectrum>::new(&tex1, &tex2, &amount)));
}
