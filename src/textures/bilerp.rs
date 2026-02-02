use crate::core::prelude::*;

use std::ops::*;
use std::sync::Arc;

pub struct BilerpTexture<T> {
    mapping: TextureMapping2D,
    v00: T,
    v01: T,
    v10: T,
    v11: T,
}

impl<T: Copy> BilerpTexture<T> {
    pub fn new(mapping: TextureMapping2D, v00: &T, v01: &T, v10: &T, v11: &T) -> Self {
        return BilerpTexture::<T> {
            mapping,
            v00: *v00,
            v01: *v01,
            v10: *v10,
            v11: *v11,
        };
    }
}

impl<T: Copy + Add<T, Output = T> + Mul<Float, Output = T>> Texture<T> for BilerpTexture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let (st, _dstdx, _dstdy) = self.mapping.map(si);
        let a = (1.0 - st[0]) * (1.0 - st[1]);
        let b = (1.0 - st[0]) * (st[1]);
        let c = (st[0]) * (1.0 - st[1]);
        let d = (st[0]) * (st[1]);
        return self.v00 * a + self.v01 * b + self.v10 * c + self.v11 * d;
    }
}

pub fn create_bilerp_float_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let map = create_texture_mapping2d(tex2world, tp)?;
    let v00 = tp.find_float("v00", 0.0);
    let v01 = tp.find_float("v01", 1.0);
    let v10 = tp.find_float("v10", 0.0);
    let v11 = tp.find_float("v11", 1.0);
    return Ok(Arc::new(BilerpTexture::<Float>::new(
        map, &v00, &v01, &v10, &v11,
    )));
}

pub fn create_bilerp_spectrum_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let map = create_texture_mapping2d(tex2world, tp)?;
    let v00 = tp.find_spectrum("v00", &Spectrum::zero());
    let v01 = tp.find_spectrum("v01", &Spectrum::one());
    let v10 = tp.find_spectrum("v10", &Spectrum::zero());
    let v11 = tp.find_spectrum("v11", &Spectrum::one());
    return Ok(Arc::new(BilerpTexture::<Spectrum>::new(
        map, &v00, &v01, &v10, &v11,
    )));
}
