use crate::core::prelude::*;

use std::sync::Arc;

struct WindyTexture {
    mapping: TextureMapping3D,
}

impl WindyTexture {
    pub fn new(mapping: TextureMapping3D) -> Self {
        Self { mapping }
    }

    fn evaluate_f(&self, si: &SurfaceInteraction) -> Float {
        let (p, dpdx, dpdy) = self.mapping.map(si);
        let wind_strength = fbm(&(0.1 * p), &(0.1 * dpdx), &(0.1 * dpdy), 0.5, 3);
        let wave_height = fbm(&p, &dpdx, &dpdy, 0.5, 6);
        return Float::abs(wind_strength) * wave_height;
    }
}

impl Texture<Float> for WindyTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Float {
        return self.evaluate_f(si);
    }
}

impl Texture<Spectrum> for WindyTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        return Spectrum::from(self.evaluate_f(si));
    }
}

pub fn create_windy_float_texture(
    tex2world: &Transform,
    _tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let map = TextureMapping3D::Identity(IdentityMapping3D::new(tex2world));
    return Ok(Arc::new(WindyTexture::new(map)));
}

pub fn create_windy_spectrum_texture(
    tex2world: &Transform,
    _tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let map = TextureMapping3D::Identity(IdentityMapping3D::new(tex2world));
    return Ok(Arc::new(WindyTexture::new(map)));
}
