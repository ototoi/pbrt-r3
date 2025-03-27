use crate::core::prelude::*;

use std::sync::Arc;

struct WrinkledTexture {
    mapping: Box<dyn TextureMapping3D>,
    octaves: u32,
    omega: Float,
}

impl WrinkledTexture {
    pub fn new(mapping: Box<dyn TextureMapping3D>, octaves: u32, omega: Float) -> Self {
        Self {
            mapping,
            octaves,
            omega,
        }
    }

    fn evaluate_f(&self, si: &SurfaceInteraction) -> Float {
        let (p, dpdx, dpdy) = self.mapping.as_ref().map(si);
        return turbulence(&p, &dpdx, &dpdy, self.omega, self.octaves);
    }
}

impl Texture<Float> for WrinkledTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Float {
        return self.evaluate_f(si);
    }
}

impl Texture<Spectrum> for WrinkledTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        return Spectrum::from(self.evaluate_f(si));
    }
}

pub fn create_wrinkled_float_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let map: Box<dyn TextureMapping3D> = Box::new(IdentityMapping3D::new(tex2world));
    let octaves = tp.find_int("octaves", 8) as u32;
    let roughness = tp.find_float("roughness", 0.5);
    return Ok(Arc::new(WrinkledTexture::new(map, octaves, roughness)));
}

pub fn create_wrinkled_spectrum_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let map: Box<dyn TextureMapping3D> = Box::new(IdentityMapping3D::new(tex2world));
    let octaves = tp.find_int("octaves", 8) as u32;
    let roughness = tp.find_float("roughness", 0.5);
    return Ok(Arc::new(WrinkledTexture::new(map, octaves, roughness)));
}
