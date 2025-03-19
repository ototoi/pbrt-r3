use crate::core::error::*;
use crate::core::interaction::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

use std::sync::Arc;

struct FBmTexture {
    mapping: Box<dyn TextureMapping3D>,
    octaves: u32,
    omega: Float,
}

impl FBmTexture {
    pub fn new(mapping: Box<dyn TextureMapping3D>, octaves: u32, omega: Float) -> Self {
        Self {
            mapping,
            octaves,
            omega,
        }
    }

    fn evaluate_f(&self, si: &SurfaceInteraction) -> Float {
        let (p, dpdx, dpdy) = self.mapping.as_ref().map(si);
        return fbm(&p, &dpdx, &dpdy, self.omega, self.octaves);
    }
}

impl Texture<Float> for FBmTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Float {
        return self.evaluate_f(si);
    }
}

impl Texture<Spectrum> for FBmTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        return Spectrum::from(self.evaluate_f(si));
    }
}

pub fn create_fbm_float_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let map: Box<dyn TextureMapping3D> = Box::new(IdentityMapping3D::new(tex2world));
    let octaves = tp.find_int("octaves", 8) as u32;
    let roughness = tp.find_float("roughness", 0.5);
    return Ok(Arc::new(FBmTexture::new(map, octaves, roughness)));
}

pub fn create_fbm_spectrum_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let map: Box<dyn TextureMapping3D> = Box::new(IdentityMapping3D::new(tex2world));
    let octaves = tp.find_int("octaves", 8) as u32;
    let roughness = tp.find_float("roughness", 0.5);
    return Ok(Arc::new(FBmTexture::new(map, octaves, roughness)));
}
