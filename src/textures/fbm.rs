use crate::core::prelude::*;

use std::sync::Arc;

struct FBmTexture {
    mapping: TextureMapping3D,
    octaves: u32,
    omega: Float,
}

impl FBmTexture {
    pub fn new(mapping: TextureMapping3D, octaves: u32, omega: Float) -> Self {
        Self {
            mapping,
            octaves,
            omega,
        }
    }

    fn evaluate_f(&self, si: &SurfaceInteraction) -> Float {
        let (p, dpdx, dpdy) = self.mapping.map(si);
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
    let map = TextureMapping3D::Identity(IdentityMapping3D::new(tex2world));
    let octaves = tp.find_int("octaves", 8) as u32;
    let roughness = tp.find_float("roughness", 0.5);
    return Ok(Arc::new(FBmTexture::new(map, octaves, roughness)));
}

pub fn create_fbm_spectrum_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let map = TextureMapping3D::Identity(IdentityMapping3D::new(tex2world));
    let octaves = tp.find_int("octaves", 8) as u32;
    let roughness = tp.find_float("roughness", 0.5);
    return Ok(Arc::new(FBmTexture::new(map, octaves, roughness)));
}
