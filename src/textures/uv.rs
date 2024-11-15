use crate::core::pbrt::*;

use std::sync::Arc;

pub struct UVTexture {
    mapping: Box<dyn TextureMapping2D>,
}

impl UVTexture {
    pub fn new(mapping: Box<dyn TextureMapping2D>) -> Self {
        UVTexture { mapping }
    }
}

impl Texture<Spectrum> for UVTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        let (st, _dstdx, _dstdy) = self.mapping.map(si);
        let rgb = [
            st[0] - Float::floor(st[0]),
            st[1] - Float::floor(st[1]),
            0.0,
        ];
        return Spectrum::from_rgb(&rgb, SpectrumType::Reflectance);
    }
}

pub fn create_uv_spectrum_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let map = create_texture_mapping2d(tex2world, tp)?;
    return Ok(Arc::new(UVTexture::new(map)));
}
