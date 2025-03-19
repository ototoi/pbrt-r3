use crate::core::error::*;
use crate::core::interaction::*;
use crate::core::material::*;
use crate::core::memory::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::reflection::*;
use crate::core::texture::*;

use std::sync::Arc;

pub struct FourierMaterial {
    bsdf_table: Arc<FourierBSDFTable>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
}

impl FourierMaterial {
    pub fn new(
        bsdf_table: &Arc<FourierBSDFTable>,
        bumpmap: &Option<Arc<dyn Texture<Float>>>,
    ) -> Self {
        FourierMaterial {
            bsdf_table: bsdf_table.clone(),
            bumpmap: bumpmap.clone(),
        }
    }
}

impl Material for FourierMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        // Perform bump mapping with _bumpMap_, if present
        if let Some(bump) = self.bumpmap.as_ref() {
            self.bump(bump, si);
        }
        let mut b = arena.alloc_bsdf(si, 1.0);

        // Checking for zero channels works as a proxy for checking whether the
        // table was successfully read from the file.
        if self.bsdf_table.as_ref().n_channels > 0 {
            let fourier: Arc<dyn BxDF> = Arc::new(FourierBSDF::new(&self.bsdf_table, mode));
            b.add(&fourier);
        }
        si.bsdf = Some(Arc::new(b));
    }
}

pub fn create_fourier_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let filename = mp.find_filename("bsdffile", "");
    let bsdf_table = FourierBSDFTable::read(&filename)?;
    let bsdf_table = Arc::new(bsdf_table);
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    return Ok(Arc::new(FourierMaterial::new(&bsdf_table, &bumpmap)));
}
