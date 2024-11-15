use crate::core::pbrt::*;
use std::sync::Arc;
use std::sync::RwLock;

pub struct FourierMaterial {
    bsdf_table: Arc<RwLock<FourierBSDFTable>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
}

impl FourierMaterial {
    pub fn new(
        bsdf_table: &Arc<RwLock<FourierBSDFTable>>,
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
        if let Some(bump) = self.bumpmap.as_ref() {
            self.bump(bump, si);
        }
        let mut b = arena.alloc_bsdf(si, 1.0);
        {
            let fourier: Arc<dyn BxDF> = Arc::new(FourierBSDF::new(&self.bsdf_table, mode));
            b.add(&fourier);
        }
        si.bsdf = Some(Arc::new(b));
    }
}

pub fn create_fourier_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let filename = mp.find_filename("bsdffile", "");
    let bsdf_table = FourierBSDFTable::read(&filename)?;
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    return Ok(Arc::new(FourierMaterial::new(&bsdf_table, &bumpmap)));
}
