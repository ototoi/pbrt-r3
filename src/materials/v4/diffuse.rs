use super::bxdfs::*;
use crate::core::prelude::*;
use std::sync::Arc;

pub struct DiffuseMaterial {
    reflectance: Arc<dyn Texture<Spectrum>>,
    displacement: Option<Arc<dyn Texture<Float>>>,
}

impl DiffuseMaterial {
    pub fn new(
        reflectance: &Arc<dyn Texture<Spectrum>>,
        displacement: &Option<Arc<dyn Texture<Float>>>,
    ) -> Self {
        DiffuseMaterial {
            reflectance: reflectance.clone(),
            displacement: displacement.clone(),
        }
    }
}

impl Material for DiffuseMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        _: TransportMode,
        _: bool,
    ) {
        if let Some(bump) = self.displacement.as_ref() {
            self.bump(bump, si);
        }

        let mut b = arena.alloc_bsdf(si, 1.0);

        let r = self.reflectance.as_ref().evaluate(si);
        if !r.is_black() {
            let r: Arc<dyn BxDF> = Arc::new(DiffuseBxDF::new(&r));
            b.add(&r);
        }

        si.bsdf = Some(Arc::new(b));
    }
}

pub fn create_diffuse_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let reflectance = mp.get_spectrum_texture("reflectance", &Spectrum::from(0.5));
    let displacement = mp.get_float_texture_or_null("displacement");
    return Ok(Arc::new(DiffuseMaterial::new(&reflectance, &displacement)));
}
