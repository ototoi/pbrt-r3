use crate::core::prelude::*;

use std::sync::Arc;

pub struct MirrorMaterial {
    kr: Arc<dyn Texture<Spectrum>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
}

impl MirrorMaterial {
    pub fn new(kr: &Arc<dyn Texture<Spectrum>>, bumpmap: &Option<Arc<dyn Texture<Float>>>) -> Self {
        MirrorMaterial {
            kr: kr.clone(),
            bumpmap: bumpmap.clone(),
        }
    }
}

impl Material for MirrorMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        if let Some(bump) = self.bumpmap.as_ref() {
            self.bump(bump, si);
        }
        let mut b = arena.alloc_bsdf(si, 1.0);

        let r = self.kr.as_ref().evaluate(si).clamp_zero();
        if !r.is_black() {
            let fresnel: Box<dyn Fresnel> = Box::new(FresnelNoOp::new());
            let r: Arc<dyn BxDF> = Arc::new(SpecularReflection::new(&r, fresnel));
            b.add(r);
        }

        si.bsdf = Some(Arc::new(b));
    }
}

pub fn create_mirror_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let kr = mp.get_spectrum_texture("Kr", &Spectrum::from(0.9));
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    return Ok(Arc::new(MirrorMaterial::new(&kr, &bumpmap)));
}
