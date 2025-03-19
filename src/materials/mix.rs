use crate::core::error::*;
use crate::core::interaction::*;
use crate::core::material::*;
use crate::core::memory::*;
use crate::core::param_set::*;
use crate::core::refrection::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

use std::sync::Arc;

pub struct MixMaterial {
    m1: Arc<dyn Material>,
    m2: Arc<dyn Material>,
    scale: Arc<dyn Texture<Spectrum>>,
}

impl MixMaterial {
    pub fn new(
        m1: &Arc<dyn Material>,
        m2: &Arc<dyn Material>,
        scale: &Arc<dyn Texture<Spectrum>>,
    ) -> Self {
        Self {
            m1: Arc::clone(m1),
            m2: Arc::clone(m2),
            scale: Arc::clone(scale),
        }
    }
}

fn scaled_bsdf(b1: &BSDF, b2: &BSDF, s1: &Spectrum, s2: &Spectrum, b: &mut BSDF) {
    assert!(b.bxdfs.len() == 0);
    {
        //let n = b1.num_components(BSDF_ALL) as usize;
        //assert!(b1.bxdfs.len() == n);
        let n = b1.bxdfs.len();
        for i in 0..n {
            let new_bxdf: Arc<dyn BxDF> = Arc::new(ScaledBxDF::new(&b1.bxdfs[i], &s1));
            b.add(&new_bxdf);
        }
    }
    {
        //let n = b2.num_components(BSDF_ALL) as usize;
        //assert!(b2.bxdfs.len() == n);
        let n = b2.bxdfs.len();
        for i in 0..n {
            let new_bxdf: Arc<dyn BxDF> = Arc::new(ScaledBxDF::new(&b2.bxdfs[i], &s2));
            b.add(&new_bxdf);
        }
    }
    {
        let n1 = b1.bxdfs.len();
        let n2 = b2.bxdfs.len();
        let nn = b.bxdfs.len();
        assert!(nn == n1 + n2);
    }
}

impl Material for MixMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        // Compute weights and original _BxDF_s for mix material
        let s1 = self.scale.as_ref().evaluate(si).clamp_zero();
        let s2 = (Spectrum::one() - s1).clamp_zero();

        let mut si2 = si.clone();

        self.m1
            .as_ref()
            .compute_scattering_functions(si, arena, mode, allow_multiple_lobes);
        self.m2
            .as_ref()
            .compute_scattering_functions(&mut si2, arena, mode, allow_multiple_lobes);

        assert!(si.bsdf.is_some());
        assert!(si2.bsdf.is_some());

        if let Some(bsdf1) = si.bsdf.as_ref() {
            if let Some(bsdf2) = si2.bsdf.as_ref() {
                let mut b = arena.alloc_bsdf(si, bsdf1.eta);
                let bsdf1 = bsdf1.as_ref();
                let bsdf2 = bsdf2.as_ref();
                scaled_bsdf(bsdf1, bsdf2, &s1, &s2, &mut b);
                si.bsdf = Some(Arc::new(b));
            }
        }
    }
}

pub fn create_mix_material(
    mp: &TextureParams,
    m1: &Arc<dyn Material>,
    m2: &Arc<dyn Material>,
) -> Result<Arc<dyn Material>, PbrtError> {
    let scale = mp.get_spectrum_texture("amount", &Spectrum::from(0.5));
    return Ok(Arc::new(MixMaterial::new(m1, m2, &scale)));
}
