use crate::core::prelude::*;

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

fn scaled_bsdf(mut b1: BSDF, mut b2: BSDF, s1: &Spectrum, s2: &Spectrum, b: &mut BSDF) {
    assert!(b.bxdfs.len() == 0);
    {
        let n = b1.bxdfs.len();
        for bxdf in b1.bxdfs.drain(..) {
            b.add(ScaledBxDF::new(bxdf, s1));
        }
        debug_assert_eq!(n, b.bxdfs.len());
    }
    {
        let n_before = b.bxdfs.len();
        let n = b2.bxdfs.len();
        for bxdf in b2.bxdfs.drain(..) {
            b.add(ScaledBxDF::new(bxdf, s2));
        }
        debug_assert_eq!(n, b.bxdfs.len() - n_before);
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

        // Save the original interaction state so m2 can be evaluated from the
        // same input as m1 without cloning the entire SurfaceInteraction.
        let n0 = si.n;
        let wo0 = si.wo;
        let shading0 = si.shading;

        self.m1
            .as_ref()
            .compute_scattering_functions(si, arena, mode, allow_multiple_lobes);

        assert!(si.bsdf.is_some());
        let bsdf1 = Arc::try_unwrap(si.bsdf.take().unwrap())
            .expect("MixMaterial expects unique BSDF ownership for m1");

        // Preserve m1's resulting shading and optional BSSRDF; original code
        // returned m1-mutated interaction state.
        let n1 = si.n;
        let wo1 = si.wo;
        let shading1 = si.shading;
        let bssrdf1 = si.bssrdf.take();

        // Reset only fields materials may mutate before evaluating m2.
        si.n = n0;
        si.wo = wo0;
        si.shading = shading0;
        si.bsdf = None;
        si.bssrdf = None;

        self.m2
            .as_ref()
            .compute_scattering_functions(si, arena, mode, allow_multiple_lobes);

        assert!(si.bsdf.is_some());
        let bsdf2 = Arc::try_unwrap(si.bsdf.take().unwrap())
            .expect("MixMaterial expects unique BSDF ownership for m2");

        let mut b = arena.alloc_bsdf(si, bsdf1.eta);
        scaled_bsdf(bsdf1, bsdf2, &s1, &s2, &mut b);

        // Restore m1 interaction side effects to keep behavior consistent.
        si.n = n1;
        si.wo = wo1;
        si.shading = shading1;
        si.bssrdf = bssrdf1;
        si.bsdf = Some(Arc::new(b));
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
