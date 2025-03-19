use crate::core::distribution::*;
use crate::core::error::*;
use crate::core::interaction::*;
use crate::core::material::*;
use crate::core::memory::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::refrection::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

use std::sync::Arc;

pub struct SubstrateMaterial {
    kd: Arc<dyn Texture<Spectrum>>,
    ks: Arc<dyn Texture<Spectrum>>,
    nu: Arc<dyn Texture<Float>>,
    nv: Arc<dyn Texture<Float>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
    remaproughness: bool,
}

impl SubstrateMaterial {
    pub fn new(
        kd: &Arc<dyn Texture<Spectrum>>,
        ks: &Arc<dyn Texture<Spectrum>>,
        u_roughness: &Arc<dyn Texture<Float>>,
        v_roughness: &Arc<dyn Texture<Float>>,
        bumpmap: &Option<Arc<dyn Texture<Float>>>,
        remaproughness: bool,
    ) -> Self {
        SubstrateMaterial {
            kd: kd.clone(),
            ks: ks.clone(),
            nu: u_roughness.clone(),
            nv: v_roughness.clone(),
            bumpmap: bumpmap.clone(),
            remaproughness,
        }
    }
}

impl Material for SubstrateMaterial {
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

        let d = self.kd.as_ref().evaluate(si).clamp_zero();
        let s = self.ks.as_ref().evaluate(si).clamp_zero();

        let mut b = arena.alloc_bsdf(si, 1.0);

        if !d.is_black() && !s.is_black() {
            let mut u_rough = self.nu.as_ref().evaluate(si);
            let mut v_rough = self.nv.as_ref().evaluate(si);

            if self.remaproughness {
                u_rough = TrowbridgeReitzDistribution::roughness_to_alpha(u_rough);
                v_rough = TrowbridgeReitzDistribution::roughness_to_alpha(v_rough);
            }

            let distrib: Box<dyn MicrofacetDistribution> =
                Box::new(TrowbridgeReitzDistribution::new(u_rough, v_rough, true));
            let spec: Arc<dyn BxDF> = Arc::new(FresnelBlend::new(&d, &s, distrib));
            b.add(&spec);
        }

        si.bsdf = Some(Arc::new(b));
    }
}

pub fn create_substrate_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let kd = mp.get_spectrum_texture("Kd", &Spectrum::from(0.5));
    let ks = mp.get_spectrum_texture("Ks", &Spectrum::from(0.5));
    let uroughness = mp.get_float_texture("uroughness", 0.1);
    let vroughness = mp.get_float_texture("vroughness", 0.1);
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    let remaproughness = mp.find_bool("remaproughness", true);

    return Ok(Arc::new(SubstrateMaterial::new(
        &kd,
        &ks,
        &uroughness,
        &vroughness,
        &bumpmap,
        remaproughness,
    )));
}
