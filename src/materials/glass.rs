use crate::core::base::*;
use crate::core::distribution::*;
use crate::core::error::*;
use crate::core::interaction::*;
use crate::core::material::*;
use crate::core::memory::*;
use crate::core::param_set::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

use std::sync::Arc;

pub struct GlassMaterial {
    kr: Arc<dyn Texture<Spectrum>>,
    kt: Arc<dyn Texture<Spectrum>>,
    u_roughness: Arc<dyn Texture<Float>>,
    v_roughness: Arc<dyn Texture<Float>>,
    index: Arc<dyn Texture<Float>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
    remaproughness: bool,
}

impl GlassMaterial {
    pub fn new(
        kr: &Arc<dyn Texture<Spectrum>>,
        kt: &Arc<dyn Texture<Spectrum>>,
        u_roughness: &Arc<dyn Texture<Float>>,
        v_roughness: &Arc<dyn Texture<Float>>,
        index: &Arc<dyn Texture<Float>>,
        bumpmap: &Option<Arc<dyn Texture<Float>>>,
        remaproughness: bool,
    ) -> Self {
        GlassMaterial {
            kr: kr.clone(),
            kt: kt.clone(),
            u_roughness: u_roughness.clone(),
            v_roughness: v_roughness.clone(),
            index: index.clone(),
            bumpmap: bumpmap.clone(),
            remaproughness,
        }
    }
}

impl Material for GlassMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        if let Some(bump) = self.bumpmap.as_ref() {
            self.bump(bump, si);
        }

        let eta = self.index.as_ref().evaluate(si);

        let mut u_rough = self.u_roughness.as_ref().evaluate(si);
        let mut v_rough = self.v_roughness.as_ref().evaluate(si);
        let r = self.kr.as_ref().evaluate(si);
        let t = self.kt.as_ref().evaluate(si);

        if r.is_black() && t.is_black() {
            return;
        }

        let mut b = arena.alloc_bsdf(si, eta);

        let is_specular = u_rough == 0.0 && v_rough == 0.0;
        if is_specular && allow_multiple_lobes {
            let reflection: Arc<dyn BxDF> = Arc::new(FresnelSpecular::new(&r, &t, 1.0, eta, mode));
            b.add(reflection);
        } else {
            if self.remaproughness {
                u_rough = TrowbridgeReitzDistribution::roughness_to_alpha(u_rough);
                v_rough = TrowbridgeReitzDistribution::roughness_to_alpha(v_rough);
            }
            if !r.is_black() {
                let fresnel: Box<dyn Fresnel> = Box::new(FresnelDielectric::new(1.0, eta));
                if is_specular {
                    let reflection: Arc<dyn BxDF> = Arc::new(SpecularReflection::new(&r, fresnel));
                    b.add(reflection);
                } else {
                    let distrib: Box<dyn MicrofacetDistribution> =
                        Box::new(TrowbridgeReitzDistribution::new(u_rough, v_rough, true));
                    let reflection: Arc<dyn BxDF> =
                        Arc::new(MicrofacetReflection::new(&r, distrib, fresnel));
                    b.add(reflection);
                }
            }
            if !t.is_black() {
                if is_specular {
                    let trans: Arc<dyn BxDF> =
                        Arc::new(SpecularTransmission::new(&t, 1.0, eta, mode));
                    b.add(trans);
                } else {
                    let distrib: Box<dyn MicrofacetDistribution> =
                        Box::new(TrowbridgeReitzDistribution::new(u_rough, v_rough, true));
                    let trans: Arc<dyn BxDF> =
                        Arc::new(MicrofacetTransmission::new(&t, distrib, 1.0, eta, mode));
                    b.add(trans);
                }
            }
        }

        si.bsdf = Some(Arc::new(b));
    }
}

fn get_float_texture_helper(
    mp: &TextureParams,
    keys: &[&str],
    value: Float,
) -> Arc<dyn Texture<Float>> {
    for key in keys {
        if let Some(tex) = mp.get_float_texture_or_null(key) {
            return tex;
        }
    }
    return mp.get_float_texture("", value);
}

pub fn create_glass_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let kr = mp.get_spectrum_texture("Kr", &Spectrum::from(1.0));
    let kt = mp.get_spectrum_texture("Kt", &Spectrum::from(1.0));
    let uroughness = mp.get_float_texture("uroughness", 0.0);
    let vroughness = mp.get_float_texture("vroughness", 0.0);
    let eta = get_float_texture_helper(mp, &["eta", "index"], 1.5);
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    let remaproughness = mp.find_bool("remaproughness", true);

    return Ok(Arc::new(GlassMaterial::new(
        &kr,
        &kt,
        &uroughness,
        &vroughness,
        &eta,
        &bumpmap,
        remaproughness,
    )));
}
