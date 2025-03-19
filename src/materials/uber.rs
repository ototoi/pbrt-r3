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

pub struct UberMaterial {
    kd: Arc<dyn Texture<Spectrum>>,
    ks: Arc<dyn Texture<Spectrum>>,
    kr: Arc<dyn Texture<Spectrum>>,
    kt: Arc<dyn Texture<Spectrum>>,
    opacity: Arc<dyn Texture<Spectrum>>,
    roughness: Arc<dyn Texture<Float>>,
    u_roughness: Option<Arc<dyn Texture<Float>>>,
    v_roughness: Option<Arc<dyn Texture<Float>>>,
    eta: Arc<dyn Texture<Float>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
    remaproughness: bool,
}

impl UberMaterial {
    pub fn new(
        kd: &Arc<dyn Texture<Spectrum>>,
        ks: &Arc<dyn Texture<Spectrum>>,
        kr: &Arc<dyn Texture<Spectrum>>,
        kt: &Arc<dyn Texture<Spectrum>>,
        roughness: &Arc<dyn Texture<Float>>,
        u_roughness: &Option<Arc<dyn Texture<Float>>>,
        v_roughness: &Option<Arc<dyn Texture<Float>>>,
        opacity: &Arc<dyn Texture<Spectrum>>,
        eta: &Arc<dyn Texture<Float>>,
        bumpmap: &Option<Arc<dyn Texture<Float>>>,
        remaproughness: bool,
    ) -> Self {
        UberMaterial {
            kd: kd.clone(),
            ks: ks.clone(),
            kr: kr.clone(),
            kt: kt.clone(),
            opacity: opacity.clone(),
            roughness: roughness.clone(),
            u_roughness: u_roughness.clone(),
            v_roughness: v_roughness.clone(),
            eta: eta.clone(),
            bumpmap: bumpmap.clone(),
            remaproughness,
        }
    }

    fn evaluate_helper(
        a: &Option<Arc<dyn Texture<Float>>>,
        b: &Arc<dyn Texture<Float>>,
        si: &mut SurfaceInteraction,
    ) -> Float {
        if let Some(aa) = a.as_ref() {
            let tex = aa.as_ref();
            return tex.evaluate(si);
        } else {
            let tex = b.as_ref();
            return tex.evaluate(si);
        }
    }
}

impl Material for UberMaterial {
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

        let e = self.eta.as_ref().evaluate(si);
        let op = self.opacity.as_ref().evaluate(si);
        let t = (Spectrum::one() - op).clamp_zero();

        let ee = if !t.is_black() { 1.0 } else { e };

        let mut b = arena.alloc_bsdf(si, ee);

        if !t.is_black() {
            let tr: Arc<dyn BxDF> = Arc::new(SpecularTransmission::new(&t, 1.0, 1.0, mode));
            b.add(&tr);
        }

        let kd = op * self.kd.as_ref().evaluate(si).clamp_zero();
        if !kd.is_black() {
            let r: Arc<dyn BxDF> = Arc::new(LambertianReflection::new(&kd));
            b.add(&r);
        }

        let ks = op * self.ks.as_ref().evaluate(si).clamp_zero();

        if !ks.is_black() {
            let fresnel: Box<dyn Fresnel> = Box::new(FresnelDielectric::new(1.0, e));

            let mut u_rough = Self::evaluate_helper(&self.u_roughness, &self.roughness, si);
            let mut v_rough = Self::evaluate_helper(&self.v_roughness, &self.roughness, si);
            if self.remaproughness {
                u_rough = TrowbridgeReitzDistribution::roughness_to_alpha(u_rough);
                v_rough = TrowbridgeReitzDistribution::roughness_to_alpha(v_rough);
            }

            let distrib: Box<dyn MicrofacetDistribution> =
                Box::new(TrowbridgeReitzDistribution::new(u_rough, v_rough, true));
            let spec: Arc<dyn BxDF> = Arc::new(MicrofacetReflection::new(&ks, distrib, fresnel));
            b.add(&spec);
        }

        let kr = op * self.kr.as_ref().evaluate(si).clamp_zero();
        if !kr.is_black() {
            let fresnel: Box<dyn Fresnel> = Box::new(FresnelDielectric::new(1.0, e));
            let reflection: Arc<dyn BxDF> = Arc::new(SpecularReflection::new(&kr, fresnel));
            b.add(&reflection);
        }

        let kt = op * self.kt.as_ref().evaluate(si).clamp_zero();
        if !kt.is_black() {
            let trans: Arc<dyn BxDF> = Arc::new(SpecularTransmission::new(&kt, 1.0, e, mode));
            b.add(&trans);
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

pub fn create_uber_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let kd = mp.get_spectrum_texture("Kd", &Spectrum::from(0.25));
    let ks = mp.get_spectrum_texture("Ks", &Spectrum::from(0.25));
    let kr = mp.get_spectrum_texture("Kr", &Spectrum::from(0.0));
    let kt = mp.get_spectrum_texture("Kt", &Spectrum::from(0.0));
    let roughness = mp.get_float_texture("roughness", 0.1);
    let uroughness = mp.get_float_texture_or_null("uroughness");
    let vroughness = mp.get_float_texture_or_null("vroughness");
    let eta = get_float_texture_helper(mp, &["eta", "index"], 1.5);
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    let opacity = mp.get_spectrum_texture("opacity", &Spectrum::from(1.0));
    let remaproughness = mp.find_bool("remaproughness", true);

    return Ok(Arc::new(UberMaterial::new(
        &kd,
        &ks,
        &kr,
        &kt,
        &roughness,
        &uroughness,
        &vroughness,
        &opacity,
        &eta,
        &bumpmap,
        remaproughness,
    )));
}
