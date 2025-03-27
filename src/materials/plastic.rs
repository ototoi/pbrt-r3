use crate::core::prelude::*;

use std::sync::Arc;

pub struct PlasticMaterial {
    kd: Arc<dyn Texture<Spectrum>>,
    ks: Arc<dyn Texture<Spectrum>>,
    roughness: Arc<dyn Texture<Float>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
    remaproughness: bool,
}

impl PlasticMaterial {
    pub fn new(
        kd: &Arc<dyn Texture<Spectrum>>,
        ks: &Arc<dyn Texture<Spectrum>>,
        roughness: &Arc<dyn Texture<Float>>,
        bumpmap: Option<Arc<dyn Texture<Float>>>,
        remaproughness: bool,
    ) -> Self {
        PlasticMaterial {
            kd: kd.clone(),
            ks: ks.clone(),
            roughness: roughness.clone(),
            bumpmap: bumpmap.clone(),
            remaproughness,
        }
    }
}

impl Material for PlasticMaterial {
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

        let kd = self.kd.as_ref();
        // Initialize diffuse component of plastic material
        let kd = kd.evaluate(si).clamp_zero();
        if !kd.is_black() {
            let r: Arc<dyn BxDF> = Arc::new(LambertianReflection::new(&kd));
            b.add(&r);
        }

        // Initialize specular component of plastic material
        let ks = self.ks.as_ref();
        let ks = ks.evaluate(si).clamp_zero();
        if !ks.is_black() {
            let roughness = self.roughness.as_ref();
            let mut rough = roughness.evaluate(si);
            if self.remaproughness {
                rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);
            }
            let fresnel: Box<dyn Fresnel> = Box::new(FresnelDielectric::new(1.5, 1.0));
            let distrib: Box<dyn MicrofacetDistribution> =
                Box::new(TrowbridgeReitzDistribution::new(rough, rough, true));
            let spec: Arc<dyn BxDF> = Arc::new(MicrofacetReflection::new(&ks, distrib, fresnel));
            b.add(&spec);
        }

        si.bsdf = Some(Arc::new(b));
    }
}

pub fn create_plastic_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let kd = mp.get_spectrum_texture("Kd", &Spectrum::from(0.25));
    let ks = mp.get_spectrum_texture("Ks", &Spectrum::from(0.25));
    let roughness = mp.get_float_texture("roughness", 0.1);
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    let remaproughness = mp.find_bool("remaproughness", true);
    return Ok(Arc::new(PlasticMaterial::new(
        &kd,
        &ks,
        &roughness,
        bumpmap,
        remaproughness,
    )));
}
