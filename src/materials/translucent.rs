use crate::core::camera::*;
use crate::core::distribution::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::lightdistrib::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use crate::core::texture::*;

use std::sync::Arc;

pub struct TranslucentMaterial {
    kd: Arc<dyn Texture<Spectrum>>,
    ks: Arc<dyn Texture<Spectrum>>,
    reflect: Arc<dyn Texture<Spectrum>>,
    transmit: Arc<dyn Texture<Spectrum>>,
    roughness: Arc<dyn Texture<Float>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
    remaproughness: bool,
}

impl TranslucentMaterial {
    pub fn new(
        kd: &Arc<dyn Texture<Spectrum>>,
        ks: &Arc<dyn Texture<Spectrum>>,
        reflect: &Arc<dyn Texture<Spectrum>>,
        transmit: &Arc<dyn Texture<Spectrum>>,
        roughness: &Arc<dyn Texture<Float>>,
        bumpmap: &Option<Arc<dyn Texture<Float>>>,
        remaproughness: bool,
    ) -> Self {
        TranslucentMaterial {
            kd: kd.clone(),
            ks: ks.clone(),
            reflect: reflect.clone(),
            transmit: transmit.clone(),
            roughness: roughness.clone(),
            bumpmap: bumpmap.clone(),
            remaproughness,
        }
    }
}

impl Material for TranslucentMaterial {
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
        let eta = 1.5;
        let reflect = self.reflect.as_ref();
        let transmit = self.transmit.as_ref();
        let r = reflect.evaluate(si).clamp_zero();
        let t = transmit.evaluate(si).clamp_zero();

        let mut b = arena.alloc_bsdf(si, eta);

        if r.is_black() && t.is_black() {
            return;
        }

        let kd = self.kd.as_ref();
        let kd = kd.evaluate(si).clamp_zero();
        if !kd.is_black() {
            if !r.is_black() {
                let refl: Arc<dyn BxDF> = Arc::new(LambertianReflection::new(&(r * kd)));
                b.add(&refl);
            }
            if !t.is_black() {
                let refl: Arc<dyn BxDF> = Arc::new(LambertianTransmission::new(&(t * kd)));
                b.add(&refl);
            }
        }

        let ks = self.ks.as_ref();
        let ks = ks.evaluate(si).clamp_zero();
        if !ks.is_black() {
            if !r.is_black() || !t.is_black() {
                let roughness = self.roughness.as_ref();
                let mut rough = roughness.evaluate(si);
                if self.remaproughness {
                    rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);
                }

                if !r.is_black() {
                    let distrib: Box<dyn MicrofacetDistribution> =
                        Box::new(TrowbridgeReitzDistribution::new(rough, rough, true));
                    let fresnel: Box<dyn Fresnel> = Box::new(FresnelDielectric::new(1.0, eta));
                    let rs: Arc<dyn BxDF> =
                        Arc::new(MicrofacetReflection::new(&(r * ks), distrib, fresnel));
                    b.add(&rs);
                }
                if !t.is_black() {
                    let distrib: Box<dyn MicrofacetDistribution> =
                        Box::new(TrowbridgeReitzDistribution::new(rough, rough, true));
                    let rt: Arc<dyn BxDF> = Arc::new(MicrofacetTransmission::new(
                        &(t * ks),
                        distrib,
                        1.0,
                        eta,
                        mode,
                    ));
                    b.add(&rt);
                }
            }
        }

        si.bsdf = Some(Arc::new(b));
    }
}

pub fn create_translucent_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let kd = mp.get_spectrum_texture("Kd", &Spectrum::from(0.25));
    let ks = mp.get_spectrum_texture("Ks", &Spectrum::from(0.25));
    let reflect = mp.get_spectrum_texture("reflect", &Spectrum::from(0.5));
    let transmit = mp.get_spectrum_texture("transmit", &Spectrum::from(0.5));
    let roughness = mp.get_float_texture("roughness", 0.1);
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    let remaproughness = mp.find_bool("remaproughness", true);
    return Ok(Arc::new(TranslucentMaterial::new(
        &kd,
        &ks,
        &reflect,
        &transmit,
        &roughness,
        &bumpmap,
        remaproughness,
    )));
}
