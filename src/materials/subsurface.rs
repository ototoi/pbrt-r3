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
use crate::media::*;

use std::sync::Arc;

use log::*;

struct SubsurfaceMaterial {
    scale: Float,
    kr: Arc<dyn Texture<Spectrum>>,
    kt: Arc<dyn Texture<Spectrum>>,
    sigma_a: Arc<dyn Texture<Spectrum>>,
    sigma_s: Arc<dyn Texture<Spectrum>>,
    eta: Float,
    u_roughness: Arc<dyn Texture<Float>>,
    v_roughness: Arc<dyn Texture<Float>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
    remaproughness: bool,
    table: Arc<BSSRDFTable>,
}

impl SubsurfaceMaterial {
    pub fn new(
        scale: Float,
        kr: &Arc<dyn Texture<Spectrum>>,
        kt: &Arc<dyn Texture<Spectrum>>,
        sigma_a: &Arc<dyn Texture<Spectrum>>,
        sigma_s: &Arc<dyn Texture<Spectrum>>,
        g: Float,
        eta: Float,
        u_roughness: &Arc<dyn Texture<Float>>,
        v_roughness: &Arc<dyn Texture<Float>>,
        bumpmap: &Option<Arc<dyn Texture<Float>>>,
        remaproughness: bool,
    ) -> Self {
        let mut table = BSSRDFTable::new(100, 64);
        compute_beam_diffusion_bssrdf(g, eta, &mut table);

        SubsurfaceMaterial {
            scale,
            kr: kr.clone(),
            kt: kt.clone(),
            sigma_a: sigma_a.clone(),
            sigma_s: sigma_s.clone(),
            eta,
            u_roughness: u_roughness.clone(),
            v_roughness: v_roughness.clone(),
            bumpmap: bumpmap.clone(),
            remaproughness,
            table: Arc::new(table),
        }
    }
}

impl Material for SubsurfaceMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        // Perform bump mapping with _bumpMap_, if present
        if let Some(bump) = self.bumpmap.as_ref() {
            self.bump(bump, si);
        }

        // Initialize BSDF for _SubsurfaceMaterial_
        let eta = self.eta;
        let r = self.kr.as_ref().evaluate(si);
        let t = self.kt.as_ref().evaluate(si);

        {
            let mut b = arena.alloc_bsdf(si, eta);
            if !r.is_black() || !t.is_black() {
                {
                    let mut u_rough = self.u_roughness.as_ref().evaluate(si);
                    let mut v_rough = self.v_roughness.as_ref().evaluate(si);
                    let is_specular = u_rough == 0.0 && v_rough == 0.0;
                    if is_specular && allow_multiple_lobes {
                        let reflection: Arc<dyn BxDF> =
                            Arc::new(FresnelSpecular::new(&r, &t, 1.0, eta, mode));
                        b.add(&reflection);
                    } else {
                        if self.remaproughness {
                            u_rough = TrowbridgeReitzDistribution::roughness_to_alpha(u_rough);
                            v_rough = TrowbridgeReitzDistribution::roughness_to_alpha(v_rough);
                        }
                        if !r.is_black() {
                            let fresnel: Box<dyn Fresnel> =
                                Box::new(FresnelDielectric::new(1.0, eta));
                            if is_specular {
                                let reflection: Arc<dyn BxDF> =
                                    Arc::new(SpecularReflection::new(&r, fresnel));
                                b.add(&reflection);
                            } else {
                                let distrib: Box<dyn MicrofacetDistribution> = Box::new(
                                    TrowbridgeReitzDistribution::new(u_rough, v_rough, true),
                                );
                                let reflection: Arc<dyn BxDF> =
                                    Arc::new(MicrofacetReflection::new(&r, distrib, fresnel));
                                b.add(&reflection);
                            }
                        }
                        if !t.is_black() {
                            if is_specular {
                                let trans: Arc<dyn BxDF> =
                                    Arc::new(SpecularTransmission::new(&t, 1.0, eta, mode));
                                b.add(&trans);
                            } else {
                                let distrib: Box<dyn MicrofacetDistribution> = Box::new(
                                    TrowbridgeReitzDistribution::new(u_rough, v_rough, true),
                                );
                                let trans: Arc<dyn BxDF> = Arc::new(MicrofacetTransmission::new(
                                    &t, distrib, 1.0, eta, mode,
                                ));
                                b.add(&trans);
                            }
                        }
                    }
                }

                {
                    let scale = self.scale;
                    let sig_a = self.sigma_a.as_ref().evaluate(si).clamp_zero() * scale;
                    let sig_s = self.sigma_s.as_ref().evaluate(si).clamp_zero() * scale;

                    let material = self as BSSRDFMaterialRawPointer;
                    let bssrdf: Arc<dyn BSSRDF> = Arc::new(TabulatedBSSRDF::new(
                        si,
                        eta,
                        material,
                        mode,
                        &sig_a,
                        &sig_s,
                        &self.table,
                    ));
                    si.bssrdf = Some(bssrdf);
                }
            }
            si.bsdf = Some(Arc::new(b));
        }
    }
}

pub fn create_subsurface_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    const SIG_A_RGB: [Float; 3] = [0.0011, 0.0024, 0.014];
    const SIG_S_RGB: [Float; 3] = [2.55, 3.21, 3.77];
    let mut sig_a = Spectrum::from(SIG_A_RGB);
    let mut sig_s = Spectrum::from(SIG_S_RGB);

    let mut g = mp.find_float("g", 0.0);
    let name = mp.find_string("name", "");
    if !name.is_empty() {
        if let Some((a, s)) = get_medium_scattering_properties(&name) {
            sig_a = a;
            sig_s = s;
            //Enforce g=0 (the database specifies reduced scattering coefficients)
            g = 0.0;
        } else {
            warn!("Named material \"{}\" not found.  Using defaults.", name);
        }
    }

    let scale = mp.find_float("scale", 1.0);
    let eta = mp.find_float("eta", 1.33);

    let sigma_a = mp.get_spectrum_texture("sigma_a", &sig_a);
    let sigma_s = mp.get_spectrum_texture("sigma_s", &sig_s);

    let kr = mp.get_spectrum_texture("Kr", &Spectrum::from(1.0));
    let kt = mp.get_spectrum_texture("Kt", &Spectrum::from(1.0));

    let uroughness = mp.get_float_texture("uroughness", 0.0);
    let vroughness = mp.get_float_texture("vroughness", 0.0);
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    let remaproughness = mp.find_bool("remaproughness", true);

    return Ok(Arc::new(SubsurfaceMaterial::new(
        scale,
        &kr,
        &kt,
        &sigma_a,
        &sigma_s,
        g,
        eta,
        &uroughness,
        &vroughness,
        &bumpmap,
        remaproughness,
    )));
}
