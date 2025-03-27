use crate::core::prelude::*;

use std::sync::Arc;

struct KdSubsurfaceMaterial {
    scale: Float,
    kd: Arc<dyn Texture<Spectrum>>,
    kr: Arc<dyn Texture<Spectrum>>,
    kt: Arc<dyn Texture<Spectrum>>,
    mfp: Arc<dyn Texture<Spectrum>>,
    eta: Float,
    u_roughness: Arc<dyn Texture<Float>>,
    v_roughness: Arc<dyn Texture<Float>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
    remaproughness: bool,
    table: Arc<BSSRDFTable>,
}

impl KdSubsurfaceMaterial {
    pub fn new(
        scale: Float,
        kd: &Arc<dyn Texture<Spectrum>>,
        kr: &Arc<dyn Texture<Spectrum>>,
        kt: &Arc<dyn Texture<Spectrum>>,
        mfp: &Arc<dyn Texture<Spectrum>>,
        g: Float,
        eta: Float,
        u_roughness: &Arc<dyn Texture<Float>>,
        v_roughness: &Arc<dyn Texture<Float>>,
        bumpmap: &Option<Arc<dyn Texture<Float>>>,
        remaproughness: bool,
    ) -> Self {
        let mut table = BSSRDFTable::new(100, 64);
        compute_beam_diffusion_bssrdf(g, eta, &mut table);

        KdSubsurfaceMaterial {
            scale,
            kd: kd.clone(),
            kr: kr.clone(),
            kt: kt.clone(),
            mfp: mfp.clone(),
            eta,
            u_roughness: u_roughness.clone(),
            v_roughness: v_roughness.clone(),
            bumpmap: bumpmap.clone(),
            remaproughness,
            table: Arc::new(table),
        }
    }
}

impl Material for KdSubsurfaceMaterial {
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
                        let fresnel: Box<dyn Fresnel> = Box::new(FresnelDielectric::new(1.0, eta));
                        if is_specular {
                            let reflection: Arc<dyn BxDF> =
                                Arc::new(SpecularReflection::new(&r, fresnel));
                            b.add(&reflection);
                        } else {
                            let distrib: Box<dyn MicrofacetDistribution> =
                                Box::new(TrowbridgeReitzDistribution::new(u_rough, v_rough, true));
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
                            //println!("tt!");
                            let distrib: Box<dyn MicrofacetDistribution> =
                                Box::new(TrowbridgeReitzDistribution::new(u_rough, v_rough, true));
                            let trans: Arc<dyn BxDF> =
                                Arc::new(MicrofacetTransmission::new(&t, distrib, 1.0, eta, mode));
                            b.add(&trans);
                        }
                    }
                }
            }
            si.bsdf = Some(Arc::new(b));
        }

        {
            let mfp = self.mfp.as_ref().evaluate(si);
            let mfree = mfp * self.scale;
            let kd = self.kd.as_ref().evaluate(si);

            let (sig_a, sig_s) = subsurface_from_diffuse(self.table.as_ref(), &kd, &mfree);

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
}

pub fn create_kdsubsurface_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    const KD: [Float; 3] = [0.5, 0.5, 0.5];

    let g = mp.find_float("g", 0.0);

    let scale = mp.find_float("scale", 1.0);
    let eta = mp.find_float("eta", 1.33);

    let kd = mp.get_spectrum_texture("Kd", &Spectrum::from(KD));
    let mfp = mp.get_spectrum_texture("mfp", &Spectrum::from(1.0));
    let kr = mp.get_spectrum_texture("Kr", &Spectrum::from(1.0));
    let kt = mp.get_spectrum_texture("Kt", &Spectrum::from(1.0));

    let uroughness = mp.get_float_texture("uroughness", 0.0);
    let vroughness = mp.get_float_texture("vroughness", 0.0);
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    let remaproughness = mp.find_bool("remaproughness", true);

    return Ok(Arc::new(KdSubsurfaceMaterial::new(
        scale,
        &kd,
        &kr,
        &kt,
        &mfp,
        g,
        eta,
        &uroughness,
        &vroughness,
        &bumpmap,
        remaproughness,
    )));
}
