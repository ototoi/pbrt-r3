use crate::core::pbrt::*;
use std::sync::Arc;

#[inline]
fn sqr(x: Float) -> Float {
    x * x
}

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
// The Schlick Fresnel approximation is:
//
// R = R(0) + (1 - R(0)) (1 - cos theta)^5,
//
// where R(0) is the reflectance at normal indicence.
#[inline]
fn schlick_weight(cos_theta: Float) -> Float {
    let m = (1.0 - cos_theta).clamp(0.0, 1.0);
    return m * m * m * m * m; // (1 - cos_theta)^5
}

#[inline]
fn fr_schlick(r0: Float, cos_theta: Float) -> Float {
    return lerp(schlick_weight(cos_theta), r0, 1.0);
}

impl Spectrum {
    pub fn fr_schlick(&self, cos_theta: Float) -> Spectrum {
        let t = schlick_weight(cos_theta);
        return Spectrum::lerp(t, self, &Spectrum::one());
    }
}

// For a dielectric, R(0) = (eta - 1)^2 / (eta + 1)^2, assuming we're
// coming from air..
fn schlick_r0_from_eta(eta: Float) -> Float {
    return sqr(eta - 1.0) / sqr(eta + 1.0);
}

///////////////////////////////////////////////////////////////////////////
// DisneyDiffuse

pub struct DisneyDiffuse {
    r: Spectrum,
}
impl DisneyDiffuse {
    pub fn new(r: Spectrum) -> Self {
        Self { r }
    }
}
impl BxDF for DisneyDiffuse {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let fo = schlick_weight(abs_cos_theta(wo));
        let fi = schlick_weight(abs_cos_theta(wi));

        // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing.
        // Burley 2015, eq (4).
        return self.r * INV_PI * (1.0 - fo / 2.0) * (1.0 - fi / 2.0);
    }
    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        self.sample_f_default(wo, u)
    }
    fn rho(&self, _wo: &Vector3f, _samples: &[Point2f]) -> Spectrum {
        self.r
    }
    fn rho2(&self, _samples: &[(Point2f, Point2f)]) -> Spectrum {
        self.r
    }
    fn pdf(&self, _wo: &Vector3f, _wi: &Vector3f) -> Float {
        0.0
    }
    fn get_type(&self) -> BxDFType {
        BSDF_REFLECTION | BSDF_DIFFUSE
    }
}

///////////////////////////////////////////////////////////////////////////
// DisneyFakeSS

pub struct DisneyFakeSS {
    r: Spectrum,
    roughness: Float,
}
impl DisneyFakeSS {
    pub fn new(r: Spectrum, roughness: Float) -> Self {
        Self { r, roughness }
    }
}
impl BxDF for DisneyFakeSS {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let wh = *wo + *wi;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::zero();
        }
        let wh = wh.normalize();
        let cos_theta_d = wi.dot(&wh);

        // Fss90 used to "flatten" retroreflection based on roughness
        let fss90 = cos_theta_d * cos_theta_d * self.roughness;
        let fo = schlick_weight(abs_cos_theta(wo));
        let fi = schlick_weight(abs_cos_theta(wi));
        let fss = lerp(fo, 1.0, fss90) * lerp(fi, 1.0, fss90);
        let ss = 1.25 * fss * (1.0 / (abs_cos_theta(wo) + abs_cos_theta(wi)) - 0.5);
        return self.r * INV_PI * ss;
    }
    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        self.sample_f_default(wo, u)
    }
    fn rho(&self, _wo: &Vector3f, _samples: &[Point2f]) -> Spectrum {
        self.r
    }
    fn rho2(&self, _samples: &[(Point2f, Point2f)]) -> Spectrum {
        self.r
    }
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        return self.pdf_default(wo, wi);
    }
    fn get_type(&self) -> BxDFType {
        BSDF_REFLECTION | BSDF_DIFFUSE
    }
}

///////////////////////////////////////////////////////////////////////////
// DisneyRetro

pub struct DisneyRetro {
    r: Spectrum,
    roughness: Float,
}
impl DisneyRetro {
    pub fn new(r: Spectrum, roughness: Float) -> Self {
        Self { r, roughness }
    }
}
impl BxDF for DisneyRetro {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        unimplemented!()
    }
    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        self.sample_f_default(wo, u)
    }
    fn rho(&self, _wo: &Vector3f, _samples: &[Point2f]) -> Spectrum {
        self.r
    }
    fn rho2(&self, _samples: &[(Point2f, Point2f)]) -> Spectrum {
        self.r
    }
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        return self.pdf_default(wo, wi);
    }
    fn get_type(&self) -> BxDFType {
        BSDF_REFLECTION | BSDF_DIFFUSE
    }
}

///////////////////////////////////////////////////////////////////////////
// DisneySheen

///////////////////////////////////////////////////////////////////////////
// DisneyClearcoat

///////////////////////////////////////////////////////////////////////////
// DisneyFresnel

///////////////////////////////////////////////////////////////////////////
// DisneyMicrofacetDistribution

///////////////////////////////////////////////////////////////////////////
// DisneyBSSRDF

pub struct DisneyMaterial {
    color: Arc<dyn Texture<Spectrum>>,
    metallic: Arc<dyn Texture<Float>>,
    eta: Arc<dyn Texture<Float>>,
    roughness: Arc<dyn Texture<Float>>,
    specular_tint: Arc<dyn Texture<Float>>,
    anisotropic: Arc<dyn Texture<Float>>,
    sheen: Arc<dyn Texture<Float>>,
    sheen_tint: Arc<dyn Texture<Float>>,
    clearcoat: Arc<dyn Texture<Float>>,
    clearcoat_gloss: Arc<dyn Texture<Float>>,
    spectrans: Arc<dyn Texture<Float>>,
    scatter_distance: Arc<dyn Texture<Spectrum>>,
    thin: bool,
    flatness: Arc<dyn Texture<Float>>,
    diff_trans: Arc<dyn Texture<Float>>,
    bumpmap: Option<Arc<dyn Texture<Float>>>,
}

impl DisneyMaterial {
    pub fn new(
        color: &Arc<dyn Texture<Spectrum>>,
        metallic: &Arc<dyn Texture<Float>>,
        eta: &Arc<dyn Texture<Float>>,
        roughness: &Arc<dyn Texture<Float>>,
        specular_tint: &Arc<dyn Texture<Float>>,
        anisotropic: &Arc<dyn Texture<Float>>,
        sheen: &Arc<dyn Texture<Float>>,
        sheen_tint: &Arc<dyn Texture<Float>>,
        clearcoat: &Arc<dyn Texture<Float>>,
        clearcoat_gloss: &Arc<dyn Texture<Float>>,
        spectrans: &Arc<dyn Texture<Float>>,
        scatter_distance: &Arc<dyn Texture<Spectrum>>,
        thin: bool,
        flatness: &Arc<dyn Texture<Float>>,
        diff_trans: &Arc<dyn Texture<Float>>,
        bumpmap: &Option<Arc<dyn Texture<Float>>>,
    ) -> Self {
        DisneyMaterial {
            color: color.clone(),
            metallic: metallic.clone(),
            eta: eta.clone(),
            roughness: roughness.clone(),
            specular_tint: specular_tint.clone(),
            anisotropic: anisotropic.clone(),
            sheen: sheen.clone(),
            sheen_tint: sheen_tint.clone(),
            clearcoat: clearcoat.clone(),
            clearcoat_gloss: clearcoat_gloss.clone(),
            spectrans: spectrans.clone(),
            scatter_distance: scatter_distance.clone(),
            thin,
            flatness: flatness.clone(),
            diff_trans: diff_trans.clone(),
            bumpmap: bumpmap.clone(),
        }
    }
}

impl Material for DisneyMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        _: TransportMode,
        _: bool,
    ) {
        // Perform bump mapping with _bumpMap_, if present
        if let Some(bump) = self.bumpmap.as_ref() {
            self.bump(bump, si);
        }

        // Evaluate textures for _DisneyMaterial_ material and allocate BRDF
        let mut b = arena.alloc_bsdf(si, 1.0);

        // Diffuse
        let c = self.color.evaluate(si).clamp_zero();
        let metallic_weight = self.metallic.evaluate(si);
        let e = self.eta.evaluate(si);
        let strans = self.spectrans.evaluate(si);
        let diffuse_weight = (1.0 - 0.5 * metallic_weight) * (1.0 - strans);
        let dt = self.diff_trans.evaluate(si) * (1.0 / 2.0); // 0: all diffuse is reflected -> 1, transmitted
        let rough = self.roughness.evaluate(si);
        let lum = c.y();
        // normalize lum. to isolate hue+sat
        let c_tint = if lum > 0.0 {
            c * (1.0 / lum)
        } else {
            Spectrum::from(1.0)
        };

        // Sheen
        let sheen_weight = self.sheen.evaluate(si);
        let c_sheen = if sheen_weight > 0.0 {
            let stint = self.sheen_tint.evaluate(si);
            Spectrum::lerp(stint, &Spectrum::from(1.0), &c_tint)
        } else {
            Spectrum::from(0.0)
        };

        if diffuse_weight > 0.0 {
            if self.thin {
                let flat = self.flatness.evaluate(si);
                // Blend between DisneyDiffuse and fake subsurface based on
                // flatness.  Additionally, weight using diffTrans.
                // DisneyDiffuse
                // DisneyFakeSS
            } else {
                let sd = self.scatter_distance.evaluate(si);
                if sd.is_black() {
                    // DisneyDiffuse
                } else {
                    // SpecularTransmission
                    // DisneyBSSRDF
                }
            }

            // Retro-reflection.

            // Sheen (if enabled)
        }

        // Create the microfacet distribution for metallic and/or specular
        // transmission.
        let aspect = Float::sqrt(1.0 - self.anisotropic.evaluate(si) * 0.9);
        let ax = Float::max(0.001, rough * rough / aspect);
        let ay = Float::max(0.001, rough * rough * aspect);
        //let distrib = MicrofacetDistribution::new("ggx", ax, ay, true);

        // Specular is Trowbridge-Reitz with a modified Fresnel function.
        let spec_tint = self.specular_tint.evaluate(si);

        // Clearcoat

        // BTDF
        if strans > 0.0 {
            // SpecularTransmission
        }
        if self.thin {
            // Lambertian, weighted by (1 - diffTrans)
        }
        si.bsdf = Some(Arc::new(b));
    }
}

pub fn create_disney_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let color = mp.get_spectrum_texture("color", &Spectrum::from(0.5));
    let metallic = mp.get_float_texture("metallic", 0.0);
    let eta = mp.get_float_texture("eta", 1.5);
    let roughness = mp.get_float_texture("roughness", 0.5);
    let specular_tint = mp.get_float_texture("speculartint", 0.0);
    let anisotropic = mp.get_float_texture("anisotropic", 0.0);
    let sheen = mp.get_float_texture("sheen", 0.0);
    let sheen_tint = mp.get_float_texture("sheentint", 0.5);
    let clearcoat = mp.get_float_texture("clearcoat", 0.0);
    let clearcoat_gloss = mp.get_float_texture("clearcoatgloss", 1.0);
    let spectrans = mp.get_float_texture("spectrans", 0.0);
    let scatter_distance = mp.get_spectrum_texture("scatterdistance", &Spectrum::from(0.0));
    let thin = mp.find_bool("thin", false);
    let flatness = mp.get_float_texture("flatness", 0.0);
    let diff_trans = mp.get_float_texture("difftrans", 1.0);
    let bumpmap = mp.get_float_texture_or_null("bumpmap");
    return Ok(Arc::new(DisneyMaterial::new(
        &color,
        &metallic,
        &eta,
        &roughness,
        &specular_tint,
        &anisotropic,
        &sheen,
        &sheen_tint,
        &clearcoat,
        &clearcoat_gloss,
        &spectrans,
        &scatter_distance,
        thin,
        &flatness,
        &diff_trans,
        &bumpmap,
    )));
}
