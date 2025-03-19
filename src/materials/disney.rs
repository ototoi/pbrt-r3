use crate::core::distribution::*;
use crate::core::error::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::material::*;
use crate::core::memory::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::refrection::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use crate::core::texture::*;

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
        let wh = *wo + *wi;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::zero();
        }
        let wh = wh.normalize();
        let cos_theta_d = wi.dot(&wh);

        let fo = schlick_weight(abs_cos_theta(wo));
        let fi = schlick_weight(abs_cos_theta(wi));
        let rr = 2.0 * self.roughness * cos_theta_d * cos_theta_d;

        // Burley 2015, eq (4).
        return self.r * INV_PI * rr * (fo + fi + fo * fi * (rr - 1.0));
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

pub struct DisneySheen {
    r: Spectrum,
}
impl DisneySheen {
    pub fn new(r: Spectrum) -> Self {
        Self { r }
    }
}
impl BxDF for DisneySheen {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let wh = *wo + *wi;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::zero();
        }
        let wh = wh.normalize();
        let cos_theta_d = wi.dot(&wh);

        return self.r * schlick_weight(cos_theta_d);
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
// DisneyClearcoat

pub struct DisneyClearcoat {
    weight: Float,
    gloss: Float,
}
impl DisneyClearcoat {
    pub fn new(weight: Float, gloss: Float) -> Self {
        Self { weight, gloss }
    }
}

#[inline]
fn gtr1(cos_theta: Float, alpha: Float) -> Float {
    let alpha2 = alpha * alpha;
    let cos_theta2 = cos_theta * cos_theta;
    return (alpha2 - 1.0) / (PI * Float::ln(alpha2) * (1.0 + (alpha2 - 1.0) * cos_theta2));
}

// Smith masking/shadowing term.
#[inline]
fn smith_g_ggx(cos_theta: Float, alpha: Float) -> Float {
    let alpha2 = alpha * alpha;
    let cos_theta2 = cos_theta * cos_theta;
    return 1.0 / (cos_theta + Float::sqrt(alpha2 + cos_theta2 - alpha2 * cos_theta2));
}

impl BxDF for DisneyClearcoat {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let wh = *wo + *wi;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::zero();
        }
        let wh = wh.normalize();

        // Clearcoat has ior = 1.5 hardcoded -> F0 = 0.04. It then uses the
        // GTR1 distribution, which has even fatter tails than Trowbridge-Reitz
        // (which is GTR2).
        let dr = gtr1(abs_cos_theta(&wh), self.gloss);
        let fr = fr_schlick(0.04, Vector3f::dot(wo, &wh));
        // The geometric term always based on alpha = 0.25.
        let gr = smith_g_ggx(abs_cos_theta(wo), 0.25) * smith_g_ggx(abs_cos_theta(wi), 0.25);
        return (self.weight * gr * fr * dr / 4.0).into();
    }
    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        // TODO: double check all this: there still seem to be some very
        // occasional fireflies with clearcoat; presumably there is a bug
        // somewhere.
        if wo.z == 0.0 {
            return None;
        }

        let alpha2 = self.gloss * self.gloss;
        let cos_theta = Float::sqrt(Float::max(
            0.0,
            1.0 - Float::powf(alpha2, 1.0 - u[0]) / (1.0 - alpha2),
        ));
        let sin_theta = Float::sqrt(Float::max(0.0, 1.0 - cos_theta * cos_theta));
        let phi = 2.0 * PI * u[1];
        let mut wh = spherical_direction(sin_theta, cos_theta, phi);
        if !same_hemisphere(wo, &wh) {
            wh = -wh;
        }
        let wi = reflect(wo, &wh);
        if !same_hemisphere(wo, &wi) {
            return None;
        }
        let pdf = self.pdf(wo, &wi);
        let f = self.f(wo, &wi);
        let t = self.get_type();
        return Some((f, wi, pdf, t));
    }
    // fn rho(&self, wo: &Vector3f, samples: &[Point2f]) -> Spectrum;
    // fn rho2(&self, _samples: &[(Point2f, Point2f)]) -> Spectrum;
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !same_hemisphere(wo, wi) {
            return 0.0;
        }
        let wh = *wo + *wi;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return 0.0;
        }
        let wh = wh.normalize();

        // The sampling routine samples wh exactly from the GTR1 distribution.
        // Thus, the final value of the PDF is just the value of the
        // distribution for wh converted to a mesure with respect to the
        // surface normal.
        let dr = gtr1(abs_cos_theta(&wh), self.gloss);
        return dr + abs_cos_theta(&wh) / (4.0 * abs_cos_theta(wo));
    }
    fn get_type(&self) -> BxDFType {
        BSDF_REFLECTION | BSDF_GLOSSY
    }
}

///////////////////////////////////////////////////////////////////////////
// DisneyFresnel

// Specialized Fresnel function used for the specular component, based on
// a mixture between dielectric and the Schlick Fresnel approximation.
pub struct DisneyFresnel {
    r0: Spectrum,
    metallic: Float,
    eta: Float,
}
impl DisneyFresnel {
    pub fn new(r0: Spectrum, metallic: Float, eta: Float) -> Self {
        Self { r0, metallic, eta }
    }
}
impl Fresnel for DisneyFresnel {
    fn evaluate(&self, cos_i: Float) -> Spectrum {
        let r0 = lerp(self.metallic, schlick_r0_from_eta(self.eta), self.r0.y());
        return self.r0 + (Spectrum::from(r0) - self.r0) * schlick_weight(cos_i);
    }
}

///////////////////////////////////////////////////////////////////////////
// DisneyMicrofacetDistribution

#[derive(Clone, Debug)]
pub struct DisneyMicrofacetDistribution {
    base: TrowbridgeReitzDistribution,
}
impl DisneyMicrofacetDistribution {
    pub fn new(alpha_x: Float, alpha_y: Float) -> Self {
        Self {
            base: TrowbridgeReitzDistribution::new(alpha_x, alpha_y, true),
        }
    }
}
impl MicrofacetDistribution for DisneyMicrofacetDistribution {
    fn d(&self, wh: &Vector3f) -> Float {
        self.base.d(wh)
    }
    fn lambda(&self, w: &Vector3f) -> Float {
        self.base.lambda(w)
    }
    fn g(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        // Disney uses the separable masking-shadowing model.
        return self.base.g1(wo) * self.base.g1(wi);
    }
    fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f {
        self.base.sample_wh(wo, u)
    }
    fn pdf(&self, wo: &Vector3f, wh: &Vector3f) -> Float {
        self.base.pdf(wo, wh)
    }
}

///////////////////////////////////////////////////////////////////////////
// DisneyBSSRDF

// Implementation of the empirical BSSRDF described in "Extending the
// Disney BRDF to a BSDF with integrated subsurface scattering" (Brent
// Burley) and "Approximate Reflectance Profiles for Efficient Subsurface
// Scattering (Christensen and Burley).
#[derive(Clone, Debug)]
pub struct DisneyBSSRDF {
    base: Arc<BaseSeparableBSSRDF>,
    r: Spectrum,
    d: Spectrum,
}
impl DisneyBSSRDF {
    pub fn new(
        r: Spectrum,
        d: Spectrum,
        po: &SurfaceInteraction,
        eta: Float,
        material: BSSRDFMaterialRawPointer,
        mode: TransportMode,
    ) -> Self {
        // 0.2 factor comes from personal communication from Brent Burley
        // and Matt Chiang.
        let d = d * 0.2;
        let base = Arc::new(BaseSeparableBSSRDF::new(po, eta, material, mode));
        Self { base, r, d }
    }

    fn sr(&self, r: Float) -> Spectrum {
        let r = r.max(1e-6); // Avoid singularity at r == 0.
        let d = self.d;
        return self.r * Spectrum::exp(&(Spectrum::from(-r) / d))
            + Spectrum::exp(&(Spectrum::from(-r) / (3.0 * d))) / (8.0 * PI * d * r);
    }

    fn pdf_sr(&self, ch: usize, r: Float) -> Float {
        let r = r.max(1e-6); // Avoid singularity at r == 0.
                             // Weight the two individual PDFs as per the sampling frequency in
                             // Sample_Sr().
        let d = self.d;
        return 0.25 * Float::exp(-r / d[ch]) / (2.0 * PI * d[ch] * r)
            + 0.75 * Float::exp(-r / (3.0 * d[ch])) / (6.0 * PI * d[ch] * r);
    }
}

impl BSSRDF for DisneyBSSRDF {
    fn s(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum {
        // Fade based on relative orientations of the two surface normals to
        // better handle surface cavities. (Details via personal communication
        // from Brent Burley; these details aren't published in the course
        // notes.)
        //
        // TODO: test
        // TODO: explain
        let po = &self.base.base;
        let a = (pi.p - po.p).normalize();
        let mut fade = 1.0;
        let n = self.base.ns;
        let cos_theta = a.dot(&n);
        if cos_theta > 0.0 {
            // Point on or above surface plane
            let sin_theta = Float::sqrt(Float::max(0.0, 1.0 - cos_theta * cos_theta));
            let a2 = n * sin_theta - (a - n * cos_theta) * (cos_theta / sin_theta);
            fade = a2.dot(&pi.shading.n).max(0.0);
        }

        let fo = schlick_weight(abs_cos_theta(&po.wo));
        let fi = schlick_weight(abs_cos_theta(wi));
        return fade * (1.0 - fo / 2.0) * (1.0 - fi / 2.0) / PI * self.sp(pi);
    }

    fn sample_s(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        arena: &mut MemoryArena,
    ) -> Option<(Spectrum, SurfaceInteraction, Float)> {
        if let Some((sp, mut si, pdf)) = self.sample_sp(scene, u1, u2, arena) {
            if !sp.is_black() {
                let mut b = arena.alloc_bsdf(&si, 1.0);
                let adapter: Arc<dyn BxDF> = Arc::new(SeparableBSSRDFAdapter::new(
                    self.base.base.eta,
                    self.base.mode,
                ));
                b.add(&adapter);
                let bsdf = Arc::new(b);
                si.bsdf = Some(bsdf);
                si.wo = si.shading.n;
            }
            return Some((sp, si, pdf));
        }
        return None;
    }
}

impl SeparableBSSRDF for DisneyBSSRDF {
    fn projection_axis(&self, u1: Float) -> (Vector3f, Vector3f, Vector3f, Float) {
        self.base.projection_axis(u1)
    }

    fn get_p(&self) -> Point3f {
        self.base.base.p
    }

    fn get_time(&self) -> Float {
        self.base.base.time
    }

    fn get_material(&self) -> BSSRDFMaterialRawPointer {
        self.base.material
    }

    fn sp(&self, pi: &SurfaceInteraction) -> Spectrum {
        let p = self.base.base.p;
        return self.sr(Vector3f::distance(&p, &pi.p));
    }

    fn sample_sr(&self, ch: usize, u: Float) -> Float {
        // The good news is that diffusion profile implemented in Sr is
        // normalized---integrating in polar coordinates, we have:
        //
        // int_0^2pi int_0^Infinity Sr(r) r dr dphi == 1.
        //
        // The CDF can be found in closed-form. It is:
        //
        // 1 - e^(-x/d) / 4 - (3 / 4) e^(-x / (3d)).
        //
        // Unfortunately, inverting the CDF requires solving a cubic, which
        // would be nice to sidestep. Therefore, following Christensen and
        // Burley's suggestion (section 6), we will sample from each of the two
        // exponential terms individually (which can be done directly) and then
        // compute an overall PDF using MIS.  There are a few details to work
        // through...
        //
        // For the first exponential term, we can find:
        // normalized PDF: e^(-r/d) / (2 Pi d r)
        // CDF: 1 - e^(-r/d)
        // sampling recipe: r = d log(1 / (1 - u))
        //
        // For the second:
        // PDF: e^(-r/(3d)) / (6 Pi d r)
        // CDF: 1 - e^(-r/(3d))
        // sampling: r = 3 d log(1 / (1 - u))
        //
        // The last question is what fraction of samples to use for each
        // technique.  The second exponential has 3x the contribution to the
        // final value as the first does, so therefore we'll take three samples
        // from that for every one sample we take from the first.
        if u <= 0.25 {
            // Sample the first exponential
            let u = Float::min(u * 4.0, ONE_MINUS_EPSILON);
            return self.d[ch] * Float::ln(1.0 / (1.0 - u));
        } else {
            // Second exponenital
            let u = Float::min((u - 0.25) / 0.75, ONE_MINUS_EPSILON);
            return 3.0 * self.d[ch] * Float::ln(1.0 / (1.0 - u));
        }
    }

    fn pdf_sp(&self, pi: &SurfaceInteraction) -> Float {
        // Express $\pti-\pto$ and $\bold{n}_i$ with respect to local coordinates at
        // $\pto$
        let po = &self.base.base;
        let d = po.p - pi.p;

        assert!(pi.n.length() > 0.0);

        let d_local = Vector3f::new(
            Vector3f::dot(&self.base.ss, &d),
            Vector3f::dot(&self.base.ts, &d),
            Vector3f::dot(&self.base.ns, &d),
        );
        let n_local = Vector3f::new(
            Vector3f::dot(&self.base.ss, &pi.n),
            Vector3f::dot(&self.base.ts, &pi.n),
            Vector3f::dot(&self.base.ns, &pi.n),
        );
        // Compute BSSRDF profile radius under projection along each axis
        let r_proj = [
            Float::sqrt(d_local.y * d_local.y + d_local.z * d_local.z),
            Float::sqrt(d_local.z * d_local.z + d_local.x * d_local.x),
            Float::sqrt(d_local.x * d_local.x + d_local.y * d_local.y),
        ];

        // Return combined probability from all BSSRDF sampling strategies
        let mut pdf = 0.0;
        let axis_prob = [0.25, 0.25, 0.5];
        let n_samples = Spectrum::N_SAMPLES;
        let ch_prob = 1.0 / (n_samples as Float);
        for axis in 0..3 {
            for ch in 0..n_samples {
                pdf += self.pdf_sr(ch, r_proj[axis])
                    * Float::abs(n_local[axis])
                    * ch_prob
                    * axis_prob[axis];
            }
        }
        return pdf;
    }
}

///////////////////////////////////////////////////////////////////////////
// DisneyMaterial Declarations
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
        mode: TransportMode,
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
        let thin = self.thin;
        // normalize lum. to isolate hue+sat
        let c_tint = if lum > 0.0 {
            c / lum
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
            if thin {
                let flat = self.flatness.evaluate(si);
                // Blend between DisneyDiffuse and fake subsurface based on
                // flatness.  Additionally, weight using diffTrans.
                let r: Arc<dyn BxDF> = Arc::new(DisneyDiffuse::new(
                    diffuse_weight * (1.0 - flat) * (1.0 - dt) * c,
                ));
                b.add(&r);
                let r: Arc<dyn BxDF> = Arc::new(DisneyFakeSS::new(
                    diffuse_weight * flat * (1.0 - dt) * c,
                    rough,
                ));
                b.add(&r);
            } else {
                let sd = self.scatter_distance.evaluate(si);
                if sd.is_black() {
                    // No subsurface scattering; use regular (Fresnel modified)
                    // diffuse.
                    let r: Arc<dyn BxDF> = Arc::new(DisneyDiffuse::new(diffuse_weight * c));
                    b.add(&r);
                } else {
                    // Use a BSSRDF instead.
                    let r: Arc<dyn BxDF> =
                        Arc::new(SpecularTransmission::new(&Spectrum::one(), 1.0, e, mode));
                    b.add(&r);
                    todo!();
                    //let r: Arc<dyn BxDF> = Arc::new(DisneyBSSRDF::new(diffuse_weight * c, sd, si, e, BSSRDFMaterialRawPointer::Disney, TransportMode::Radiance));//todo
                    // DisneyBSSRDF
                }
            }

            // Retro-reflection.
            {
                let r: Arc<dyn BxDF> = Arc::new(DisneyRetro::new(diffuse_weight * c, rough));
                b.add(&r);
            }

            // Sheen (if enabled)
            if sheen_weight > 0.0 {
                let r: Arc<dyn BxDF> =
                    Arc::new(DisneySheen::new(diffuse_weight * sheen_weight * c_sheen));
                b.add(&r);
            }
        }

        // Create the microfacet distribution for metallic and/or specular
        // transmission.
        let aspect = Float::sqrt(1.0 - self.anisotropic.evaluate(si) * 0.9);
        let ax = Float::max(0.001, rough * rough / aspect);
        let ay = Float::max(0.001, rough * rough * aspect);

        let distrib = DisneyMicrofacetDistribution::new(ax, ay);

        // Specular is Trowbridge-Reitz with a modified Fresnel function.
        {
            let distrib: Box<dyn MicrofacetDistribution> = Box::new(distrib.clone());

            let spec_tint = self.specular_tint.evaluate(si);
            let cspec0 = Spectrum::lerp(
                metallic_weight,
                &(schlick_r0_from_eta(e) * Spectrum::lerp(spec_tint, &Spectrum::one(), &c_tint)),
                &c,
            );
            let fresnel: Box<dyn Fresnel> =
                Box::new(DisneyFresnel::new(cspec0, metallic_weight, e));
            let r: Arc<dyn BxDF> = Arc::new(MicrofacetReflection::new(
                &Spectrum::from(1.0),
                distrib,
                fresnel,
            ));
            b.add(&r);
        }

        // Clearcoat
        let cc = self.clearcoat.evaluate(si);
        if cc > 0.0 {
            let gloss = lerp(self.clearcoat_gloss.evaluate(si), 0.1, 0.001);
            let r: Arc<dyn BxDF> = Arc::new(DisneyClearcoat::new(cc, gloss));
            b.add(&r);
        }

        // BTDF
        if strans > 0.0 {
            // Walter et al's model, with the provided transmissive term scaled
            // by sqrt(color), so that after two refractions, we're back to the
            // provided color.
            let t = strans * c.sqrt();
            if thin {
                // Scale roughness based on IOR (Burley 2015, Figure 15).
                let rscaled = (0.65 * e - 0.35) * rough;
                let ax = Float::max(0.001, sqr(rscaled) / aspect);
                let ay = Float::max(0.001, sqr(rscaled) * aspect);
                let scaled_distrib: Box<dyn MicrofacetDistribution> =
                    Box::new(TrowbridgeReitzDistribution::new(ax, ay, true));
                let r: Arc<dyn BxDF> = Arc::new(MicrofacetTransmission::new(
                    &t,
                    scaled_distrib,
                    1.0,
                    e,
                    mode,
                ));
                b.add(&r);
            } else {
                let distrib: Box<dyn MicrofacetDistribution> = Box::new(distrib.clone());

                let r: Arc<dyn BxDF> =
                    Arc::new(MicrofacetTransmission::new(&t, distrib, 1.0, e, mode));
                b.add(&r);
            }
        }
        if thin {
            // Lambertian, weighted by (1 - diffTrans)
            let r: Arc<dyn BxDF> = Arc::new(LambertianTransmission::new(&(dt * c)));
            b.add(&r);
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
