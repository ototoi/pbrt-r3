use crate::core::pbrt::*;
use std::sync::Arc;

use log;

const P_MAX: usize = 3;
const SQRT_PI_OVER8: Float = 0.626657069;

// Hair Local Declarations
#[inline]
fn mp(
    cos_theta_i: Float,
    cos_theta_o: Float,
    sin_theta_i: Float,
    sin_theta_o: Float,
    v: Float,
) -> Float {
    const LN2: Float = std::f32::consts::LN_2 as Float;

    let a = cos_theta_i * cos_theta_o / v;
    let b = sin_theta_i * sin_theta_o / v;
    assert!(v.is_finite() && v > 0.0);
    assert!(cos_theta_i.is_finite());
    assert!(cos_theta_o.is_finite());
    assert!(a.is_finite());
    assert!(b.is_finite());

    let mp = if v <= 0.1 {
        Float::exp(log_i0(a) - b - 1.0 / v + LN2 + Float::ln(1.0 / (2.0 * v)))
    } else {
        Float::exp(-b) * i0(a) / (Float::sinh(1.0 / v) * 2.0 * v)
    };
    assert!(!mp.is_infinite());
    assert!(!mp.is_nan());

    // pbrt-r3:
    let mp = mp.max(0.0);
    // pbrt-r3:

    return mp;
}

#[inline]
fn i0(x: Float) -> Float {
    let mut val = 0.0;
    let mut x2i = 1.0;
    let mut ifact = 1;
    let mut i4 = 1;
    // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
    for i in 0..10 {
        if i > 1 {
            ifact *= i;
        }
        val += x2i / (i4 as Float * sqr(ifact as Float));
        x2i *= x * x;
        i4 *= 4;
    }
    val
}

#[inline]
fn log_i0(x: Float) -> Float {
    if x > 12.0 {
        return x + 0.5 * (-Float::ln(2.0 * PI) + Float::ln(1.0 / x) + 1.0 / (8.0 * x));
    } else {
        return Float::ln(i0(x));
    }
}

#[inline]
fn ap(cos_theta_o: Float, eta: Float, h: Float, t: &Spectrum) -> [Spectrum; P_MAX + 1] {
    let mut ap = [Spectrum::zero(); P_MAX + 1];
    // Compute $p=0$ attenuation at initial cylinder intersection
    let cos_gamma_o = safe_sqrt(1.0 - h * h);
    let cos_theta = cos_theta_o * cos_gamma_o;
    let f = fr_dielectric(cos_theta, 1.0, eta);
    let tf = *t * f;

    ap[0] = Spectrum::from(f);
    // Compute $p=1$ attenuation term
    ap[1] = *t * sqr(1.0 - f);
    // Compute attenuation terms up to $p=_pMax_$
    for p in 2..P_MAX {
        ap[p] = ap[p - 1] * tf;
    }

    // Compute attenuation term accounting for remaining orders of scattering
    ap[P_MAX] = (ap[P_MAX - 1] * tf) / (Spectrum::one() - tf);

    return ap;
}

#[inline]
fn phi(p: usize, gamma_o: Float, gamma_t: Float) -> Float {
    let p = p as Float;
    return 2.0 * p * gamma_t - 2.0 * gamma_o + p * PI;
}

#[inline]
fn logistic(x: Float, s: Float) -> Float {
    let x = Float::abs(x);
    return Float::exp(-x / s) / (s * sqr(1.0 + Float::exp(-x / s)));
}

#[inline]
fn logistic_cdf(x: Float, s: Float) -> Float {
    return 1.0 / (1.0 + Float::exp(-x / s));
}

#[inline]
fn trimmed_logistic_cdf(x: Float, s: Float, a: Float, b: Float) -> Float {
    assert!(a < b);
    return logistic(x, s) / (logistic_cdf(b, s) - logistic_cdf(a, s));
}

#[inline]
fn np(phi_: Float, p: usize, s: Float, gamma_o: Float, gamma_t: Float) -> Float {
    let mut dphi = phi_ - phi(p, gamma_o, gamma_t);
    // Remap _dphi_ to $[-\pi,\pi]$
    while dphi > PI {
        dphi -= 2.0 * PI;
    }
    while dphi < -PI {
        dphi += 2.0 * PI;
    }
    let np = trimmed_logistic_cdf(dphi, s, -PI, PI);
    // pbrt-r3:
    let np = np.max(0.0);
    // pbrt-r3:
    return np;
}

#[inline]
fn sample_trimmed_logistic(u: Float, s: Float, a: Float, b: Float) -> Float {
    assert!(a < b);
    let k = logistic_cdf(b, s) - logistic_cdf(a, s);
    let x = -s * Float::ln(1.0 / (u * k + logistic_cdf(a, s)) - 1.0);
    assert!(!Float::is_nan(x));
    return Float::clamp(x, a, b);
}

#[inline]
fn safe_asin(x: Float) -> Float {
    assert!(x >= -1.0001 && x <= 1.0001);
    return Float::asin(Float::clamp(x, -1.0, 1.0));
}

#[inline]
fn safe_sqrt(x: Float) -> Float {
    assert!(x >= -1e-4);
    return Float::sqrt(Float::max(0.0, x));
}

#[inline]
fn sqr(x: Float) -> Float {
    return x * x;
}

#[inline]
fn radians(x: Float) -> Float {
    return x * (PI / 180.0);
}

// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
#[inline]
fn compact1_by1(mut x: u32) -> u32 {
    // TODO: as of Haswell, the PEXT instruction could do all this in a
    // single instruction.
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x &= 0x55555555;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 1)) & 0x33333333;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff;
    return x;
}

#[inline]
fn to_u32(x: u64) -> u32 {
    return (x & 0xFFFFFFFF) as u32;
}

#[inline]
fn demux_float(f: Float) -> Point2f {
    assert!(f >= 0.0 && f < 1.0);
    let v = (f as f64 * ((1u64 << 32) as f64)) as u64;
    assert!(v < 0x100000000);
    let bits = [compact1_by1(to_u32(v)), compact1_by1(to_u32(v >> 1))];
    return Point2f::new(
        bits[0] as Float / (1 << 16) as Float,
        bits[1] as Float / (1 << 16) as Float,
    );
}

pub struct HairMaterial {
    sigma_a: Option<Arc<dyn Texture<Spectrum>>>,
    color: Option<Arc<dyn Texture<Spectrum>>>,
    eumelanin: Option<Arc<dyn Texture<Float>>>,
    pheomelanin: Option<Arc<dyn Texture<Float>>>,
    eta: Arc<dyn Texture<Float>>,
    beta_m: Arc<dyn Texture<Float>>,
    beta_n: Arc<dyn Texture<Float>>,
    alpha: Arc<dyn Texture<Float>>,
}

impl HairMaterial {
    pub fn new(
        sigma_a: &Option<Arc<dyn Texture<Spectrum>>>,
        color: &Option<Arc<dyn Texture<Spectrum>>>,
        eumelanin: &Option<Arc<dyn Texture<Float>>>,
        pheomelanin: &Option<Arc<dyn Texture<Float>>>,
        eta: &Arc<dyn Texture<Float>>,
        beta_m: &Arc<dyn Texture<Float>>,
        beta_n: &Arc<dyn Texture<Float>>,
        alpha: &Arc<dyn Texture<Float>>,
    ) -> Self {
        HairMaterial {
            sigma_a: sigma_a.clone(),
            color: color.clone(),
            eumelanin: eumelanin.clone(),
            pheomelanin: pheomelanin.clone(),
            eta: eta.clone(),
            beta_m: beta_m.clone(),
            beta_n: beta_n.clone(),
            alpha: alpha.clone(),
        }
    }
}

impl Material for HairMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        let bm = self.beta_m.as_ref().evaluate(si);
        let bn = self.beta_n.as_ref().evaluate(si);
        let a = self.alpha.as_ref().evaluate(si);
        let e = self.eta.as_ref().evaluate(si);

        let mut b = arena.alloc_bsdf(si, e);

        let sig_a = if let Some(sigma_a) = self.sigma_a.as_ref() {
            sigma_a.as_ref().evaluate(si).clamp_zero()
        } else if let Some(color) = self.color.as_ref() {
            let c = color.as_ref().evaluate(si).clamp_zero();
            HairBSDF::sigma_a_from_reflectance(&c, bn)
        } else {
            let eumelanin = self.eumelanin.as_ref();
            let pheomelanin = self.pheomelanin.as_ref();
            assert!(eumelanin.is_some() || pheomelanin.is_some());
            HairBSDF::sigma_a_from_concentration(
                if let Some(e) = eumelanin {
                    e.as_ref().evaluate(si).max(0.0)
                } else {
                    0.0
                },
                if let Some(p) = pheomelanin {
                    p.as_ref().evaluate(si).max(0.0)
                } else {
                    0.0
                },
            )
        };
        {
            // Offset along width
            let h = -1.0 + 2.0 * si.uv[1];
            let hair: Arc<dyn BxDF> = Arc::new(HairBSDF::new(h, e, sig_a, bm, bn, a));
            b.add(&hair);
        }

        si.bsdf = Some(Arc::new(b));
    }
}

pub struct HairBSDF {
    h: Float,
    gamma_o: Float,
    eta: Float,
    sigma_a: Spectrum,
    //beta_m: Float,
    //beta_n: Float,
    v: [Float; P_MAX + 1],
    s: Float,
    sin2k_alpha: [Float; 3],
    cos2k_alpha: [Float; 3],
}

// HairBSDF Method Definitions
impl HairBSDF {
    pub fn new(
        h: Float,
        eta: Float,
        sigma_a: Spectrum,
        beta_m: Float,
        beta_n: Float,
        alpha: Float,
    ) -> Self {
        let gamma_o = safe_asin(h);
        assert!(h >= -1.0 && h <= 1.0);
        assert!(beta_m >= 0.0 && beta_m <= 1.0);
        assert!(beta_n >= 0.0 && beta_n <= 1.0);

        assert!(
            P_MAX >= 3,
            "Longitudinal variance code must be updated to handle low pMax"
        );

        let mut v = [0.0; P_MAX + 1];
        let v0 = sqr(0.726 * beta_m + 0.812 * sqr(beta_m) + 3.7 * Float::powf(beta_m, 20.0));
        assert!(v0.is_finite() && v0 > 0.0);
        v[0] = v0;
        v[1] = 0.25 * v0;
        v[2] = 4.0 * v0;
        for i in 3..v.len() {
            v[i] = v[2]; // TODO: is there anything better here?
        }
        // Compute azimuthal logistic scale factor from $\beta_n$
        let s = SQRT_PI_OVER8
            * (0.265 * beta_n + 1.194 * sqr(beta_n) + 5.37 * Float::powf(beta_n, 22.0));

        // Compute $\alpha$ terms for hair scales
        let sin0 = Float::sin(radians(alpha));
        let cos0 = safe_sqrt(1.0 - sqr(sin0));
        let mut sin2k_alpha = [sin0, 0.0, 0.0];
        let mut cos2k_alpha = [cos0, 0.0, 0.0];
        for i in 1..3 {
            sin2k_alpha[i] = 2.0 * cos2k_alpha[i - 1] * sin2k_alpha[i - 1];
            cos2k_alpha[i] = sqr(cos2k_alpha[i - 1]) - sqr(sin2k_alpha[i - 1]);
        }

        HairBSDF {
            h,
            gamma_o,
            eta,
            sigma_a,
            //beta_m,
            //beta_n,
            v,
            s,
            sin2k_alpha,
            cos2k_alpha,
        }
    }

    fn compute_ap_pdf(&self, cos_theta_o: Float) -> [Float; P_MAX + 1] {
        let h = self.h;
        let sigma_a = self.sigma_a;
        let eta = self.eta;
        // Compute array of $A_p$ values for _cosThetaO_
        let sin_theta_o = safe_sqrt(1.0 - cos_theta_o * cos_theta_o);

        // Compute $\cos \thetat$ for refracted ray
        let sin_theta_t = sin_theta_o / eta;
        let cos_theta_t = safe_sqrt(1.0 - sqr(sin_theta_t));

        // Compute $\gammat$ for refracted ray
        let etap = Float::sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
        let sin_gamma_t = h / etap;
        let cos_gamma_t = safe_sqrt(1.0 - sqr(sin_gamma_t));

        // Compute the transmittance _T_ of a single path through the cylinder
        let t = Spectrum::exp(&(-sigma_a * (2.0 * cos_gamma_t / cos_theta_t)));
        let ap = ap(cos_theta_o, eta, h, &t);
        let ap_pdf: Vec<_> = ap.iter().map(|s| s.y()).collect();
        let sum_y = ap_pdf.iter().sum::<Float>();
        let ap_pdf: Vec<_> = ap_pdf.iter().map(|c| c / sum_y).collect();
        let ap_pdf = ap_pdf.iter().map(|x| x.max(0.0)).collect::<Vec<_>>(); //pbrt-r3
        let ap_pdf: [Float; P_MAX + 1] = ap_pdf.try_into().unwrap();
        return ap_pdf;
    }

    fn sigma_a_from_concentration(ce: Float, cp: Float) -> Spectrum {
        let mut sigma_a = [0.0; 3];
        let eumelanin_sigma_a = [0.419, 0.697, 1.37];
        let pheomelanin_sigma_a = [0.187, 0.4, 1.05];
        for i in 0..3 {
            sigma_a[i] = ce * eumelanin_sigma_a[i] + cp * pheomelanin_sigma_a[i];
        }
        return Spectrum::from(sigma_a);
    }

    fn sigma_a_from_reflectance(c: &Spectrum, beta_n: Float) -> Spectrum {
        let mut sigma_a = Spectrum::default();
        for i in 0..sigma_a.len() {
            sigma_a[i] = sqr(Float::ln(c[i])
                / (5.969 - 0.215 * beta_n + 2.532 * sqr(beta_n)
                    - 10.73 * Float::powf(beta_n, 3.0)
                    + 5.574 * Float::powf(beta_n, 4.0)
                    + 0.245 * Float::powf(beta_n, 5.0)));
        }
        return sigma_a;
    }

    fn compute_theta_op(&self, p: usize, sin_theta_o: Float, cos_theta_o: Float) -> (Float, Float) {
        // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
        let sin2k_alpha = &self.sin2k_alpha;
        let cos2k_alpha = &self.cos2k_alpha;
        return match p {
            0 => (
                sin_theta_o * cos2k_alpha[1] - cos_theta_o * sin2k_alpha[1],
                cos_theta_o * cos2k_alpha[1] + sin_theta_o * sin2k_alpha[1],
            ),
            1 => (
                sin_theta_o * cos2k_alpha[0] + cos_theta_o * sin2k_alpha[0],
                cos_theta_o * cos2k_alpha[0] - sin_theta_o * sin2k_alpha[0],
            ),
            2 => (
                sin_theta_o * cos2k_alpha[2] + cos_theta_o * sin2k_alpha[2],
                cos_theta_o * cos2k_alpha[2] - sin_theta_o * sin2k_alpha[2],
            ),
            _ => (sin_theta_o, cos_theta_o),
        };
    }
}

impl BxDF for HairBSDF {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let h = self.h;
        let eta = self.eta;
        let v = &self.v;
        let gamma_o = self.gamma_o;
        let s = self.s;
        let sigma_a = self.sigma_a;

        // Compute hair coordinate system terms related to _wo_
        let sin_theta_o = wo.x;
        let cos_theta_o = safe_sqrt(1.0 - sqr(sin_theta_o));
        let phi_o = Float::atan2(wo.z, wo.y);

        // Compute hair coordinate system terms related to _wi_
        let sin_theta_i = wi.x;
        let cos_theta_i = safe_sqrt(1.0 - sqr(sin_theta_i));
        let phi_i = Float::atan2(wi.z, wi.y);

        // Compute $\cos \thetat$ for refracted ray
        let sin_theta_t = sin_theta_o / eta;
        let cos_theta_t = safe_sqrt(1.0 - sqr(sin_theta_t));

        // Compute $\gammat$ for refracted ray
        let etap = Float::sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
        let sin_gamma_t = h / etap;
        let cos_gamma_t = safe_sqrt(1.0 - sqr(sin_gamma_t));
        let gamma_t = safe_asin(sin_gamma_t);

        // Compute the transmittance _T_ of a single path through the cylinder
        let t = Spectrum::exp(&(-sigma_a * (2.0 * cos_gamma_t / cos_theta_t)));

        // Compute PDF for $A_p$ terms
        let ap = ap(cos_theta_o, eta, h, &t);

        // Evaluate hair BSDF
        let phi = phi_i - phi_o;
        let mut fsum = Spectrum::zero();
        for p in 0..P_MAX {
            // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
            let (sin_theta_op, cos_theta_op) = self.compute_theta_op(p, sin_theta_o, cos_theta_o);
            // Handle out-of-range $\cos \thetao$ from scale adjustment
            let cos_theta_op = Float::abs(cos_theta_op);
            fsum += mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p])
                * ap[p]
                * np(phi, p, s, gamma_o, gamma_t);
        }

        // Compute contribution of remaining terms after _pMax_
        fsum += mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[P_MAX])
            * ap[P_MAX]
            * (1.0 / (2.0 * PI));
        let abs_cos_theta_wi = abs_cos_theta(wi);
        if abs_cos_theta_wi > 0.0 {
            fsum *= 1.0 / abs_cos_theta_wi;
        }
        assert!(fsum.y().is_finite());
        return fsum;
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        u2: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        let h = self.h;
        let eta = self.eta;
        let v = &self.v;
        let gamma_o = self.gamma_o;
        let s = self.s;

        // Compute hair coordinate system terms related to _wo_
        let sin_theta_o = wo.x;
        let cos_theta_o = safe_sqrt(1.0 - sqr(sin_theta_o));
        let phi_o = Float::atan2(wo.z, wo.y);

        // Derive four random samples from _u2_
        let mut u = [demux_float(u2[0]), demux_float(u2[1])];

        // Determine which term $p$ to sample for hair scattering
        let ap_pdf = self.compute_ap_pdf(cos_theta_o);
        let mut p = P_MAX;
        for i in 0..(P_MAX + 1) {
            if u[0][0] < ap_pdf[i] {
                p = i;
                break;
            }
            u[0][0] -= ap_pdf[i];
        }

        // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
        let (sin_theta_op, cos_theta_op) = self.compute_theta_op(p, sin_theta_o, cos_theta_o);

        // Sample $M_p$ to compute $\thetai$
        u[1][0] = Float::max(u[1][0], 1e-5);
        let cos_theta = 1.0 + v[p] * Float::ln(u[1][0] + (1.0 - u[1][0]) * Float::exp(-2.0 / v[p]));
        let sin_theta = safe_sqrt(1.0 - sqr(cos_theta));
        let cos_phi = Float::cos(2.0 * PI * u[1][1]);
        let sin_theta_i = -cos_theta * sin_theta_op + sin_theta * cos_phi * cos_theta_op;
        let cos_theta_i = safe_sqrt(1.0 - sqr(sin_theta_i));

        // Sample $N_p$ to compute $\Delta\phi$

        // Compute $\gammat$ for refracted ray
        let etap = Float::sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
        let sin_gamma_t = h / etap;
        let gamma_t = safe_asin(sin_gamma_t);
        let dphi = if p < P_MAX {
            phi(p, gamma_o, gamma_t) + sample_trimmed_logistic(u[0][1], s, -PI, PI)
        } else {
            2.0 * PI * u[0][1]
        };

        // Compute _wi_ from sampled hair scattering angles
        let phi_i = phi_o + dphi;
        let wi = Vector3f::new(
            sin_theta_i,
            cos_theta_i * Float::cos(phi_i),
            cos_theta_i * Float::sin(phi_i),
        );
        let mut pdf = 0.0;
        for p in 0..P_MAX {
            let (sin_theta_op, cos_theta_op) = self.compute_theta_op(p, sin_theta_o, cos_theta_o);
            // Handle out-of-range $\cos \thetao$ from scale adjustment
            let cos_theta_op = Float::abs(cos_theta_op);
            pdf += mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p])
                * ap_pdf[p]
                * np(dphi, p, s, gamma_o, gamma_t);
        }
        pdf += mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[P_MAX])
            * ap_pdf[P_MAX]
            * (1.0 / (2.0 * PI));
        if pdf > 0.0 {
            let spc = self.f(wo, &wi);
            return Some((spc, wi, pdf, 0));
        } else {
            return None;
        }
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        let h = self.h;
        let eta = self.eta;
        let v = &self.v;
        let gamma_o = self.gamma_o;
        let s = self.s;
        // Compute hair coordinate system terms related to _wo_
        let sin_theta_o = wo.x;
        let cos_theta_o = safe_sqrt(1.0 - sqr(sin_theta_o));
        let phi_o = Float::atan2(wo.z, wo.y);

        // Compute hair coordinate system terms related to _wi_
        let sin_theta_i = wi.x;
        let cos_theta_i = safe_sqrt(1.0 - sqr(sin_theta_i));
        let phi_i = Float::atan2(wi.z, wi.y);

        // Compute $\gammat$ for refracted ray
        let etap = Float::sqrt(eta * eta - sqr(sin_theta_o)) / cos_theta_o;
        let sin_gamma_t = h / etap;
        let gamma_t = safe_asin(sin_gamma_t);

        // Compute PDF for $A_p$ terms
        let ap_pdf = self.compute_ap_pdf(cos_theta_o);
        // Compute PDF sum for hair scattering events
        let phi = phi_i - phi_o;
        let mut pdf = 0.0;
        for p in 0..P_MAX {
            // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
            let (sin_theta_op, cos_theta_op) = self.compute_theta_op(p, sin_theta_o, cos_theta_o);
            // Handle out-of-range $\cos \thetao$ from scale adjustment
            let cos_theta_op = Float::abs(cos_theta_op);
            pdf += mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p])
                * ap_pdf[p]
                * np(phi, p, s, gamma_o, gamma_t);
        }
        pdf += mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[P_MAX])
            * ap_pdf[P_MAX]
            * (1.0 / (2.0 * PI));
        return pdf;
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION;
    }

    fn to_string(&self) -> String {
        return format!("[ HairBSDF ]");
    }
}

pub fn create_hair_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let mut sigma_a = mp.get_spectrum_texture_or_null("sigma_a");
    let color = mp.get_spectrum_texture_or_null("color");
    let eumelanin = mp.get_float_texture_or_null("eumelanin");
    let pheomelanin = mp.get_float_texture_or_null("pheomelanin");

    if sigma_a.is_some() {
        if color.is_some() {
            log::warn!("Ignoring \"color\" parameter since \"sigma_a\" was provided.");
        }
        if eumelanin.is_some() {
            log::warn!("Ignoring \"eumelanin\" parameter since \"sigma_a\" was provided.");
        }
        if pheomelanin.is_some() {
            log::warn!("Ignoring \"pheomelanin\" parameter since \"sigma_a\" was provided.")
        }
    } else if color.is_some() {
        if eumelanin.is_some() {
            log::warn!("Ignoring \"eumelanin\" parameter since \"color\" was provided.");
        }
        if pheomelanin.is_some() {
            log::warn!("Ignoring \"pheomelanin\" parameter since \"color\" was provided.")
        }
    } else if eumelanin.is_some() || pheomelanin.is_some() {
        if sigma_a.is_some() {
            log::warn!("Ignoring \"sigma_a\" parameter since \"eumelanin\" or \"pheomelanin\" was provided.");
        }
        if color.is_some() {
            log::warn!(
                "Ignoring \"color\" parameter since \"eumelanin\" or \"pheomelanin\" was provided."
            );
        }
    } else {
        // Default: brown-ish hair.
        let spc = Spectrum::one();
        sigma_a = Some(Arc::new(ConstantTexture::<Spectrum>::new(&spc)));
    }

    let eta = mp.get_float_texture("eta", 1.55);
    let beta_m: Arc<dyn Texture<Float>> = mp.get_float_texture("beta_m", 0.3);
    let beta_n: Arc<dyn Texture<Float>> = mp.get_float_texture("beta_n", 0.3);
    let alpha = mp.get_float_texture("alpha", 2.0);

    return Ok(Arc::new(HairMaterial::new(
        &sigma_a,
        &color,
        &eumelanin,
        &pheomelanin,
        &eta,
        &beta_m,
        &beta_n,
        &alpha,
    )));
}
