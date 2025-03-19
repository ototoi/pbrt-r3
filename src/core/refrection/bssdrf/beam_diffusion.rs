use super::functions::*;
use crate::core::interpolation::*;
use crate::core::medium::*;
use crate::core::pbrt::*;
use crate::core::refrection::*;

#[derive(Debug)]
pub struct BSSRDFTable {
    pub rho_samples: Vec<Float>,
    pub radius_samples: Vec<Float>,
    pub profile: Vec<Float>,
    pub rho_eff: Vec<Float>,
    pub profile_cdf: Vec<Float>,
}

impl BSSRDFTable {
    pub fn new(n_rho_samples: usize, n_radius_samples: usize) -> Self {
        let rho_samples = vec![0.0; n_rho_samples];
        let radius_samples = vec![0.0; n_radius_samples];
        let profile = vec![0.0; n_rho_samples * n_radius_samples];
        let rho_eff = vec![0.0; n_rho_samples];
        let profile_cdf = vec![0.0; n_rho_samples * n_radius_samples];

        BSSRDFTable {
            rho_samples,
            radius_samples,
            profile,
            rho_eff,
            profile_cdf,
        }
    }

    pub fn eval_profile(&self, rho_index: usize, radius_index: usize) -> Float {
        return self.profile[rho_index * self.radius_samples.len() + radius_index];
    }
}

const N_SAMPLES: usize = 100;
pub fn beam_diffusion_ms(sigma_s: Float, sigma_a: Float, g: Float, eta: Float, r: Float) -> Float {
    let mut ed = 0.0;
    // Precompute information for dipole integrand

    // Compute reduced scattering coefficients $\sigmaps, \sigmapt$ and albedo
    // $\rhop$
    let sigmap_s = sigma_s * (1.0 - g);
    let sigmap_t = sigma_a + sigmap_s;
    let rhop = sigmap_s / sigmap_t;

    // Compute non-classical diffusion coefficient $D_\roman{G}$ using
    // Equation (15.24)
    let d_g = (2.0 * sigma_a + sigmap_s) / (3.0 * sigmap_t * sigmap_t);

    // Compute effective transport coefficient $\sigmatr$ based on $D_\roman{G}$
    let sigma_tr = Float::sqrt(sigma_a / d_g);

    // Determine linear extrapolation distance $\depthextrapolation$ using
    // Equation (15.28)
    let fm1 = fresnel_moment1(eta);
    let fm2 = fresnel_moment2(eta);
    let ze = -2.0 * d_g * (1.0 + 3.0 * fm2) / (1.0 - 2.0 * fm1);

    // Determine exitance scale factors using Equations (15.31) and (15.32)
    let c_phi = 0.25 * (1.0 - 2.0 * fm1);
    let c_e = 0.5 * (1.0 - 3.0 * fm2);
    for i in 0..N_SAMPLES {
        // Sample real point source depth $\depthreal$
        let zr = -Float::ln(1.0 - (i as Float + 0.5) / (N_SAMPLES as Float)) / sigmap_t;

        // Evaluate dipole integrand $E_{\roman{d}}$ at $\depthreal$ and add to
        // _Ed_
        let zv = -zr + 2.0 * ze;
        let dr = Float::sqrt(r * r + zr * zr);
        let dv = Float::sqrt(r * r + zv * zv);

        // Compute dipole fluence rate $\dipole(r)$ using Equation (15.27)
        let phi_d =
            INV_4_PI / d_g * (Float::exp(-sigma_tr * dr) / dr - Float::exp(-sigma_tr * dv) / dv);

        // Compute dipole vector irradiance $-\N{}\cdot\dipoleE(r)$ using
        // Equation (15.27)
        let edn = INV_4_PI
            * (zr * (1.0 + sigma_tr * dr) * Float::exp(-sigma_tr * dr) / (dr * dr * dr)
                - zv * (1.0 + sigma_tr * dv) * Float::exp(-sigma_tr * dv) / (dv * dv * dv));

        // Add contribution from dipole for depth $\depthreal$ to _Ed_
        let e = phi_d * c_phi + edn * c_e;
        let kappa = 1.0 - Float::exp(-2.0 * sigmap_t * (dr + zr));
        ed += kappa * rhop * rhop * e;
    }

    return ed / N_SAMPLES as Float;
}

pub fn beam_diffusion_ss(sigma_s: Float, sigma_a: Float, g: Float, eta: Float, r: Float) -> Float {
    // Compute material parameters and minimum $t$ below the critical angle
    let sigma_t = sigma_a + sigma_s;
    let rho = sigma_s / sigma_t;
    let t_crit = r * Float::sqrt(eta * eta - 1.0);
    let mut ess = 0.0;

    for i in 0..N_SAMPLES {
        // Evaluate single scattering integrand and add to _Ess_
        let ti = t_crit - Float::ln(1.0 - (i as Float + 0.5) / (N_SAMPLES as Float)) / sigma_t;

        // Determine length $d$ of connecting segment and $\cos\theta_\roman{o}$
        let d = Float::sqrt(r * r + ti * ti);
        let cos_theta_o = ti / d;

        // Add contribution of single scattering at depth $t$
        ess += rho * Float::exp(-sigma_t * (d + t_crit)) / (d * d)
            * phase_hg(cos_theta_o, g)
            * (1.0 - fr_dielectric(-cos_theta_o, 1.0, eta))
            * Float::abs(cos_theta_o);
    }
    return ess / N_SAMPLES as Float;
}

pub fn compute_beam_diffusion_bssrdf(g: Float, eta: Float, t: &mut BSSRDFTable) {
    // Choose radius values of the diffusion profile discretization
    let n_radius_samples = t.radius_samples.len();
    t.radius_samples[0] = 0.0;
    t.radius_samples[1] = 2.5e-3;
    for i in 2..n_radius_samples {
        t.radius_samples[i] = t.radius_samples[i - 1] * 1.2;
    }

    // Choose albedo values of the diffusion profile discretization
    let n_rho_samples = t.rho_samples.len();
    for i in 0..n_rho_samples {
        let a = 1.0 - Float::exp(-8.0 * (i as Float) / (n_rho_samples - 1) as Float);
        let b = 1.0 - Float::exp(-8.0);
        t.rho_samples[i] = a / b;
    }

    for i in 0..n_rho_samples {
        // Compute the diffusion profile for the _i_th albedo sample

        // Compute scattering profile for chosen albedo $\rho$
        for j in 0..n_radius_samples {
            let rho = t.rho_samples[i];
            let r = t.radius_samples[j];

            let ss = beam_diffusion_ss(rho, 1.0 - rho, g, eta, r);
            let ms = beam_diffusion_ms(rho, 1.0 - rho, g, eta, r);
            t.profile[i * n_radius_samples + j] = 2.0 * PI * r * (ss + ms);
        }

        // Compute effective albedo $\rho_{\roman{eff}}$ and CDF for importance
        // sampling
        let s0 = i * n_radius_samples;
        let s1 = s0 + n_radius_samples;
        let cdf = integrate_catmull_rom(&t.radius_samples, &t.profile[s0..s1]);
        t.rho_eff[i] = cdf[n_radius_samples - 1]; //The last element of cdf is sum
        t.profile_cdf[s0..s1].copy_from_slice(&cdf);
    }
}
