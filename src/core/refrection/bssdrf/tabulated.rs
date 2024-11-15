use super::beam_diffusion::BSSRDFTable;
use super::separable::BSSRDFMaterialRawPointer;
use super::separable::*;
use crate::core::interpolation::catmull_rom_weights;
use crate::core::interpolation::sample_catmull_rom_2d;
use crate::core::pbrt::*;

use std::sync::Arc;

#[derive(Debug)]
pub struct TabulatedBSSRDF {
    base: BaseSeparableBSSRDF,
    table: Arc<BSSRDFTable>,
    sigma_t: Spectrum,
    rho: Spectrum,
}

impl TabulatedBSSRDF {
    pub fn new(
        po: &SurfaceInteraction,
        eta: Float,
        material: BSSRDFMaterialRawPointer,
        mode: TransportMode,
        sigma_a: &Spectrum,
        sigma_s: &Spectrum,
        table: &Arc<BSSRDFTable>,
    ) -> Self {
        let sigma_t = *sigma_a + *sigma_s;
        let rho = Self::get_rho(sigma_s, &sigma_t);
        TabulatedBSSRDF {
            base: BaseSeparableBSSRDF::new(po, eta, material, mode),
            table: table.clone(),
            sigma_t,
            rho,
        }
    }

    fn get_rho(sigma_s: &Spectrum, sigma_t: &Spectrum) -> Spectrum {
        let sv = sigma_s.to_vec();
        let tv = sigma_t.to_vec();
        let rv: Vec<Float> = sv
            .iter()
            .zip(tv)
            .map(|(s, t)| if t != 0.0 { s / t } else { 0.0 })
            .collect();
        return Spectrum::from(rv);
    }

    fn sr_core(&self, ch: usize, r: Float) -> (Float, Float) {
        // Convert $r$ into unitless optical radius $r_{\roman{optical}}$
        let rho = self.rho[ch];
        let r_optical = r * self.sigma_t[ch];

        let table = self.table.as_ref();
        // Compute spline weights to interpolate BSSRDF density on channel _ch_
        if let Some((rho_offset, rho_weights)) = catmull_rom_weights(&table.rho_samples, rho) {
            if let Some((radius_offset, radius_weights)) =
                catmull_rom_weights(&table.radius_samples, r_optical)
            {
                // Return BSSRDF profile density for channel _ch_
                let mut sr = 0.0;
                let mut rho_eff = 0.0;
                for i in 0..4 {
                    if rho_weights[i as usize] == 0.0 {
                        continue;
                    }
                    rho_eff += table.rho_eff[(rho_offset + i) as usize] * rho_weights[i as usize];
                    for j in 0..4 {
                        if radius_weights[j as usize] == 0.0 {
                            continue;
                        }
                        sr += table
                            .eval_profile((rho_offset + i) as usize, (radius_offset + j) as usize)
                            * rho_weights[i as usize]
                            * radius_weights[j as usize];
                    }
                }

                // Cancel marginal PDF factor from tabulated BSSRDF profile
                if r_optical != 0.0 {
                    sr /= 2.0 * PI * r_optical;
                }
                let sigma_t = self.sigma_t[ch];
                let sigma_2 = sigma_t * sigma_t;
                return (sr * sigma_2, rho_eff);
            }
        }
        return (0.0, 1.0);
    }
}

impl BSSRDF for TabulatedBSSRDF {
    fn s(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum {
        return self.s_as_separable(pi, wi);
    }

    fn sample_s(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        arena: &mut MemoryArena,
    ) -> Option<(Spectrum, SurfaceInteraction, Float)> {
        return self.sample_s_as_separable(scene, u1, u2, arena);
    }
}

impl SeparableBSSRDF for TabulatedBSSRDF {
    fn sr(&self, r: Float) -> Spectrum {
        let n_samples = Spectrum::N_SAMPLES;
        let mut values = Spectrum::zero().to_rgb();
        for ch in 0..n_samples {
            let (sr, _) = self.sr_core(ch, r);
            values[ch] = Float::max(0.0, sr);
        }
        return Spectrum::from(values).clamp_zero();
    }

    fn sample_sr(&self, ch: usize, u: Float) -> Float {
        if self.sigma_t[ch] == 0.0 {
            return -1.0;
        }
        let rho = self.rho[ch];
        let table = self.table.as_ref();
        if let Some((value, _, _)) = sample_catmull_rom_2d(
            &table.rho_samples,
            &table.radius_samples,
            &table.profile,
            &table.profile_cdf,
            rho,
            u,
        ) {
            return value / self.sigma_t[ch];
        } else {
            return -1.0;
        }
    }

    fn pdf_sr(&self, ch: usize, r: Float) -> Float {
        let (sr, eff) = self.sr_core(ch, r);
        return Float::max(0.0, sr / eff);
    }

    fn as_separable(&self) -> &BaseSeparableBSSRDF {
        return &self.base;
    }
}
