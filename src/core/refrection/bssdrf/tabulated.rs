use super::beam_diffusion::BSSRDFTable;
use super::separable::BSSRDFMaterialRawPointer;
use super::separable::*;
use crate::core::interpolation::catmull_rom_weights;
use crate::core::interpolation::sample_catmull_rom_2d;
use crate::core::pbrt::*;

use std::sync::Arc;

#[derive(Debug)]
pub struct TabulatedBSSRDF {
    base: Arc<BaseSeparableBSSRDF>,
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
        let base = Arc::new(BaseSeparableBSSRDF::new(po, eta, material, mode));
        TabulatedBSSRDF {
            base,
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

    fn sw(&self, w: &Vector3f) -> Spectrum {
        return self.base.sw(w);
    }

    fn sp(&self, pi: &SurfaceInteraction) -> Spectrum {
        let p = self.base.base.p;
        return self.sr(Vector3f::distance(&p, &pi.p));
    }

    fn sample_sp(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        _arena: &mut MemoryArena,
    ) -> Option<(Spectrum, SurfaceInteraction, Float)> {
        let (vx, vy, vz, u1) = self.base.projection_axis(u1);
        let n_samples = Spectrum::N_SAMPLES;
        // Choose spectral channel for BSSRDF sampling
        let ch = usize::clamp((u1 * n_samples as Float) as usize, 0, n_samples - 1);
        let u1 = u1 * n_samples as Float - ch as Float;

        // Sample BSSRDF profile in polar coordinates
        let r = self.sample_sr(ch, u2[0]);
        if r < 0.0 {
            return None;
        }
        let phi = 2.0 * PI * u2[1];

        // Compute BSSRDF profile bounds and intersection height
        let r_max = self.sample_sr(ch, 0.999);
        if r_max < 0.0 {
            return None;
        }
        if r >= r_max {
            return None;
        }
        let l = 2.0 * Float::sqrt(r_max * r_max - r * r);

        assert!(l > 1e-6);
        //println!("l:{:?}", l);

        // Compute BSSRDF sampling ray segment
        let p = self.base.base.p + r * (vx * Float::cos(phi) + vy * Float::sin(phi)) - l * vz * 0.5;
        let time = self.base.base.time;
        let mut base = SurfaceInteraction {
            p: p,
            time: time,
            ..Default::default()
        };
        let p_target = p + l * vz;

        // Intersect BSSRDF sampling ray against the scene geometry

        let mut chain = Vec::new();
        loop {
            let distance = Vector3f::distance(&base.p, &p_target);
            if distance <= MACHINE_EPSILON {
                break;
            }
            let r = base.spawn_ray_to_point(&p_target);

            if let Some(si) = scene.intersect(&r) {
                base = si.clone();
                if let Some(primitive) = si.primitive.as_ref() {
                    let primitive = primitive.upgrade().unwrap();
                    let primitive = primitive.as_ref();
                    if let Some(material) = primitive.get_material() {
                        let material = material.as_ref();
                        let ptr = material as BSSRDFMaterialRawPointer;
                        // Append admissible intersection to _IntersectionChain_
                        #[allow(clippy::vtable_address_comparisons)]
                        if std::ptr::eq(ptr, self.base.material) {
                            chain.push(si);
                        }
                    }
                }
            } else {
                break;
            }
        }

        // Randomly choose one of several intersections during BSSRDF sampling
        if chain.is_empty() {
            return None;
        }
        let n_found = chain.len();
        //println!("n_found:{:?}", n_found);
        let selected = usize::clamp((u1 * n_found as Float) as usize, 0, n_found - 1);
        let pi = chain[selected].clone();

        // Compute sample PDF and return the spatial BSSRDF term $\Sp$
        let pdf = self.pdf_sp(&pi) / (n_found as Float);
        let f = self.sp(&pi);
        return Some((f, pi, pdf));
    }

    pub fn pdf_sp(&self, pi: &SurfaceInteraction) -> Float {
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
}

impl BSSRDF for TabulatedBSSRDF {
    fn s(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum {
        let eta = self.base.base.eta;
        let wo = self.base.base.wo;
        let ft = fr_dielectric(cos_theta(&wo), 1.0, eta);
        return (self.sp(pi) * self.sw(wi)) * (1.0 - ft);
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
                let adapter: Arc<dyn BxDF> =
                    Arc::new(SeparableBSSRDFAdapter::new(self.base.clone()));
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
