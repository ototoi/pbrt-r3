use crate::core::interpolation::{catmull_rom_weights, sample_catmull_rom_2d};
use crate::core::interpolation::{evaluate_fourier, sample_fourier};
use crate::core::pbrt::*;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::ops::Deref;
use std::path::Path;
use std::sync::Arc;
use std::sync::RwLock;

const HEADER_EXP: [char; 8] = ['S', 'C', 'A', 'T', 'F', 'U', 'N', '\x01'];

fn read_i32(reader: &mut BufReader<File>) -> Result<i32, PbrtError> {
    let mut buffer: [u8; 4] = [0; 4];
    reader.read_exact(&mut buffer)?;
    let f = i32::from_le_bytes(buffer);
    return Ok(f);
}

fn read_f32(reader: &mut BufReader<File>) -> Result<f32, PbrtError> {
    let mut buffer: [u8; 4] = [0; 4];
    reader.read_exact(&mut buffer)?;
    let f = f32::from_le_bytes(buffer);
    return Ok(f);
}

pub struct FourierBSDFTable {
    pub n_channels: u32,
    pub n_bases: u32,
    pub eta: Float,
    pub mu: Vec<Float>,
    pub m: Vec<i32>,
    pub a_offset: Vec<i32>,
    pub a: Vec<Float>,
    pub a0: Vec<Float>,
    pub cdf: Vec<Float>,
    pub recip: Vec<Float>,
}

impl FourierBSDFTable {
    fn error(filename: &str) -> PbrtError {
        let msg = format!(
            "Tabulated BSDF file \"{}\" has an incompatible file format or version.",
            filename
        );
        return PbrtError::error(&msg);
    }

    pub fn read(filename: &str) -> Result<FourierBSDFTable, PbrtError> {
        let path = Path::new(filename);
        let fp = File::open(path)?;
        let mut reader = BufReader::new(fp);
        let mut header: [u8; 8] = [0; 8];

        reader.read_exact(&mut header)?;
        for i in 0..8 {
            if header[i] as char != HEADER_EXP[i] {
                return Err(Self::error(filename));
            }
        }
        let flags: i32 = read_i32(&mut reader)?;
        let n_mu: usize = read_i32(&mut reader)? as usize;
        let n_coeffs: usize = read_i32(&mut reader)? as usize;
        let m_max: usize = read_i32(&mut reader)? as usize;
        let n_channels: u32 = read_i32(&mut reader)? as u32;
        let n_bases: u32 = read_i32(&mut reader)? as u32;
        for _ in 0..3 {
            let _: i32 = read_i32(&mut reader)?;
        }
        let eta = read_f32(&mut reader)?;
        for _ in 0..4 {
            let _: i32 = read_i32(&mut reader)?;
        }

        if flags != 1 || (n_channels != 1 && n_channels != 3) || n_bases != 1 {
            return Err(Self::error(filename));
        }

        let mut mu: Vec<f32> = vec![0.0; n_mu];
        let mut cdf: Vec<f32> = vec![0.0; n_mu * n_mu];
        let mut a0: Vec<f32> = vec![0.0; n_mu * n_mu];
        let mut offset_and_length: Vec<i32> = vec![0; (n_mu * n_mu * 2) as usize];
        let mut a_offset: Vec<i32> = vec![0; n_mu * n_mu];
        let mut m: Vec<i32> = vec![0; n_mu * n_mu];
        let mut a: Vec<f32> = vec![0.0; n_coeffs];

        for i in 0..n_mu {
            mu[i] = read_f32(&mut reader)?;
        }
        for i in 0..(n_mu * n_mu) {
            cdf[i] = read_f32(&mut reader)?;
        }
        for i in 0..(n_mu * n_mu * 2) {
            offset_and_length[i] = read_i32(&mut reader)?;
        }
        for i in 0..n_coeffs {
            a[i] = read_f32(&mut reader)?;
        }

        for i in 0..(n_mu * n_mu) {
            let offset = offset_and_length[2 * i];
            let length = offset_and_length[2 * i + 1];
            a_offset[i] = offset;
            m[i] = length;

            a0[i] = if length > 0 { a[offset as usize] } else { 0.0 }
        }

        let mut recip = vec![0.0; m_max];
        for i in 0..m_max {
            recip[i] = if i == 0 {
                Float::INFINITY
            } else {
                1.0 / (i as Float)
            };
        }

        let table = FourierBSDFTable {
            n_channels: n_channels as u32,
            n_bases: n_bases as u32,
            eta,
            mu,
            m,
            a_offset,
            a,
            a0,
            cdf,
            recip,
        };
        return Ok(table);
    }

    pub fn get_weights_and_offset(&self, cos_theta: Float) -> Option<(i32, [Float; 4])> {
        return catmull_rom_weights(&self.mu, cos_theta);
    }

    pub fn get_ak(&self, offset_i: i32, offset_o: i32) -> (usize, usize) {
        let n_mu = self.mu.len() as i32;
        let offset_i = i32::clamp(offset_i, 0, n_mu - 1);
        let offset_o = i32::clamp(offset_o, 0, n_mu - 1);
        let offset = (offset_o * n_mu + offset_i) as usize;
        let m = self.m[offset] as usize;
        let a_offset = self.a_offset[offset] as usize;
        return (m, a_offset);
    }
}

pub struct FourierBSDF {
    bsdf_table: Arc<FourierBSDFTable>,
    mode: TransportMode,
}

pub type FourierOffsetsAndWeights = ((i32, [Float; 4]), (i32, [Float; 4]));

impl FourierBSDF {
    pub fn new(bsdf_table: &Arc<FourierBSDFTable>, mode: TransportMode) -> Self {
        FourierBSDF {
            bsdf_table: bsdf_table.clone(),
            mode,
        }
    }

    fn get_offsets_and_weights(
        bsdf_table: &FourierBSDFTable,
        mu_i: Float,
        mu_o: Float,
    ) -> Option<FourierOffsetsAndWeights> {
        if let Some((offset_i, weights_i)) = bsdf_table.get_weights_and_offset(mu_i) {
            if let Some((offset_o, weights_o)) = bsdf_table.get_weights_and_offset(mu_o) {
                return Some(((offset_i, weights_i), (offset_o, weights_o)));
            }
        }
        return None;
    }

    fn create_ak(
        bsdf_table: &FourierBSDFTable,
        offset_i: i32,
        weights_i: &[Float],
        offset_o: i32,
        weights_o: &[Float],
        n_channels: usize,
    ) -> Vec<Float> {
        let mut m_max = 1;

        let mut weights = [0.0; 16];
        let mut m_and_offset = [(0, 0); 16];
        for b in 0..4 {
            for a in 0..4 {
                // Add contribution of _(a, b)_ to $a_k$ values
                let weight = weights_i[a as usize] * weights_o[b as usize];
                if weight != 0.0 {
                    let (m, a_offset) = bsdf_table.get_ak(offset_i + a, offset_o + b);
                    m_max = usize::max(m_max, m);
                    m_and_offset[(4 * b + a) as usize] = (m, a_offset);
                }
                weights[(4 * b + a) as usize] = weight;
            }
        }

        let mut ak = vec![0.0; m_max * n_channels];
        for b in 0..4 {
            for a in 0..4 {
                // Add contribution of _(a, b)_ to $a_k$ values
                let weight = weights[(4 * b + a) as usize];
                if weight != 0.0 {
                    let (m, a_offset) = m_and_offset[(4 * b + a) as usize];
                    for c in 0..n_channels {
                        for k in 0..m {
                            ak[c * m_max + k] += weight * bsdf_table.a[a_offset + c * m + k];
                        }
                    }
                }
            }
        }
        return ak;
    }
}

impl BxDF for FourierBSDF {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        // Find the zenith angle cosines and azimuth difference angle
        let mu_i = cos_theta(&-*wi);
        let mu_o = cos_theta(wo);
        let cos_phi = cos_d_phi(&-*wi, wo);

        // Compute Fourier coefficients $a_k$ for $(\mui, \muo)$

        // Determine offsets and weights for $\mui$ and $\muo$
        let bsdf_table = self.bsdf_table.as_ref();
        let n_channels = bsdf_table.n_channels as usize;
        if let Some(((offset_i, weights_i), (offset_o, weights_o))) =
            Self::get_offsets_and_weights(bsdf_table, mu_i, mu_o)
        {
            /*
            let i_sums:Float = weights_i.iter().sum();
            assert_eq!(i_sums, 1.0);
            let o_sums:Float = weights_o.iter().sum();
            assert_eq!(o_sums, 1.0);
            */

            let ak = Self::create_ak(
                bsdf_table, offset_i, &weights_i, offset_o, &weights_o, n_channels,
            );
            assert!(!ak.is_empty());
            let m_max = ak.len() / n_channels;

            let offset_y = 0 * m_max;
            let offset_r = 1 * m_max;
            let offset_b = 2 * m_max;

            // Evaluate Fourier expansion for angle $\phi$
            let y = Float::max(
                0.0,
                evaluate_fourier(&ak[offset_y..(offset_y + m_max)], cos_phi),
            );
            let mut scale = if mu_i != 0.0 {
                1.0 / Float::abs(mu_i)
            } else {
                0.0
            };

            if self.mode == TransportMode::Radiance && mu_i * mu_o > 0.0 {
                let eta = if mu_i > 0.0 {
                    1.0 / bsdf_table.eta
                } else {
                    bsdf_table.eta
                };
                scale *= eta * eta;
            }

            if n_channels == 1 {
                return Spectrum::from(y * scale);
            } else {
                let r = evaluate_fourier(&ak[offset_r..(offset_r + m_max)], cos_phi);
                let b = evaluate_fourier(&ak[offset_b..(offset_b + m_max)], cos_phi);
                let g = 1.39829 * y - 0.100913 * b - 0.297375 * r;
                let rgb = [r * scale, g * scale, b * scale];
                return Spectrum::from(rgb).clamp_zero();
            }
        }
        return Spectrum::zero();
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        u: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        let bsdf_table = self.bsdf_table.as_ref();
        // Sample zenith angle component for _FourierBSDF_
        let mu_o = cos_theta(wo);
        let (mu_i, _, pdf_mu) = sample_catmull_rom_2d(
            &bsdf_table.mu,
            &bsdf_table.mu,
            &bsdf_table.a0,
            &bsdf_table.cdf,
            mu_o,
            u[1],
        )?;
        //

        // Compute Fourier coefficients $a_k$ for $(\mui, \muo)$

        // Determine offsets and weights for $\mui$ and $\muo$

        let n_channels = bsdf_table.n_channels as usize;
        let ((offset_i, weights_i), (offset_o, weights_o)) =
            Self::get_offsets_and_weights(bsdf_table, mu_i, mu_o)?;

        let ak = Self::create_ak(
            bsdf_table, offset_i, &weights_i, offset_o, &weights_o, n_channels,
        );
        assert!(!ak.is_empty());
        let m_max = ak.len() / n_channels;

        let offset_y = 0 * m_max;
        let offset_r = 1 * m_max;
        let offset_b = 2 * m_max;

        // Importance sample the luminance Fourier expansion
        let (y, pdf_phi, phi) =
            sample_fourier(&ak[offset_y..(offset_y + m_max)], &bsdf_table.recip, u[0]);
        //let y = 1.0;
        //let phi = 1.0;
        //let pdf_phi = 1.0;
        let pdf = Float::max(0.0, pdf_phi * pdf_mu);
        if pdf <= 0.0 {
            return None;
        }

        // Compute the scattered direction for _FourierBSDF_
        let sin2_theta_i = Float::max(0.0, 1.0 - mu_i * mu_i);
        let mut norm = Float::sqrt(sin2_theta_i / sin_2_theta(wo));
        if Float::is_infinite(norm) {
            norm = 0.0;
        }
        let sin_phi = Float::sin(phi);
        let cos_phi = Float::cos(phi);
        let wi = -Vector3f::new(
            norm * (cos_phi * wo.x - sin_phi * wo.y),
            norm * (sin_phi * wo.x + cos_phi * wo.y),
            mu_i,
        )
        .normalize();

        // Mathematically, wi will be normalized (if wo was). However, in
        // practice, floating-point rounding error can cause some error to
        // accumulate in the computed value of wi here. This can be
        // catastrophic: if the ray intersects an object with the FourierBSDF
        // again and the wo (based on such a wi) is nearly perpendicular to the
        // surface, then the wi computed at the next intersection can end up
        // being substantially (like 4x) longer than normalized, which leads to
        // all sorts of errors, including negative spectral values. Therefore,
        // we normalize again here.
        // wi = wi.normalize();

        // Evaluate remaining Fourier expansions for angle $\phi$
        let mut scale = if mu_i != 0.0 {
            1.0 / Float::abs(mu_i)
        } else {
            0.0
        };
        if self.mode == TransportMode::Radiance && mu_i * mu_o > 0.0 {
            let eta = if mu_i > 0.0 {
                1.0 / bsdf_table.eta
            } else {
                bsdf_table.eta
            };
            scale *= eta * eta;
        }

        if n_channels == 1 {
            let f = Spectrum::from(y * scale);
            return Some((f, wi, pdf, 0));
        } else {
            let r = evaluate_fourier(&ak[offset_r..(offset_r + m_max)], cos_phi);
            let b = evaluate_fourier(&ak[offset_b..(offset_b + m_max)], cos_phi);
            let g = 1.39829 * y - 0.100913 * b - 0.297375 * r;
            let rgb = [r * scale, g * scale, b * scale];
            let f = Spectrum::from(rgb).clamp_zero();
            //let f = Spectrum::one();
            return Some((f, wi, pdf, 0));
        }
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        // Find the zenith angle cosines and azimuth difference angle
        let mu_i = cos_theta(&-*wi);
        let mu_o = cos_theta(wo);
        let cos_phi = cos_d_phi(&-*wi, wo);

        // Compute Fourier coefficients $a_k$ for $(\mui, \muo)$

        // Determine offsets and weights for $\mui$ and $\muo$
        let bsdf_table = self.bsdf_table.as_ref();
        let n_channels = 1;
        if let Some(((offset_i, weights_i), (offset_o, weights_o))) =
            Self::get_offsets_and_weights(bsdf_table, mu_i, mu_o)
        {
            let ak = Self::create_ak(
                bsdf_table, offset_i, &weights_i, offset_o, &weights_o, n_channels,
            );
            assert!(!ak.is_empty());
            let m_max = ak.len() / n_channels;

            // Evaluate probability of sampling _wi_
            let n_mu = bsdf_table.mu.len() as i32;
            let mut rho = 0.0;
            for o in 0..4 {
                if weights_o[o as usize] == 0.0 {
                    continue;
                }
                let index = i32::clamp(offset_o + o, 0, n_mu - 1);

                rho += weights_o[o as usize]
                    * bsdf_table.cdf[(index * n_mu + n_mu - 1) as usize]
                    * (2.0 * PI);
            }

            let offset_y = 0 * m_max;
            let y = evaluate_fourier(&ak[offset_y..(offset_y + m_max)], cos_phi);
            return if rho > 0.0 && y > 0.0 { y / rho } else { 0.0 };
        }
        return 0.0;
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY;
    }

    fn to_string(&self) -> String {
        return format!("FourierBSDF {:?}", self.get_type());
    }
}
