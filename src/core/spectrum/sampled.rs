use super::blackbody::*;
use super::config::*;
use super::constants;
use super::convert::*;
use super::load::*;
use super::rgb::RGBSpectrum;
use super::utils::*;
use crate::core::pbrt::lerp;
use crate::core::pbrt::*;

use std::ops;

#[derive(Copy, Clone, Debug)]
pub enum SpectrumType {
    Reflectance,
    Illuminant,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct SampledSpectrum {
    pub c: [Float; SPECTRAL_SAMPLES],
}

impl SampledSpectrum {
    pub const N_SAMPLES: usize = SPECTRAL_SAMPLES;

    pub fn zero() -> Self {
        SampledSpectrum {
            c: [0.0; SPECTRAL_SAMPLES],
        }
    }

    pub fn one() -> Self {
        SampledSpectrum {
            c: [1.0; SPECTRAL_SAMPLES],
        }
    }

    pub fn clamp(&self, low: Float, hi: Float) -> Self {
        let c = &self.c;
        let mut r = [0.0; SPECTRAL_SAMPLES];
        for i in 0..c.len() {
            r[i] = Float::clamp(c[i], low, hi);
        }
        return SampledSpectrum::from(r);
    }

    pub fn clamp_zero(&self) -> Self {
        return self.clamp(0.0, Float::INFINITY);
    }

    pub fn to_xyz(&self) -> [Float; 3] {
        let c = &self.c;
        let mut xyz: [Float; 3] = [0.0; 3];
        for i in 0..c.len() {
            xyz[0] += data::xyz::ARRAY_CIE_X[i] * c[i];
            xyz[1] += data::xyz::ARRAY_CIE_Y[i] * c[i];
            xyz[2] += data::xyz::ARRAY_CIE_Z[i] * c[i];
        }
        let scale = (config::SAMPLED_LAMBDA_END - config::SAMPLED_LAMBDA_START) as Float
            / (data::CIE_Y_INTEGRAL * config::SPECTRAL_SAMPLES as Float);
        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;
        return xyz;
    }

    pub fn y(&self) -> Float {
        let c = &self.c;
        let mut yy = 0.0;
        for i in 0..c.len() {
            yy += data::ARRAY_CIE_Y[i] * c[i];
        }
        let scale = (config::SAMPLED_LAMBDA_END - config::SAMPLED_LAMBDA_START) as Float
            / (data::CIE_Y_INTEGRAL * config::SPECTRAL_SAMPLES as Float);
        return yy * scale;
    }

    pub fn max_component_value(&self) -> Float {
        let c = &self.c;
        let mut m = c[0];
        for i in 1..c.len() {
            m = Float::max(m, c[i]);
        }
        return m;
    }

    pub fn to_rgb(&self) -> [Float; 3] {
        let xyz = self.to_xyz();
        return xyz_to_rgb(&xyz);
    }

    pub fn mul_scalar(&mut self, s: Float) {
        let c = &mut self.c;
        for i in 0..c.len() {
            c[i] *= s;
        }
    }

    pub fn div_scalar(&mut self, s: Float) {
        let c = &mut self.c;
        for i in 0..c.len() {
            c[i] /= s;
        }
    }

    pub fn to_vec(&self) -> Vec<Float> {
        let c = &self.c;
        return c.to_vec();
    }

    pub fn set_vec(&mut self, v: &[Float]) {
        let c = &mut self.c;
        c.copy_from_slice(v);
        //for i in 0..c.len() {
        //    c[i] = v[i];
        //}
    }

    #[inline]
    pub fn is_rgb(&self) -> bool {
        return false;
    }

    pub fn is_black(&self) -> bool {
        let c = &self.c;
        return c.iter().all(|x| -> bool { *x <= 0.0 });
    }

    pub fn is_valid(&self) -> bool {
        let c = &self.c;
        return c.iter().all(|x| -> bool { x.is_finite() });
    }

    pub fn sampled_from_sampled(lambda: &[Float], vals: &[Float]) -> SampledSpectrum {
        if !spectrum_samples_sorted(lambda, vals) {
            let mut slambda = Vec::new();
            let mut sv = Vec::new();
            slambda.extend_from_slice(lambda);
            sv.extend_from_slice(vals);
            sort_spectrum_samples(&mut slambda, &mut sv);
            return Self::sampled_from_sampled(&slambda, &sv);
        }
        return SampledSpectrum {
            c: sample_spectrum(lambda, vals),
        };
    }

    pub fn sampled_from_rgb(rgb: &[Float], t: SpectrumType) -> SampledSpectrum {
        let mut r: SampledSpectrum = SampledSpectrum {
            c: [0.0; SPECTRAL_SAMPLES],
        };
        match t {
            SpectrumType::Reflectance => {
                if rgb[0] <= rgb[1] && rgb[0] <= rgb[2] {
                    r += constants::RGBREFL2SPECT_WHITE * rgb[0];
                    if rgb[1] <= rgb[2] {
                        r += constants::RGBREFL2SPECT_CYAN * (rgb[1] - rgb[0]);
                        r += constants::RGBREFL2SPECT_BLUE * (rgb[2] - rgb[1]);
                    } else {
                        r += constants::RGBREFL2SPECT_CYAN * (rgb[2] - rgb[0]);
                        r += constants::RGBREFL2SPECT_GREEN * (rgb[1] - rgb[2]);
                    }
                } else if rgb[1] <= rgb[0] && rgb[1] <= rgb[2] {
                    r += constants::RGBREFL2SPECT_WHITE * rgb[1];
                    if rgb[0] <= rgb[2] {
                        r += constants::RGBREFL2SPECT_MAGENTA * (rgb[0] - rgb[1]);
                        r += constants::RGBREFL2SPECT_BLUE * (rgb[2] - rgb[0]);
                    } else {
                        r += constants::RGBREFL2SPECT_MAGENTA * (rgb[2] - rgb[1]);
                        r += constants::RGBREFL2SPECT_RED * (rgb[0] - rgb[2]);
                    }
                } else {
                    r += constants::RGBREFL2SPECT_WHITE * rgb[2];
                    if rgb[0] <= rgb[2] {
                        r += constants::RGBREFL2SPECT_YELLOW * (rgb[0] - rgb[2]);
                        r += constants::RGBREFL2SPECT_GREEN * (rgb[1] - rgb[0]);
                    } else {
                        r += constants::RGBREFL2SPECT_YELLOW * (rgb[1] - rgb[2]);
                        r += constants::RGBREFL2SPECT_RED * (rgb[0] - rgb[1]);
                    }
                }
                r *= 0.94;
            }
            SpectrumType::Illuminant => {
                if rgb[0] <= rgb[1] && rgb[0] <= rgb[2] {
                    r += constants::RGBILLUM2SPECT_WHITE * rgb[0];
                    if rgb[1] <= rgb[2] {
                        r += constants::RGBILLUM2SPECT_CYAN * (rgb[1] - rgb[0]);
                        r += constants::RGBILLUM2SPECT_BLUE * (rgb[2] - rgb[1]);
                    } else {
                        r += constants::RGBILLUM2SPECT_CYAN * (rgb[2] - rgb[0]);
                        r += constants::RGBILLUM2SPECT_GREEN * (rgb[1] - rgb[2]);
                    }
                } else if rgb[1] <= rgb[0] && rgb[1] <= rgb[2] {
                    r += constants::RGBILLUM2SPECT_WHITE * rgb[1];
                    if rgb[0] <= rgb[2] {
                        r += constants::RGBILLUM2SPECT_MAGENTA * (rgb[0] - rgb[1]);
                        r += constants::RGBILLUM2SPECT_BLUE * (rgb[2] - rgb[0]);
                    } else {
                        r += constants::RGBILLUM2SPECT_MAGENTA * (rgb[2] - rgb[1]);
                        r += constants::RGBILLUM2SPECT_RED * (rgb[0] - rgb[2]);
                    }
                } else {
                    r += constants::RGBILLUM2SPECT_WHITE * rgb[2];
                    if rgb[0] <= rgb[2] {
                        r += constants::RGBILLUM2SPECT_YELLOW * (rgb[0] - rgb[2]);
                        r += constants::RGBILLUM2SPECT_GREEN * (rgb[1] - rgb[0]);
                    } else {
                        r += constants::RGBILLUM2SPECT_YELLOW * (rgb[1] - rgb[2]);
                        r += constants::RGBILLUM2SPECT_RED * (rgb[0] - rgb[1]);
                    }
                }
                r *= 0.86445;
            }
        }

        return r.clamp(0.0, Float::INFINITY);
    }

    pub fn sampled_from_xyz(xyz: &[Float], t: SpectrumType) -> SampledSpectrum {
        let rgb = convert::xyz_to_rgb(xyz);
        return Self::sampled_from_rgb(&rgb, t);
    }

    pub fn load_sampled_spectrum_file(path: &str) -> Result<SampledSpectrum, PbrtError> {
        return load_from_file(path);
    }

    pub fn from_rgb(rgb: &[Float], t: SpectrumType) -> SampledSpectrum {
        return Self::sampled_from_rgb(rgb, t);
    }

    pub fn from_sampled(lambda: &[Float], vals: &[Float]) -> SampledSpectrum {
        return Self::sampled_from_sampled(lambda, vals);
    }

    //-----------------------------------------------
    // ref. ParamSet::AddBlackbodySpectrum
    pub fn from_blackbody(values: &[Float]) -> SampledSpectrum {
        let n_values = values.len();
        debug_assert_eq!(n_values % 2, 0);
        let n_values = n_values / 2;
        let mut s = SampledSpectrum::zero();
        for i in 0..n_values {
            let v = blackbody_normalized(&CIE_LAMBDA, values[2 * i + 0]);
            s += Self::sampled_from_sampled(&CIE_LAMBDA, &v) * values[2 * i + 1];
        }
        return s;
    }

    pub fn near_equal(a: &SampledSpectrum, b: &SampledSpectrum, eps: Float) -> bool {
        let aa = a.to_vec();
        let bb = b.to_vec();
        if aa.len() != bb.len() {
            return false;
        }
        let s = aa
            .iter()
            .zip(bb)
            .map(|(x, y)| -> Float { Float::abs(x - y) })
            .sum::<Float>()
            / (aa.len() as Float);
        return s < eps;
    }

    pub fn sqrt(&self) -> Self {
        let v: Vec<_> = self
            .c
            .iter()
            .map(|x| -> Float { Float::sqrt(*x) })
            .collect();
        let a: [Float; 60] = v.try_into().unwrap();
        return SampledSpectrum { c: a };
    }

    pub fn exp(&self) -> Self {
        let v: Vec<_> = self.c.iter().map(|x| -> Float { Float::exp(*x) }).collect();
        let a: [Float; 60] = v.try_into().unwrap();
        return SampledSpectrum { c: a };
    }

    pub fn len(&self) -> usize {
        return SPECTRAL_SAMPLES;
    }
    pub fn is_empty(&self) -> bool {
        return false;
    }

    pub fn lerp(t: Float, s1: &Self, s2: &Self) -> Self {
        let a = &s1.c;
        let b = &s2.c;
        let c: Vec<_> = a
            .iter()
            .zip(b.iter())
            .map(|(aa, bb)| lerp(t, *aa, *bb))
            .collect();
        return SampledSpectrum::from(c);
    }
}

impl ops::Index<usize> for SampledSpectrum {
    type Output = Float;
    #[inline]
    fn index(&self, i: usize) -> &Self::Output {
        return &self.c[i];
    }
}

impl ops::IndexMut<usize> for SampledSpectrum {
    #[inline]
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        return &mut self.c[i];
    }
}

impl ops::Mul<Float> for SampledSpectrum {
    type Output = SampledSpectrum;
    #[inline]
    fn mul(self, s: Float) -> SampledSpectrum {
        let mut spc = self;
        spc.mul_scalar(s);
        return spc;
    }
}

impl ops::Div<Float> for SampledSpectrum {
    type Output = SampledSpectrum;
    #[inline]
    fn div(self, s: Float) -> SampledSpectrum {
        let mut spc = self;
        spc.div_scalar(s);
        return spc;
    }
}

impl ops::Mul<SampledSpectrum> for Float {
    type Output = SampledSpectrum;
    #[inline]
    fn mul(self, rhs: SampledSpectrum) -> SampledSpectrum {
        let mut spc = rhs;
        spc.mul_scalar(self);
        return spc;
    }
}

impl ops::Add<SampledSpectrum> for SampledSpectrum {
    type Output = SampledSpectrum;
    #[inline]
    fn add(self, s: SampledSpectrum) -> SampledSpectrum {
        let a = &self.c;
        let b = &s.c;
        let c: Vec<_> = a.iter().zip(b.iter()).map(|(aa, bb)| aa + bb).collect();
        return SampledSpectrum::from(c);
    }
}

impl ops::Sub<SampledSpectrum> for SampledSpectrum {
    type Output = SampledSpectrum;
    #[inline]
    fn sub(self, s: SampledSpectrum) -> SampledSpectrum {
        let a = &self.c;
        let b = &s.c;
        let c: Vec<_> = a.iter().zip(b.iter()).map(|(aa, bb)| aa - bb).collect();
        return SampledSpectrum::from(c);
    }
}

impl ops::Mul<SampledSpectrum> for SampledSpectrum {
    type Output = SampledSpectrum;
    #[inline]
    fn mul(self, s: SampledSpectrum) -> SampledSpectrum {
        let a = &self.c;
        let b = &s.c;
        let c: Vec<_> = a.iter().zip(b.iter()).map(|(aa, bb)| aa * bb).collect();
        return SampledSpectrum::from(c);
    }
}

impl ops::Div<SampledSpectrum> for SampledSpectrum {
    type Output = SampledSpectrum;
    #[inline]
    fn div(self, s: SampledSpectrum) -> SampledSpectrum {
        let a = &self.c;
        let b = &s.c;
        let c: Vec<_> = a.iter().zip(b.iter()).map(|(aa, bb)| aa / bb).collect();
        return SampledSpectrum::from(c);
    }
}

impl ops::AddAssign<SampledSpectrum> for SampledSpectrum {
    #[inline]
    fn add_assign(&mut self, s: SampledSpectrum) {
        let mut a = self.to_vec();
        let b = s.to_vec();
        for i in 0..a.len() {
            a[i] += b[i];
        }
        self.set_vec(&a);
    }
}

impl ops::SubAssign<SampledSpectrum> for SampledSpectrum {
    #[inline]
    fn sub_assign(&mut self, s: SampledSpectrum) {
        let mut a = self.to_vec();
        let b = s.to_vec();

        for i in 0..a.len() {
            a[i] -= b[i];
        }
        self.set_vec(&a);
    }
}

impl ops::MulAssign<Float> for SampledSpectrum {
    #[inline]
    fn mul_assign(&mut self, s: Float) {
        self.mul_scalar(s);
    }
}

impl ops::MulAssign<SampledSpectrum> for SampledSpectrum {
    #[inline]
    fn mul_assign(&mut self, s: SampledSpectrum) {
        let mut a = self.to_vec();
        let b = s.to_vec();

        for i in 0..a.len() {
            a[i] *= b[i];
        }
        self.set_vec(&a);
    }
}

impl ops::DivAssign<Float> for SampledSpectrum {
    #[inline]
    fn div_assign(&mut self, s: Float) {
        self.div_scalar(s);
    }
}

impl ops::Neg for SampledSpectrum {
    type Output = SampledSpectrum;
    #[inline]
    fn neg(self) -> Self::Output {
        let c: Vec<Float> = self.c.iter().map(|v| -v).collect();
        return SampledSpectrum::from(c);
    }
}

impl Default for SampledSpectrum {
    #[inline]
    fn default() -> Self {
        SampledSpectrum::zero()
    }
}

impl From<Float> for SampledSpectrum {
    #[inline]
    fn from(value: Float) -> Self {
        let mut c = [0.0; SPECTRAL_SAMPLES];
        for i in 0..SPECTRAL_SAMPLES {
            c[i] = value;
        }
        SampledSpectrum { c }
    }
}

impl From<&[Float; SPECTRAL_SAMPLES]> for SampledSpectrum {
    #[inline]
    fn from(value: &[Float; SPECTRAL_SAMPLES]) -> Self {
        SampledSpectrum { c: *value }
    }
}

impl From<[Float; SPECTRAL_SAMPLES]> for SampledSpectrum {
    #[inline]
    fn from(value: [Float; SPECTRAL_SAMPLES]) -> Self {
        SampledSpectrum { c: value }
    }
}

impl From<[Float; 3]> for SampledSpectrum {
    #[inline]
    fn from(value: [Float; 3]) -> Self {
        return SampledSpectrum::from_rgb(&value, SpectrumType::Reflectance);
    }
}

impl From<&[Float; 3]> for SampledSpectrum {
    #[inline]
    fn from(value: &[Float; 3]) -> Self {
        return SampledSpectrum::from_rgb(value, SpectrumType::Reflectance);
    }
}

impl From<Vec<Float>> for SampledSpectrum {
    #[inline]
    fn from(value: Vec<Float>) -> Self {
        SampledSpectrum {
            c: value.try_into().unwrap(),
        }
    }
}

impl From<&SampledSpectrum> for SampledSpectrum {
    fn from(value: &SampledSpectrum) -> Self {
        return value.clone();
    }
}

impl From<&RGBSpectrum> for SampledSpectrum {
    fn from(value: &RGBSpectrum) -> Self {
        return SampledSpectrum::from_rgb(&value.to_rgb(), SpectrumType::Reflectance);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_004() {
        let v1 = SampledSpectrum::from([1.0; SPECTRAL_SAMPLES]);
        let v2 = SampledSpectrum::from([2.0; SPECTRAL_SAMPLES]);
        let v3 = SampledSpectrum::from([3.0; SPECTRAL_SAMPLES]);
        assert_eq!(v1 + v2, v3);
    }

    #[test]
    fn test_005() {
        let v1 = SampledSpectrum::from([1.0; SPECTRAL_SAMPLES]);
        let v2 = SampledSpectrum::from([2.0; SPECTRAL_SAMPLES]);
        let v3 = SampledSpectrum::from([1.0; SPECTRAL_SAMPLES]);
        assert_eq!(v2 - v1, v3);
    }

    #[test]
    fn test_006() {
        let v1 = SampledSpectrum::from([2.0; SPECTRAL_SAMPLES]);
        let v2 = SampledSpectrum::from([3.0; SPECTRAL_SAMPLES]);
        let v3 = SampledSpectrum::from([6.0; SPECTRAL_SAMPLES]);
        assert_eq!(v2 * v1, v3);
    }

    #[test]
    fn test_007() {
        let v1 = SampledSpectrum::sampled_from_rgb(&[1.0, 2.0, 3.0], SpectrumType::Reflectance);
        let v2 = SampledSpectrum::sampled_from_rgb(&[4.0, 5.0, 6.0], SpectrumType::Reflectance);
        let v3 = SampledSpectrum::sampled_from_rgb(&[5.0, 7.0, 9.0], SpectrumType::Reflectance);
        assert!(SampledSpectrum::near_equal(&(v1 + v2), &v3, 1e-6));
    }

    #[test]
    fn test_008() {
        let v1 = SampledSpectrum::sampled_from_rgb(&[1.0, 2.0, 3.0], SpectrumType::Reflectance);
        let v2 = SampledSpectrum::sampled_from_rgb(&[4.0, 5.0, 6.0], SpectrumType::Reflectance);
        let v3 = SampledSpectrum::sampled_from_rgb(&[3.0, 3.0, 3.0], SpectrumType::Reflectance);
        assert!(SampledSpectrum::near_equal(&(v2 - v1), &v3, 1e-6));
    }
}
