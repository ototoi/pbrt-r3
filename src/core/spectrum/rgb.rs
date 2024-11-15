use super::convert::*;
use super::sampled::SampledSpectrum;
use super::utils::*;
use crate::core::pbrt::*;
use std::ops;

const YWEIGHT: [f32; 3] = [0.212671, 0.715160, 0.072169];

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct RGBSpectrum {
    c: [f32; 3],
}

impl RGBSpectrum {
    pub const N_SAMPLES: usize = 3;

    #[inline]
    pub fn new(r: Float, g: Float, b: Float) -> Self {
        RGBSpectrum { c: [r, g, b] }
    }

    #[inline]
    pub fn zero() -> Self {
        RGBSpectrum { c: [0.0, 0.0, 0.0] }
    }

    #[inline]
    pub fn one() -> Self {
        RGBSpectrum { c: [1.0, 1.0, 1.0] }
    }

    pub fn clamp(&self, low: Float, hi: Float) -> Self {
        let c = &self.c;
        let mut r: [Float; 3] = [0.0; 3];
        for i in 0..c.len() {
            r[i] = Float::clamp(c[i], low, hi);
        }
        return RGBSpectrum::from(r);
    }

    pub fn clamp_zero(&self) -> Self {
        return self.clamp(0.0, Float::INFINITY);
    }

    pub fn max_component_value(&self) -> Float {
        let c = &self.c;
        let mut m = c[0];
        for i in 1..c.len() {
            m = Float::max(m, c[i]);
        }
        return m;
    }

    pub fn to_xyz(&self) -> [Float; 3] {
        let c = &self.c;
        return rgb_to_xyz(c);
    }

    pub fn y(&self) -> Float {
        let c = &self.c;
        return YWEIGHT[0] * c[0] + YWEIGHT[1] * c[1] + YWEIGHT[2] * c[2];
    }

    pub fn to_rgb(&self) -> [Float; 3] {
        let c = &self.c;
        return *c;
    }

    #[inline]
    pub fn mul_scalar(&mut self, s: Float) {
        let c = &mut self.c;
        c[0] *= s;
        c[1] *= s;
        c[2] *= s;
    }

    pub fn to_vec(&self) -> Vec<f32> {
        let c = &self.c;
        return c.to_vec();
    }

    pub fn set_vec(&mut self, v: &[f32]) {
        let c = &mut self.c;
        c.copy_from_slice(v);
        //for i in 0..c.len() {
        //    c[i] = v[i];
        //}
    }

    #[inline]
    pub fn is_rgb(&self) -> bool {
        return true;
    }

    #[inline]
    pub fn is_black(&self) -> bool {
        let c = &self.c;
        return c.iter().all(|x| -> bool { x.abs() <= 0.0 });
    }

    pub fn is_valid(&self) -> bool {
        let c = &self.c;
        return c.iter().all(|x| -> bool { x.is_finite() && !x.is_nan() });
    }

    pub fn rgb_from_xyz(xyz: &[Float]) -> RGBSpectrum {
        RGBSpectrum {
            c: convert::xyz_to_rgb(xyz),
        }
    }

    pub fn rgb_from_sampled(lambda: &[Float], vals: &[Float]) -> RGBSpectrum {
        // Sort samples if unordered, use sorted for returned spectrum
        if !spectrum_samples_sorted(lambda, vals) {
            let mut slambda = Vec::new();
            let mut sv = Vec::new();
            slambda.extend_from_slice(lambda);
            sv.extend_from_slice(vals);
            sort_spectrum_samples(&mut slambda, &mut sv);
            return Self::rgb_from_sampled(&slambda, &sv);
        }
        let mut xyz: [Float; 3] = [0.0; 3];
        for i in 0..CIE_SAMPLES {
            let val: f32 = interpolate_spectrum_samples(lambda, vals, CIE_LAMBDA[i]);
            xyz[0] += val * CIE_X[i];
            xyz[1] += val * CIE_Y[i];
            xyz[2] += val * CIE_Z[i];
        }

        let scale = (CIE_LAMBDA[CIE_SAMPLES - 1] - CIE_LAMBDA[0]) as Float
            / (CIE_Y_INTEGRAL * CIE_SAMPLES as Float);
        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;

        return Self::rgb_from_xyz(&xyz);
    }

    pub fn rgb_from_rgb(rgb: &[Float]) -> RGBSpectrum {
        RGBSpectrum {
            c: [rgb[0], rgb[1], rgb[2]],
        }
    }

    pub fn from_rgb(rgb: &[Float], _t: SpectrumType) -> RGBSpectrum {
        return Self::rgb_from_rgb(rgb);
    }

    pub fn from_sampled(lambda: &[Float], vals: &[Float]) -> RGBSpectrum {
        return Self::rgb_from_sampled(lambda, vals);
    }

    //-----------------------------------------------

    pub fn near_equal(a: &RGBSpectrum, b: &RGBSpectrum, eps: f32) -> bool {
        let aa = a.to_vec();
        let bb = b.to_vec();
        if aa.len() != bb.len() {
            return false;
        }
        let s = aa
            .iter()
            .zip(bb)
            .map(|(x, y)| -> f32 { f32::abs(x - y) })
            .sum::<f32>()
            / (aa.len() as f32);
        return s < eps;
    }

    pub fn sqrt(&self) -> Self {
        let c = &self.c;
        return RGBSpectrum {
            c: [Float::sqrt(c[0]), Float::sqrt(c[1]), Float::sqrt(c[2])],
        };
    }

    pub fn exp(&self) -> Self {
        let c = &self.c;
        return RGBSpectrum {
            c: [Float::exp(c[0]), Float::exp(c[1]), Float::exp(c[2])],
        };
    }

    pub fn len(&self) -> usize {
        return 3;
    }
    pub fn is_empty(&self) -> bool {
        return false;
    }
}

impl ops::Index<usize> for RGBSpectrum {
    type Output = Float;
    #[inline]
    fn index(&self, i: usize) -> &Self::Output {
        return &self.c[i];
    }
}

impl ops::IndexMut<usize> for RGBSpectrum {
    #[inline]
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        return &mut self.c[i];
    }
}

impl ops::Mul<Float> for RGBSpectrum {
    type Output = RGBSpectrum;
    #[inline]
    fn mul(self, s: f32) -> RGBSpectrum {
        return RGBSpectrum::from([self[0] * s, self[1] * s, self[2] * s]);
    }
}

impl ops::Mul<RGBSpectrum> for Float {
    type Output = RGBSpectrum;
    #[inline]
    fn mul(self, rhs: RGBSpectrum) -> RGBSpectrum {
        return RGBSpectrum::from([self * rhs[0], self * rhs[1], self * rhs[2]]);
    }
}

impl ops::Add<RGBSpectrum> for RGBSpectrum {
    type Output = RGBSpectrum;
    #[inline]
    fn add(self, s: RGBSpectrum) -> RGBSpectrum {
        let a = &self;
        let b = &s;
        return RGBSpectrum::from([a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
    }
}

impl ops::Sub<RGBSpectrum> for RGBSpectrum {
    type Output = RGBSpectrum;
    #[inline]
    fn sub(self, s: RGBSpectrum) -> RGBSpectrum {
        let a = &self;
        let b = &s;
        return RGBSpectrum::from([a[0] - b[0], a[1] - b[1], a[2] - b[2]]);
    }
}

impl ops::Mul<RGBSpectrum> for RGBSpectrum {
    type Output = RGBSpectrum;
    #[inline]
    fn mul(self, s: RGBSpectrum) -> RGBSpectrum {
        let a = &self;
        let b = &s;
        return RGBSpectrum::from([a[0] * b[0], a[1] * b[1], a[2] * b[2]]);
    }
}

impl ops::Div<RGBSpectrum> for RGBSpectrum {
    type Output = RGBSpectrum;
    #[inline]
    fn div(self, s: RGBSpectrum) -> RGBSpectrum {
        let a = &self;
        let b = &s;
        return RGBSpectrum::from([a[0] / b[0], a[1] / b[1], a[2] / b[2]]);
    }
}

impl ops::AddAssign<RGBSpectrum> for RGBSpectrum {
    #[inline]
    fn add_assign(&mut self, s: RGBSpectrum) {
        let a = self;
        let b = &s;
        a[0] += b[0];
        a[1] += b[1];
        a[2] += b[2];
    }
}

impl ops::SubAssign<RGBSpectrum> for RGBSpectrum {
    #[inline]
    fn sub_assign(&mut self, s: RGBSpectrum) {
        let a = self;
        let b = &s;
        a[0] -= b[0];
        a[1] -= b[1];
        a[2] -= b[2];
    }
}

impl ops::MulAssign<Float> for RGBSpectrum {
    #[inline]
    fn mul_assign(&mut self, s: Float) {
        let a = self;
        a[0] *= s;
        a[1] *= s;
        a[2] *= s;
    }
}

impl ops::MulAssign<RGBSpectrum> for RGBSpectrum {
    #[inline]
    fn mul_assign(&mut self, s: RGBSpectrum) {
        let a = self;
        let b = &s;
        a[0] *= b[0];
        a[1] *= b[1];
        a[2] *= b[2];
    }
}

impl ops::Neg for RGBSpectrum {
    type Output = RGBSpectrum;
    #[inline]
    fn neg(self) -> Self::Output {
        let c = [-self.c[0], -self.c[1], -self.c[2]];
        RGBSpectrum { c }
    }
}

impl Default for RGBSpectrum {
    #[inline]
    fn default() -> Self {
        RGBSpectrum::zero()
    }
}

impl From<f32> for RGBSpectrum {
    #[inline]
    fn from(value: f32) -> Self {
        RGBSpectrum {
            c: [value, value, value],
        }
    }
}

impl From<(f32, f32, f32)> for RGBSpectrum {
    #[inline]
    fn from(value: (f32, f32, f32)) -> Self {
        RGBSpectrum {
            c: [value.0, value.1, value.2],
        }
    }
}

impl From<&[f32; 3]> for RGBSpectrum {
    #[inline]
    fn from(value: &[f32; 3]) -> Self {
        RGBSpectrum {
            c: [value[0], value[1], value[2]],
        }
    }
}

impl From<[f32; 3]> for RGBSpectrum {
    #[inline]
    fn from(value: [f32; 3]) -> Self {
        RGBSpectrum {
            c: [value[0], value[1], value[2]],
        }
    }
}

impl From<Vec<f32>> for RGBSpectrum {
    #[inline]
    fn from(value: Vec<f32>) -> Self {
        RGBSpectrum {
            c: value.try_into().unwrap(),
        }
    }
}

impl From<&RGBSpectrum> for RGBSpectrum {
    fn from(value: &RGBSpectrum) -> Self {
        return *value;
    }
}

impl From<&SampledSpectrum> for RGBSpectrum {
    fn from(value: &SampledSpectrum) -> Self {
        return RGBSpectrum::from(&value.to_rgb());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let v1 = RGBSpectrum::from([1.0, 2.0, 3.0]);
        let v2 = RGBSpectrum::from([4.0, 5.0, 6.0]);
        let v3 = RGBSpectrum::from([5.0, 7.0, 9.0]);
        assert_eq!(v1 + v2, v3);
    }

    #[test]
    fn test_002() {
        let v1 = RGBSpectrum::from([1.0, 2.0, 3.0]);
        let v2 = RGBSpectrum::from([4.0, 5.0, 6.0]);
        let v3 = RGBSpectrum::from([3.0, 3.0, 3.0]);
        assert_eq!(v2 - v1, v3);
    }

    #[test]
    fn test_003() {
        let v1 = RGBSpectrum::from([1.0, 2.0, 3.0]);
        let v2 = RGBSpectrum::from([4.0, 5.0, 6.0]);
        let v3 = RGBSpectrum::from([4.0, 10.0, 18.0]);
        assert_eq!(v2 * v1, v3);
    }
}
