use crate::core::misc::*;

use std::ops;

#[derive(Debug, PartialEq, Default, Copy, Clone)]
pub struct EFloat {
    pub v: f32,
    pub low: f32,
    pub high: f32,
}

impl EFloat {
    pub fn from_f32(v: f32, err: f32) -> Self {
        if err == 0.0 {
            EFloat { v, low: v, high: v }
        } else {
            EFloat {
                v,
                low: next_float_down_32(v - err),
                high: next_float_up_32(v + err),
            }
        }
    }

    pub fn from_f64(v: f64, err: f64) -> Self {
        return Self::from_f32(v as f32, err as f32);
    }

    pub fn upper_bound(&self) -> f32 {
        self.high
    }

    pub fn lower_bound(&self) -> f32 {
        self.low
    }

    pub fn get_absolute_error(&self) -> f32 {
        return next_float_up_32(f32::max(
            f32::abs(self.high - self.v),
            f32::abs(self.v - self.low),
        ));
    }

    pub fn sqrt(&self) -> Self {
        let v = f32::sqrt(self.v);
        let low = next_float_down_32(f32::sqrt(self.low));
        let high = next_float_up_32(f32::sqrt(self.high));
        EFloat { v, low, high }
    }

    pub fn abs(&self) -> Self {
        if self.low >= 0.0 {
            return self.clone();
        } else if self.high <= 0.0 {
            let v = -self.v;
            let low = -self.high;
            let high = -self.low;
            EFloat { v, low, high }
        } else {
            let v = f32::abs(self.v);
            let low = 0.0;
            let high = f32::max(-self.low, self.high);
            EFloat { v, low, high }
        }
    }

    pub fn quadratic(a: EFloat, b: EFloat, c: EFloat) -> Option<(EFloat, EFloat)> {
        const EPS: f64 = f64::EPSILON;
        // Find quadratic discriminant
        let av = a.v as f64;
        let bv = b.v as f64;
        let cv = c.v as f64;
        let discrim = bv * bv - 4.0 * av * cv;
        if discrim < 0.0 {
            return None;
        }
        let root_discrim = discrim.sqrt();
        let float_root_discrim = EFloat::from_f64(root_discrim, EPS * root_discrim);

        // Compute quadratic _t_ values
        let q = if b.v < 0.0 {
            (b - float_root_discrim) * -0.5
        } else {
            (b + float_root_discrim) * -0.5
        };
        let t0 = q / a;
        let t1 = c / q;
        if t0.v <= t1.v {
            return Some((t0, t1));
        } else {
            return Some((t1, t0));
        }
    }
}

impl ops::Add<EFloat> for EFloat {
    type Output = EFloat;
    fn add(self, ef: EFloat) -> EFloat {
        return EFloat {
            v: self.v + ef.v,
            low: next_float_down_32(self.low + ef.low),
            high: next_float_up_32(self.high + ef.high),
        };
    }
}

impl ops::Sub<EFloat> for EFloat {
    type Output = EFloat;
    fn sub(self, ef: EFloat) -> EFloat {
        return EFloat {
            v: self.v - ef.v,
            low: next_float_down_32(self.low - ef.high),
            high: next_float_up_32(self.high - ef.low),
        };
    }
}

impl ops::Mul<EFloat> for EFloat {
    type Output = EFloat;
    fn mul(self, ef: EFloat) -> EFloat {
        let v = self.v * ef.v;
        let prod = [
            self.low * ef.low,
            self.high * ef.low,
            self.low * ef.high,
            self.high * ef.high,
        ];
        let low = next_float_down_32(f32::min(
            f32::min(prod[0], prod[1]),
            f32::min(prod[2], prod[3]),
        ));
        let high = next_float_up_32(f32::max(
            f32::max(prod[0], prod[1]),
            f32::max(prod[2], prod[3]),
        ));
        return EFloat { v, low, high };
    }
}

impl ops::Div<EFloat> for EFloat {
    type Output = EFloat;
    fn div(self, ef: EFloat) -> EFloat {
        let v = self.v / ef.v;
        if ef.low < 0.0 && ef.high > 0.0 {
            let low = -f32::INFINITY;
            let high = f32::INFINITY;
            return EFloat { v, low, high };
        } else {
            let div = [
                self.lower_bound() / ef.lower_bound(),
                self.upper_bound() / ef.lower_bound(),
                self.lower_bound() / ef.upper_bound(),
                self.upper_bound() / ef.upper_bound(),
            ];
            let low =
                next_float_down_32(f32::min(f32::min(div[0], div[1]), f32::min(div[2], div[3])));
            let high =
                next_float_up_32(f32::max(f32::max(div[0], div[1]), f32::max(div[2], div[3])));
            return EFloat { v, low, high };
        }
    }
}

impl ops::Mul<f32> for EFloat {
    type Output = EFloat;
    fn mul(self, f: f32) -> EFloat {
        return self * EFloat::from_f32(f, 0.0);
    }
}

impl ops::Mul<EFloat> for f32 {
    type Output = EFloat;
    fn mul(self, rhs: EFloat) -> EFloat {
        return EFloat::from_f32(self, 0.0) * rhs;
    }
}

impl From<f32> for EFloat {
    fn from(value: f32) -> Self {
        EFloat::from_f32(value, 0.0)
    }
}

impl From<(f32, f32)> for EFloat {
    fn from(value: (f32, f32)) -> Self {
        EFloat::from_f32(value.0, value.1)
    }
}

impl From<f64> for EFloat {
    fn from(value: f64) -> Self {
        EFloat::from_f64(value, 0.0)
    }
}

impl From<(f64, f64)> for EFloat {
    fn from(value: (f64, f64)) -> Self {
        EFloat::from_f64(value.0, value.1)
    }
}

impl From<EFloat> for f32 {
    fn from(value: EFloat) -> f32 {
        value.v
    }
}

impl From<EFloat> for f64 {
    fn from(value: EFloat) -> f64 {
        value.v as f64
    }
}
