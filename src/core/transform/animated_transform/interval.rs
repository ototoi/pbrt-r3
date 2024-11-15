use crate::core::pbrt::*;

use std::ops;

#[derive(Debug, Clone, Copy)]
pub struct Interval {
    pub low: Float,
    pub high: Float,
}

impl Interval {
    pub fn new(low: Float, high: Float) -> Self {
        Self { low, high }
    }

    pub fn sin(i: Interval) -> Self {
        assert!(0.0 <= i.low);
        assert!(i.high <= 2.0001 * PI);

        let mut sin_low = Float::sin(i.low);
        let mut sin_high = Float::sin(i.high);
        if sin_low > sin_high {
            std::mem::swap(&mut sin_low, &mut sin_high);
        }
        if i.low < PI / 2.0 && i.high > PI / 2.0 {
            sin_high = 1.;
        }
        if i.low < (3.0 / 2.0) * PI && i.high > (3.0 / 2.0) * PI {
            sin_low = -1.0;
        }
        return Interval::new(sin_low, sin_high);
    }

    pub fn cos(i: Interval) -> Self {
        assert!(0.0 <= i.low);
        assert!(i.high <= 2.0001 * PI);

        let mut cos_low = Float::cos(i.low);
        let mut cos_high = Float::cos(i.high);
        if cos_low > cos_high {
            std::mem::swap(&mut cos_low, &mut cos_high);
        }
        if i.low < PI && i.high > PI {
            cos_low = -1.0;
        }
        return Interval::new(cos_low, cos_high);
    }
}

impl From<Float> for Interval {
    fn from(f: Float) -> Self {
        return Interval::new(f, f);
    }
}

impl ops::Add<Interval> for Interval {
    type Output = Interval;
    fn add(self, rhs: Interval) -> Self::Output {
        return Interval::new(self.low + rhs.low, self.high + rhs.high);
    }
}

impl ops::Sub<Interval> for Interval {
    type Output = Interval;
    fn sub(self, rhs: Interval) -> Self::Output {
        return Interval::new(self.low - rhs.low, self.high - rhs.high);
    }
}

impl ops::Mul<Interval> for Interval {
    type Output = Interval;
    fn mul(self, rhs: Interval) -> Self::Output {
        let albl = self.low * rhs.low;
        let ahbl = self.high * rhs.low;
        let albh = self.low * rhs.high;
        let ahbh = self.high * rhs.high;
        let min = Float::min(Float::min(albl, ahbl), Float::min(albh, ahbh));
        let max = Float::max(Float::max(albl, ahbl), Float::max(albh, ahbh));
        return Interval::new(min, max);
    }
}

pub fn interval_find_zeros_(
    c1: Float,
    c2: Float,
    c3: Float,
    c4: Float,
    c5: Float,
    theta: Float,
    t: Interval,
    depth: u32,
    zeros: &mut Vec<Float>,
) {
    // Evaluate motion derivative in interval form, return if no zeros
    let range = Interval::from(c1)
        + (Interval::from(c2) + Interval::from(c3) * t)
            * Interval::cos(Interval::from(2.0 * theta) * t)
        + (Interval::from(c4) + Interval::from(c5) * t)
            * Interval::sin(Interval::from(2.0 * theta) * t);
    if range.low > 0.0 || range.high < 0.0 || range.low == range.high {
        return;
    }
    if depth > 0 {
        // Split the interval and recurse
        let tm = 0.5 * (t.low + t.high);
        interval_find_zeros_(
            c1,
            c2,
            c3,
            c4,
            c5,
            theta,
            Interval::new(t.low, tm),
            depth - 1,
            zeros,
        );
        interval_find_zeros_(
            c1,
            c2,
            c3,
            c4,
            c5,
            theta,
            Interval::new(tm, t.high),
            depth - 1,
            zeros,
        );
    } else {
        // Use Newton's method to refine zero
        let mut t_newton = (t.low + t.high) * 0.5;
        for _i in 0..4 {
            let f_newton = c1
                + (c2 + c3 * t_newton) * Float::cos(2.0 * theta * t_newton)
                + (c4 + c5 * t_newton) * Float::sin(2.0 * theta * t_newton);
            let f_prime_newton = (c3 + 2.0 * (c4 + c5 * t_newton) * theta)
                * Float::cos(2.0 * theta * t_newton)
                + (c5 - 2.0 * (c2 + c3 * t_newton) * theta) * Float::sin(2.0 * theta * t_newton);
            if f_newton == 0.0 || f_prime_newton == 0.0 {
                break;
            }
            t_newton -= f_newton / f_prime_newton;
        }
        if t_newton >= (t.low - 1e-3) && t_newton < (t.high + 1e-3) {
            zeros.push(t_newton);
        }
    }
}

pub fn interval_find_zeros(
    c1: Float,
    c2: Float,
    c3: Float,
    c4: Float,
    c5: Float,
    theta: Float,
    t: Interval,
    depth: u32,
) -> Vec<Float> {
    let mut zeros = Vec::new();
    interval_find_zeros_(c1, c2, c3, c4, c5, theta, t, depth, &mut zeros);
    return zeros;
}
