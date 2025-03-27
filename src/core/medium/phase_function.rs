use crate::core::base::*;

pub trait PhaseFunction: std::fmt::Display {
    fn p(&self, wo: &Vector3f, wi: &Vector3f) -> Float;
    fn sample_p(&self, wo: &Vector3f, u: &Point2f) -> (Float, Vector3f);
}

// Media Inline Functions
#[inline]
pub fn phase_hg(cos_theta: Float, g: Float) -> Float {
    let denom = 1.0 + g * g + 2.0 * g * cos_theta;
    return INV_4_PI * (1.0 - g * g) / (denom * Float::sqrt(denom));
}
