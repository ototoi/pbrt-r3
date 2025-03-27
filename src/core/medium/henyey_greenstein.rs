use super::phase_function::*;
use crate::core::base::*;
use crate::core::geometry::*;
use crate::core::profile::*;

pub struct HenyeyGreenstein {
    g: Float,
}

impl HenyeyGreenstein {
    pub fn new(g: Float) -> Self {
        HenyeyGreenstein { g }
    }
}

impl PhaseFunction for HenyeyGreenstein {
    fn p(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        let _p = ProfilePhase::new(Prof::PhaseFuncEvaluation);
        return phase_hg(Vector3f::dot(wo, wi), self.g);
    }
    fn sample_p(&self, wo: &Vector3f, u: &Point2f) -> (Float, Vector3f) {
        let _p = ProfilePhase::new(Prof::PhaseFuncSampling);
        // Compute $\cos \theta$ for Henyey--Greenstein sample
        let g = self.g;
        let cos_theta = if g < 1e-3 {
            1.0 - 2.0 * u[0]
        } else {
            let sqr_term = (1.0 - g * g) / (1.0 + g - 2.0 * g * u[0]);
            -(1.0 + g * g - sqr_term * sqr_term) / (2.0 * g)
        };

        // Compute direction _wi_ for Henyey--Greenstein sample
        let sin_theta = Float::max(0.0, 1.0 - cos_theta * cos_theta).sqrt();
        let phi = 2.0 * PI * u[1];
        let (v1, v2) = coordinate_system(&wo);
        let wi = spherical_direction_axes(sin_theta, cos_theta, phi, &v1, &v2, &wo);
        return (phase_hg(cos_theta, g), wi);
    }
}

impl std::fmt::Display for HenyeyGreenstein {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        return write!(f, "[ HenyeyGreenstein g: {} ]", &self.g);
    }
}
