use crate::core::camera::*;
use crate::core::distribution::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::lightdistrib::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use crate::core::texture::*;
use crate::core::transform::*;

use std::sync::Arc;

/*
    const Spectrum sigma_a, sigma_s, sigma_t;
    const Float g;
*/
#[derive(Debug, Clone)]
pub struct HomogeneousMedium {
    pub sigma_a: Spectrum,
    pub sigma_s: Spectrum,
    pub sigma_t: Spectrum,
    pub g: Float,
}

impl HomogeneousMedium {
    pub fn new(sigma_a: &Spectrum, sigma_s: &Spectrum, g: Float) -> Self {
        HomogeneousMedium {
            sigma_a: *sigma_a,
            sigma_s: *sigma_s,
            sigma_t: *sigma_a + *sigma_s,
            g,
        }
    }
}

impl Medium for HomogeneousMedium {
    fn tr(&self, ray: &Ray, _sampler: &mut dyn Sampler) -> Spectrum {
        let _p = ProfilePhase::new(Prof::MediumTr);
        let sigma_t = self.sigma_t;
        let ray_length = ray.d.length();
        let tmax = ray.t_max.get();
        return (-sigma_t * Float::min(ray_length * tmax, std::f32::MAX as Float)).exp();
    }
    fn sample(
        &self,
        ray: &Ray,
        sampler: &mut dyn Sampler,
        _arena: &mut MemoryArena,
    ) -> (Spectrum, Option<MediumInteraction>) {
        let _p = ProfilePhase::new(Prof::MediumSample);
        // Sample a channel and distance along the ray
        let channel = usize::min(
            (sampler.get_1d() * Spectrum::N_SAMPLES as Float) as usize,
            Spectrum::N_SAMPLES - 1,
        );
        let dist = -Float::ln(1.0 - sampler.get_1d()) / self.sigma_t[channel];
        let t_max = ray.t_max.get();
        let t = Float::min(dist / ray.d.length(), t_max);
        let sampled_medium = t < t_max;

        let mi = if sampled_medium {
            let phase: Arc<dyn PhaseFunction> = Arc::new(HenyeyGreenstein::new(self.g));
            Some(MediumInteraction::new(
                &ray.position(t),
                &(-ray.d),
                ray.time,
                &phase,
            ))
        } else {
            None
        };

        // Compute the transmittance and sampling density
        let tr = (-self.sigma_t * Float::min(t, std::f32::MAX as Float) * ray.d.length()).exp();
        // Return weighting factor for scattering from homogeneous medium
        let density = if sampled_medium {
            self.sigma_t * tr
        } else {
            tr
        };
        let mut pdf = 0.0;
        for i in 0..Spectrum::N_SAMPLES {
            pdf += density[i];
        }
        pdf /= Spectrum::N_SAMPLES as Float;

        if pdf == 0.0 {
            pdf = 1.0;
        }

        let spec = if sampled_medium {
            tr * self.sigma_s / pdf
        } else {
            tr / pdf
        };
        return (spec, mi);
    }
}
