use super::lightdistrib::*;
use crate::core::pbrt::*;
use crate::core::sampling::*;
use crate::core::scene::*;

use std::sync::Arc;

// Ported from integrator.cpp
pub fn compute_light_power_distribution(scene: &Scene) -> Arc<Distribution1D> {
    assert!(!scene.lights.is_empty());
    let mut light_power = Vec::new();
    for light in scene.lights.iter() {
        let light = light.as_ref();
        light_power.push(light.power().y());
    }
    return Arc::new(Distribution1D::new(&light_power));
}

pub struct PowerLightDistribution {
    distrib: Arc<Distribution1D>,
}

impl PowerLightDistribution {
    pub fn new(scene: &Scene) -> Self {
        PowerLightDistribution {
            distrib: compute_light_power_distribution(scene),
        }
    }
}

impl LightDistribution for PowerLightDistribution {
    fn lookup(&self, _p: &Point3f) -> Arc<Distribution1D> {
        return self.distrib.clone();
    }
}
