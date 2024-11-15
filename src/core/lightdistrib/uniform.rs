use crate::core::pbrt::*;
use std::sync::Arc;
use std::sync::RwLock;

pub struct UniformLightDistribution {
    distrib: Arc<Distribution1D>,
}

impl UniformLightDistribution {
    pub fn new(scene: &Scene) -> Self {
        let prob = vec![1.0; scene.lights.len()];
        UniformLightDistribution {
            distrib: Arc::new(Distribution1D::new(&prob)),
        }
    }
}

impl LightDistribution for UniformLightDistribution {
    fn lookup(&self, _p: &Point3f) -> Arc<Distribution1D> {
        return self.distrib.clone();
    }
}
