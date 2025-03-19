use super::lightdistrib::*;
use super::power::*;
use super::spatial::*;
use super::uniform::*;
use crate::core::error::PbrtError;
use crate::core::scene::*;

use log::*;
use std::sync::Arc;

pub fn create_light_sample_distribution(
    name: &str,
    scene: &Scene,
) -> Result<Arc<dyn LightDistribution>, PbrtError> {
    if scene.lights.is_empty() {
        let msg = format!(
            "Light sample distribution type \"{}\" cannot create since no light.",
            name
        );
        return Err(PbrtError::error(&msg));
    }
    match name {
        "uniform" => {
            if scene.lights.len() != 1 {
                return create_light_sample_distribution("spatial", scene);
            } else {
                return Ok(Arc::new(UniformLightDistribution::new(scene)));
            }
        }
        "power" => {
            return Ok(Arc::new(PowerLightDistribution::new(scene)));
        }
        "spatial" => {
            let max_voxels = 64;
            return Ok(Arc::new(SpatialLightDistribution::new(scene, max_voxels)));
        }
        s => {
            warn!(
                "Light sample distribution type \"{}\" unknown. Using \"spatial\".",
                s
            );
            return create_light_sample_distribution("spatial", scene);
        }
    }
}
