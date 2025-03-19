use crate::core::medium::*;
use crate::core::transform::*;

#[derive(Clone, Default, Debug)]
pub struct BaseLight {
    pub flags: u32,
    pub n_samples: u32,
    pub medium_interface: MediumInterface,
    pub light_to_world: Transform,
    pub world_to_light: Transform,
}

impl BaseLight {
    pub fn new(
        flags: u32,
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        n_samples: u32,
    ) -> Self {
        BaseLight {
            flags,
            n_samples,
            medium_interface: medium_interface.clone(),
            light_to_world: light_to_world.clone(),
            world_to_light: Transform::inverse(light_to_world),
        }
    }
}
