use crate::core::pbrt::*;

#[derive(Clone, Debug)]
pub struct BaseInteraction {
    pub p: Point3f,
    pub p_error: Vector3f,
    pub wo: Vector3f,
    pub n: Normal3f,
    pub time: Float,
    pub medium_interface: MediumInterface,
}

impl Default for BaseInteraction {
    fn default() -> Self {
        BaseInteraction {
            p: Point3f::default(),
            p_error: Vector3f::default(),
            wo: Vector3f::default(),
            n: Normal3f::new(0.0, 0.0, 1.0),
            time: 0.0,
            medium_interface: MediumInterface::default(),
        }
    }
}

impl From<&Ray> for BaseInteraction {
    fn from(ray: &Ray) -> Self {
        BaseInteraction {
            p: ray.o,
            time: ray.time,
            p_error: Vector3f::zero(),
            n: Vector3f::zero(),
            wo: Vector3f::zero(),
            medium_interface: MediumInterface::from(&ray.medium),
        }
    }
}
