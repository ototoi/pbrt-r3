use crate::core::base::*;
use crate::core::geometry::*;
use crate::core::medium::*;
use crate::core::interaction::*;
use std::sync::Arc;

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
            p: Point3f::zero(),
            p_error: Vector3f::zero(),
            wo: Vector3f::zero(),
            n: Normal3f::zero(),
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
            n: Normal3f::zero(),
            wo: Vector3f::zero(),
            medium_interface: MediumInterface::from(&ray.medium),
        }
    }
}

impl BaseInteraction {
    pub fn get_medium(&self, w: &Vector3f) -> Option<Arc<dyn Medium>> {
        if Vector3f::dot(w, &self.n) > 0.0 {
            self.medium_interface.get_outside()
        } else {
            self.medium_interface.get_inside()
        }
    }

    pub fn spawn_ray_to(&self, it: &BaseInteraction) -> Ray {
        let origin = offset_ray_origin(&self.p, &self.p_error, &self.n, &(it.p - self.p));
        let target = offset_ray_origin(&it.p, &it.p_error, &it.n, &(origin - it.p));
        let d = target - origin;
        let mut r = Ray::new(&origin, &d, 1.0 - SHADOW_EPSILON, self.time);
        r.medium = self.get_medium(&d);
        r
    }
}

impl From<&SurfaceInteraction> for BaseInteraction {
    fn from(inter: &SurfaceInteraction) -> Self {
        BaseInteraction {
            p: inter.p,
            p_error: inter.p_error,
            wo: inter.wo,
            n: inter.n,
            time: inter.time,
            medium_interface: inter.medium_interface.clone(),
        }
    }
}

impl From<&MediumInteraction> for BaseInteraction {
    fn from(inter: &MediumInteraction) -> Self {
        BaseInteraction {
            p: inter.p,
            p_error: inter.p_error,
            wo: inter.wo,
            n: inter.n,
            time: inter.time,
            medium_interface: inter.medium_interface.clone(),
        }
    }
}

impl From<&Interaction> for BaseInteraction {
    fn from(inter: &Interaction) -> Self {
        match inter {
            Interaction::Base(base) => base.clone(),
            Interaction::Surface(surface) => BaseInteraction::from(surface),
            Interaction::Medium(medium) => BaseInteraction::from(medium),
        }
    }
}
