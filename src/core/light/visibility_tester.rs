use crate::core::base::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::medium::medium::Medium;
use crate::core::medium::medium_interface::MediumInterface;
use crate::core::geometry::ray::Ray;
use crate::core::sampler::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use std::sync::Arc;

#[derive(Clone)]
pub struct VisibilityEndpoint {
    pub p: Point3f,
    pub p_error: Vector3f,
    pub n: Normal3f,
    pub time: Float,
    pub medium_interface: MediumInterface,
}

impl Default for VisibilityEndpoint {
    fn default() -> Self {
        Self {
            p: Point3f::zero(),
            p_error: Vector3f::zero(),
            n: Normal3f::zero(),
            time: 0.0,
            medium_interface: MediumInterface::default(),
        }
    }
}

impl VisibilityEndpoint {
    #[inline]
    pub fn get_medium(&self, w: &Vector3f) -> Option<Arc<dyn Medium>> {
        if Vector3f::dot(w, &self.n) > 0.0 {
            self.medium_interface.get_outside()
        } else {
            self.medium_interface.get_inside()
        }
    }

    #[inline]
    pub fn spawn_ray_to(&self, it: &VisibilityEndpoint) -> Ray {
        let origin = offset_ray_origin(&self.p, &self.p_error, &self.n, &(it.p - self.p));
        let target = offset_ray_origin(&it.p, &it.p_error, &it.n, &(origin - it.p));
        let d = target - origin;
        let mut r = Ray::new(&origin, &d, 1.0 - SHADOW_EPSILON, self.time);
        r.medium = self.get_medium(&d);
        r
    }

    #[inline]
    pub fn to_interaction(&self) -> Interaction {
        Interaction::Base(BaseInteraction {
            p: self.p,
            p_error: self.p_error,
            wo: Vector3f::zero(),
            n: self.n,
            time: self.time,
            medium_interface: self.medium_interface.clone(),
        })
    }
}

impl From<&BaseInteraction> for VisibilityEndpoint {
    fn from(inter: &BaseInteraction) -> Self {
        Self {
            p: inter.p,
            p_error: inter.p_error,
            n: inter.n,
            time: inter.time,
            medium_interface: inter.medium_interface.clone(),
        }
    }
}

impl From<&Interaction> for VisibilityEndpoint {
    fn from(inter: &Interaction) -> Self {
        let base = BaseInteraction::from(inter);
        Self::from(&base)
    }
}

#[derive(Clone)]
pub struct VisibilityTester {
    pub p0: VisibilityEndpoint,
    pub p1: VisibilityEndpoint,
}

impl VisibilityTester {
    pub fn new() -> Self {
        VisibilityTester {
            p0: VisibilityEndpoint::default(),
            p1: VisibilityEndpoint::default(),
        }
    }

    pub fn unoccluded(&self, scene: &Scene) -> bool {
        return !scene.intersect_p(&self.p0.spawn_ray_to(&self.p1));
    }

    pub fn tr(&self, scene: &Scene, sampler: &mut dyn Sampler) -> Spectrum {
        let mut ray = self.p0.spawn_ray_to(&self.p1);

        let mut tr = Spectrum::one();

        // Fast path: when the segment starts in vacuum, skip medium work as long
        // as we don't encounter a medium transition boundary.
        if ray.medium.is_none() {
            loop {
                let Some(isect) = scene.intersect(&ray) else {
                    return tr;
                };
                if let Some(primitive) = isect.get_primitive() {
                    if primitive.get_material().is_some() {
                        return Spectrum::zero();
                    }
                }
                let tisect = BaseInteraction::from(&isect);
                let next_ray = VisibilityTester::spawn_ray_from_base_to_endpoint(&tisect, &self.p1);
                if next_ray.medium.is_some() {
                    ray = next_ray;
                    break;
                }
                ray = next_ray;
            }
        }

        loop {
            if let Some(isect) = scene.intersect(&ray) {
                if let Some(primitive) = isect.get_primitive() {
                    if primitive.get_material().is_some() {
                        return Spectrum::zero();
                    }
                }

                if let Some(medium) = ray.medium.as_ref() {
                    tr *= medium.as_ref().tr(&ray, sampler);
                }

                let tisect = BaseInteraction::from(&isect);
                ray = VisibilityTester::spawn_ray_from_base_to_endpoint(&tisect, &self.p1);
            } else {
                if let Some(medium) = ray.medium.as_ref() {
                    tr *= medium.as_ref().tr(&ray, sampler);
                }
                break;
            }
        }
        return tr;
    }

    #[inline]
    fn spawn_ray_from_base_to_endpoint(from: &BaseInteraction, to: &VisibilityEndpoint) -> Ray {
        let origin = offset_ray_origin(&from.p, &from.p_error, &from.n, &(to.p - from.p));
        let target = offset_ray_origin(&to.p, &to.p_error, &to.n, &(origin - to.p));
        let d = target - origin;
        let mut r = Ray::new(&origin, &d, 1.0 - SHADOW_EPSILON, from.time);
        r.medium = from.get_medium(&d);
        r
    }
}

impl From<(&Interaction, &Interaction)> for VisibilityTester {
    fn from(value: (&Interaction, &Interaction)) -> Self {
        VisibilityTester {
            p0: VisibilityEndpoint::from(value.0),
            p1: VisibilityEndpoint::from(value.1),
        }
    }
}

impl From<(Interaction, Interaction)> for VisibilityTester {
    fn from(value: (Interaction, Interaction)) -> Self {
        VisibilityTester {
            p0: VisibilityEndpoint::from(&value.0),
            p1: VisibilityEndpoint::from(&value.1),
        }
    }
}

impl From<(BaseInteraction, BaseInteraction)> for VisibilityTester {
    fn from(value: (BaseInteraction, BaseInteraction)) -> Self {
        VisibilityTester {
            p0: VisibilityEndpoint::from(&value.0),
            p1: VisibilityEndpoint::from(&value.1),
        }
    }
}
