use std::sync::Arc;

use super::base_interaction::BaseInteraction;
use super::medium_interaction::MediumInteraction;
use super::surface_interaction::SurfaceInteraction;
use crate::core::pbrt::*;

#[derive(Clone, Debug)]
pub enum Interaction {
    Base(BaseInteraction),
    Surface(SurfaceInteraction),
    Medium(MediumInteraction),
}

impl Interaction {
    pub fn zero() -> Self {
        Interaction::Base(BaseInteraction::default())
    }

    pub fn is_surface_interaction(&self) -> bool {
        match self {
            Self::Base(inter) => {
                return inter.n.length_squared() > 0.0;
            }
            Self::Surface(inter) => {
                return inter.n.length_squared() > 0.0;
            }
            _ => {
                return false;
            }
        }
    }

    pub fn as_surface_interaction(&self) -> Option<&SurfaceInteraction> {
        match self {
            Self::Base(_) => None,
            Self::Surface(inter) => {
                if inter.n.length_squared() > 0.0 {
                    Some(inter)
                } else {
                    return None;
                }
            }
            Self::Medium(_) => None,
        }
    }

    pub fn as_medium_interaction(&self) -> Option<&MediumInteraction> {
        match self {
            Self::Base(_) => None,
            Self::Surface(_) => None,
            Self::Medium(inter) => Some(inter),
        }
    }

    pub fn as_surface_interaction_mut(&mut self) -> Option<&mut SurfaceInteraction> {
        match self {
            Self::Base(_) => None,
            Self::Surface(inter) => {
                if inter.n.length_squared() > 0.0 {
                    Some(inter)
                } else {
                    return None;
                }
            }
            Self::Medium(_) => None,
        }
    }

    pub fn as_medium_interaction_mut(&mut self) -> Option<&mut MediumInteraction> {
        match self {
            Self::Base(_) => None,
            Self::Surface(_) => None,
            Self::Medium(inter) => Some(inter),
        }
    }

    pub fn get_p(&self) -> Point3f {
        match self {
            Self::Base(inter) => inter.p,
            Self::Surface(inter) => inter.p,
            Self::Medium(inter) => inter.p,
        }
    }

    pub fn get_p_error(&self) -> Vector3f {
        match self {
            Self::Base(inter) => inter.p_error,
            Self::Surface(inter) => inter.p_error,
            Self::Medium(inter) => inter.p_error,
        }
    }

    pub fn get_n(&self) -> Normal3f {
        match self {
            Self::Base(inter) => inter.n,
            Self::Surface(inter) => inter.n,
            Self::Medium(inter) => inter.n,
        }
    }

    pub fn get_time(&self) -> Float {
        match self {
            Self::Base(inter) => inter.time,
            Self::Surface(inter) => inter.time,
            Self::Medium(inter) => inter.time,
        }
    }

    pub fn get_base_tuple(&self) -> (Point3f, Vector3f, Normal3f, Float) {
        match self {
            Self::Base(inter) => (inter.p, inter.p_error, inter.n, inter.time),
            Self::Surface(inter) => (inter.p, inter.p_error, inter.n, inter.time),
            Self::Medium(inter) => (inter.p, inter.p_error, inter.n, inter.time),
        }
    }

    pub fn spawn_ray(&self, d: &Vector3f) -> Ray {
        let (p, p_error, n, time) = self.get_base_tuple();
        let o = offset_ray_origin(&p, &p_error, &n, d);
        let mut r = Ray::new(&o, d, Float::INFINITY, time);
        r.medium = self.get_medium(d);
        return r;
    }

    pub fn spawn_ray_to(&self, it: &Interaction) -> Ray {
        let (ap, ap_error, an, atime) = self.get_base_tuple();
        let (bp, bp_error, bn, _btime) = it.get_base_tuple();
        let origin = offset_ray_origin(&ap, &ap_error, &an, &(bp - ap));
        let target = offset_ray_origin(&bp, &bp_error, &bn, &(origin - bp));
        let d = target - origin;
        let mut r = Ray::new(&origin, &d, 1.0 - SHADOW_EPSILON, atime);
        r.medium = self.get_medium(&d);
        return r;
    }

    pub fn get_medium(&self, w: &Vector3f) -> Option<Arc<dyn Medium>> {
        match self {
            Self::Base(_inter) => None,
            Self::Surface(inter) => inter.get_medium(w),
            Self::Medium(inter) => inter.get_medium(w),
        }
    }
}

impl Default for Interaction {
    fn default() -> Self {
        Interaction::zero()
    }
}

impl From<BaseInteraction> for Interaction {
    fn from(value: BaseInteraction) -> Self {
        Self::Base(value)
    }
}

impl From<SurfaceInteraction> for Interaction {
    fn from(value: SurfaceInteraction) -> Self {
        Self::Surface(value)
    }
}

impl From<&SurfaceInteraction> for Interaction {
    fn from(value: &SurfaceInteraction) -> Self {
        Self::Surface(value.clone())
    }
}

impl From<&mut SurfaceInteraction> for Interaction {
    fn from(value: &mut SurfaceInteraction) -> Self {
        Self::Surface(value.clone())
    }
}

impl From<MediumInteraction> for Interaction {
    fn from(value: MediumInteraction) -> Self {
        Self::Medium(value)
    }
}

impl From<&MediumInteraction> for Interaction {
    fn from(value: &MediumInteraction) -> Self {
        Self::Medium(value.clone())
    }
}

impl From<(Point3f, Float)> for Interaction {
    fn from(value: (Point3f, Float)) -> Self {
        Self::Base(BaseInteraction {
            p: value.0,
            time: value.1,
            p_error: Vector3f::zero(),
            n: Vector3f::new(0.0, 0.0, 1.0),
            wo: Vector3f::zero(),
            medium_interface: MediumInterface::new(),
        })
    }
}

impl From<(Point3f, Float, MediumInterface)> for Interaction {
    fn from(value: (Point3f, Float, MediumInterface)) -> Self {
        Self::Base(BaseInteraction {
            p: value.0,
            time: value.1,
            p_error: Vector3f::zero(),
            n: Vector3f::new(0.0, 0.0, 1.0),
            wo: Vector3f::zero(),
            medium_interface: value.2,
        })
    }
}

impl From<(Point3f, Vector3f, Normal3f, Float)> for Interaction {
    fn from(value: (Point3f, Vector3f, Normal3f, Float)) -> Self {
        Self::Base(BaseInteraction {
            p: value.0,
            time: value.3,
            p_error: value.1,
            n: value.2,
            wo: Vector3f::zero(),
            medium_interface: MediumInterface::new(),
        })
    }
}

/*
impl
    From<(
        Point3f,
        Vector3f,
        Normal3f,
        Vector3f,
        Float,
        MediumInterface,
    )> for Interaction
{
    fn from(
        value: (
            Point3f,
            Vector3f,
            Normal3f,
            Vector3f,
            Float,
            MediumInterface,
        ),
    ) -> Self {
        Self::Base(BaseInteraction {
            p: value.0,
            time: value.4,
            p_error: value.1,
            n: value.2,
            wo: value.3,
            medium_interface: value.5,
        })
    }
}
*/