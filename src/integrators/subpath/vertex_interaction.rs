use crate::core::prelude::*;

use std::sync::Arc;

pub type CameraInteraction = (Interaction, Option<Arc<dyn Camera>>);
pub type LightInteraction = (Interaction, Option<Arc<dyn Light>>);

#[derive(Clone)]
pub enum EndpointInteraction {
    Camera(CameraInteraction),
    Light(LightInteraction),
}

impl EndpointInteraction {
    pub fn is_camera(&self) -> bool {
        match self {
            Self::Camera(_) => true,
            _ => false,
        }
    }

    pub fn is_light(&self) -> bool {
        match self {
            Self::Light(_) => true,
            _ => false,
        }
    }

    pub fn is_delta_light(&self) -> bool {
        match self {
            Self::Light(inter) => {
                if let Some(light) = inter.1.as_ref() {
                    return light.is_delta();
                }
                return false;
            }
            _ => false,
        }
    }

    // Vertex::IsInfiniteLight
    pub fn is_infinite_light(&self) -> bool {
        match self {
            Self::Light(inter) => {
                if let Some(light) = inter.1.as_ref() {
                    let flags = light.get_light_flags();
                    return ((flags & LightFlags::Infinite as u32) != 0)
                        || ((flags & LightFlags::DeltaDirection as u32) != 0);
                } else {
                    return true;
                }
            }
            _ => false,
        }
    }

    pub fn as_camera(&self) -> Option<&CameraInteraction> {
        match self {
            Self::Camera(inter) => Some(inter),
            _ => None,
        }
    }

    pub fn as_light(&self) -> Option<&LightInteraction> {
        match self {
            Self::Light(inter) => Some(inter),
            _ => None,
        }
    }

    pub fn get_camera(&self) -> Option<Arc<dyn Camera>> {
        match self {
            Self::Camera(inter) => {
                return inter.1.clone();
            }
            _ => None,
        }
    }

    pub fn get_light(&self) -> Option<Arc<dyn Light>> {
        match self {
            Self::Light(inter) => {
                return inter.1.clone();
            }
            _ => None,
        }
    }

    pub fn from_camera_ray(camera: &Arc<dyn Camera>, ray: &Ray) -> Self {
        let mut inter = BaseInteraction::default();
        inter.p = ray.o;
        inter.time = ray.time;
        inter.medium_interface = MediumInterface::from(&ray.medium);
        let inter = Interaction::Base(inter);
        return Self::Camera((inter, Some(camera.clone())));
    }

    pub fn from_camera_interaction(camera: &Arc<dyn Camera>, inter: &Interaction) -> Self {
        return Self::Camera((inter.clone(), Some(camera.clone())));
    }

    pub fn from_light_ray(light: &Arc<dyn Light>, ray: &Ray, n_light: &Normal3f) -> Self {
        let mut inter = BaseInteraction::default();
        inter.p = ray.o;
        inter.time = ray.time;
        inter.n = *n_light;
        inter.medium_interface = MediumInterface::from(&ray.medium);
        let inter = Interaction::Base(inter);
        return Self::Light((inter, Some(light.clone())));
    }

    pub fn from_light_interaction(light: &Arc<dyn Light>, inter: &Interaction) -> Self {
        return Self::Light((inter.clone(), Some(light.clone())));
    }

    pub fn from_ray(ray: &Ray) -> Self {
        let p = ray.position(1.0);
        let n = -ray.d;
        let mut inter = BaseInteraction::default();
        inter.p = p;
        inter.time = ray.time;
        inter.n = n;
        inter.medium_interface = MediumInterface::from(&ray.medium);
        let inter = Interaction::Base(inter);
        return Self::Light((inter, None));
    }

    pub fn get_p(&self) -> Point3f {
        match self {
            Self::Camera(inter) => inter.0.get_p(),
            Self::Light(inter) => inter.0.get_p(),
        }
    }

    pub fn get_time(&self) -> Float {
        match self {
            Self::Camera(inter) => inter.0.get_time(),
            Self::Light(inter) => inter.0.get_time(),
        }
    }

    pub fn get_n(&self) -> Normal3f {
        match self {
            Self::Camera(inter) => inter.0.get_n(),
            Self::Light(inter) => inter.0.get_n(),
        }
    }
}

impl Default for EndpointInteraction {
    fn default() -> Self {
        return Self::Camera((Interaction::default(), None));
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum VertexType {
    Camera,
    Light,
    Surface,
    Medium,
}

#[derive(Clone)]
pub enum VertexInteraction {
    EndPoint(EndpointInteraction),
    Medium(MediumInteraction),
    Surface(SurfaceInteraction),
}

impl Default for VertexInteraction {
    fn default() -> Self {
        return Self::EndPoint(EndpointInteraction::default());
    }
}

impl VertexInteraction {
    pub fn from_camera_ray(camera: &Arc<dyn Camera>, ray: &Ray) -> Self {
        return Self::EndPoint(EndpointInteraction::from_camera_ray(camera, ray));
    }

    pub fn from_camera_interaction(camera: &Arc<dyn Camera>, inter: &Interaction) -> Self {
        return Self::EndPoint(EndpointInteraction::from_camera_interaction(camera, inter));
    }

    pub fn from_light_ray(light: &Arc<dyn Light>, ray: &Ray, n_light: &Normal3f) -> Self {
        return Self::EndPoint(EndpointInteraction::from_light_ray(light, ray, n_light));
    }

    pub fn get_type(&self) -> VertexType {
        match self {
            Self::EndPoint(inter) => {
                if inter.is_camera() {
                    return VertexType::Camera;
                } else {
                    return VertexType::Light;
                }
            }
            Self::Medium(_) => return VertexType::Medium,
            Self::Surface(_) => return VertexType::Surface,
        }
    }

    pub fn as_camera(&self) -> Option<&CameraInteraction> {
        match self {
            Self::EndPoint(inter) => inter.as_camera(),
            _ => None,
        }
    }

    pub fn as_light(&self) -> Option<&LightInteraction> {
        match self {
            Self::EndPoint(inter) => inter.as_light(),
            _ => None,
        }
    }

    pub fn as_surface(&self) -> Option<&SurfaceInteraction> {
        match self {
            Self::Surface(inter) => Some(inter),
            _ => None,
        }
    }

    pub fn as_medium(&self) -> Option<&MediumInteraction> {
        match self {
            Self::Medium(inter) => Some(inter),
            _ => None,
        }
    }

    pub fn get_p(&self) -> Point3f {
        match self {
            Self::EndPoint(inter) => inter.get_p(),
            Self::Medium(inter) => inter.p,
            Self::Surface(inter) => inter.p,
        }
    }

    pub fn get_time(&self) -> Float {
        match self {
            Self::EndPoint(inter) => inter.get_time(),
            Self::Medium(inter) => inter.time,
            Self::Surface(inter) => inter.time,
        }
    }

    pub fn get_ng(&self) -> Normal3f {
        match self {
            Self::EndPoint(inter) => inter.get_n(),
            Self::Medium(inter) => inter.n,
            Self::Surface(inter) => inter.n,
        }
    }

    pub fn get_ns(&self) -> Normal3f {
        match self {
            Self::EndPoint(inter) => inter.get_n(),
            Self::Medium(inter) => inter.n,
            Self::Surface(inter) => inter.shading.n,
        }
    }

    pub fn get_camera(&self) -> Option<Arc<dyn Camera>> {
        match self {
            Self::EndPoint(inter) => inter.get_camera(),
            _ => None,
        }
    }

    pub fn get_light(&self) -> Option<Arc<dyn Light>> {
        match self {
            Self::EndPoint(inter) => inter.get_light(),
            Self::Surface(si) => {
                if let Some(primitive) = &si.primitive {
                    let primitive = primitive.upgrade().unwrap();
                    return primitive.get_area_light();
                }
                return None;
            }
            _ => None,
        }
    }

    pub fn get_interaction(&self) -> Interaction {
        match self {
            Self::EndPoint(inter) => match inter {
                EndpointInteraction::Camera(inter) => inter.0.clone(),
                EndpointInteraction::Light(inter) => inter.0.clone(),
            },
            Self::Medium(inter) => Interaction::Medium(inter.clone()),
            Self::Surface(inter) => Interaction::Surface(inter.clone()),
        }
    }
}
