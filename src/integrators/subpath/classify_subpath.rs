use super::vertex::*;
use super::vertex_interaction::*;
use crate::core::pbrt::*;
use std::ops::Deref;
use std::sync::Arc;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum VertexClass {
    Any,
    Camera,
    Light,
    InfiniteLight,
    DistantLight,
    Surface,
    TransmissionSurface,
    SpecularSurface,
}

fn is_vertex_class(vc: VertexClass, vertex: &Vertex) -> bool {
    match vc {
        VertexClass::Any => true,
        VertexClass::Camera => vertex.get_camera().is_some(),
        VertexClass::Light => vertex.get_light().is_some(),
        VertexClass::InfiniteLight => {
            if let Some(light) = vertex.get_light() {
                let flags = light.get_light_flags();
                return (flags & LightFlags::Infinite as u32) != 0;
            }
            return false;
        }
        VertexClass::DistantLight => {
            if let Some(light) = vertex.get_light() {
                let flags = light.get_light_flags();
                return (flags & LightFlags::DeltaDirection as u32) != 0;
            }
            return false;
        }
        VertexClass::Surface => {
            let interaction = vertex.interaction.value.read().unwrap();
            if let VertexInteraction::Surface(_si) = interaction.deref() {
                return true;
            }
            return false;
        }
        VertexClass::TransmissionSurface => {
            let interaction = vertex.interaction.value.read().unwrap();
            if let VertexInteraction::Surface(si) = interaction.deref() {
                if let Some(bsdf) = si.bsdf.as_ref() {
                    return bsdf.has_components(BSDF_TRANSMISSION);
                }
            }
            return false;
        }
        VertexClass::SpecularSurface => {
            let interaction = vertex.interaction.value.read().unwrap();
            if let VertexInteraction::Surface(si) = interaction.deref() {
                if let Some(bsdf) = si.bsdf.as_ref() {
                    return bsdf.has_components(BSDF_SPECULAR);
                }
            }
            return false;
        }
    }
}

pub fn classify_subpath(
    camera_classes: &[VertexClass],
    light_classes: &[VertexClass],
    camera_vertices: &[Arc<Vertex>],
    light_vertices: &[Arc<Vertex>],
) -> bool {
    if camera_vertices.len() < camera_classes.len() {
        return false;
    }
    if light_vertices.len() < light_classes.len() {
        return false;
    }
    for i in 0..camera_classes.len() {
        if !is_vertex_class(camera_classes[i], &*camera_vertices[i]) {
            return false;
        }
    }
    for i in 0..light_classes.len() {
        if !is_vertex_class(light_classes[i], &*light_vertices[i]) {
            return false;
        }
    }
    return true;
}
