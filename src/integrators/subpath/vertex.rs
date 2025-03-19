use super::vertex_interaction::*;
use crate::core::camera::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::material::*;
use crate::core::pbrt::*;
use crate::core::refrection::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::spectrum::*;

use std::collections::HashMap;
use std::ops::Deref;
use std::sync::Arc;
use std::sync::RwLock;

//pub type LightIndexMap<'a> = HashMap<&'a dyn Light, usize>;
//pub type LightIndexMap = HashMap<*const dyn Light, usize>;
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct LightKeyType {
    ptr: *const dyn Light,
}
impl From<*const dyn Light> for LightKeyType {
    fn from(ptr: *const dyn Light) -> Self {
        LightKeyType { ptr }
    }
}
unsafe impl Send for LightKeyType {}
unsafe impl Sync for LightKeyType {}

pub type LightIndexMap = HashMap<LightKeyType, usize>;

/// Forward declaration (correction term for adjoint BSDF with shading normals)
pub fn correct_shading_normal(
    isect: &SurfaceInteraction,
    wo: &Vector3f,
    wi: &Vector3f,
    mode: TransportMode,
) -> Float {
    if mode == TransportMode::Importance {
        let num = Vector3f::abs_dot(&wo, &isect.shading.n) * Vector3f::abs_dot(&wi, &isect.n);
        let denom = Vector3f::abs_dot(&wo, &isect.n) * Vector3f::abs_dot(&wi, &isect.shading.n);
        if denom == 0.0 {
            return 0.0;
        }
        return num / denom;
    } else {
        return 1.0;
    }
}

pub fn infinite_light_density(
    scene: &Scene,
    light_distr: &Distribution1D,
    light_to_index: &LightIndexMap,
    w: &Vector3f,
) -> Float {
    let mut pdf = 0.0;
    for light in scene.infinite_lights.iter() {
        let light = light.as_ref();
        let light_ptr = light as *const dyn Light;
        let key = LightKeyType::from(light_ptr);
        let index = light_to_index.get(&key).unwrap();
        let index = *index;
        assert!(index < light_distr.func.len());
        let pdf_i = light.pdf_li(&Interaction::default(), &-*w) * light_distr.func[index];
        assert!(pdf_i >= 0.0);
        pdf += pdf_i;
    }
    return pdf / (light_distr.func_int * light_distr.count() as Float);
}

#[derive(Clone, Debug)]
pub struct VertexValue<T: Clone> {
    pub value: Arc<RwLock<T>>,
}

impl<T: Clone> VertexValue<T> {
    pub fn new(value: T) -> Self {
        VertexValue {
            value: Arc::new(RwLock::new(value)),
        }
    }

    pub fn get(&self) -> T {
        return self.value.read().unwrap().clone();
    }

    pub fn set(&self, value: T) {
        *self.value.write().unwrap() = value;
    }
}

#[derive(Clone)]
pub struct Vertex {
    pub beta: VertexValue<Spectrum>,
    pub interaction: VertexValue<VertexInteraction>,
    pub delta: VertexValue<bool>,
    pub pdf_fwd: VertexValue<Float>,
    pub pdf_rev: VertexValue<Float>,
}

impl Vertex {
    pub fn new(
        beta: Spectrum,
        interaction: VertexInteraction,
        delta: bool,
        pdf_fwd: Float,
        pdf_rev: Float,
    ) -> Self {
        Vertex {
            beta: VertexValue::new(beta),
            interaction: VertexValue::new(interaction),
            delta: VertexValue::new(delta),
            pdf_fwd: VertexValue::new(pdf_fwd),
            pdf_rev: VertexValue::new(pdf_rev),
        }
    }

    pub fn as_tuple(
        &self,
    ) -> (
        Arc<RwLock<Spectrum>>,
        Arc<RwLock<VertexInteraction>>,
        Arc<RwLock<bool>>,
        Arc<RwLock<Float>>,
        Arc<RwLock<Float>>,
    ) {
        return (
            self.beta.value.clone(),
            self.interaction.value.clone(),
            self.delta.value.clone(),
            self.pdf_fwd.value.clone(),
            self.pdf_rev.value.clone(),
        );
    }

    pub fn get_type(&self) -> VertexType {
        let interaction = self.interaction.value.read().unwrap();
        return interaction.get_type();
    }

    pub fn convert_density(&self, pdf: Float, next: &Vertex) -> Float {
        // Return solid angle density if _next_ is an infinite area light
        if next.is_infinite_light() {
            return pdf;
        }
        let w = next.get_p() - self.get_p();
        if w.length_squared() == 0.0 {
            return 0.0;
        }
        let inv_dist2 = 1.0 / w.length_squared();
        let mut pdf = pdf;
        if next.is_on_surface() {
            pdf *= Vector3f::abs_dot(&next.get_ng(), &(w * Float::sqrt(inv_dist2)));
        }
        return pdf * inv_dist2;
    }

    pub fn pdf(&self, scene: &Scene, prev: &Option<Arc<Vertex>>, next: &Vertex) -> Float {
        let t = self.get_type();
        if t == VertexType::Light {
            return self.pdf_light(scene, next);
        }

        let wn = next.get_p() - self.get_p();
        if wn.length_squared() == 0.0 {
            return 0.0;
        }
        let wn = wn.normalize();
        assert!(t != VertexType::Light);
        // Compute directional density depending on the vertex types
        let interaction = self.interaction.value.read().unwrap();
        match interaction.deref() {
            VertexInteraction::EndPoint(ei) => {
                if let EndpointInteraction::Camera(ei) = ei {
                    let ray = ei.0.spawn_ray(&wn);
                    let camera = ei.1.as_ref().unwrap();
                    if let Some((_pdf_pos, pdf_dir)) = camera.pdf_we(&ray) {
                        let pdf = self.convert_density(pdf_dir, next);
                        return pdf;
                    }
                }
                return 0.0;
            }
            VertexInteraction::Surface(si) => {
                let prev = prev.as_ref().unwrap();
                let wp = prev.get_p() - self.get_p();
                if wp.length_squared() == 0.0 {
                    return 0.0;
                }
                let wp = wp.normalize();

                let bsdf = si.bsdf.as_ref().unwrap();
                let pdf = bsdf.pdf(&wp, &wn, BSDF_ALL);
                let pdf = self.convert_density(pdf, next);
                return pdf;
            }
            VertexInteraction::Medium(mi) => {
                let prev = prev.as_ref().unwrap();
                let wp = prev.get_p() - self.get_p();
                if wp.length_squared() == 0.0 {
                    return 0.0;
                }
                let wp = wp.normalize();

                let pdf = mi.phase.p(&wp, &wn);
                let pdf = self.convert_density(pdf, next);
                return pdf;
            }
        }
    }

    //PdfLight
    pub fn pdf_light(&self, scene: &Scene, v: &Vertex) -> Float {
        let w = v.get_p() - self.get_p();
        let dist2 = w.length_squared();
        if dist2 == 0.0 {
            return 0.0;
        }
        let inv_dist2 = 1.0 / dist2;
        let w = w * Float::sqrt(inv_dist2);
        let mut pdf = 0.0;
        if self.is_infinite_light() {
            let (_world_center, world_radius) = scene.world_bound().bounding_sphere();
            pdf = 1.0 / (PI * world_radius * world_radius);
        } else {
            if let Some(light) = self.get_light() {
                let ray = Ray::new(&self.get_p(), &w, Float::INFINITY, self.get_time());
                if let Some((_pdf_pos, pdf_dir)) = light.pdf_le(&ray, &self.get_ng()) {
                    assert!(pdf_dir >= 0.0);
                    pdf = pdf_dir * inv_dist2;
                }
            }
        }
        if v.is_on_surface() {
            pdf *= Vector3f::abs_dot(&v.get_ng(), &w);
        }
        return pdf;
    }

    //PdfLightOrigin
    pub fn pdf_light_origin(
        &self,
        scene: &Scene,
        v: &Vertex,
        light_distr: &Distribution1D,
        light_to_index: &LightIndexMap,
    ) -> Float {
        let w = v.get_p() - self.get_p();
        if w.length_squared() == 0.0 {
            return 0.0;
        }
        let w = w.normalize();
        if self.is_infinite_light() {
            // Return solid angle density for infinite light sources
            let pdf = infinite_light_density(scene, light_distr, light_to_index, &w);
            assert!(pdf >= 0.0);
            return pdf;
        } else {
            // Return solid angle density for non-infinite light sources
            // Get pointer _light_ to the light source at the vertex
            if let Some(light) = self.get_light() {
                let light = light.as_ref();
                let light_ptr = light as *const dyn Light;
                let key = LightKeyType::from(light_ptr);
                let index = light_to_index.get(&LightKeyType::from(key)).unwrap();
                let index = *index;
                assert!(index < light_distr.func.len());
                let pdf_choice = light_distr.discrete_pdf(index);
                let ray = Ray::new(&self.get_p(), &w, Float::INFINITY, self.get_time());
                if let Some((pdf_pos, _pdf_dir)) = light.pdf_le(&ray, &self.get_ng()) {
                    return pdf_choice * pdf_pos;
                }
            }
            return 0.0;
        }
    }

    pub fn is_light(&self) -> bool {
        let interaction = self.interaction.value.read().unwrap();
        match interaction.deref() {
            VertexInteraction::EndPoint(ei) => {
                return ei.is_light();
            }
            VertexInteraction::Surface(si) => {
                let premitive = si.primitive.clone();
                if let Some(premitive) = premitive.as_ref() {
                    let premitive = premitive.upgrade().unwrap();
                    return premitive.get_area_light().is_some();
                }
                return false;
            }
            _ => {
                return false;
            }
        }
    }

    // Vertex::IsDeltaLight
    pub fn is_delta_light(&self) -> bool {
        let interaction = self.interaction.value.read().unwrap();
        match interaction.deref() {
            VertexInteraction::EndPoint(ei) => {
                return ei.is_delta_light();
            }
            _ => {
                return false;
            }
        }
    }

    // Vertex::IsInfiniteLight
    pub fn is_infinite_light(&self) -> bool {
        let interaction = self.interaction.value.read().unwrap();
        match interaction.deref() {
            VertexInteraction::EndPoint(ei) => {
                return ei.is_infinite_light();
            }
            _ => {
                return false;
            }
        }
    }

    pub fn is_connectible(&self) -> bool {
        let t = self.get_type();
        match t {
            VertexType::Medium => {
                return true;
            }
            VertexType::Light => {
                let interaction = self.interaction.value.read().unwrap();
                if let Some(light) = interaction.get_light() {
                    return !light.is_delta_direction();
                }
                return false;
            }
            VertexType::Camera => {
                return true;
            }
            VertexType::Surface => {
                let interaction = self.interaction.value.read().unwrap();
                let si = interaction.as_surface().unwrap();
                let bsdf = si.bsdf.as_ref().unwrap();
                return bsdf.num_components(
                    BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION,
                ) > 0;
            }
        }
    }

    pub fn is_on_surface(&self) -> bool {
        let ng = self.get_ng();
        return ng.length_squared() != 0.0;
    }

    pub fn f(&self, next: &Vertex, mode: TransportMode) -> Spectrum {
        let wi = next.get_p() - self.get_p();
        if wi.length_squared() == 0.0 {
            return Spectrum::zero();
        }
        let wi = wi.normalize();
        let interaction = self.interaction.value.read().unwrap();
        match interaction.deref() {
            VertexInteraction::Surface(si) => {
                let bsdf = si.bsdf.as_ref().unwrap();
                let f = bsdf.f(&si.wo, &wi, BSDF_ALL);
                let ns = correct_shading_normal(si, &si.wo, &wi, mode);
                return f * ns;
            }
            VertexInteraction::Medium(mi) => {
                let p = mi.phase.p(&mi.wo, &wi);
                return Spectrum::from(p);
            }
            _ => {
                return Spectrum::zero();
            }
        }
    }

    pub fn le(&self, scene: &Scene, v: &Vertex) -> Spectrum {
        if !self.is_light() {
            return Spectrum::zero();
        }
        let w = v.get_p() - self.get_p();
        if w.length_squared() == 0.0 {
            return Spectrum::zero();
        }
        let w = w.normalize();

        let interaction = self.interaction.value.read().unwrap();
        match interaction.deref() {
            VertexInteraction::EndPoint(ei) => {
                let mut le = Spectrum::zero();
                if ei.is_infinite_light() {
                    let ray =
                        RayDifferential::from(Ray::new(&self.get_p(), &-w, Float::INFINITY, 0.0));
                    // Return emitted radiance for infinite light sources
                    for light in scene.infinite_lights.iter() {
                        le += light.le(&ray);
                    }
                }
                return le;
            }
            VertexInteraction::Surface(si) => {
                let premitive = si.primitive.as_ref().unwrap().upgrade().unwrap();
                let area_light = premitive.get_area_light().unwrap();
                let area_light = area_light.as_area_light().unwrap();
                return area_light.l(&Interaction::from(si), &w);
            }
            _ => {
                return Spectrum::zero();
            }
        }
    }

    pub fn create_camera_from_ray(camera: &Arc<dyn Camera>, ray: &Ray, beta: &Spectrum) -> Self {
        let v = Vertex::new(
            beta.clone(),
            VertexInteraction::from_camera_ray(camera, ray),
            false,
            0.0,
            0.0,
        );
        assert!(v.get_type() == VertexType::Camera);
        return v;
    }

    pub fn create_camera_from_interaction(
        camera: &Arc<dyn Camera>,
        it: &Interaction,
        beta: &Spectrum,
    ) -> Self {
        let v = Vertex::new(
            beta.clone(),
            VertexInteraction::from_camera_interaction(camera, it),
            false,
            0.0,
            0.0,
        );
        assert!(v.get_type() == VertexType::Camera);
        return v;
    }

    pub fn create_light_from_ray(
        light: &Arc<dyn Light>,
        ray: &Ray,
        n_light: &Normal3f,
        le: &Spectrum,
        pdf: Float,
    ) -> Self {
        let v = Vertex::new(
            le.clone(),
            VertexInteraction::from_light_ray(light, ray, n_light),
            false,
            pdf,
            0.0,
        );
        assert!(v.get_type() == VertexType::Light);
        return v;
    }

    pub fn create_light_from_endpoint(
        ei: &EndpointInteraction,
        beta: &Spectrum,
        pdf: Float,
    ) -> Self {
        let v = Vertex::new(
            beta.clone(),
            VertexInteraction::EndPoint(ei.clone()),
            false,
            pdf,
            0.0,
        );
        assert!(v.get_type() == VertexType::Light);
        return v;
    }

    pub fn create_surface(
        si: &SurfaceInteraction,
        beta: &Spectrum,
        pdf: Float,
        prev: &Vertex,
    ) -> Self {
        let v = Vertex::new(
            beta.clone(),
            VertexInteraction::Surface(si.clone()),
            false,
            0.0,
            0.0,
        );
        v.pdf_fwd.set(prev.convert_density(pdf, &v));
        assert!(v.get_type() == VertexType::Surface);
        return v;
    }

    pub fn create_medium(
        mi: &MediumInteraction,
        beta: &Spectrum,
        pdf: Float,
        prev: &Vertex,
    ) -> Self {
        let v = Vertex::new(
            beta.clone(),
            VertexInteraction::Medium(mi.clone()),
            false,
            0.0,
            0.0,
        );
        v.pdf_fwd.set(prev.convert_density(pdf, &v));
        assert!(v.get_type() == VertexType::Medium);
        return v;
    }

    pub fn get_interaction(&self) -> Interaction {
        let interaction = self.interaction.value.read().unwrap();
        return interaction.get_interaction();
    }

    pub fn get_p(&self) -> Point3f {
        let interaction = self.interaction.value.read().unwrap();
        return interaction.get_p();
    }

    pub fn get_time(&self) -> Float {
        let interaction = self.interaction.value.read().unwrap();
        return interaction.get_time();
    }

    pub fn get_ng(&self) -> Normal3f {
        let interaction = self.interaction.value.read().unwrap();
        return interaction.get_ng();
    }

    pub fn get_ns(&self) -> Normal3f {
        let interaction = self.interaction.value.read().unwrap();
        return interaction.get_ns();
    }

    pub fn get_light(&self) -> Option<Arc<dyn Light>> {
        let interaction = self.interaction.value.read().unwrap();
        return interaction.get_light();
    }

    pub fn get_camera(&self) -> Option<Arc<dyn Camera>> {
        let interaction = self.interaction.value.read().unwrap();
        return interaction.get_camera();
    }
}
