use super::vertex_interaction::*;
use crate::core::prelude::*;

use std::collections::HashMap;
use std::sync::Arc;

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

pub fn correct_shading_normal(
    isect: &SurfaceInteraction,
    wo: &Vector3f,
    wi: &Vector3f,
    mode: TransportMode,
) -> Float {
    if mode == TransportMode::Importance {
        let num = Vector3f::abs_dot(wo, &isect.shading.n) * Vector3f::abs_dot(wi, &isect.n);
        let denom =
            Vector3f::abs_dot(wo, &isect.n) * Vector3f::abs_dot(wi, &isect.shading.n);
        if denom == 0.0 {
            return 0.0;
        }
        num / denom
    } else {
        1.0
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
        let key = LightKeyType::from(light as *const dyn Light);
        let index = *light_to_index.get(&key).unwrap();
        assert!(index < light_distr.func.len());
        let pdf_i = light.pdf_li(&Interaction::default(), &-*w) * light_distr.func[index];
        assert!(pdf_i >= 0.0);
        pdf += pdf_i;
    }
    pdf / (light_distr.func_int * light_distr.count() as Float)
}

#[derive(Clone)]
pub struct Vertex {
    pub beta: Spectrum,
    pub interaction: VertexInteraction,
    pub delta: bool,
    pub pdf_fwd: Float,
    pub pdf_rev: Float,
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
            beta,
            interaction,
            delta,
            pdf_fwd,
            pdf_rev,
        }
    }

    pub fn get_type(&self) -> VertexType {
        self.interaction.get_type()
    }

    pub fn convert_density(&self, pdf: Float, next: &Vertex) -> Float {
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
        pdf * inv_dist2
    }

    pub fn pdf(&self, scene: &Scene, prev: Option<&Vertex>, next: &Vertex) -> Float {
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

        match &self.interaction {
            VertexInteraction::EndPoint(ei) => {
                if let EndpointInteraction::Camera(ei) = ei {
                    let ray = ei.0.spawn_ray(&wn);
                    let camera = ei.1.as_ref().unwrap();
                    if let Some((_pdf_pos, pdf_dir)) = camera.pdf_we(&ray) {
                        return self.convert_density(pdf_dir, next);
                    }
                }
                0.0
            }
            VertexInteraction::Surface(si) => {
                let prev = prev.unwrap();
                let wp = prev.get_p() - self.get_p();
                if wp.length_squared() == 0.0 {
                    return 0.0;
                }
                let wp = wp.normalize();

                let bsdf = si.bsdf.as_ref().unwrap();
                let pdf = bsdf.pdf(&wp, &wn, BSDF_ALL);
                self.convert_density(pdf, next)
            }
            VertexInteraction::Medium(mi) => {
                let prev = prev.unwrap();
                let wp = prev.get_p() - self.get_p();
                if wp.length_squared() == 0.0 {
                    return 0.0;
                }
                let wp = wp.normalize();

                let pdf = mi.phase.p(&wp, &wn);
                self.convert_density(pdf, next)
            }
        }
    }

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
        } else if let Some(light) = self.get_light() {
            let ray = Ray::new(&self.get_p(), &w, Float::INFINITY, self.get_time());
            if let Some((_pdf_pos, pdf_dir)) = light.pdf_le(&ray, &self.get_ng()) {
                assert!(pdf_dir >= 0.0);
                pdf = pdf_dir * inv_dist2;
            }
        }
        if v.is_on_surface() {
            pdf *= Vector3f::abs_dot(&v.get_ng(), &w);
        }
        pdf
    }

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
            let pdf = infinite_light_density(scene, light_distr, light_to_index, &w);
            assert!(pdf >= 0.0);
            pdf
        } else if let Some(light) = self.get_light() {
            let light_ref = light.as_ref();
            let key = LightKeyType::from(light_ref as *const dyn Light);
            let index = *light_to_index.get(&key).unwrap();
            assert!(index < light_distr.func.len());
            let pdf_choice = light_distr.discrete_pdf(index);
            let ray = Ray::new(&self.get_p(), &w, Float::INFINITY, self.get_time());
            if let Some((pdf_pos, _pdf_dir)) = light_ref.pdf_le(&ray, &self.get_ng()) {
                pdf_choice * pdf_pos
            } else {
                0.0
            }
        } else {
            0.0
        }
    }

    pub fn is_light(&self) -> bool {
        match &self.interaction {
            VertexInteraction::EndPoint(ei) => ei.is_light(),
            VertexInteraction::Surface(si) => {
                if let Some(primitive) = &si.primitive {
                    let primitive = primitive.upgrade().unwrap();
                    primitive.get_area_light().is_some()
                } else {
                    false
                }
            }
            _ => false,
        }
    }

    pub fn is_delta_light(&self) -> bool {
        match &self.interaction {
            VertexInteraction::EndPoint(ei) => ei.is_delta_light(),
            _ => false,
        }
    }

    pub fn is_infinite_light(&self) -> bool {
        match &self.interaction {
            VertexInteraction::EndPoint(ei) => ei.is_infinite_light(),
            _ => false,
        }
    }

    pub fn is_connectible(&self) -> bool {
        match self.get_type() {
            VertexType::Medium => true,
            VertexType::Light => {
                if let Some(light) = self.interaction.get_light() {
                    !light.is_delta_direction()
                } else {
                    false
                }
            }
            VertexType::Camera => true,
            VertexType::Surface => {
                let si = self.interaction.as_surface().unwrap();
                let bsdf = si.bsdf.as_ref().unwrap();
                bsdf.num_components(BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)
                    > 0
            }
        }
    }

    pub fn is_on_surface(&self) -> bool {
        self.get_ng().length_squared() != 0.0
    }

    pub fn f(&self, next: &Vertex, mode: TransportMode) -> Spectrum {
        let wi = next.get_p() - self.get_p();
        if wi.length_squared() == 0.0 {
            return Spectrum::zero();
        }
        let wi = wi.normalize();
        match &self.interaction {
            VertexInteraction::Surface(si) => {
                let bsdf = si.bsdf.as_ref().unwrap();
                let f = bsdf.f(&si.wo, &wi, BSDF_ALL);
                let ns = correct_shading_normal(si, &si.wo, &wi, mode);
                f * ns
            }
            VertexInteraction::Medium(mi) => Spectrum::from(mi.phase.p(&mi.wo, &wi)),
            _ => Spectrum::zero(),
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

        match &self.interaction {
            VertexInteraction::EndPoint(ei) => {
                let mut le = Spectrum::zero();
                if ei.is_infinite_light() {
                    let ray = RayDifferential::from(Ray::new(&self.get_p(), &-w, Float::INFINITY, 0.0));
                    for light in scene.infinite_lights.iter() {
                        le += light.le(&ray);
                    }
                }
                le
            }
            VertexInteraction::Surface(si) => {
                let primitive = si.primitive.as_ref().unwrap().upgrade().unwrap();
                let area_light = primitive.get_area_light().unwrap();
                let area_light = area_light.as_area_light().unwrap();
                area_light.l(&Interaction::from(si), &w)
            }
            _ => Spectrum::zero(),
        }
    }

    pub fn create_camera_from_ray(camera: &Arc<dyn Camera>, ray: &Ray, beta: &Spectrum) -> Self {
        let v = Vertex::new(
            *beta,
            VertexInteraction::from_camera_ray(camera, ray),
            false,
            0.0,
            0.0,
        );
        assert!(v.get_type() == VertexType::Camera);
        v
    }

    pub fn create_camera_from_interaction(
        camera: &Arc<dyn Camera>,
        it: &Interaction,
        beta: &Spectrum,
    ) -> Self {
        let v = Vertex::new(
            *beta,
            VertexInteraction::from_camera_interaction(camera, it),
            false,
            0.0,
            0.0,
        );
        assert!(v.get_type() == VertexType::Camera);
        v
    }

    pub fn create_light_from_ray(
        light: &Arc<dyn Light>,
        ray: &Ray,
        n_light: &Normal3f,
        le: &Spectrum,
        pdf: Float,
    ) -> Self {
        let v = Vertex::new(
            *le,
            VertexInteraction::from_light_ray(light, ray, n_light),
            false,
            pdf,
            0.0,
        );
        assert!(v.get_type() == VertexType::Light);
        v
    }

    pub fn create_light_from_endpoint(
        ei: &EndpointInteraction,
        beta: &Spectrum,
        pdf: Float,
    ) -> Self {
        let v = Vertex::new(*beta, VertexInteraction::EndPoint(ei.clone()), false, pdf, 0.0);
        assert!(v.get_type() == VertexType::Light);
        v
    }

    pub fn create_surface(
        si: &SurfaceInteraction,
        beta: &Spectrum,
        pdf: Float,
        prev: &Vertex,
    ) -> Self {
        let mut v = Vertex::new(*beta, VertexInteraction::Surface(si.clone()), false, 0.0, 0.0);
        v.pdf_fwd = prev.convert_density(pdf, &v);
        assert!(v.get_type() == VertexType::Surface);
        v
    }

    pub fn create_medium(
        mi: &MediumInteraction,
        beta: &Spectrum,
        pdf: Float,
        prev: &Vertex,
    ) -> Self {
        let mut v = Vertex::new(*beta, VertexInteraction::Medium(mi.clone()), false, 0.0, 0.0);
        v.pdf_fwd = prev.convert_density(pdf, &v);
        assert!(v.get_type() == VertexType::Medium);
        v
    }

    pub fn get_interaction(&self) -> Interaction {
        self.interaction.get_interaction()
    }

    pub fn get_p(&self) -> Point3f {
        self.interaction.get_p()
    }

    pub fn get_time(&self) -> Float {
        self.interaction.get_time()
    }

    pub fn get_ng(&self) -> Normal3f {
        self.interaction.get_ng()
    }

    pub fn get_ns(&self) -> Normal3f {
        self.interaction.get_ns()
    }

    pub fn get_light(&self) -> Option<Arc<dyn Light>> {
        self.interaction.get_light()
    }

    pub fn get_camera(&self) -> Option<Arc<dyn Camera>> {
        self.interaction.get_camera()
    }
}
