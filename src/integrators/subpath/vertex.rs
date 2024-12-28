use super::vertex_interaction::*;
use crate::core::pbrt::*;
use std::collections::HashMap;
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
        if let Some(index) = light_to_index.get(&key) {
            let index = *index;
            if index < light_distr.func.len() {
                let pdf_i = light.pdf_li(&Interaction::default(), &-*w) * light_distr.func[index];
                assert!(pdf_i >= 0.0);
                pdf += pdf_i;
            }
        }
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

#[derive(Clone, Debug)]
pub struct Vertex {
    pub beta: VertexValue<Spectrum>,
    pub interaction: VertexValue<VertexInteraction>,
    pub delta: VertexValue<bool>,
    pub pdf_fwd: VertexValue<Float>,
    pub pdf_rev: VertexValue<Float>,
}

impl Vertex {
    pub fn new(beta: Spectrum, interaction: VertexInteraction, delta: bool, pdf_fwd: Float, pdf_rev: Float) -> Self {
        Vertex {
            beta: VertexValue::new(beta),
            interaction: VertexValue::new(interaction),
            delta: VertexValue::new(delta),
            pdf_fwd: VertexValue::new(pdf_fwd),
            pdf_rev: VertexValue::new(pdf_rev),
        }
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

    pub fn pdf(&self, scene: &Scene, prev: &Option<Arc<RwLock<Vertex>>>, next: &Vertex) -> Float {
        let t = self.get_type();
        if t == VertexType::Light {
            let pdf = self.pdf_light(scene, next);
            assert!(pdf >= 0.0);
            return pdf;
        }
        let wn = next.get_p() - self.get_p();
        if wn.length_squared() == 0.0 {
            return 0.0;
        }
        let wn = wn.normalize();

        assert!(t != VertexType::Light);
        // Compute directional density depending on the vertex types
        if t == VertexType::Camera {
            let interaction = self.interaction.value.read().unwrap();
            let ei = interaction.as_camera().unwrap();
            let ray = ei.0.spawn_ray(&wn);
            let camera = ei.1.as_ref().unwrap().upgrade().unwrap();
            if let Some((_pdf_pos, pdf_dir)) = camera.pdf_we(&ray) {
                assert!(pdf_dir >= 0.0);
                let pdf = self.convert_density(pdf_dir, next);
                assert!(pdf >= 0.0);
                return pdf;
            } else {
                return 0.0;
            }
        } else {
            assert!(prev.is_some());
            let prev = prev.as_ref().unwrap();
            let prev = prev.read().unwrap();
            let wp = prev.get_p() - self.get_p();
            if wp.length_squared() == 0.0 {
                return 0.0;
            }
            let wp = wp.normalize();
            if t == VertexType::Surface {
                let interaction = self.interaction.value.read().unwrap();
                let si = interaction.as_surface().unwrap();
                let bsdf = si.bsdf.as_ref().unwrap().clone();
                let pdf = bsdf.pdf(&wp, &wn, BSDF_ALL);
                assert!(pdf >= 0.0);
                let pdf = self.convert_density(pdf, next);
                assert!(pdf >= 0.0);
                return pdf;
            } else {
                assert!(t == VertexType::Medium);
                let interaction = self.interaction.value.read().unwrap();
                let mi = interaction.as_medium().unwrap();
                let phase = mi.phase.as_ref().unwrap();
                let pdf = phase.p(&wp, &wn);
                assert!(pdf >= 0.0);
                let pdf = self.convert_density(pdf, next);
                assert!(pdf >= 0.0);
                return pdf;
            }
        }
    }

    pub fn pdf_light(&self, scene: &Scene, v: &Vertex) -> Float {
        let w = v.get_p() - self.get_p();
        let inv_dist2 = 1.0 / w.length_squared();
        if w.length_squared() == 0.0 {
            return 0.0;
        }
        let w = w * Float::sqrt(inv_dist2);
        let mut pdf = 0.0;
        if self.is_infinite_light() {
            let (_world_center, world_radius) = scene.world_bound().bounding_sphere();
            pdf = 1.0 / (PI * world_radius * world_radius);
        } else {
            if let Some(light) = self.get_light() {
                let light = light.as_ref();
                let ray = Ray::new(&self.get_p(), &w, Float::INFINITY, 0.0);
                if let Some((_pdf_pos, pdf_dir)) = light.pdf_le(&ray, &self.get_ng()) {
                    assert!(pdf_dir >= 0.0);
                    pdf = pdf_dir * inv_dist2;
                } else {
                    pdf = 0.0;
                }
            }
        }
        if v.is_on_surface() {
            pdf *= Vector3f::abs_dot(&v.get_ng(), &w);
        }
        return pdf;
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
                if let Some(index) = light_to_index.get(&LightKeyType::from(key)) {
                    let index = *index;
                    let pdf_choice = light_distr.discrete_pdf(index);
                    let ray = Ray::new(&self.get_p(), &w, Float::INFINITY, 0.0);
                    if let Some((pdf_pos, _pdf_dir)) = light.pdf_le(&ray, &self.get_ng()) {
                        return pdf_choice * pdf_pos;
                    }
                }
            }
            return 0.0;
        }
    }

    pub fn is_light(&self) -> bool {
        let t = self.get_type();
        match t {
            VertexType::Light => {
                return true;
            }
            VertexType::Surface => {
                let interaction = self.interaction.value.read().unwrap();
                let si = interaction.as_surface().unwrap();
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

    pub fn is_delta_light(&self) -> bool {
        let t = self.get_type();
        if t == VertexType::Light {
            let interaction = self.interaction.value.read().unwrap();
            if let Some(light) = interaction.get_light() {
                return light.is_delta();
            }
        }
        return false;
    }

    pub fn is_infinite_light(&self) -> bool {
        let t = self.get_type();
        if t == VertexType::Light {
            let interaction = self.interaction.value.read().unwrap();
            if let Some(light) = interaction.get_light() {
                let flags = light.get_light_flags();
                return ((flags & LightFlags::Infinite as u32) != 0)
                    || ((flags & LightFlags::DeltaDirection as u32) != 0);
            } else {
                return true;
            }
        }
        return false;
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
                assert!(false);
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
        let t = self.get_type();
        match t {
            VertexType::Surface => {
                let interaction = self.interaction.value.read().unwrap();
                let si = interaction.as_surface().unwrap();
                let bsdf = si.bsdf.as_ref().unwrap();
                let f = bsdf.f(&si.wo, &wi, BSDF_ALL);
                let ns = correct_shading_normal(si, &si.wo, &wi, mode);
                return f * ns;
            }
            VertexType::Medium => {
                let interaction = self.interaction.value.read().unwrap();
                let mi = interaction.as_medium().unwrap();
                let phase = mi.phase.as_ref().unwrap();
                let p = phase.p(&mi.wo, &wi);
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
        if self.is_infinite_light() {
            // Return emitted radiance for infinite light sources
            let mut le = Spectrum::zero();
            for light in scene.infinite_lights.iter() {
                let ray = RayDifferential::from(Ray::new(&self.get_p(), &-w, Float::INFINITY, 0.0));
                le += light.le(&ray);
            }
            return le;
        } else {
            let interaction = self.interaction.value.read().unwrap();
            let si = interaction.as_surface().unwrap();
            let premitive = si.primitive.as_ref().unwrap().upgrade().unwrap();
            let area_light = premitive.get_area_light().unwrap();
            let area_light = area_light.as_area_light().unwrap();
            return area_light.l(&Interaction::from(si), &w);
        }
    }

    pub fn create_camera_from_ray(camera: &Arc<dyn Camera>, ray: &Ray, beta: &Spectrum) -> Self {
        let v = Vertex::new(beta.clone(), VertexInteraction::from_camera_ray(camera, ray), false, 0.0, 0.0);
        assert!(v.get_type() == VertexType::Camera);
        return v;
    }

    pub fn create_camera_from_interaction(
        camera: &Arc<dyn Camera>,
        it: &Interaction,
        beta: &Spectrum,
    ) -> Self {
        let v = Vertex::new(beta.clone(), VertexInteraction::from_camera_interaction(camera, it), false, 0.0, 0.0);
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
        let v = Vertex::new(le.clone(), VertexInteraction::from_light_ray(light, ray, n_light), false, pdf, 0.0);
        assert!(v.get_type() == VertexType::Light);
        return v;
    }

    pub fn create_light_from_endpoint(
        ei: &EndpointInteraction,
        beta: &Spectrum,
        pdf: Float,
    ) -> Self {
        let v = Vertex::new(beta.clone(), VertexInteraction::EndPoint(ei.clone()), false, pdf, 0.0);
        assert!(v.get_type() == VertexType::Light);
        return v;
    }

    pub fn create_surface(
        si: &SurfaceInteraction,
        beta: &Spectrum,
        pdf: Float,
        prev: &Vertex,
    ) -> Self {
        let v = Vertex::new(beta.clone(), VertexInteraction::Surface(si.clone()), false, 0.0, 0.0);
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
        let v = Vertex::new(beta.clone(), VertexInteraction::Medium(mi.clone()), false, 0.0, 0.0);
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
}
