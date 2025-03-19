use super::primitive::Primitive;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::material::*;
use crate::core::memory::*;
use crate::core::pbrt::*;
use crate::core::transform::*;

use std::sync::Arc;

pub struct TransformedPrimitive {
    primitive: Arc<dyn Primitive>,
    primitive_to_world: AnimatedTransform,
}

impl TransformedPrimitive {
    pub fn new(primitive: &Arc<dyn Primitive>, primitive_to_world: &AnimatedTransform) -> Self {
        TransformedPrimitive {
            primitive: Arc::clone(primitive),
            primitive_to_world: primitive_to_world.clone(),
        }
    }
}

impl Primitive for TransformedPrimitive {
    fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction> {
        let m = self.primitive_to_world.interpolate(r.time);
        let (ray, _, _) = m.inverse().transform_ray(r);
        let primitive = self.primitive.as_ref();
        if let Some(si) = primitive.intersect(&ray) {
            r.t_max.set(ray.t_max.get());
            let mut si = m.transform_surface_interaction(&si);
            if self.primitive.is_geometric() {
                si.primitive = Some(Arc::downgrade(&self.primitive));
            }
            return Some(si);
        }
        return None;
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        let m = self.primitive_to_world.interpolate(r.time);
        let (ray, _, _) = m.inverse().transform_ray(r);
        let primitive = self.primitive.as_ref();
        return primitive.intersect_p(&ray);
    }
    fn get_area_light(&self) -> Option<Arc<dyn Light>> {
        return None;
    }
    fn get_material(&self) -> Option<Arc<dyn Material>> {
        return None;
    }
    fn compute_scattering_functions(
        &self,
        _si: &mut SurfaceInteraction,
        _arena: &mut MemoryArena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        unimplemented!("TransformedPrimitive::compute_scattering_functions() shouldn't be called");
        //let primitive = self.primitive.as_ref().read().unwrap();
        //primitive.compute_scattering_functions(si, arena, mode, _allow_multiple_lobes);
    }
    fn world_bound(&self) -> Bounds3f {
        let b = self.primitive.as_ref().world_bound();
        return self.primitive_to_world.motion_bounds(&b);
    }
}
