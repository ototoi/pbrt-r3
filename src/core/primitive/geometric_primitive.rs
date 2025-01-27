use super::primitive::Primitive;
use crate::core::pbrt::*;

use std::sync::Arc;

pub struct GeometricPrimitive {
    pub shape: Arc<dyn Shape>,
    pub material: Option<Arc<dyn Material>>,
    pub area_light: Option<Arc<dyn Light>>,
    pub mi: MediumInterface,
}

impl GeometricPrimitive {
    pub fn new(
        shape: &Arc<dyn Shape>,
        material: &Option<Arc<dyn Material>>,
        area_light: &Option<Arc<dyn Light>>,
        mi: &MediumInterface,
    ) -> Self {
        GeometricPrimitive {
            shape: Arc::clone(shape),
            material: material.clone(),
            area_light: area_light.clone(),
            mi: mi.clone(),
        }
    }
}

impl Primitive for GeometricPrimitive {
    fn world_bound(&self) -> Bounds3f {
        let s = self.shape.as_ref();
        return s.world_bound();
    }
    fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction> {
        let _p = ProfilePhase::new(Prof::GeometricPrimitiveIntersect);

        let s = self.shape.as_ref();
        if let Some((t_hit, mut isect)) = s.intersect(r) {
            r.t_max.set(t_hit);
            isect.set_shape(&self.shape);

            if self.mi.is_medium_transition() {
                isect.medium_interface = self.mi.clone();
            } else {
                isect.medium_interface = MediumInterface::from(&r.medium);
            }

            return Some(isect);
        }
        return None;
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        let _p = ProfilePhase::new(Prof::GeometricPrimitiveIntersectP);

        let s = self.shape.as_ref();
        return s.intersect_p(r);
    }

    fn get_area_light(&self) -> Option<Arc<dyn Light>> {
        return self.area_light.clone();
    }

    fn get_material(&self) -> Option<Arc<dyn Material>> {
        return self.material.clone();
    }

    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        let _p = ProfilePhase::new(Prof::ComputeScatteringFuncs);
        if let Some(mat) = self.material.as_ref() {
            mat.compute_scattering_functions(si, arena, mode, allow_multiple_lobes);
        }
    }

    fn is_geometric(&self) -> bool {
        return true;
    }
}
