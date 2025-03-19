use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::material::*;
use crate::core::memory::*;

use std::sync::Arc;

pub trait Primitive {
    fn world_bound(&self) -> Bounds3f;
    fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction>;
    fn intersect_p(&self, r: &Ray) -> bool;
    fn get_area_light(&self) -> Option<Arc<dyn Light>> {
        None
    }
    fn get_material(&self) -> Option<Arc<dyn Material>> {
        None
    }
    fn compute_scattering_functions(
        &self,
        _si: &mut SurfaceInteraction,
        _arena: &mut MemoryArena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
    }

    fn is_geometric(&self) -> bool {
        return false;
    }
}
