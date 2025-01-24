use crate::core::pbrt::*;

use std::sync::Arc;

pub struct KDTreeAccel {
    pub prims: Vec<Arc<dyn Primitive>>,
    pub bounds: Bounds3f,
}

impl KDTreeAccel {
    pub fn new(prims: &[Arc<dyn Primitive>], _: &ParamSet) -> Self {
        let bounds = Bounds3f::new(&Vector3f::new(0.0, 0.0, 0.0), &Vector3f::new(1.0, 1.0, 1.0));
        let prims = prims.to_vec();
        KDTreeAccel { prims, bounds }
    }
}

impl Primitive for KDTreeAccel {
    fn world_bound(&self) -> Bounds3f {
        return self.bounds;
    }
    fn intersect(&self, _r: &Ray) -> Option<SurfaceInteraction> {
        None
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        for it in self.prims.iter() {
            let prim = it.as_ref();
            if prim.intersect_p(r) {
                return true;
            }
        }
        return false;
    }
}

impl Aggregate for KDTreeAccel {}

pub fn create_kdtree_accelerator(
    prims: &[Arc<dyn Primitive>],
    params: &ParamSet,
) -> Result<Arc<dyn Primitive>, PbrtError> {
    return Ok(Arc::new(KDTreeAccel::new(prims, params)));
}
