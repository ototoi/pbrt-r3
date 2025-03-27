use crate::core::base::*;
use crate::core::error::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::param_set::*;
use crate::core::primitive::*;

use std::sync::Arc;

fn get_bounds(prims: &[Arc<dyn Primitive>]) -> Bounds3f {
    let bounds: Vec<Bounds3f> = prims
        .iter()
        .map(|p| {
            return p.as_ref().world_bound();
        })
        .collect();
    let min = bounds
        .iter()
        .map(|b| -> Vector3f {
            return b.min;
        })
        .reduce(|a, b| -> Vector3f {
            return Vector3f::new(
                Float::min(a[0], b[0]),
                Float::min(a[1], b[1]),
                Float::min(a[2], b[2]),
            );
        })
        .unwrap();
    let max = bounds
        .iter()
        .map(|b| -> Vector3f {
            return b.max;
        })
        .reduce(|a, b| -> Vector3f {
            return Vector3f::new(
                Float::max(a[0], b[0]),
                Float::max(a[1], b[1]),
                Float::max(a[2], b[2]),
            );
        })
        .unwrap();
    return Bounds3f::from(((min[0], min[1], min[2]), (max[0], max[1], max[2])));
}

pub struct ExhaustiveAccel {
    pub prims: Vec<Arc<dyn Primitive>>,
    pub bounds: Bounds3f,
}

impl ExhaustiveAccel {
    pub fn new(prims: &[Arc<dyn Primitive>]) -> Self {
        assert!(!prims.is_empty());
        let bounds = get_bounds(prims);
        let prims = prims.to_vec();
        ExhaustiveAccel { prims, bounds }
    }
}

impl Primitive for ExhaustiveAccel {
    fn world_bound(&self) -> Bounds3f {
        return self.bounds;
    }
    fn intersect(&self, r: &Ray) -> Option<SurfaceInteraction> {
        //let (b, _, _) = self.bounds.intersect_p(r);
        //if b {
        if let Some(_) = self.bounds.intersect_p(r) {
            let mut opt_isect = None;
            for it in self.prims.iter() {
                let prim = it.as_ref();
                if let Some(mut isect) = prim.intersect(r) {
                    if prim.is_geometric() {
                        isect.primitive = Some(Arc::downgrade(it));
                    }
                    opt_isect = Some(isect);
                }
            }
            return opt_isect;
        }
        return None;
    }

    fn intersect_p(&self, r: &Ray) -> bool {
        //let (b, _, _) = self.bounds.intersect_p(r);
        //if b {
        if let Some(_) = self.bounds.intersect_p(r) {
            for it in self.prims.iter() {
                let prim = it.as_ref();
                if prim.intersect_p(r) {
                    return true;
                }
            }
        }
        return false;
    }
}

impl Aggregate for ExhaustiveAccel {}

pub fn create_exhaustive_accelerator(
    prims: &[Arc<dyn Primitive>],
    _: &ParamSet,
) -> Result<Arc<dyn Primitive>, PbrtError> {
    return Ok(Arc::new(ExhaustiveAccel::new(prims)));
}
