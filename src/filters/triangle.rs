use crate::core::error::PbrtError;
use crate::core::filter::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;

use std::sync::Arc;

pub struct TriangleFilter {
    base: BaseFilter,
}

impl TriangleFilter {
    pub fn new(radius: &Vector2f) -> Self {
        TriangleFilter {
            base: BaseFilter::new(radius),
        }
    }
}

impl Filter for TriangleFilter {
    fn evaluate(&self, p: &Point2f) -> Float {
        return Float::max(0.0, self.base.radius.x - Float::abs(p.x))
            * Float::max(0.0, self.base.radius.y - Float::abs(p.y));
    }
    fn get_radius(&self) -> Vector2f {
        self.base.get_radius()
    }
}

pub fn create_triangle_filter(params: &ParamSet) -> Result<Arc<dyn Filter>, PbrtError> {
    let xw = params.find_one_float("xwidth", 2.0);
    let yw = params.find_one_float("ywidth", 2.0);
    return Ok(Arc::new(TriangleFilter::new(&Vector2f::new(xw, yw))));
}
