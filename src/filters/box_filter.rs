use crate::core::error::*;
use crate::core::filter::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;

use std::sync::Arc;

pub struct BoxFilter {
    base: BaseFilter,
}

impl BoxFilter {
    pub fn new(radius: &Vector2f) -> Self {
        BoxFilter {
            base: BaseFilter::new(radius),
        }
    }
}

impl Filter for BoxFilter {
    fn evaluate(&self, _p: &Point2f) -> Float {
        return 1.0;
    }
    fn get_radius(&self) -> Vector2f {
        self.base.get_radius()
    }
}

pub fn create_box_filter(params: &ParamSet) -> Result<Arc<dyn Filter>, PbrtError> {
    let xw = params.find_one_float("xwidth", 0.5);
    let yw = params.find_one_float("ywidth", 0.5);
    return Ok(Arc::new(BoxFilter::new(&Vector2f::new(xw, yw))));
}
