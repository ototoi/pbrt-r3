use crate::core::prelude::*;

use std::sync::Arc;

pub struct MitchellFilter {
    base: BaseFilter,
    b: Float,
    c: Float,
}

impl MitchellFilter {
    pub fn new(radius: &Vector2f, b: Float, c: Float) -> Self {
        MitchellFilter {
            base: BaseFilter::new(radius),
            b,
            c,
        }
    }

    fn mitchell_1d(&self, x: Float) -> Float {
        let b = self.b;
        let c = self.c;
        let x = Float::abs(2.0 * x);
        if x > 1.0 {
            return ((-b - 6.0 * c) * x * x * x
                + (6.0 * b + 30.0 * c) * x * x
                + (-12.0 * b - 48.0 * c) * x
                + (8.0 * b + 24.0 * c))
                * (1.0 / 6.0);
        } else {
            return ((12.0 - 9.0 * b - 6.0 * c) * x * x * x
                + (-18.0 + 12.0 * b + 6.0 * c) * x * x
                + (6.0 - 2.0 * b))
                * (1.0 / 6.0);
        }
    }
}

impl Filter for MitchellFilter {
    fn evaluate(&self, p: &Point2f) -> Float {
        return self.mitchell_1d(p.x * self.base.inv_radius.x)
            * self.mitchell_1d(p.y * self.base.inv_radius.y);
    }
    fn get_radius(&self) -> Vector2f {
        self.base.get_radius()
    }
}

pub fn create_mitchell_filter(params: &ParamSet) -> Result<Arc<dyn Filter>, PbrtError> {
    let xw = params.find_one_float("xwidth", 2.0);
    let yw = params.find_one_float("ywidth", 2.0);
    let b = params.find_one_float("B", 1.0 / 3.0);
    let c = params.find_one_float("C", 1.0 / 3.0);
    return Ok(Arc::new(MitchellFilter::new(&Vector2f::new(xw, yw), b, c)));
}
