use crate::core::pbrt::*;

use std::sync::Arc;

pub struct GaussianFilter {
    base: BaseFilter,
    alpha: Float,
    exp_x: Float,
    exp_y: Float,
}

impl GaussianFilter {
    pub fn new(radius: &Vector2f, alpha: Float) -> Self {
        GaussianFilter {
            base: BaseFilter::new(radius),
            alpha,
            exp_x: Float::exp(-alpha * radius.x * radius.x),
            exp_y: Float::exp(-alpha * radius.y * radius.y),
        }
    }

    fn gaussian(&self, d: Float, expv: Float) -> Float {
        return Float::max(0.0, Float::exp(-self.alpha * d * d) - expv);
    }
}

impl Filter for GaussianFilter {
    fn evaluate(&self, p: &Point2f) -> Float {
        return self.gaussian(p.x, self.exp_x) * self.gaussian(p.y, self.exp_y);
    }

    fn get_radius(&self) -> Vector2f {
        self.base.get_radius()
    }
}

pub fn create_gaussian_filter(params: &ParamSet) -> Result<Arc<dyn Filter>, PbrtError> {
    let xw = params.find_one_float("xwidth", 2.0);
    let yw = params.find_one_float("ywidth", 2.0);
    let alpha = params.find_one_float("alpha", 2.0);
    return Ok(Arc::new(GaussianFilter::new(&Vector2f::new(xw, yw), alpha)));
}
