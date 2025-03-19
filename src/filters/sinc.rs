use crate::core::error::PbrtError;
use crate::core::filter::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;

use std::sync::Arc;

pub struct SincFilter {
    base: BaseFilter,
    tau: Float,
}

impl SincFilter {
    pub fn new(radius: &Vector2f, tau: Float) -> Self {
        SincFilter {
            base: BaseFilter::new(radius),
            tau,
        }
    }
    fn sinc(x: Float) -> Float {
        let x = Float::abs(x);
        if x < 1e-5 {
            return 1.0;
        } else {
            return Float::sin(PI * x) / (PI * x);
        }
    }
    fn windowed_sinc(x: Float, radius: Float, tau: Float) -> Float {
        let x = Float::abs(x);
        if x > radius {
            return 0.0;
        } else {
            let lanczos = Self::sinc(x / tau);
            return Self::sinc(x) * lanczos;
        }
    }
}

impl Filter for SincFilter {
    fn evaluate(&self, p: &Point2f) -> Float {
        return Self::windowed_sinc(p.x, self.base.radius.x, self.tau)
            * Self::windowed_sinc(p.y, self.base.radius.y, self.tau);
    }
    fn get_radius(&self) -> Vector2f {
        self.base.get_radius()
    }
}

pub fn create_sinc_filter(params: &ParamSet) -> Result<Arc<dyn Filter>, PbrtError> {
    let xw = params.find_one_float("xwidth", 4.0);
    let yw = params.find_one_float("ywidth", 4.0);
    let tau = params.find_one_float("tau", 3.0);
    return Ok(Arc::new(SincFilter::new(&Vector2f::new(xw, yw), tau)));
}
