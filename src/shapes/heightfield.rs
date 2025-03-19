use crate::core::camera::*;
use crate::core::distribution::*;
use crate::core::efloat::EFloat;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::lightdistrib::*;
use crate::core::lowdiscrepancy::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::rng::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use crate::core::texture::*;
use crate::core::transform::*;

use std::sync::Arc;

pub fn create_heightfield(
    _o2w: &Transform,
    _w2o: &Transform,
    _reverse_orientation: bool,
    _params: &ParamSet,
) -> Result<Vec<Arc<dyn Shape>>, PbrtError> {
    return Err(PbrtError::error("Not implemented yet."));
}
