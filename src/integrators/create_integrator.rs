use super::ao::*;
use super::aov::*;
use super::bdpt::*;
use super::directlighting::*;
use super::mlt::*;
use super::path::*;
use super::sppm::*;
use super::volpath::*;
use super::whitted::*;
use crate::core::camera::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::lightdistrib::*;
use crate::core::material::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use crate::core::stats::*;

use std::sync::Arc;
use std::sync::RwLock;

pub fn create_integrator(
    name: &str,
    params: &ParamSet,
    sampler: &Arc<RwLock<dyn Sampler>>,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    match name {
        "whitted" => create_whitted_integrator(params, sampler, camera),
        "directlighting" => create_direct_lighting_integrator(params, sampler, camera),
        "path" => create_path_integrator(params, sampler, camera),
        "volpath" => create_volpath_integrator(params, sampler, camera),
        "bdpt" => create_bdpt_integrator(params, sampler, camera),
        "mlt" => create_mlt_integrator(params, camera),
        "ambientocclusion" => create_ao_integrator(params, sampler, camera),
        "sppm" => create_sppm_integrator(params, camera),
        "aov" => create_aov_integrator(params, sampler, camera),
        _ => {
            let msg = format!("Integrator \"{}\" unknown.", name);
            return Err(PbrtError::warning(&msg));
        }
    }
}
