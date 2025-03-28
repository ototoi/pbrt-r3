use super::halton::create_halton_sampler;
use super::maxmin::create_maxmindist_sampler;
use super::random::create_random_sampler;
use super::sobol::create_sobol_sampler;
use super::stratified::create_stratified_sampler;
use super::zerotwosequence::create_zerotwosequence_sampler;
use crate::core::prelude::*;

use std::sync::Arc;
use std::sync::RwLock;

pub fn create_sampler(
    name: &str,
    params: &ParamSet,
    film: &Arc<RwLock<Film>>,
) -> Result<Arc<RwLock<dyn Sampler>>, PbrtError> {
    let pixel_bounds = film.read().unwrap().get_sample_bounds();
    match name {
        "lowdiscrepancy" => create_zerotwosequence_sampler(params),
        "02sequence" => create_zerotwosequence_sampler(params),
        "maxmindist" => create_maxmindist_sampler(params),
        "halton" => create_halton_sampler(params, &pixel_bounds),
        "sobol" => create_sobol_sampler(params, &pixel_bounds),
        "random" => create_random_sampler(params),
        "stratified" => create_stratified_sampler(params),
        _ => {
            return Err(PbrtError::warning(&format!(
                "Sampler \"{}\" unknown.",
                name
            )));
        }
    }
}
