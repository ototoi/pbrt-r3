use super::sampled::*;

use crate::core::error::*;
use crate::core::misc::read_float_file;
use crate::core::pbrt::*;
use log::*;

//AddSampledSpectrumFiles
pub fn load_from_file(path: &str) -> Result<SampledSpectrum, PbrtError> {
    match read_float_file(path) {
        Ok(vals) => {
            if vals.len() % 2 != 0 {
                warn!(
                    "Extra value found in spectrum file \"{}\". Ignoring it.",
                    path
                );
            }
            let mut wls = Vec::new();
            let mut v = Vec::new();
            for j in 0..(vals.len() / 2) {
                wls.push(vals[2 * j] as Float);
                v.push(vals[2 * j + 1] as Float);
            }
            return Ok(SampledSpectrum::sampled_from_sampled(&wls, &v));
        }
        Err(e) => {
            return Err(e);
        }
    }
}
