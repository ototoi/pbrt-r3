pub mod blackbody;
mod config;
mod constants;
mod convert;
mod data;
mod load;
pub mod rgb;
pub mod sampled;
pub mod utils;

pub use convert::*;
//pub use data::*;
pub use rgb::*;
pub use sampled::*;

#[cfg(not(feature = "sampled-spectrum"))]
pub type Spectrum = RGBSpectrum;

#[cfg(feature = "sampled-spectrum")]
pub type Spectrum = SampledSpectrum;
