#[allow(unused_imports)]
use super::rgb::RGBSpectrum;
#[allow(unused_imports)]
use super::sampled::SampledSpectrum;

#[cfg(not(feature = "sampled-spectrum"))]
pub type Spectrum = RGBSpectrum;

#[cfg(feature = "sampled-spectrum")]
pub type Spectrum = SampledSpectrum;
