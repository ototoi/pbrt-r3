#![allow(unused)]

use super::data::rgb_refl::*;
use super::data::xyz::*;
use super::sampled::SampledSpectrum;

pub const X: SampledSpectrum = SampledSpectrum { c: ARRAY_CIE_X };
pub const Y: SampledSpectrum = SampledSpectrum { c: ARRAY_CIE_Y };
pub const Z: SampledSpectrum = SampledSpectrum { c: ARRAY_CIE_Z };

pub const RGBREFL2SPECT_WHITE: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBREFL2SPECT_WHITE,
};
pub const RGBREFL2SPECT_CYAN: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBREFL2SPECT_CYAN,
};
pub const RGBREFL2SPECT_MAGENTA: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBREFL2SPECT_MAGENTA,
};
pub const RGBREFL2SPECT_YELLOW: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBREFL2SPECT_YELLOW,
};
pub const RGBREFL2SPECT_RED: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBREFL2SPECT_RED,
};
pub const RGBREFL2SPECT_GREEN: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBREFL2SPECT_GREEN,
};
pub const RGBREFL2SPECT_BLUE: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBREFL2SPECT_BLUE,
};

pub const RGBILLUM2SPECT_WHITE: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBILLUM2SPECT_WHITE,
};
pub const RGBILLUM2SPECT_CYAN: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBILLUM2SPECT_CYAN,
};
pub const RGBILLUM2SPECT_MAGENTA: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBILLUM2SPECT_MAGENTA,
};
pub const RGBILLUM2SPECT_YELLOW: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBILLUM2SPECT_YELLOW,
};
pub const RGBILLUM2SPECT_RED: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBILLUM2SPECT_RED,
};
pub const RGBILLUM2SPECT_GREEN: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBILLUM2SPECT_GREEN,
};
pub const RGBILLUM2SPECT_BLUE: SampledSpectrum = SampledSpectrum {
    c: ARRAY_RGBILLUM2SPECT_BLUE,
};
