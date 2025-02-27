pub use super::spectrum_config::*;

#[cfg(not(feature = "float-as-double"))]
pub type Float = f32;
#[cfg(feature = "float-as-double")]
pub type Float = f64;