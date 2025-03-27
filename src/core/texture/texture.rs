use crate::core::base::*;
use crate::core::interaction::*;

use serde::{Deserialize, Serialize};
use std::fmt::Display;

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ImageWrap {
    Repeat,
    Black,
    Clamp,
}

pub trait Texture<T: Copy> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T;
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TexInfo {
    pub filename: String,
    pub trilinear: bool,
    pub max_aniso: Float,
    pub swrap_mode: ImageWrap,
    pub twrap_mode: ImageWrap,
    pub scale: Float,
    pub gamma: bool,
    pub flip_y: bool,
}

impl Display for TexInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = serde_json::to_string(&self).unwrap();
        return write!(f, "{}", s);
    }
}
