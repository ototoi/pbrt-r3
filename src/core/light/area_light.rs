use crate::core::interaction::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;

pub trait AreaLight {
    fn l(&self, inter: &Interaction, w: &Vector3f) -> Spectrum;
}
