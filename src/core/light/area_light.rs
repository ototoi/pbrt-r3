use crate::core::pbrt::*;

pub trait AreaLight {
    fn l(&self, inter: &Interaction, w: &Vector3f) -> Spectrum;
}
