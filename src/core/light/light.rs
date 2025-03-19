use super::area_light::AreaLight;
use super::visibility_tester::VisibilityTester;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::pbrt::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use crate::core::transform::*;

pub enum LightFlags {
    DeltaPosition = 1,
    DeltaDirection = 2,
    Area = 4,
    Infinite = 8,
}
/*
pub type LightFlags = u32;
pub const LIGHT_DELTA_POSITION:LightFlags = 1;
pub const LIGHT_DELTA_DIRECTION:LightFlags = 2;
pub const LIGHT_AREA:LightFlags = 4;
pub const LIGHT_INFINITE:LightFlags = 8;
*/

pub trait Light: Sync + Send {
    fn sample_li(
        &self,
        inter: &Interaction,
        u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, VisibilityTester)>;
    fn power(&self) -> Spectrum {
        return Spectrum::zero();
    }
    fn preprocess(&self, _scene: &Scene) {}
    fn le(&self, _r: &RayDifferential) -> Spectrum {
        return Spectrum::zero();
    }
    fn sample_le(
        &self,
        _u1: &Point2f,
        _u2: &Point2f,
        _time: Float,
    ) -> Option<(Spectrum, Ray, Normal3f, Float, Float)> {
        return None;
    }
    fn pdf_li(&self, _inter: &Interaction, _wi: &Vector3f) -> Float {
        return 0.0;
    }
    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f) -> Option<(Float, Float)> {
        return None;
    }

    fn is_infinite(&self) -> bool {
        let flags = self.get_light_flags();
        return (flags & LightFlags::Infinite as u32) != 0;
    }
    //---------------------------------
    fn is_delta(&self) -> bool {
        let flags = self.get_light_flags();
        return ((flags & LightFlags::DeltaPosition as u32) != 0)
            || ((flags & LightFlags::DeltaDirection as u32) != 0);
    }

    fn is_delta_direction(&self) -> bool {
        let flags = self.get_light_flags();
        return (flags & LightFlags::DeltaDirection as u32) != 0;
    }

    fn is_area(&self) -> bool {
        let flags = self.get_light_flags();
        return (flags & LightFlags::Area as u32) != 0;
    }

    fn get_light_flags(&self) -> u32 {
        return 0;
    }

    fn get_sample_count(&self) -> u32 {
        return 1;
    }

    fn as_area_light(&self) -> Option<&dyn AreaLight> {
        return None;
    }
}
