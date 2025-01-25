use crate::core::pbrt::*;

use log::*;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, Clone, Copy)]
pub struct CameraSample {
    pub p_film: Point2f,
    pub p_lens: Point2f,
    pub time: Float,
}

pub trait Camera: Sync + Send {
    fn generate_ray(&self, sample: &CameraSample) -> Option<(Float, Ray)>;
    fn generate_ray_differential(&self, sample: &CameraSample) -> Option<(Float, RayDifferential)>;
    fn we(&self, _ray: &Ray) -> (Spectrum, Point2f) {
        error!("we() method is not implemented for this camera!"); //maybe panic?
        return (Spectrum::default(), Point2f::default());
    }
    fn pdf_we(&self, _ray: &Ray) -> (Float, Float) {
        error!("pdf_we() method is not implemented for this camera!");
        return (0.0, 0.0);
    }
    fn sample_wi(
        &self,
        _inter: &Interaction,
        _u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, Point2f, VisibilityTester)> {
        error!("sample_wi() method is not implemented for this camera!");
        None
    }
    fn get_film(&self) -> Arc<RwLock<Film>>;
    fn get_medium(&self) -> Option<Arc<dyn Medium>>;

    fn get_shutter(&self) -> (Float, Float);
}
