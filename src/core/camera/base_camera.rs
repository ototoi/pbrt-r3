use crate::core::pbrt::*;

use log::*;
use std::sync::Arc;
use std::sync::RwLock;

/*
    AnimatedTransform CameraToWorld;
    const Float shutterOpen, shutterClose;
    Film *film;
    const Medium *medium;
*/

#[derive(Clone)]
pub struct BaseCamera {
    pub camera_to_world: AnimatedTransform,
    pub shutter_open: Float,
    pub shutter_close: Float,
    pub film: Arc<RwLock<Film>>,
    pub medium: Option<Arc<dyn Medium>>,
}

impl BaseCamera {
    pub fn new(
        camera_to_world: &AnimatedTransform,
        shutter_open: Float,
        shutter_close: Float,
        film: &Arc<RwLock<Film>>,
        medium: &Option<Arc<dyn Medium>>,
    ) -> Self {
        if camera_to_world.has_scale() {
            warn!(
                "Scaling detected in world-to-camera transformation!\nThe system has numerous assumptions, implicit and explicit,\nthat this transform will have no scale factors in it.\nProceed at your own risk; your image may have errors or\nthe system may crash as a result of this."
            );
        }

        BaseCamera {
            camera_to_world: camera_to_world.clone(),
            shutter_open,
            shutter_close,
            film: film.clone(),
            medium: medium.clone(),
        }
    }

    pub fn get_film(&self) -> Arc<RwLock<Film>> {
        return Arc::clone(&self.film);
    }

    pub fn get_medium(&self) -> Option<Arc<dyn Medium>> {
        return self.medium.clone();
    }

    pub fn get_shutter(&self) -> (Float, Float) {
        return (self.shutter_open, self.shutter_close);
    }

    //ray diff helper
    pub fn generate_ray_differential(
        camera: &dyn Camera,
        sample: &CameraSample,
    ) -> Option<(Float, RayDifferential)> {
        let (wt, rd) = camera.generate_ray(sample)?;

        // Find camera ray after shifting a fraction of a pixel in the $x$ direction
        let mut rx_origin = Vector3f::zero();
        let mut rx_direction = Vector3f::zero();
        let mut wtx = 0.0;
        for eps in [0.05, -0.05] {
            let mut sshift = *sample;
            sshift.p_film.x += eps;
            if let Some((wwtx, rx)) = camera.generate_ray(&sshift) {
                wtx = wwtx;
                rx_origin = rd.o + (rx.o - rd.o) / eps;
                rx_direction = rd.d + (rx.d - rd.d) / eps;
                break;
            }
        }
        if wtx == 0.0 {
            return None;
        }

        // Find camera ray after shifting a fraction of a pixel in the $x$ direction
        let mut ry_origin = Vector3f::zero();
        let mut ry_direction = Vector3f::zero();
        let mut wty = 0.0;
        for eps in [0.05, -0.05] {
            let mut sshift = *sample;
            sshift.p_film.y += eps;
            if let Some((wwty, ry)) = camera.generate_ray(&sshift) {
                wty = wwty;
                ry_origin = rd.o + (ry.o - rd.o) / eps;
                ry_direction = rd.d + (ry.d - rd.d) / eps;
                break;
            }
        }
        if wty == 0.0 {
            return None;
        }
        let rd = RayDifferential {
            ray: rd,
            has_differentials: true,
            rx_origin,
            rx_direction,
            ry_origin,
            ry_direction,
        };

        return Some((wt, rd));
    }
}
