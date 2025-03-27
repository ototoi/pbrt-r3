use crate::core::base::*;
use crate::core::camera::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::medium::*;
use crate::core::param_set::*;
use crate::core::profile::*;
use crate::core::transform::*;

use log::*;
use std::sync::Arc;
use std::sync::RwLock;

pub struct EnvironmentCamera {
    base: BaseCamera,
}

impl EnvironmentCamera {
    pub fn new(
        camera_to_world: &AnimatedTransform,
        shutter_open: Float,
        shutter_close: Float,
        film: &Arc<RwLock<Film>>,
        medium: &Option<Arc<dyn Medium>>,
    ) -> Self {
        let base = BaseCamera::new(camera_to_world, shutter_open, shutter_close, film, medium);
        EnvironmentCamera { base }
    }
}

impl Camera for EnvironmentCamera {
    fn generate_ray(&self, sample: &CameraSample) -> Option<(Float, Ray)> {
        let _p = ProfilePhase::new(Prof::GenerateCameraRay);

        // Compute environment camera ray direction
        let (theta, phi) = {
            let film = self.base.get_film();
            let film = film.read().unwrap();
            let full_resolution = film.full_resolution;
            let theta = PI * sample.p_film.y / full_resolution.y as Float;
            let phi = 2.0 * PI * sample.p_film.x / full_resolution.x as Float;
            (theta, phi)
        };
        let dir = Vector3f::new(
            theta.sin() * phi.cos(),
            theta.cos(),
            theta.sin() * phi.sin(),
        );
        let mut ray = Ray::new(
            &Point3f::zero(),
            &dir,
            Float::INFINITY,
            lerp(sample.time, self.base.shutter_open, self.base.shutter_close),
        );
        ray.medium = self.base.get_medium();
        let (ray, _, _) = self.base.camera_to_world.transform_ray(&ray);
        return Some((1.0, ray));
    }

    fn generate_ray_differential(&self, sample: &CameraSample) -> Option<(Float, RayDifferential)> {
        return BaseCamera::generate_ray_differential(self, sample);
    }

    /*
    fn we(&self, _ray: &Ray) -> (Spectrum, Point2f) {
    }

    fn pdf_we(&self, _ray: &Ray) -> (Float, Float) {
    }

    fn sample_wi(
        &self,
        _inter: &Interaction,
        _u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, Point2f, VisibilityTester)> {
    }
    */

    fn get_film(&self) -> Arc<RwLock<Film>> {
        self.base.get_film()
    }

    fn get_medium(&self) -> Option<Arc<dyn Medium>> {
        self.base.get_medium()
    }

    fn get_shutter(&self) -> (Float, Float) {
        self.base.get_shutter()
    }
}

pub fn create_environment_camera(
    params: &ParamSet,
    cam2world: &AnimatedTransform,
    film: &Arc<RwLock<Film>>,
    medium: &Option<Arc<dyn Medium>>,
) -> Result<Arc<dyn Camera>, PbrtError> {
    let mut shutteropen = params.find_one_float("shutteropen", 0.0);
    let mut shutterclose = params.find_one_float("shutterclose", 1.0);
    if shutterclose < shutteropen {
        warn!(
            "Shutter close time [{}] < shutter open [{}].  Swapping them.",
            shutterclose, shutteropen
        );
        std::mem::swap(&mut shutteropen, &mut shutterclose);
    }
    let _lensradius = params.find_one_float("lensradius", 0.0);
    let _focaldistance = params.find_one_float("focaldistance", 1e6);

    let frame = {
        let film = film.read().unwrap();
        film.full_resolution.x as Float / film.full_resolution.y as Float
    };
    let frame = params.find_one_float("frameaspectratio", frame);
    let mut screen = if frame > 1.0 {
        Bounds2f {
            min: Point2f { x: -frame, y: -1.0 },
            max: Point2f { x: frame, y: 1.0 },
        }
    } else {
        Bounds2f {
            min: Point2f {
                x: -1.0,
                y: -1.0 / frame,
            },
            max: Point2f {
                x: 1.0,
                y: 1.0 / frame,
            },
        }
    };
    if let Some(sw) = params.get_floats_ref("screenwindow") {
        if sw.len() == 4 {
            screen.min.x = sw[0];
            screen.max.x = sw[1];
            screen.min.y = sw[2];
            screen.max.y = sw[3];
        } else {
            error!("Screen window should have four values. Using default.");
        }
    }

    return Ok(Arc::new(EnvironmentCamera::new(
        cam2world,
        shutteropen,
        shutterclose,
        film,
        medium,
    )));
}
