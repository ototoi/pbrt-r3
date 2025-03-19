use crate::core::camera::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::medium::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::sampling::*;
use crate::core::transform::*;

use log::*;
use std::sync::Arc;
use std::sync::RwLock;

pub struct OrthographicCamera {
    base: ProjectiveCamera,
    dx_camera: Vector3f,
    dy_camera: Vector3f,
}

impl OrthographicCamera {
    pub fn new(
        camera_to_world: &AnimatedTransform,
        screen_window: &Bounds2f,
        shutter_open: Float,
        shutter_close: Float,
        lens_radius: Float,
        focal_distance: Float,
        film: &Arc<RwLock<Film>>,
        medium: &Option<Arc<dyn Medium>>,
    ) -> Self {
        let camera_to_screen = Transform::orthographic(0.0, 1.0);
        let base = ProjectiveCamera::new(
            camera_to_world,
            &camera_to_screen,
            screen_window,
            shutter_open,
            shutter_close,
            lens_radius,
            focal_distance,
            film,
            medium,
        );
        let dx_camera = base
            .raster_to_camera
            .transform_vector(&Vector3f::new(1.0, 0.0, 0.0));
        let dy_camera = base
            .raster_to_camera
            .transform_vector(&Vector3f::new(0.0, 1.0, 0.0));
        OrthographicCamera {
            base,
            dx_camera,
            dy_camera,
        }
    }
}

impl Camera for OrthographicCamera {
    fn generate_ray(&self, sample: &CameraSample) -> Option<(Float, Ray)> {
        let _p = ProfilePhase::new(Prof::GenerateCameraRay);

        // Compute raster and camera sample positions
        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);
        let p_camera = self.base.raster_to_camera.transform_point(&p_film);
        let mut ray = Ray::new(
            &p_camera,
            &Vector3f::new(0.0, 0.0, 1.0),
            Float::INFINITY,
            sample.time,
        );
        // Modify ray for depth of field
        if self.base.lens_radius > 0.0 {
            // Sample point on lens
            let p_lens = self.base.lens_radius * concentric_sample_disk(&sample.p_lens);

            // Compute point on plane of focus
            let ft = self.base.focal_distance / ray.d.z;
            let p_focus = ray.position(ft);

            // Update ray for effect of lens
            ray.o = Point3f::new(p_lens.x, p_lens.y, 0.0);
            ray.d = (p_focus - ray.o).normalize();
        }
        ray.time = lerp(
            sample.time,
            self.base.base.shutter_open,
            self.base.base.shutter_close,
        );
        ray.medium = self.base.base.get_medium();
        let (ray, _, _) = self.base.base.camera_to_world.transform_ray(&ray);
        return Some((1.0, ray));
    }

    fn generate_ray_differential(&self, sample: &CameraSample) -> Option<(Float, RayDifferential)> {
        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);
        let p_camera = self.base.raster_to_camera.transform_point(&p_film);

        let mut ray = RayDifferential::new(
            &p_camera,
            &Vector3f::new(0.0, 0.0, 1.0),
            Float::INFINITY,
            sample.time,
        );
        // Modify ray for depth of field
        if self.base.lens_radius > 0.0 {
            // Sample point on lens
            let p_lens = concentric_sample_disk(&sample.p_lens) * self.base.lens_radius;

            // Compute point on plane of focus
            let ft = self.base.focal_distance / ray.ray.d.z;
            let p_focus = ray.ray.position(ft);

            // Update ray for effect of lens
            ray.ray.o = Point3f::new(p_lens.x, p_lens.y, 0.0);
            ray.ray.d = (p_focus - ray.ray.o).normalize();
        }

        // Compute offset rays for _PerspectiveCamera_ ray differentials
        if self.base.lens_radius > 0.0 {
            // Compute _PerspectiveCamera_ ray differentials accounting for lens

            // Sample point on lens
            let p_lens = concentric_sample_disk(&sample.p_lens) * self.base.lens_radius;
            {
                let dx = (p_camera + self.dx_camera).normalize();
                let ft = self.base.focal_distance / dx.z;
                let p_focus = p_camera + self.dx_camera + (ft * Vector3f::new(0.0, 0.0, 1.0));
                ray.rx_origin = Point3f::new(p_lens.x, p_lens.y, 0.0);
                ray.rx_direction = (p_focus - ray.rx_origin).normalize();
            }
            {
                let dy = (p_camera + self.dy_camera).normalize();
                let ft = self.base.focal_distance / dy.z;
                let p_focus = p_camera + self.dy_camera + (ft * Vector3f::new(0.0, 0.0, 1.0));
                ray.ry_origin = Point3f::new(p_lens.x, p_lens.y, 0.0);
                ray.ry_direction = (p_focus - ray.ry_origin).normalize();
            }
        } else {
            ray.rx_origin = ray.ray.o + self.dx_camera;
            ray.ry_origin = ray.ray.o + self.dy_camera;
            ray.rx_direction = ray.ray.d;
            ray.ry_direction = ray.ray.d;
        }

        ray.ray.time = lerp(
            sample.time,
            self.base.base.shutter_open,
            self.base.base.shutter_close,
        );
        ray.ray.medium = self.base.base.get_medium();

        let (mut ray, _, _) = self
            .base
            .base
            .camera_to_world
            .transform_ray_differential(&ray);
        ray.has_differentials = true;

        return Some((1.0, ray));
    }

    /*
    fn we(&self, _ray: &Ray) -> (Spectrum, Point2f) {
        return (Spectrum::zero(), Point2f::new(0.0, 0.0));
    }
    fn pdf_we(&self, _ray: &Ray) -> (Float, Float) {
        (1.0, 1.0)
    }
    fn sample_wi(
        &self,
        _inter: &Interaction,
        _u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, Point2f, VisibilityTester)> {
        None
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

pub fn create_orthographic_camera(
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
    let lensradius = params.find_one_float("lensradius", 0.0);
    let focaldistance = params.find_one_float("focaldistance", 1e6);

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

    return Ok(Arc::new(OrthographicCamera::new(
        cam2world,
        &screen,
        shutteropen,
        shutterclose,
        lensradius,
        focaldistance,
        film,
        &medium,
    )));
}
