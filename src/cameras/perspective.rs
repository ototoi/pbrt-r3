use crate::core::pbrt::*;
use log::*;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Clone)]
pub struct PerspectiveCamera {
    base: ProjectiveCamera,
    dx_camera: Vector3f,
    dy_camera: Vector3f,
    a: Float,
}

impl PerspectiveCamera {
    pub fn new(
        camera_to_world: &AnimatedTransform,
        screen_window: &Bounds2f,
        shutter_open: Float,
        shutter_close: Float,
        lens_radius: Float,
        focal_distance: Float,
        fov: Float,
        film: &Arc<RwLock<Film>>,
        medium: &Option<Arc<dyn Medium>>,
    ) -> Self {
        let camera_to_screen = Transform::perspective(fov, 1e-2, 1000.0);
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
            .transform_point(&Point3f::new(1.0, 0.0, 0.0))
            - base
                .raster_to_camera
                .transform_point(&Point3f::new(0.0, 0.0, 0.0));
        let dy_camera = base
            .raster_to_camera
            .transform_point(&Point3f::new(0.0, 1.0, 0.0))
            - base
                .raster_to_camera
                .transform_point(&Point3f::new(0.0, 0.0, 0.0));
        let full_resolution = film.read().unwrap().full_resolution;
        let mut min = base
            .raster_to_camera
            .transform_point(&Point3f::new(0.0, 0.0, 0.0));
        let mut max = base.raster_to_camera.transform_point(&Point3f::new(
            full_resolution.x as Float,
            full_resolution.y as Float,
            0.0,
        ));
        min *= 1.0 / min.z;
        max *= 1.0 / max.z;
        let a = Float::abs((max.x - min.x) * (max.y - min.y));
        PerspectiveCamera {
            base,
            dx_camera,
            dy_camera,
            a,
        }
    }
}

unsafe impl Sync for PerspectiveCamera {}

impl Camera for PerspectiveCamera {
    fn generate_ray(&self, sample: &CameraSample) -> Option<(Float, Ray)> {
        let _p = ProfilePhase::new(Prof::GenerateCameraRay);

        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);
        let p_camera = self.base.raster_to_camera.transform_point(&p_film);

        let mut ray = Ray::new(
            &Point3f::zero(),
            &p_camera.normalize(),
            Float::INFINITY,
            sample.time,
        );
        if self.base.lens_radius > 0.0 {
            // Sample point on lens
            let p_lens = concentric_sample_disk(&sample.p_lens) * self.base.lens_radius;

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
        ray.medium = self.base.base.medium.clone();
        let (ray, _, _) = self.base.base.camera_to_world.transform_ray(&ray);
        return Some((1.0, ray));
    }

    fn generate_ray_differential(&self, sample: &CameraSample) -> Option<(Float, RayDifferential)> {
        let p_film = Point3f::new(sample.p_film.x, sample.p_film.y, 0.0);
        let p_camera = self.base.raster_to_camera.transform_point(&p_film);
        let dir = Vector3f::from(p_camera).normalize();

        let mut ray = RayDifferential::new(&Point3f::zero(), &dir, Float::INFINITY, sample.time);
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
                let p_focus = Point3f::new(0.0, 0.0, 0.0) + (ft * dx);
                ray.rx_origin = Point3f::new(p_lens.x, p_lens.y, 0.0);
                ray.rx_direction = (p_focus - ray.rx_origin).normalize();
            }
            {
                let dy = (p_camera + self.dy_camera).normalize();
                let ft = self.base.focal_distance / dy.z;
                let p_focus = Point3f::new(0.0, 0.0, 0.0) + (ft * dy);
                ray.ry_origin = Point3f::new(p_lens.x, p_lens.y, 0.0);
                ray.ry_direction = (p_focus - ray.ry_origin).normalize();
            }
        } else {
            ray.rx_origin = ray.ray.o;
            ray.ry_origin = ray.ray.o;
            ray.rx_direction = (p_camera + self.dx_camera).normalize();
            ray.ry_direction = (p_camera + self.dy_camera).normalize();
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

    fn we(&self, ray: &Ray) -> Option<(Spectrum, Point2f)> {
        // Interpolate camera matrix and check if $\w{}$ is forward-facing
        let c2w = self.base.base.camera_to_world.interpolate(ray.time);
        let w2c = c2w.inverse();
        let cos_theta = ray
            .d
            .dot(&c2w.transform_vector(&Vector3f::new(0.0, 0.0, 1.0)));
        if cos_theta <= 0.0 {
            return None;
        }

        // Map ray $(\w{}, t)$ onto the film plane
        let lens_radius = self.base.lens_radius;
        let focal_distance = self.base.focal_distance;
        let distance = if lens_radius > 0.0 {
            focal_distance
        } else {
            1.0
        } / cos_theta;
        let p_focus = ray.position(distance);
        let camera_to_raster = self.base.raster_to_camera.inverse();
        let p_raster = camera_to_raster.transform_point(&w2c.transform_point(&p_focus)); //world -> camera -> raster
        let p_raster = Vector2f::new(p_raster.x, p_raster.y);

        // Return zero importance for out of bounds points
        let film = self.get_film();
        {
            let film = film.read().unwrap();
            let sample_bounds = film.get_sample_bounds();
            if p_raster.x < sample_bounds.min.x as Float
                || p_raster.x >= sample_bounds.max.x as Float
                || p_raster.y < sample_bounds.min.y as Float
                || p_raster.y >= sample_bounds.max.y as Float
            {
                return None;
            }
        }

        // Compute lens area of perspective camera
        let lens_area = if lens_radius != 0.0 {
            PI * lens_radius * lens_radius
        } else {
            1.0
        };

        // Return importance for point on image plane
        let cos2_theta = cos_theta * cos_theta;
        let cos4_theta = cos2_theta * cos2_theta;
        let a = self.a;
        return Some((Spectrum::from(1.0 / (a * lens_area * cos4_theta)), p_raster));
    }

    fn pdf_we(&self, ray: &Ray) -> Option<(Float, Float)> {
        // Interpolate camera matrix and check if $\w{}$ is forward-facing
        let c2w = self.base.base.camera_to_world.interpolate(ray.time);
        let w2c = c2w.inverse();
        let cos_theta = ray
            .d
            .dot(&c2w.transform_vector(&Vector3f::new(0.0, 0.0, 1.0)));
        if cos_theta <= 0.0 {
            return None;
        }

        // Map ray $(\w{}, t)$ onto the film plane
        let lens_radius = self.base.lens_radius;
        let focal_distance = self.base.focal_distance;
        let distance = if lens_radius > 0.0 {
            focal_distance
        } else {
            1.0
        } / cos_theta;
        let p_focus = ray.position(distance);
        let camera_to_raster = self.base.raster_to_camera.inverse();
        let p_raster = camera_to_raster.transform_point(&w2c.transform_point(&p_focus)); //world -> camera -> raster
                                                                                         //let p_raster2 = Vector2f::new(p_raster.x, p_raster.y);

        // Return zero importance for out of bounds points
        let film = self.get_film();
        {
            let film = film.read().unwrap();
            let sample_bounds = film.get_sample_bounds();
            if p_raster.x < sample_bounds.min.x as Float
                || p_raster.x >= sample_bounds.max.x as Float
                || p_raster.y < sample_bounds.min.y as Float
                || p_raster.y >= sample_bounds.max.y as Float
            {
                return None;
            }
        }
        // Compute lens area of perspective camera
        let lens_area = if lens_radius != 0.0 {
            PI * lens_radius * lens_radius
        } else {
            1.0
        };
        let a = self.a;
        let pdf_pos = 1.0 / lens_area;
        let pdf_dir = 1.0 / (a * cos_theta * cos_theta * cos_theta);
        return Some((pdf_pos, pdf_dir));
    }

    fn sample_wi(
        &self,
        inter: &Interaction,
        u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, Point2f, VisibilityTester)> {
        let lens_radius = self.base.lens_radius;
        let time = inter.get_time();
        let c2w = self.base.base.camera_to_world.interpolate(time);

        // Uniformly sample a lens interaction _lensIntr_
        let p_lens = lens_radius * concentric_sample_disk(&u);
        let p_lens_world = c2w.transform_point(&Point3f::new(p_lens.x, p_lens.y, 0.0));
        let medium = self.base.base.medium.clone();
        let mi = MediumInterface::from(&medium);
        let n = Normal3f::from(c2w.transform_vector(&Vector3f::new(0.0, 0.0, 1.0)));
        let lens_intr = Interaction::Base(BaseInteraction {
            p: p_lens_world,
            time: time,
            medium_interface: mi,
            n: n,
            ..Default::default()
        });

        // Populate arguments and compute the importance value
        let vis = VisibilityTester::from((inter.clone(), lens_intr.clone()));
        let wi = p_lens_world - inter.get_p();
        let dist = wi.length();
        let wi = wi / dist;

        // Compute PDF for importance arriving at _ref_

        // Compute lens area of perspective camera
        let lens_area = if lens_radius != 0.0 {
            PI * lens_radius * lens_radius
        } else {
            1.0
        };
        let pdf = (dist * dist) / (n.abs_dot(&wi) * lens_area);
        let ray = lens_intr.spawn_ray(&-wi);
        if let Some((spec, p_raster)) = self.we(&ray) {
            return Some((spec, wi, pdf, p_raster, vis));
        } else {
            return None;
        }
    }

    fn get_film(&self) -> Arc<RwLock<Film>> {
        return self.base.get_film();
    }

    fn get_medium(&self) -> Option<Arc<dyn Medium>> {
        return self.base.get_medium();
    }

    fn get_shutter(&self) -> (Float, Float) {
        return self.base.get_shutter();
    }
}

pub fn create_perspective_camera(
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
    let full_resolution = film.read().unwrap().full_resolution;
    let aspect = (full_resolution.x as Float) / (full_resolution.y as Float);
    let lensradius = params.find_one_float("lensradius", 0.0);
    let focaldistance = params.find_one_float("focaldistance", 1e6);
    let frame = params.find_one_float("frameaspectratio", aspect);
    let screen = if let Some(swi) = params.get_floats_ref("screenwindow") {
        if swi.len() >= 4 {
            Bounds2f::from(((swi[0], swi[2]), (swi[1], swi[3])))
        } else {
            let msg = format!("\"screenwindow\" should have four values");
            return Err(PbrtError::error(&msg));
        }
    } else if frame > 1.0 {
        Bounds2f::from(((-frame, -1.0), (frame, 1.0)))
    } else {
        Bounds2f::from(((-1.0, -1.0 / frame), (1.0, 1.0 / frame)))
    };
    let mut fov = params.find_one_float("fov", 90.0);
    let halffov = params.find_one_float("halffov", -1.0);
    if halffov > 0.0 {
        // hack for structure synth, which exports half of the full fov
        fov = 2.0 * halffov;
    }
    return Ok(Arc::new(PerspectiveCamera::new(
        cam2world,
        &screen,
        shutteropen,
        shutterclose,
        lensradius,
        focaldistance,
        fov,
        film,
        medium,
    )));
}
