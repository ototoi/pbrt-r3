use crate::core::pbrt::*;

use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, Default, Copy, Clone)]
struct WorldBound {
    pub center: Point3f,
    pub radius: Float,
}

pub struct DistantLight {
    base: BaseLight,
    l: Spectrum,
    w_light: Vector3f,
    bound: RwLock<WorldBound>,
}

impl DistantLight {
    pub fn new(light_to_world: &Transform, l: &Spectrum, w_light: &Vector3f) -> Self {
        let base = BaseLight::new(
            LightFlags::DeltaDirection as u32,
            light_to_world,
            &MediumInterface::new(),
            1,
        );
        let l = *l;
        let w_light = light_to_world.transform_vector(w_light).normalize();
        DistantLight {
            base,
            l,
            w_light,
            bound: RwLock::new(WorldBound {
                center: Point3f::new(0.0, 0.0, 0.0),
                radius: Float::INFINITY,
            }),
        }
    }

    fn get_bound(&self) -> WorldBound {
        let bound = self.bound.read().unwrap();
        return *bound;
    }
}

impl Light for DistantLight {
    fn preprocess(&self, scene: &Scene) {
        let (center, radius) = scene.world_bound.bounding_sphere();
        //println!("center:{:?}, radius:{:?}", center, radius);
        let mut bound = self.bound.write().unwrap();
        bound.center = center;
        bound.radius = radius;
    }

    fn sample_li(
        &self,
        inter: &Interaction,
        _u: &Point2f,
    ) -> Option<(Spectrum, Vector3f, Float, VisibilityTester)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let bound = self.bound.read().unwrap();
        let world_radius = bound.radius;
        let wi = self.w_light;
        let pdf = 1.0;
        let p_outside = inter.get_p() + self.w_light * (2.0 * world_radius);

        let inter_light = Interaction::from_light_sample(
            &p_outside,
            inter.get_time(),
            &self.base.medium_interface,
        );
        let vis = VisibilityTester::from((inter, &inter_light));
        return Some((self.l, wi, pdf, vis));
    }

    fn power(&self) -> Spectrum {
        let bound = self.bound.read().unwrap();
        let world_radius = bound.radius;
        return self.l * (PI * world_radius * world_radius);
    }

    fn pdf_li(&self, _inter: &Interaction, _wi: &Vector3f) -> Float {
        return 0.0;
    }

    fn sample_le(
        &self,
        u1: &Point2f,
        _u2: &Point2f,
        time: Float,
    ) -> Option<(Spectrum, Ray, Normal3f, Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightSample);

        let w_light = self.w_light;
        let bound = self.get_bound();
        let world_center = bound.center;
        let world_radius = bound.radius;
        let (v1, v2) = coordinate_system(&w_light);
        let cd = concentric_sample_disk(u1);
        let p_disk = world_center + world_radius * (cd.x * v1 + cd.y * v2);

        // Set ray origin and direction for infinite light ray
        let l = self.l;
        let ray = Ray::new(
            &(p_disk + world_radius * w_light),
            &-w_light,
            Float::INFINITY,
            time,
        );
        let n_light = -w_light;
        let pdf_pos = 1.0 / (PI * world_radius * world_radius);
        let pdf_dir = 1.0;
        return Some((l, ray, n_light, pdf_pos, pdf_dir));
    }

    fn pdf_le(&self, _ray: &Ray, _n: &Normal3f) -> Option<(Float, Float)> {
        let _p = ProfilePhase::new(Prof::LightPdf);

        let bound = self.bound.read().unwrap();
        let world_radius = bound.radius;
        let pdf_pos = 1.0 / (PI * world_radius * world_radius);
        let pdf_dir = 0.0;
        return Some((pdf_pos, pdf_dir));
    }

    fn get_light_flags(&self) -> u32 {
        return self.base.flags;
    }

    fn get_sample_count(&self) -> u32 {
        return self.base.n_samples;
    }
}

pub fn create_distant_light(
    light2world: &Transform,
    params: &ParamSet,
) -> Result<Arc<dyn Light>, PbrtError> {
    let l = params.find_one_spectrum("L", &Spectrum::one());
    let scale = params.find_one_spectrum("scale", &Spectrum::one());
    let from = params.find_one_point3f("from", &Point3f::new(0.0, 0.0, 0.0));
    let to = params.find_one_point3f("to", &Point3f::new(0.0, 0.0, 1.0));
    let dir = from - to;
    return Ok(Arc::new(DistantLight::new(light2world, &(l * scale), &dir)));
}
