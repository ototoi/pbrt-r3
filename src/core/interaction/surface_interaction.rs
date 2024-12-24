use crate::core::pbrt::*;
use std::sync::{Arc, Weak};

#[derive(Default, Debug, Clone, Copy)]
pub struct SurfaceInteractionShading {
    pub n: Normal3f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
}

#[derive(Default, Clone, Debug)]
pub struct SurfaceInteraction {
    pub p: Point3f,
    pub p_error: Vector3f,
    pub n: Normal3f,
    pub time: Float,
    //-------------
    pub wo: Vector3f,
    pub medium_interface: MediumInterface,
    //-------------
    pub uv: Point2f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
    //-------------
    pub shape: Option<Weak<dyn Shape>>,
    pub primitive: Option<Weak<dyn Primitive>>,
    pub shading: SurfaceInteractionShading,

    pub bsdf: Option<Arc<BSDF>>,
    pub bssrdf: Option<Arc<dyn BSSRDF>>,

    //pub material: Option<Weak<RefCell<dyn Material>>>,
    pub dpdx: Vector3f,
    pub dpdy: Vector3f,
    pub dudx: Float,
    pub dvdx: Float,
    pub dudy: Float,
    pub dvdy: Float,
    //-------------
    pub face_index: u32,
    //-------------
}

impl SurfaceInteraction {
    /*
    pub fn from_interaction(inter: &Interaction) -> Self {
        SurfaceInteraction {
            p: inter.p,
            time: inter.time,
            p_error: inter.p_error,
            wo: inter.wo,
            n: inter.n,
            medium_interface: inter.medium_interface.clone(),
            uv: Vector2f::zero(),
            dpdu: Vector3f::zero(),
            dpdv: Vector3f::zero(),
            dndu: Vector3f::zero(),
            dndv: Vector3f::zero(),
            shape: None,
            primitive: None,
            shading: SurfaceInteractionShading::default(),
            bsdf: None,
            bssrdf: None,
            //material: None,
            dpdx: Vector3f::zero(),
            dpdy: Vector3f::zero(),
            dudx: 0.0,
            dudy: 0.0,
            dvdx: 0.0,
            dvdy: 0.0,
            face_index: 0,
        }
    }
    */

    pub fn new(
        p: &Point3f,
        p_error: &Point3f,
        uv: &Point2f,
        wo: &Vector3f,
        n: &Normal3f,
        dpdu: &Vector3f,
        dpdv: &Vector3f,
        dndu: &Vector3f,
        dndv: &Vector3f,
        time: Float,
        face_index: u32,
    ) -> Self {
        //let n = Vector3f::cross(dpdu, dpdv).normalize();
        let shading = SurfaceInteractionShading {
            n: *n,
            dpdu: *dpdu,
            dpdv: *dpdv,
            dndu: *dndu,
            dndv: *dndv,
        };
        SurfaceInteraction {
            p: *p,
            p_error: *p_error,
            uv: *uv,
            wo: *wo,
            n: *n,
            medium_interface: MediumInterface::new(),
            dpdu: *dpdu,
            dpdv: *dpdv,
            dndu: *dndu,
            dndv: *dndv,
            time,
            shape: None,
            primitive: None,
            shading,
            bsdf: None,
            bssrdf: None,
            //material: None,
            dpdx: Vector3f::zero(),
            dpdy: Vector3f::zero(),
            dudx: 0.0,
            dudy: 0.0,
            dvdx: 0.0,
            dvdy: 0.0,
            face_index,
        }
    }

    pub fn set_shading_geometry(
        &mut self,
        dpdus: &Vector3f,
        dpdvs: &Vector3f,
        dndus: &Normal3f,
        dndvs: &Normal3f,
        orientation_is_authoritative: bool,
    ) {
        // Compute _shading.n_ for _SurfaceInteraction_
        self.shading.n = Vector3f::cross(dpdus, dpdvs).normalize();
        if orientation_is_authoritative {
            self.n = face_forward(&self.n, &self.shading.n);
        } else {
            self.shading.n = face_forward(&self.shading.n, &self.n);
        }

        // Initialize _shading_ partial derivative values
        self.shading.dpdu = *dpdus;
        self.shading.dpdv = *dpdvs;
        self.shading.dndu = *dndus;
        self.shading.dndv = *dndvs;
    }

    pub fn set_shape(&mut self, shape: &Arc<dyn Shape>) {
        //{
        //    let shape = shape.as_ref().read().unwrap();
        //    if shape.
        //}
        self.shape = Some(Arc::downgrade(shape));
    }

    pub fn get_shape(&self) -> Option<Arc<dyn Shape>> {
        if let Some(wp) = self.shape.as_ref() {
            return wp.upgrade();
        }
        return None;
    }

    pub fn get_primitive(&self) -> Option<Arc<dyn Primitive>> {
        if let Some(wp) = self.primitive.as_ref() {
            return wp.upgrade();
        }
        return None;
    }

    pub fn get_bsdf(&self) -> Option<Arc<BSDF>> {
        return self.bsdf.clone();
    }

    /*
       Point3f o = OffsetRayOrigin(p, p_error, n, d);
       return Ray(o, d, Infinity, time, GetMedium(d));
    */

    pub fn spawn_ray(&self, d: &Vector3f) -> Ray {
        let o = offset_ray_origin(&self.p, &self.p_error, &self.n, d);
        let mut r = Ray::new(&o, d, Float::INFINITY, self.time);
        //if self.medium_interface.inside.is_some() {
        //    print!("spawn_ray: inside\n");
        //}
        r.medium = self.get_medium(d);
        return r;
    }

    pub fn spawn_ray_to_point(&self, p2: &Point3f) -> Ray {
        let d = *p2 - self.p;
        let o = offset_ray_origin(&self.p, &self.p_error, &self.n, &d);
        let mut r = Ray::new(&o, &d, 1.0 - SHADOW_EPSILON, self.time);
        r.medium = self.get_medium(&d);
        return r;
    }

    fn compute_differentials_fail(&mut self) {
        self.dudx = 0.0;
        self.dvdx = 0.0;
        self.dudy = 0.0;
        self.dvdy = 0.0;
        self.dpdx = Vector3f::zero();
        self.dpdy = Vector3f::zero();
    }

    pub fn compute_differentials(&mut self, ray: &RayDifferential) {
        if ray.has_differentials {
            // println!("compute_differentials: has_differentials");
            // Estimate screen space change in $\pt{}$ and $(u,v)$

            let p = self.p;
            let n = self.n;
            // Compute auxiliary intersection points with plane
            let d = Vector3f::dot(&n, &Vector3f::from(p));
            let tx =
                -(Vector3f::dot(&n, &ray.rx_origin) - d) / Vector3f::dot(&n, &ray.rx_direction);
            if !tx.is_finite() {
                self.compute_differentials_fail();
                return;
            }
            let px = ray.rx_origin + tx * ray.rx_direction;
            let ty =
                -(Vector3f::dot(&n, &ray.ry_origin) - d) / Vector3f::dot(&n, &ray.ry_direction);
            if !ty.is_finite() {
                self.compute_differentials_fail();
                return;
            }
            let py = ray.ry_origin + ty * ray.ry_direction;
            self.dpdx = px - p;
            self.dpdy = py - p;

            // Compute $(u,v)$ offsets at auxiliary points

            // Choose two dimensions to use for ray offset computation
            let dim = if n.x.abs() > n.y.abs() && n.x.abs() > n.z.abs() {
                (1, 2)
            } else if n.y.abs() > n.z.abs() {
                (0, 2)
            } else {
                (0, 1)
            };

            // Initialize _A_, _Bx_, and _By_ matrices for offset computation
            let a = [
                [self.dpdu[dim.0], self.dpdv[dim.0]],
                [self.dpdu[dim.1], self.dpdv[dim.1]],
            ];
            let bx = [px[dim.0] - p[dim.0], px[dim.1] - p[dim.1]];
            let by = [py[dim.0] - p[dim.0], py[dim.1] - p[dim.1]];
            if let Some((dudx, dvdx)) = solve_linear_system_2x2(&a, &bx) {
                self.dudx = dudx;
                self.dvdx = dvdx;
            } else {
                self.dudx = 0.0;
                self.dvdx = 0.0;
            }
            if let Some((dudy, dvdy)) = solve_linear_system_2x2(&a, &by) {
                self.dudy = dudy;
                self.dvdy = dvdy;
            } else {
                self.dudy = 0.0;
                self.dvdy = 0.0;
            }
        } else {
            self.compute_differentials_fail();
        }
    }

    pub fn compute_scattering_functions(
        &mut self,
        ray: &RayDifferential,
        arena: &mut MemoryArena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        if let Some(primitive) = self.get_primitive() {
            self.compute_differentials(ray);
            primitive.compute_scattering_functions(self, arena, mode, allow_multiple_lobes);
        }
    }

    pub fn le(&self, w: &Vector3f) -> Spectrum {
        if let Some(primitive) = self.get_primitive() {
            if let Some(light) = primitive.get_area_light() {
                if let Some(area_light) = light.as_ref().as_area_light() {
                    return area_light.l(&Interaction::from(self), w);
                }
            }
        }
        return Spectrum::zero();
    }

    pub fn get_medium(&self, w: &Vector3f) -> Option<Arc<dyn Medium>> {
        if Vector3f::dot(w, &self.n) > 0.0 {
            //print!("get_medium: outside {}\n", self.medium_interface.outside.is_some());
            return self.medium_interface.outside.clone();
        } else {
            //print!("get_medium: inside\n");
            return self.medium_interface.inside.clone();
        }
    }
}

/*
impl From<Interaction> for SurfaceInteraction {
    fn from(from: Interaction) -> SurfaceInteraction {
        return SurfaceInteraction::from_interaction(&from);
    }
}
*/
