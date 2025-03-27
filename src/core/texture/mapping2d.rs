use crate::core::base::*;
use crate::core::error::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::param_set::*;
use crate::core::transform::*;

pub trait TextureMapping2D {
    fn map(&self, si: &SurfaceInteraction) -> (Point2f, Vector2f, Vector2f);
}

pub struct UVMapping2D {
    pub su: Float,
    pub sv: Float,
    pub du: Float,
    pub dv: Float,
}

impl UVMapping2D {
    pub fn new(su: Float, sv: Float, du: Float, dv: Float) -> Self {
        UVMapping2D { su, sv, du, dv }
    }
}

impl TextureMapping2D for UVMapping2D {
    fn map(&self, si: &SurfaceInteraction) -> (Point2f, Vector2f, Vector2f) {
        // Compute texture differentials for 2D identity mapping
        let dstdx = Vector2f::new(self.su * si.dudx, self.sv * si.dvdx);
        let dstdy = Vector2f::new(self.su * si.dudy, self.sv * si.dvdy);
        let st = Point2f::new(self.su * si.uv[0] + self.du, self.sv * si.uv[1] + self.dv);
        return (st, dstdx, dstdy);
    }
}

pub struct SphericalMapping2D {
    world_to_texture: Transform,
}

impl SphericalMapping2D {
    pub fn new(world_to_texture: &Transform) -> Self {
        SphericalMapping2D {
            world_to_texture: *world_to_texture,
        }
    }
    fn sphere(&self, p: &Point3f) -> Point2f {
        let vec =
            (self.world_to_texture.transform_point(p) - Point3f::new(0.0, 0.0, 0.0)).normalize();
        let theta = spherical_theta(&vec);
        let phi = spherical_phi(&vec);
        return Point2f::new(theta * INV_PI, phi * INV_2_PI);
    }
}

impl TextureMapping2D for SphericalMapping2D {
    fn map(&self, si: &SurfaceInteraction) -> (Point2f, Vector2f, Vector2f) {
        let st = self.sphere(&si.p);
        // Compute texture coordinate differentials for sphere $(u,v)$ mapping
        const DELTA: Float = 0.1;
        let st_delta_x = self.sphere(&(si.p + DELTA * si.dpdx));
        let mut dstdx = (st_delta_x - st) * (1.0 / DELTA);
        let st_delta_y = self.sphere(&(si.p + DELTA * si.dpdy));
        let mut dstdy = (st_delta_y - st) * (1.0 / DELTA);

        // Handle sphere mapping discontinuity for coordinate differentials
        if dstdx[1] > 0.5 {
            dstdx[1] = 1.0 - dstdx[1];
        } else if dstdx[1] < -0.5 {
            dstdx[1] = -(dstdx[1] + 1.0);
        }

        if dstdy[1] > 0.5 {
            dstdy[1] = 1.0 - dstdy[1];
        } else if dstdy[1] < -0.5 {
            dstdy[1] = -(dstdy[1] + 1.0);
        }

        return (st, dstdx, dstdy);
    }
}

pub struct CylindricalMapping2D {
    world_to_texture: Transform,
}

impl CylindricalMapping2D {
    pub fn new(world_to_texture: &Transform) -> Self {
        CylindricalMapping2D {
            world_to_texture: *world_to_texture,
        }
    }
    fn cylinder(&self, p: &Point3f) -> Point2f {
        let vec =
            (self.world_to_texture.transform_point(p) - Point3f::new(0.0, 0.0, 0.0)).normalize();
        return Point2f::new((PI + Float::atan2(vec.y, vec.x)) * INV_2_PI, vec.z);
    }
}

impl TextureMapping2D for CylindricalMapping2D {
    fn map(&self, si: &SurfaceInteraction) -> (Point2f, Vector2f, Vector2f) {
        let st = self.cylinder(&si.p);
        // Compute texture coordinate differentials for sphere $(u,v)$ mapping
        const DELTA: Float = 0.1;
        let st_delta_x = self.cylinder(&(si.p + DELTA * si.dpdx));
        let mut dstdx = (st_delta_x - st) * (1.0 / DELTA);
        let st_delta_y = self.cylinder(&(si.p + DELTA * si.dpdy));
        let mut dstdy = (st_delta_y - st) * (1.0 / DELTA);

        // Handle sphere mapping discontinuity for coordinate differentials
        if dstdx[1] > 0.5 {
            dstdx[1] = 1.0 - dstdx[1];
        } else if dstdx[1] < -0.5 {
            dstdx[1] = -(dstdx[1] + 1.0);
        }

        if dstdy[1] > 0.5 {
            dstdy[1] = 1.0 - dstdy[1];
        } else if dstdy[1] < -0.5 {
            dstdy[1] = -(dstdy[1] + 1.0);
        }

        return (st, dstdx, dstdy);
    }
}

pub struct PlanarMapping2D {
    vs: Vector3f,
    vt: Vector3f,
    ds: Float,
    dt: Float,
}

impl PlanarMapping2D {
    pub fn new(vs: &Vector3f, vt: &Vector3f, ds: Float, dt: Float) -> Self {
        PlanarMapping2D {
            vs: *vs,
            vt: *vt,
            ds,
            dt,
        }
    }
}

impl TextureMapping2D for PlanarMapping2D {
    fn map(&self, si: &SurfaceInteraction) -> (Point2f, Vector2f, Vector2f) {
        let vec = si.p;

        let st = Vector2f::new(
            self.ds + Vector3f::dot(&vec, &self.vs),
            self.dt + Vector3f::dot(&vec, &self.vt),
        );

        let dstdx = Vector2f::new(
            Vector3f::dot(&si.dpdx, &self.vs),
            Vector3f::dot(&si.dpdx, &self.vt),
        );

        let dstdy = Vector2f::new(
            Vector3f::dot(&si.dpdy, &self.vs),
            Vector3f::dot(&si.dpdy, &self.vt),
        );

        return (st, dstdx, dstdy);
    }
}

pub fn create_texture_mapping2d(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Box<dyn TextureMapping2D>, PbrtError> {
    let mapping = tp.find_string("mapping", "uv");
    match mapping.as_ref() {
        "uv" => {
            let su = tp.find_float("uscale", 1.0);
            let sv = tp.find_float("vscale", 1.0);
            let du = tp.find_float("udelta", 0.0);
            let dv = tp.find_float("vdelta", 0.0);
            return Ok(Box::new(UVMapping2D::new(su, sv, du, dv)));
        }
        "spherical" => {
            let it = tex2world.inverse();
            return Ok(Box::new(SphericalMapping2D::new(&it)));
        }
        "cylindrical" => {
            let it = tex2world.inverse();
            return Ok(Box::new(CylindricalMapping2D::new(&it)));
        }
        "planar" => {
            let v1 = tp.find_vector3f("v1", &Vector3f::new(1.0, 0.0, 0.0));
            let v2 = tp.find_vector3f("v2", &Vector3f::new(0.0, 1.0, 0.0));
            let du = tp.find_float("udelta", 0.0);
            let dv = tp.find_float("vdelta", 0.0);
            return Ok(Box::new(PlanarMapping2D::new(&v1, &v2, du, dv)));
        }
        _ => {
            let msg = format!("2D texture mapping \"{}\" unknown", mapping);
            return Err(PbrtError::error(&msg));
        }
    }
}
