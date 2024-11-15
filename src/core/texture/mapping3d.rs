use crate::core::pbrt::*;

pub trait TextureMapping3D {
    fn map(&self, si: &SurfaceInteraction) -> (Point3f, Vector3f, Vector3f);
}

pub struct IdentityMapping3D {
    world_to_tex: Transform,
}

impl IdentityMapping3D {
    pub fn new(world_to_tex: &Transform) -> Self {
        IdentityMapping3D {
            world_to_tex: *world_to_tex,
        }
    }
}

impl TextureMapping3D for IdentityMapping3D {
    fn map(&self, si: &SurfaceInteraction) -> (Point3f, Vector3f, Vector3f) {
        let p = self.world_to_tex.transform_point(&si.p);
        let dpdx = self.world_to_tex.transform_vector(&si.dpdx);
        let dpdy = self.world_to_tex.transform_vector(&si.dpdx);
        return (p, dpdx, dpdy);
    }
}
