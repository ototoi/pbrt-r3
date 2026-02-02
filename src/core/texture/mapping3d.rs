use crate::core::base::*;
use crate::core::interaction::*;
use crate::core::transform::*;

#[derive(Debug, Clone, Copy)]
pub enum TextureMapping3D {
    Identity(IdentityMapping3D),
}

impl TextureMapping3D {
    pub fn map(&self, si: &SurfaceInteraction) -> (Point3f, Vector3f, Vector3f) {
        match self {
            TextureMapping3D::Identity(m) => m.map(si),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct IdentityMapping3D {
    world_to_tex: Transform,
}

impl IdentityMapping3D {
    pub fn new(world_to_tex: &Transform) -> Self {
        IdentityMapping3D {
            world_to_tex: *world_to_tex,
        }
    }

    fn map(&self, si: &SurfaceInteraction) -> (Point3f, Vector3f, Vector3f) {
        let p = self.world_to_tex.transform_point(&si.p);
        let dpdx = self.world_to_tex.transform_vector(&si.dpdx);
        let dpdy = self.world_to_tex.transform_vector(&si.dpdx);
        return (p, dpdx, dpdy);
    }
}
