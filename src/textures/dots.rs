use crate::core::camera::*;
use crate::core::distribution::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::lightdistrib::*;
use crate::core::lowdiscrepancy::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::refrection::*;
use crate::core::rng::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use crate::core::texture::*;
use crate::core::transform::*;

use std::sync::Arc;

pub struct DotsTexture<T> {
    mapping: Box<dyn TextureMapping2D>,
    outside_dot: Arc<dyn Texture<T>>,
    inside_dot: Arc<dyn Texture<T>>,
}

impl<T> DotsTexture<T> {
    pub fn new(
        mapping: Box<dyn TextureMapping2D>,
        tex1: &Arc<dyn Texture<T>>,
        tex2: &Arc<dyn Texture<T>>,
    ) -> Self {
        return DotsTexture::<T> {
            mapping,
            outside_dot: Arc::clone(tex1),
            inside_dot: Arc::clone(tex2),
        };
    }
}

impl<T: Copy> Texture<T> for DotsTexture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let (st, _dstdx, _dstdy) = self.mapping.map(si);
        let s_cell = Float::floor(st[0] + 0.5);
        let t_cell = Float::floor(st[1] + 0.5);
        if noise(s_cell + 0.5, t_cell + 0.5, 0.0) > 0.0 {
            let radius = 0.35;
            let max_shift = 0.5 - radius;
            let s_center = s_cell + max_shift * noise(s_cell + 1.5, t_cell + 2.8, 0.0);
            let t_center = t_cell + max_shift * noise(s_cell + 4.5, t_cell + 9.8, 0.0);
            let dst = st - Point2f::new(s_center, t_center);
            if dst.length_squared() < radius * radius {
                self.inside_dot.as_ref().evaluate(si);
            }
        }
        return self.outside_dot.as_ref().evaluate(si);
    }
}

pub fn create_dots_float_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let map = create_texture_mapping2d(tex2world, tp)?;
    let tex1 = tp.get_float_texture("tex1", 1.0);
    let tex2 = tp.get_float_texture("tex2", 0.0);
    return Ok(Arc::new(DotsTexture::<Float>::new(map, &tex1, &tex2)));
}

pub fn create_dots_spectrum_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let map = create_texture_mapping2d(tex2world, tp)?;
    let tex1 = tp.get_spectrum_texture("tex1", &Spectrum::from(1.0));
    let tex2 = tp.get_spectrum_texture("tex2", &Spectrum::from(0.0));
    return Ok(Arc::new(DotsTexture::<Spectrum>::new(map, &tex1, &tex2)));
}
