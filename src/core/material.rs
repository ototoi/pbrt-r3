use crate::core::interaction::SurfaceInteraction;
use crate::core::memory::MemoryArena;
use crate::core::pbrt::*;
use crate::core::texture::Texture;

use std::sync::Arc;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum TransportMode {
    Radiance,
    Importance,
}

pub fn material_bump(d: &Arc<dyn Texture<Float>>, si: &mut SurfaceInteraction) {
    // Compute offset positions and evaluate displacement texture
    let mut si_eval: SurfaceInteraction = si.clone();

    // Shift _siEval_ _du_ in the $u$ direction
    let mut du = 0.5 * (Float::abs(si.dudx) + Float::abs(si.dudy));
    // The most common reason for du to be zero is for ray that start from
    // light sources, where no differentials are available. In this case,
    // we try to choose a small enough du so that we still get a decently
    // accurate bump value.
    if du == 0.0 {
        du = 0.0005;
    }
    si_eval.p = si.p + du * si.shading.dpdu;
    si_eval.uv = si.uv + Vector2f::new(du, 0.0);
    si_eval.n = (Vector3f::cross(&si.shading.dpdu, &si.shading.dpdv) + du * si.dndu).normalize();
    let u_displace = d.as_ref().evaluate(&si_eval);

    // Shift _siEval_ _dv_ in the $v$ direction
    let mut dv = 0.5 * (Float::abs(si.dvdx) + Float::abs(si.dvdy));
    if dv == 0.0 {
        dv = 0.0005;
    }
    si_eval.p = si.p + dv * si.shading.dpdv;
    si_eval.uv = si.uv + Vector2f::new(0.0, dv);
    si_eval.n = (Vector3f::cross(&si.shading.dpdu, &si.shading.dpdv) + dv * si.dndv).normalize();
    let v_displace = d.as_ref().evaluate(&si_eval);
    let displace = d.as_ref().evaluate(si);

    // Compute bump-mapped differential geometry
    {
        let dpdu = si.shading.dpdu
            + (u_displace - displace) / du * si.shading.n
            + displace * si.shading.dndu;
        let dpdv = si.shading.dpdv
            + (v_displace - displace) / dv * si.shading.n
            + displace * si.shading.dndv;
        let dndu = si.shading.dndu;
        let dndv = si.shading.dndv;
        si.set_shading_geometry(&dpdu, &dpdv, &dndu, &dndv, false);
    }
}

pub trait Material {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        mode: TransportMode,
        _allow_multiple_lobes: bool,
    );

    fn bump(&self, d: &Arc<dyn Texture<Float>>, si: &mut SurfaceInteraction) {
        material_bump(d, si);
    }
}
