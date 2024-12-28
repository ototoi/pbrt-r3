use crate::core::pbrt::*;

use log::*;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, Clone, Copy)]
pub enum AOVTarget {
    Distance,
    Depth,
    N,
    NS,
    UV,
    RDXC,
    RDYC,
    DRODX,
    DRDDX,
    DPDX,
    DPDY,
    DPDU,
    DPDV,
    DUVDX,
    DUVDY,
    DPDUS,
    DPDVS,
}

fn get_aov_target(name: &str) -> Result<AOVTarget, PbrtError> {
    match name {
        "distance" => Ok(AOVTarget::Distance),
        "depth" => Ok(AOVTarget::Depth),
        "n" => Ok(AOVTarget::N),
        "ng" => Ok(AOVTarget::N),
        "ns" => Ok(AOVTarget::NS),
        "uv" => Ok(AOVTarget::UV),
        "rdxc" => Ok(AOVTarget::RDXC),
        "rdyc" => Ok(AOVTarget::RDYC),
        "drodx" => Ok(AOVTarget::DRODX),
        "drddx" => Ok(AOVTarget::DRDDX),
        "dpdx" => Ok(AOVTarget::DPDX),
        "dpdy" => Ok(AOVTarget::DPDY),
        "dpdu" => Ok(AOVTarget::DPDU),
        "dpdv" => Ok(AOVTarget::DPDV),
        "dstdx" => Ok(AOVTarget::DUVDX),
        "dstdy" => Ok(AOVTarget::DUVDY),
        "duvdx" => Ok(AOVTarget::DUVDX),
        "dpdus" => Ok(AOVTarget::DPDUS),
        "dpdvs" => Ok(AOVTarget::DPDVS),
        "shading.n" => Ok(AOVTarget::NS),
        "shading.dpdu" => Ok(AOVTarget::DPDUS),
        "shading.dpdv" => Ok(AOVTarget::DPDVS),
        _ => {
            let msg = format!("AOV target \"{}\" unknown.", name);
            return Err(PbrtError::error(&msg));
        }
    }
}

fn v2c(v: &Vector3f) -> Spectrum {
    let v = 0.5 * *v + Vector3f::from(0.5);
    return Spectrum::new(v[0], v[1], v[2]).clamp(0.0, 1.0);
}

pub struct AOVIntegrator {
    base: BaseSamplerIntegrator,
    target: AOVTarget,
    scale: Float,
}

impl AOVIntegrator {
    pub fn new(
        camera: &Arc<dyn Camera>,
        sampler: &Arc<RwLock<dyn Sampler>>,
        pixel_bounds: &Bounds2i,
        target: AOVTarget,
        scale: Float,
    ) -> Self {
        AOVIntegrator {
            base: BaseSamplerIntegrator::new(camera, sampler, pixel_bounds),
            target,
            scale,
        }
    }
}

impl Integrator for AOVIntegrator {
    fn render(&mut self, scene: &Scene) {
        BaseSamplerIntegrator::render(self, scene);
    }

    fn get_camera(&self) -> Arc<dyn Camera> {
        return self.base.camera.clone();
    }
}

impl SamplerIntegrator for AOVIntegrator {
    fn preprocess(&mut self, scene: &Scene, sampler: &mut dyn Sampler) {
        //
    }

    fn li(
        &self,
        r: &RayDifferential,
        scene: &Scene,
        _sampler: &mut dyn Sampler,
        arena: &mut MemoryArena,
        _depth: i32,
    ) -> Spectrum {
        let _p = ProfilePhase::new(Prof::SamplerIntegratorLi);

        let scale = self.scale;
        if let Some(mut si) = scene.intersect(&r.ray) {
            si.compute_scattering_functions(r, arena, TransportMode::Radiance, true);
            match self.target {
                AOVTarget::Distance => {
                    let dist = r.ray.t_max.get() / r.ray.d.length();
                    return Spectrum::new(dist, dist, dist) * scale;
                }
                AOVTarget::Depth => {
                    let dist = si.p - r.ray.o;
                    return Spectrum::new(dist.length(), dist.length(), dist.length()) * scale;
                }
                AOVTarget::N => {
                    return v2c(&si.n) * scale;
                }
                AOVTarget::NS => {
                    return v2c(&si.shading.n) * scale;
                }
                AOVTarget::UV => {
                    let v = si.uv;
                    return Spectrum::new(v[0], v[1], 0.0) * scale;
                }
                AOVTarget::RDXC => {
                    let v = r.rx_origin;
                    return v2c(&v) * scale;
                }
                AOVTarget::DRODX => {
                    let v = r.rx_origin;
                    return v2c(&v) * scale;
                }
                AOVTarget::DRDDX => {
                    let v = r.rx_direction;
                    return v2c(&v) * scale;
                }
                AOVTarget::DPDX => {
                    let v = si.dpdx;
                    return v2c(&v) * scale;
                }
                AOVTarget::DPDY => {
                    let v = si.dpdy;
                    return v2c(&v) * scale;
                }
                AOVTarget::DPDU => {
                    let v = si.dpdu;
                    return v2c(&v) * scale;
                }
                AOVTarget::DPDV => {
                    let v = si.dpdv;
                    return v2c(&v) * scale;
                }
                AOVTarget::DUVDX => {
                    let v = Vector3::new(Float::abs(si.dudx), Float::abs(si.dvdx), 0.0);
                    return v2c(&v) * scale;
                }
                AOVTarget::DUVDY => {
                    let v = Vector3::new(Float::abs(si.dudy), Float::abs(si.dvdy), 0.0);
                    return v2c(&v) * scale;
                }
                AOVTarget::DPDUS => {
                    let v = si.shading.dpdu;
                    return v2c(&v) * scale;
                }
                AOVTarget::DPDVS => {
                    let v = si.shading.dpdv;
                    return v2c(&v) * scale;
                }
                _ => {}
            }
        }
        return Spectrum::zero();
    }

    fn get_sampler(&self) -> Arc<RwLock<dyn Sampler>> {
        return Arc::clone(&self.base.sampler);
    }

    fn get_pixel_bounds(&self) -> Bounds2i {
        return self.base.pixel_bounds;
    }
}

unsafe impl Sync for AOVIntegrator {}

pub fn create_aov_integrator(
    params: &ParamSet,
    sampler: &Arc<RwLock<dyn Sampler>>,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    let target = params.find_one_string("target", "uv");
    let scale = params.find_one_float("scale", 1.0);

    let pixel_bounds = camera.get_film().read().unwrap().get_sample_bounds();

    let target = get_aov_target(&target)?;

    return Ok(Arc::new(RwLock::new(AOVIntegrator::new(
        &camera,
        &sampler,
        &pixel_bounds,
        target,
        scale,
    ))));
}
