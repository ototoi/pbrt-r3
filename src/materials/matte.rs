use crate::core::pbrt::*;
use std::sync::Arc;

pub struct MatteMaterial {
    kd: Arc<dyn Texture<Spectrum>>,
    sigma: Arc<dyn Texture<Float>>,
}

impl MatteMaterial {
    pub fn new(kd: &Arc<dyn Texture<Spectrum>>, sigma: &Arc<dyn Texture<Float>>) -> Self {
        MatteMaterial {
            kd: kd.clone(),
            sigma: sigma.clone(),
        }
    }
}

impl Material for MatteMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        arena: &mut MemoryArena,
        _: TransportMode,
        _: bool,
    ) {
        let mut b = arena.alloc_bsdf(si, 1.0);

        let r = self.kd.as_ref().evaluate(si);
        let sig = Float::clamp(self.sigma.as_ref().evaluate(si), 0.0, 90.0);
        if !r.is_black() {
            if sig == 0.0 {
                let r: Arc<dyn BxDF> = Arc::new(LambertianReflection::new(&r));
                b.add(&r);
            } else {
                let r: Arc<dyn BxDF> = Arc::new(OrenNayar::new(&r, sig));
                b.add(&r);
            }
        }

        si.bsdf = Some(Arc::new(b));
    }
}

pub fn create_matte_material(mp: &TextureParams) -> Result<Arc<dyn Material>, PbrtError> {
    let kd = mp.get_spectrum_texture("Kd", &Spectrum::from(0.5));
    let sigma = mp.get_float_texture("sigma", 0.0);
    return Ok(Arc::new(MatteMaterial::new(&kd, &sigma)));
}