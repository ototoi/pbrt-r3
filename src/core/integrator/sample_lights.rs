use crate::core::base::*;
use crate::core::geometry::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::memory::*;
use crate::core::profile::*;
use crate::core::reflection::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::spectrum::*;

fn surface_as_interaction(it: &SurfaceInteraction) -> Interaction {
    Interaction::from(BaseInteraction {
        p: it.p,
        p_error: it.p_error,
        wo: it.wo,
        n: it.n,
        time: it.time,
        medium_interface: it.medium_interface.clone(),
    })
}

pub fn uniform_sample_all_lights(
    it: &Interaction,
    scene: &Scene,
    arena: &mut MemoryArena,
    sampler: &mut dyn Sampler,
    n_light_samples: &[u32],
    handle_media: bool,
) -> Spectrum {
    let _p = ProfilePhase::new(Prof::DirectLighting);

    let mut l = Spectrum::zero();
    for j in 0..n_light_samples.len() {
        let light = scene.lights[j].as_ref();
        let n_samples = n_light_samples[j];
        let u_light_array = sampler.get_2d_array(n_samples);
        let u_scattering_array = sampler.get_2d_array(n_samples);
        match (u_light_array.as_ref(), u_scattering_array.as_ref()) {
            (Some(u_light_array), Some(u_scattering_array)) => {
                // Estimate direct lighting using sample arrays
                let mut ld = Spectrum::zero();
                for k in 0..n_samples as usize {
                    ld += estimate_direct(
                        it,
                        &u_scattering_array[k],
                        light,
                        &u_light_array[k],
                        scene,
                        sampler,
                        arena,
                        handle_media,
                        false,
                    );
                }
                l += ld / (n_samples as Float);
            }
            _ => {
                // Use a single sample for illumination from _light_
                let u_light = sampler.get_2d();
                let u_scattering = sampler.get_2d();
                l += estimate_direct(
                    it,
                    &u_scattering,
                    light,
                    &u_light,
                    scene,
                    sampler,
                    arena,
                    handle_media,
                    false,
                );
            }
        }
    }
    return l;
}

pub fn uniform_sample_one_light(
    it: &Interaction,
    scene: &Scene,
    arena: &mut MemoryArena,
    sampler: &mut dyn Sampler,
    handle_media: bool,
    distrib: Option<&Distribution1D>,
) -> Spectrum {
    let _p = ProfilePhase::new(Prof::DirectLighting);

    let n_lights = scene.lights.len();
    if n_lights == 0 {
        return Spectrum::zero();
    }

    let light_num;
    let light_pdf;
    if let Some(distrib) = distrib {
        let (offset, pdf, _) = distrib.sample_discrete(sampler.get_1d()); //offset, pdf, remapped
        if pdf <= 0.0 {
            return Spectrum::zero();
        }
        light_num = offset;
        light_pdf = pdf;
    } else {
        light_num = usize::min(
            (sampler.get_1d() * n_lights as Float) as usize,
            n_lights - 1,
        );
        light_pdf = 1.0 / n_lights as Float;
    }
    assert!(light_pdf >= 0.0);

    let u_light = sampler.get_2d();
    let u_scattering = sampler.get_2d();
    let light = scene.lights[light_num].as_ref();
    return estimate_direct(
        it,
        &u_scattering,
        light,
        &u_light,
        scene,
        sampler,
        arena,
        handle_media,
        false,
    ) / light_pdf;
}

pub fn uniform_sample_one_light_surface(
    it: &SurfaceInteraction,
    scene: &Scene,
    arena: &mut MemoryArena,
    sampler: &mut dyn Sampler,
    handle_media: bool,
    distrib: Option<&Distribution1D>,
) -> Spectrum {
    let _p = ProfilePhase::new(Prof::DirectLighting);

    let n_lights = scene.lights.len();
    if n_lights == 0 {
        return Spectrum::zero();
    }

    let light_num;
    let light_pdf;
    if let Some(distrib) = distrib {
        let (offset, pdf, _) = distrib.sample_discrete(sampler.get_1d());
        if pdf <= 0.0 {
            return Spectrum::zero();
        }
        light_num = offset;
        light_pdf = pdf;
    } else {
        light_num = usize::min(
            (sampler.get_1d() * n_lights as Float) as usize,
            n_lights - 1,
        );
        light_pdf = 1.0 / n_lights as Float;
    }
    assert!(light_pdf >= 0.0);

    let u_light = sampler.get_2d();
    let u_scattering = sampler.get_2d();
    let light = scene.lights[light_num].as_ref();
    return estimate_direct_surface(
        it,
        &u_scattering,
        light,
        &u_light,
        scene,
        sampler,
        arena,
        handle_media,
        false,
    ) / light_pdf;
}

pub fn estimate_direct(
    it: &Interaction,
    u_scattering: &Point2f,
    light: &dyn Light,
    u_light: &Point2f,
    scene: &Scene,
    sampler: &mut dyn Sampler,
    _arena: &mut MemoryArena,
    handle_media: bool,
    specular: bool,
) -> Spectrum {
    let _p = ProfilePhase::new(Prof::EstimateDirect);

    let bsdf_flags = if specular {
        BSDF_ALL
    } else {
        BSDF_ALL & !BSDF_SPECULAR
    };
    let mut ld = Spectrum::zero();
    //println!("{:?}", bsdf_flags);
    // Sample light source with multiple importance sampling
    {
        let _p = ProfilePhase::new(Prof::SampleLightImportance);

        if let Some((mut li, wi, light_pdf, visibilty)) = light.sample_li(it, u_light) {
            if light_pdf > 0.0 && !li.is_black() {
                // Compute BSDF or phase function's value for light sample
                let mut f = Spectrum::zero();
                let mut scattering_pdf = 0.0;
                {
                    let _p = ProfilePhase::new(Prof::ComputeBSDF);

                    if let Some(isect) = it.as_surface_interaction() {
                        if let Some(bsdf) = isect.bsdf.as_ref() {
                            //assert!(isect.n.length() > 0.0);
                            //assert!(isect.wo.length() > 0.0);
                            //assert!(isect.shading.n.length() > 0.0);
                            f = bsdf.f(&isect.wo, &wi, bsdf_flags)
                                * Vector3f::abs_dot(&wi, &isect.shading.n);
                            scattering_pdf = bsdf.pdf(&isect.wo, &wi, bsdf_flags);

                            //assert!(f.y() >= 0.0);
                            //assert!(scattering_pdf >= 0.0);
                        }
                    } else if let Some(mi) = it.as_medium_interaction() {
                        let p = mi.phase.p(&mi.wo, &wi);
                        f = Spectrum::from(p);
                        scattering_pdf = p;
                    } else {
                        panic!("Interaction is neither SurfaceInteraction nor MediumInteraction");
                    }
                }

                {
                    let _p = ProfilePhase::new(Prof::ComputeLightImportance);

                    if !f.is_black() {
                        if handle_media {
                            li *= visibilty.tr(scene, sampler);
                        } else {
                            if !visibilty.unoccluded(scene) {
                                li = Spectrum::zero();
                            }
                        }

                        if !li.is_black() {
                            if light.is_delta() {
                                ld += f * li / light_pdf;
                            } else {
                                let weight = power_heuristic(1, light_pdf, 1, scattering_pdf);
                                //println!("{:?}: {:?} / {:?}, {:?}", (weight / light_pdf), weight, light_pdf, scattering_pdf);
                                ld += f * li * (weight / light_pdf);
                            }
                        }
                    }
                }
            }
        }
    }

    // Sample BSDF with multiple importance sampling
    {
        let _p = ProfilePhase::new(Prof::SampleBSDFImportance);

        if !light.is_delta() {
            let mut f = Spectrum::zero();
            let mut wi = Vector3f::zero();
            let mut sampled_specular = false;
            let mut scattering_pdf = 0.0;
            if let Some(isect) = it.as_surface_interaction() {
                // Sample scattered direction for surface interactions
                if let Some(bsdf) = isect.bsdf.as_ref() {
                    if let Some((f2, wi2, scattering_pdf2, sampled_type)) =
                        bsdf.sample_f(&isect.wo, u_scattering, bsdf_flags)
                    {
                        assert!(isect.shading.n.length() > 0.0);

                        f = f2 * Vector3f::abs_dot(&wi2, &isect.shading.n);
                        wi = wi2;
                        scattering_pdf = scattering_pdf2;
                        sampled_specular = (sampled_type & BSDF_SPECULAR) != 0;
                    }
                }
            } else if let Some(mi) = it.as_medium_interaction() {
                let (p, wi2) = mi.phase.sample_p(&mi.wo, u_scattering);
                wi = wi2;
                f = Spectrum::from(p);
                scattering_pdf = p;
            }

            if !f.is_black() && scattering_pdf > 0.0 {
                let mut weight = 1.0;
                if !sampled_specular {
                    let light_pdf = light.pdf_li(it, &wi);
                    if light_pdf == 0.0 {
                        return ld;
                    }
                    weight = power_heuristic(1, scattering_pdf, 1, light_pdf);
                }

                let ray = it.spawn_ray(&wi);
                let (opt_light_isect, tr) = if handle_media {
                    scene.intersect_tr(&ray, sampler)
                } else {
                    (scene.intersect(&ray), Spectrum::one())
                };

                let mut li = Spectrum::zero();
                if let Some(light_isect) = opt_light_isect.as_ref() {
                    if let Some(primitive) = light_isect.get_primitive() {
                        if let Some(other) = primitive.get_area_light() {
                            let alight_ptr: *const dyn Light = light;
                            let blight = other.as_ref();
                            let blight_ptr: *const dyn Light = blight;
                            #[allow(clippy::vtable_address_comparisons)]
                            if std::ptr::eq(alight_ptr, blight_ptr) {
                                li = light_isect.le(&-wi);
                            }
                        }
                    }
                } else {
                    li = light.le(&RayDifferential::from(&ray));
                }
                if !li.is_black() {
                    ld += f * li * tr * (weight / scattering_pdf);
                }
            }
        }
    }
    return ld;
}

pub fn estimate_direct_surface(
    it: &SurfaceInteraction,
    u_scattering: &Point2f,
    light: &dyn Light,
    u_light: &Point2f,
    scene: &Scene,
    sampler: &mut dyn Sampler,
    _arena: &mut MemoryArena,
    handle_media: bool,
    specular: bool,
) -> Spectrum {
    let _p = ProfilePhase::new(Prof::EstimateDirect);

    let bsdf_flags = if specular {
        BSDF_ALL
    } else {
        BSDF_ALL & !BSDF_SPECULAR
    };
    let mut ld = Spectrum::zero();
    let it_light = surface_as_interaction(it);

    // Sample light source with multiple importance sampling
    {
        let _p = ProfilePhase::new(Prof::SampleLightImportance);

        if let Some((mut li, wi, light_pdf, visibilty)) = light.sample_li(&it_light, u_light) {
            if light_pdf > 0.0 && !li.is_black() {
                // Compute BSDF's value for light sample
                let mut f = Spectrum::zero();
                let mut scattering_pdf = 0.0;
                {
                    let _p = ProfilePhase::new(Prof::ComputeBSDF);

                    if let Some(bsdf) = it.bsdf.as_ref() {
                        f = bsdf.f(&it.wo, &wi, bsdf_flags) * Vector3f::abs_dot(&wi, &it.shading.n);
                        scattering_pdf = bsdf.pdf(&it.wo, &wi, bsdf_flags);
                    }
                }

                {
                    let _p = ProfilePhase::new(Prof::ComputeLightImportance);

                    if !f.is_black() {
                        if handle_media {
                            li *= visibilty.tr(scene, sampler);
                        } else if !visibilty.unoccluded(scene) {
                            li = Spectrum::zero();
                        }

                        if !li.is_black() {
                            if light.is_delta() {
                                ld += f * li / light_pdf;
                            } else {
                                let weight = power_heuristic(1, light_pdf, 1, scattering_pdf);
                                ld += f * li * (weight / light_pdf);
                            }
                        }
                    }
                }
            }
        }
    }

    // Sample BSDF with multiple importance sampling
    {
        let _p = ProfilePhase::new(Prof::SampleBSDFImportance);

        if !light.is_delta() {
            let mut f = Spectrum::zero();
            let mut wi = Vector3f::zero();
            let mut sampled_specular = false;
            let mut scattering_pdf = 0.0;
            if let Some(bsdf) = it.bsdf.as_ref() {
                if let Some((f2, wi2, scattering_pdf2, sampled_type)) =
                    bsdf.sample_f(&it.wo, u_scattering, bsdf_flags)
                {
                    f = f2 * Vector3f::abs_dot(&wi2, &it.shading.n);
                    wi = wi2;
                    scattering_pdf = scattering_pdf2;
                    sampled_specular = (sampled_type & BSDF_SPECULAR) != 0;
                }
            }

            if !f.is_black() && scattering_pdf > 0.0 {
                let mut weight = 1.0;
                if !sampled_specular {
                    let light_pdf = light.pdf_li(&it_light, &wi);
                    if light_pdf == 0.0 {
                        return ld;
                    }
                    weight = power_heuristic(1, scattering_pdf, 1, light_pdf);
                }

                let ray = it.spawn_ray(&wi);
                let (opt_light_isect, tr) = if handle_media {
                    scene.intersect_tr(&ray, sampler)
                } else {
                    (scene.intersect(&ray), Spectrum::one())
                };

                let mut li = Spectrum::zero();
                if let Some(light_isect) = opt_light_isect.as_ref() {
                    if let Some(primitive) = light_isect.get_primitive() {
                        if let Some(other) = primitive.get_area_light() {
                            let alight_ptr: *const dyn Light = light;
                            let blight = other.as_ref();
                            let blight_ptr: *const dyn Light = blight;
                            #[allow(clippy::vtable_address_comparisons)]
                            if std::ptr::eq(alight_ptr, blight_ptr) {
                                li = light_isect.le(&-wi);
                            }
                        }
                    }
                } else {
                    li = light.le(&RayDifferential::from(&ray));
                }
                if !li.is_black() {
                    ld += f * li * tr * (weight / scattering_pdf);
                }
            }
        }
    }
    ld
}
