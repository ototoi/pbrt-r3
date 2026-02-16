use super::vertex::*;
use super::vertex_interaction::*;
use crate::core::prelude::*;

use std::sync::Arc;

thread_local!(static PATHS: StatPercent = StatPercent::new("Integrator/Zero-radiance paths"));
thread_local!(static PATH_LENGTH: StatIntDistribution = StatIntDistribution::new("Integrator/Path length"));

// RandomWalk
fn random_walk(
    scene: &Scene,
    ray: &RayDifferential,
    sampler: &mut dyn Sampler,
    arena: &mut MemoryArena,
    beta: &Spectrum,
    pdf: Float,
    max_depth: usize,
    mode: TransportMode,
    path: &mut Vec<Vertex>,
) {
    assert!(!path.is_empty());
    let mut ray = ray.clone();

    // Declare variables for forward and reverse probability densities
    let mut beta = *beta;

    #[allow(unused_assignments)]
    let mut pdf_fwd = pdf;
    #[allow(unused_assignments)]
    let mut pdf_rev = 0.0;
    loop {
        // println!("bounces: {}", bounces);
        // Attempt to create the next subpath vertex in _path_
        let mut mi = None;

        // vlog(2) << "Random walk. Path: " << *this << ", bounces: " << bounces << ", beta: " << beta << ", pdfFwd: " << pdfFwd << ", pdfRev: " << pdfRev;

        // Trace a ray and sample the medium, if any
        let found_intersection = scene.intersect(&ray.ray);
        if let Some(medium) = ray.ray.medium.as_ref() {
            let (spec, m) = medium.sample(&ray.ray, sampler, arena);
            if let Some(mut m) = m {
                m.medium_interface = MediumInterface::from(medium);
                mi = Some(m);
            }
            beta *= spec;
        };
        if beta.is_black() {
            break;
        }

        let prev_index = path.len() - 1;
        let curr_index = path.len();

        if let Some(mi) = mi.as_ref() {
            // mi.is_valid() && mi.phase.is_some() {
            // Record medium interaction in _path_ and compute forward density
            let prev = &path[prev_index];
            let vertex = Vertex::create_medium(mi, &beta, pdf_fwd, prev);
            path.push(vertex);
            if path.len() >= max_depth {
                break;
            }

            // Sample direction and compute reverse density at preceding vertex
            let wo = -ray.ray.d;
            let (pdf, wi) = mi.phase.sample_p(&wo, &sampler.get_2d());
            assert!(pdf >= 0.0);
            pdf_fwd = pdf;
            pdf_rev = pdf;

            ray = mi.spawn_ray(&wi).into();
        } else {
            // Handle surface interaction for path generation
            if found_intersection.is_none() {
                // Capture escaped rays when tracing from the camera
                if mode == TransportMode::Radiance {
                    let ei = EndpointInteraction::from_ray(&ray.ray);
                    let vertex = Vertex::create_light_from_endpoint(&ei, &beta, pdf_fwd);
                    path.push(vertex);
                }
                break;
            }

            assert!(found_intersection.is_some());

            // Compute scattering functions for _mode_ and skip over medium
            // boundaries
            let mut isect = found_intersection.unwrap();
            isect.compute_scattering_functions(&ray, arena, mode, true);
            if isect.bsdf.is_none() {
                ray = isect.spawn_ray(&ray.ray.d).into();
                continue;
            }

            // Initialize _vertex_ with surface intersection information
            let prev = &path[prev_index];
            let vertex = Vertex::create_surface(&isect, &beta, pdf_fwd, prev);
            path.push(vertex);
            if path.len() >= max_depth {
                break;
            }

            // Sample BSDF at current vertex and compute reverse probability
            {
                let wo = isect.wo;
                let bsdf = isect.bsdf.as_ref().unwrap();
                // Sample BSDF at current vertex and compute reverse probability
                if let Some((f, wi, pdf, t)) = bsdf.sample_f(&wo, &sampler.get_2d(), BSDF_ALL) {
                    pdf_fwd = pdf;
                    if f.is_black() || pdf == 0.0 {
                        break;
                    }
                    beta *= f * (Vector3f::abs_dot(&wi, &isect.shading.n) / pdf_fwd);
                    pdf_rev = bsdf.pdf(&wi, &wo, BSDF_ALL);
                    assert!(pdf_rev >= 0.0);
                    if (t & BSDF_SPECULAR) != 0 {
                        path[curr_index].delta = true;
                        pdf_rev = 0.0;
                        pdf_fwd = 0.0;
                    }
                    beta *= correct_shading_normal(&isect, &wo, &wi, mode);

                    ray = isect.spawn_ray(&wi).into();
                } else {
                    break;
                };
            }
        }

        // Compute reverse area density at preceding vertex
        {
            let pdf_rev = path[curr_index].convert_density(pdf_rev, &path[prev_index]);
            assert!(pdf_rev >= 0.0);
            path[prev_index].pdf_rev = pdf_rev;
        }
    }
}

// GenerateLightSubpath

pub fn generate_camera_subpath(
    scene: &Scene,
    sampler: &mut dyn Sampler,
    arena: &mut MemoryArena,
    max_depth: usize,
    camera: &Arc<dyn Camera>,
    p_film: &Point2f,
    path: &mut Vec<Vertex>,
) -> usize {
    let _p = ProfilePhase::new(Prof::BDPTGenerateSubpath);

    if max_depth == 0 {
        return path.len();
    }

    // Sample initial ray for camera subpath
    let camera_sample = CameraSample {
        p_film: *p_film,
        time: sampler.get_1d(),
        p_lens: sampler.get_2d(),
    };

    if let Some((beta, mut ray)) = camera.generate_ray_differential(&camera_sample) {
        let beta = Spectrum::from(beta);
        ray.scale_differentials(1.0 / Float::sqrt(sampler.get_samples_per_pixel() as Float));
        if let Some((_pdf_pos, pdf_dir)) = camera.pdf_we(&ray.ray) {
            let new_vertex = Vertex::create_camera_from_ray(camera, &ray.ray, &beta);
            path.push(new_vertex);
            if max_depth > 1 {
                random_walk(
                    scene,
                    &ray,
                    sampler,
                    arena,
                    &beta,
                    pdf_dir,
                    max_depth,
                    TransportMode::Radiance,
                    path,
                );
            }
        }
    }
    return path.len();
}

pub fn generate_light_subpath(
    scene: &Scene,
    sampler: &mut dyn Sampler,
    arena: &mut MemoryArena,
    max_depth: usize,
    time: Float,
    light_distr: &Distribution1D,
    light_to_index: &LightIndexMap,
    path: &mut Vec<Vertex>,
) -> usize {
    let _p = ProfilePhase::new(Prof::BDPTGenerateSubpath);

    if max_depth == 0 {
        return path.len();
    }

    let (light_num, light_pdf, _remapped) = light_distr.sample_discrete(sampler.get_1d());
    let light = scene.lights[light_num].clone();
    if let Some((le, ray, n_light, pdf_pos, pdf_dir)) =
        light.sample_le(&sampler.get_2d(), &sampler.get_2d(), time)
    {
        if pdf_pos <= 0.0 || pdf_dir <= 0.0 || le.is_black() {
            return path.len();
        }

        // Generate first vertex on light subpath and start random walk
        let vertex = Vertex::create_light_from_ray(
            &light,
            &ray,
            &n_light,
            &le,
            pdf_pos * light_pdf,
        );
        path.push(vertex);
        if max_depth == 1 {
            return path.len();
        }

        let beta = le * (n_light.abs_dot(&ray.d) / (light_pdf * pdf_pos * pdf_dir));
        let ray = RayDifferential::from(&ray);
        assert!(max_depth > 1);

        random_walk(
            scene,
            &ray,
            sampler,
            arena,
            &beta,
            pdf_dir,
            max_depth,
            TransportMode::Importance,
            path,
        );

        // Correct subpath sampling densities for infinite area lights
        if path[0].is_infinite_light() {
            // Set spatial density of _path[1]_ for infinite area light
            if path.len() > 1 {
                let mut pdf_fwd = pdf_pos;
                if path[1].is_on_surface() {
                    pdf_fwd *= Vector3f::abs_dot(&ray.ray.d, &path[1].get_ng());
                }
                path[1].pdf_fwd = pdf_fwd;
            }
            // Set spatial density of _path[0]_ for infinite area light
            let pdf_fwd = infinite_light_density(scene, light_distr, light_to_index, &ray.ray.d);
            path[0].pdf_fwd = pdf_fwd;
        }
    }
    return path.len();
}

#[inline]
fn remap0(f: Float) -> Float {
    return if f != 0.0 { f } else { 1.0 };
}

fn g(scene: &Scene, sampler: &mut dyn Sampler, v0: &Vertex, v1: &Vertex) -> Spectrum {
    let p0 = v0.get_p();
    let p1 = v1.get_p();
    let mut d = p0 - p1;
    let mut g = 1.0 / d.length_squared();
    d *= Float::sqrt(g);
    if v0.is_on_surface() {
        g *= Vector3f::abs_dot(&v0.get_ns(), &d);
    }
    if v1.is_on_surface() {
        g *= Vector3f::abs_dot(&v1.get_ns(), &d);
    }
    let vis = VisibilityTester::from((v0.get_interaction(), v1.get_interaction()));
    return g * vis.tr(scene, sampler);
}

fn mis_weight(
    scene: &Scene,
    light_vertices: &[Vertex],
    camera_vertices: &[Vertex],
    sampled: &Option<Vertex>,
    s: i32,
    t: i32,
    light_pdf: &Distribution1D,
    light_to_index: &LightIndexMap,
) -> Float {
    if (s + t) == 2 {
        return 1.0;
    }

    // Temporarily update vertex properties for current strategy (value-local copies)
    let mut qs = if s > 0 {
        Some(light_vertices[(s - 1) as usize].clone())
    } else {
        None
    };
    let mut pt = if t > 0 {
        Some(camera_vertices[(t - 1) as usize].clone())
    } else {
        None
    };
    let mut qs_minus = if s > 1 {
        Some(light_vertices[(s - 2) as usize].clone())
    } else {
        None
    };
    let mut pt_minus = if t > 1 {
        Some(camera_vertices[(t - 2) as usize].clone())
    } else {
        None
    };

    if s == 1 {
        qs = sampled.clone();
    } else if t == 1 {
        pt = sampled.clone();
    }

    if let Some(v) = pt.as_mut() {
        v.delta = false;
    }
    if let Some(v) = qs.as_mut() {
        v.delta = false;
    }

    if let Some(v) = pt.as_mut() {
        let pdf = if s > 0 {
            qs.as_ref().unwrap()
                .pdf(scene, qs_minus.as_ref(), v)
        } else {
            v.pdf_light_origin(scene, pt_minus.as_ref().unwrap(), light_pdf, light_to_index)
        };
        assert!(pdf >= 0.0);
        v.pdf_rev = pdf;
    }

    if let Some(v) = pt_minus.as_mut() {
        let pdf = if s > 0 {
            pt.as_ref().unwrap().pdf(scene, qs.as_ref(), v)
        } else {
            pt.as_ref().unwrap().pdf_light(scene, v)
        };
        assert!(pdf >= 0.0);
        v.pdf_rev = pdf;
    }

    if let Some(v) = qs.as_mut() {
        if let Some(ptv) = pt.as_ref() {
            let pdf = ptv.pdf(scene, pt_minus.as_ref(), v);
            assert!(pdf >= 0.0);
            v.pdf_rev = pdf;
        }
    }

    if let Some(v) = qs_minus.as_mut() {
        if let Some(qsv) = qs.as_ref() {
            let pdf = qsv.pdf(scene, pt.as_ref(), v);
            assert!(pdf >= 0.0);
            v.pdf_rev = pdf;
        }
    }

    let mut sum_ri = 0.0;
    // Consider hypothetical connection strategies along the camera subpath
    {
        let mut ri = 1.0;
        let mut i = t - 1; //0..
        while i > 0 {
            let vert = if i == t - 1 {
                pt.as_ref().unwrap_or(&camera_vertices[i as usize])
            } else if i == t - 2 {
                pt_minus.as_ref().unwrap_or(&camera_vertices[i as usize])
            } else {
                &camera_vertices[i as usize]
            };
            let prev = if (i - 1) == t - 1 {
                pt.as_ref().unwrap_or(&camera_vertices[(i - 1) as usize])
            } else if (i - 1) == t - 2 {
                pt_minus.as_ref().unwrap_or(&camera_vertices[(i - 1) as usize])
            } else {
                &camera_vertices[(i - 1) as usize]
            };

            let pdf_rev = vert.pdf_rev;
            let pdf_fwd = vert.pdf_fwd;
            assert!(pdf_fwd >= 0.0);
            assert!(pdf_rev >= 0.0);

            let pdf_delta = remap0(pdf_rev) / remap0(pdf_fwd);
            assert!(pdf_delta.is_finite());
            assert!(pdf_delta >= 0.0);
            ri *= pdf_delta;

            let cur_delta = vert.delta;
            let prv_delta = prev.delta;
            if !cur_delta && !prv_delta {
                sum_ri += ri;
            }
            i -= 1;
        }
    }

    // Consider hypothetical connection strategies along the light subpath
    {
        let mut ri = 1.0;
        let mut i = s - 1; //0..
        while i >= 0 {
            let vert = if i == s - 1 {
                qs.as_ref().unwrap_or(&light_vertices[i as usize])
            } else if i == s - 2 {
                qs_minus.as_ref().unwrap_or(&light_vertices[i as usize])
            } else {
                &light_vertices[i as usize]
            };
            //let prev = &light_vertices[(i - 1) as usize];

            let pdf_rev = vert.pdf_rev;
            let pdf_fwd = vert.pdf_fwd;
            assert!(pdf_fwd >= 0.0);
            assert!(pdf_rev >= 0.0);

            let pdf_delta = remap0(pdf_rev) / remap0(pdf_fwd);
            assert!(pdf_delta.is_finite());
            assert!(pdf_delta >= 0.0);
            // pbrt-r3
            // let pdf_delta = cushion_pdf(pdf_delta);
            // pbrt-r3
            ri *= pdf_delta;

            let cur_delta = vert.delta;

            let delta_light_vertex = if i > 0 {
                let prev = if (i - 1) == s - 1 {
                    qs.as_ref().unwrap_or(&light_vertices[(i - 1) as usize])
                } else if (i - 1) == s - 2 {
                    qs_minus.as_ref().unwrap_or(&light_vertices[(i - 1) as usize])
                } else {
                    &light_vertices[(i - 1) as usize]
                };
                prev.delta
            } else {
                let v0 = &light_vertices[0];
                v0.is_delta_light()
            };
            if !cur_delta && !delta_light_vertex {
                sum_ri += ri;
            }
            i -= 1;
        }
    }

    return 1.0 / (1.0 + sum_ri);
}

// ConnectBDPT
pub fn connect_bdpt(
    scene: &Scene,
    light_vertices: &[Vertex],
    camera_vertices: &[Vertex],
    s: i32,
    t: i32,
    light_distr: &Distribution1D,
    light_to_index: &LightIndexMap,
    camera: &Arc<dyn Camera>,
    sampler: &mut dyn Sampler,
    p_raster: &Point2f,
) -> Option<(Spectrum, Float, Point2f)> {
    let _p = ProfilePhase::new(Prof::BDPTConnectSubpaths);

    let mut p_raster = *p_raster;
    // Ignore invalid connections related to infinite area lights
    if t > 1 && s != 0 {
        let t = camera_vertices[(t - 1) as usize].get_type();
        if t == VertexType::Light {
            return None;
        }
    }

    let mut l = Spectrum::zero();
    let mut sampled: Option<Vertex> = None;
                            // Perform connection and write contribution to _L_
    if s == 0 {
        assert!(t >= 1);
        // Interpret the camera subpath as a complete path
        let pt = &camera_vertices[(t - 1) as usize];
        if pt.is_light() {
            assert!(t >= 2);
            let target = &camera_vertices[(t - 2) as usize];
            l = pt.le(scene, target) * pt.beta;
        }
    } else if t == 1 {
        assert!(s >= 1);
        // Sample a point on the camera and connect it to the light subpath
        let qs = &light_vertices[(s - 1) as usize];
        if qs.is_connectible() {
            let intersection = &qs.get_interaction();
            if let Some((spec, wi, pdf, pr, vis)) =
                camera.as_ref().sample_wi(intersection, &sampler.get_2d())
            {
                p_raster = pr;
                if pdf > 0.0 && !spec.is_black() {
                    // Initialize dynamically sampled vertex and _L_ for $t=1$ case
                    let sampled_beta = spec / pdf;
                    let sampled_v =
                        Vertex::create_camera_from_interaction(camera, &vis.p1, &sampled_beta);
                    l = qs.beta * qs.f(&sampled_v, TransportMode::Importance) * sampled_beta;
                    if qs.is_on_surface() {
                        l *= Vector3f::abs_dot(&wi, &qs.get_ns());
                    }
                    assert!(l.is_valid());
                    // Only check visibility after we know that the path would
                    // make a non-zero contribution.
                    if !l.is_black() {
                        l *= vis.tr(scene, sampler);
                    }

                    sampled = Some(sampled_v);
                }
            }
        }
    } else if s == 1 {
        assert!(t >= 1);
        // Sample a point on a light and connect it to the camera subpath
        let pt = &camera_vertices[(t - 1) as usize];
        if pt.is_connectible() {
            let (light_num, light_pdf, _remapped) = light_distr.sample_discrete(sampler.get_1d());
            let light = scene.lights[light_num].clone();
            let inter = pt.get_interaction();
            if let Some((light_weight, wi, pdf, vis)) = light.sample_li(&inter, &sampler.get_2d()) {
                if pdf > 0.0 && !light_weight.is_black() {
                    let ei = EndpointInteraction::from_light_interaction(&light, &vis.p1);
                    let sampled_beta = light_weight / (pdf * light_pdf);
                    let mut sampled_v =
                        Vertex::create_light_from_endpoint(&ei, &sampled_beta, 0.0);
                    let pdf_fwd = sampled_v.pdf_light_origin(scene, pt, light_distr, light_to_index);
                    assert!(pdf_fwd >= 0.0);
                    sampled_v.pdf_fwd = pdf_fwd;
                    l = pt.beta * pt.f(&sampled_v, TransportMode::Radiance) * sampled_beta;

                    if pt.is_on_surface() {
                        l *= Vector3f::abs_dot(&wi, &pt.get_ns());
                    }
                    // Only check visibility if the path would carry radiance.
                    if !l.is_black() {
                        l *= vis.tr(scene, sampler);
                    }

                    sampled = Some(sampled_v);
                }
            }
        }
    } else {
        assert!(s >= 1);
        assert!(t >= 1);
        // Handle all other bidirectional connection cases
        let qs = &light_vertices[(s - 1) as usize];
        let pt = &camera_vertices[(t - 1) as usize];
        if qs.is_connectible() && pt.is_connectible() {
            l = qs.beta
                * qs.f(pt, TransportMode::Importance)
                * pt.beta
                * pt.f(qs, TransportMode::Radiance);
            if !l.is_black() {
                l *= g(scene, sampler, qs, pt);
            }
        }
    }

    if l.is_black() {
        PATHS.with(|stat| {
            stat.add_num(1); //zeroRadiancePathså
            stat.add_denom(1); //totalPaths
        });
    } else {
        PATHS.with(|stat| {
            stat.add_denom(1); //totalPaths
        });
    }

    PATH_LENGTH.with(|stat| {
        stat.add((s + t - 2).max(0) as u64);
    });

    // Compute MIS weight for connection strategy
    let mis_weight = if l.is_black() {
        0.0
    } else {
        mis_weight(
            scene,
            light_vertices,
            camera_vertices,
            &sampled,
            s,
            t,
            light_distr,
            light_to_index,
        )
    };

    assert!(l.is_valid());
    assert!(mis_weight.is_finite());

    l *= mis_weight;

    return Some((l, mis_weight, p_raster));
}
