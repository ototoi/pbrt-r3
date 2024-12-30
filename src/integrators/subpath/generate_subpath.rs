use super::scoped_assignment::*;
use super::vertex::*;
use super::vertex_interaction::*;
use crate::core::pbrt::*;

use std::ops::Deref;
use std::sync::Arc;
use std::sync::RwLock;

thread_local!(static PATHS: StatPercent = StatPercent::new("Integrator/Zero-radiance paths"));
thread_local!(static PATH_LENGTH: StatIntDistribution = StatIntDistribution::new("Integrator/Path length"));

fn random_walk(
    scene: &Scene,
    ray: &RayDifferential,
    sampler: &mut dyn Sampler,
    arena: &mut MemoryArena,
    beta: &Spectrum,
    pdf: Float,
    mode: TransportMode,
    start_index: usize,
    path: &mut Vec<Option<Arc<RwLock<Vertex>>>>,
) -> usize {
    let max_depth = path.len();
    assert!(!path.is_empty());
    let path_offset = start_index - 1;
    let mut ray = ray.clone();

    let mut bounces: usize = 0;
    // Declare variables for forward and reverse probability densities
    let mut beta = *beta;
    let mut pdf_fwd = pdf;
    let mut pdf_rev = 0.0;
    loop {
        // println!("bounces: {}", bounces);
        // Attempt to create the next subpath vertex in _path_
        let mut mi = MediumInteraction::default();

        // vlog(2) << "Random walk. Path: " << *this << ", bounces: " << bounces << ", beta: " << beta << ", pdfFwd: " << pdfFwd << ", pdfRev: " << pdfRev;

        // Trace a ray and sample the medium, if any
        let found_intersection = scene.intersect(&ray.ray);
        if let Some(medium) = ray.ray.medium.as_ref() {
            let (spec, m) = medium.sample(&ray.ray, sampler, arena);
            if let Some(mut m) = m {
                m.medium_interface = MediumInterface::from(medium);
                mi = m;
            }
            beta *= spec;
        };
        if beta.is_black() {
            break;
        }
        let prv_index = path_offset + bounces;
        let cur_index = prv_index + 1;
        if mi.is_valid() {
            if let Some(phase) = mi.phase.as_ref() {
                // Record medium interaction in _path_ and compute forward density
                let prev = path[prv_index].as_ref().unwrap().clone();
                let prev = prev.read().unwrap();
                let vertex = Arc::new(RwLock::new(Vertex::create_medium(
                    &mi, &beta, pdf_fwd, &prev,
                )));
                path[cur_index] = Some(vertex.clone());
                bounces += 1;
                if cur_index + 1 >= max_depth {
                    break;
                }

                // Sample direction and compute reverse density at preceding vertex
                let wo = -ray.ray.d;
                let (pdf, wi) = phase.sample_p(&wo, &sampler.get_2d());
                assert!(pdf >= 0.0);
                pdf_fwd = pdf;
                pdf_rev = pdf;

                ray = mi.spawn_ray(&wi).into();
            }
        } else {
            // Handle surface interaction for path generation
            if found_intersection.is_none() {
                // Capture escaped rays when tracing from the camera
                if mode == TransportMode::Radiance {
                    let ei = EndpointInteraction::from_ray(&ray.ray);
                    let vertex = Arc::new(RwLock::new(Vertex::create_light_from_endpoint(
                        &ei, &beta, pdf_fwd,
                    )));
                    path[cur_index] = Some(vertex);
                    bounces += 1;
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
            let prev = path[prv_index].as_ref().unwrap().clone();
            let prev = prev.read().unwrap();
            let vertex = Arc::new(RwLock::new(Vertex::create_surface(
                &isect, &beta, pdf_fwd, &prev,
            )));
            path[cur_index] = Some(vertex.clone());
            bounces += 1;
            if cur_index + 1 >= max_depth {
                break;
            }

            // Sample BSDF at current vertex and compute reverse probability
            {
                let wo = isect.wo;
                let bsdf = isect.bsdf.as_ref().unwrap().clone();
                let bsdf = bsdf.as_ref();
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
                        let vertex = vertex.as_ref().read().unwrap();
                        vertex.delta.set(true); // = true;
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
            let prev = path[prv_index].as_ref().unwrap();
            let vertex = path[cur_index].as_ref().unwrap();
            let prev = prev.read().unwrap();
            let vertex = vertex.read().unwrap();
            let pdf_rev = vertex.convert_density(pdf_rev, &prev);
            assert!(pdf_rev >= 0.0);
            prev.pdf_rev.set(pdf_rev);
        }
    }
    //println!("breaked!");
    return bounces;
}

// GenerateCameraSubpath
fn generate_camera_subpath_core(
    scene: &Scene,
    sampler: &mut dyn Sampler,
    arena: &mut MemoryArena,
    max_depth: usize,
    camera: &Arc<dyn Camera>,
    p_film: &Point2f,
) -> Vec<Option<Arc<RwLock<Vertex>>>> {
    assert!(max_depth > 0);
    let mut path = vec![None; max_depth];
    // Sample initial ray for camera subpath
    let camera_sample = CameraSample {
        p_film: *p_film,
        time: sampler.get_1d(),
        p_lens: sampler.get_2d(),
    };
    if let Some((beta, mut ray)) = camera.generate_ray_differential(&camera_sample) {
        let beta = Spectrum::from(beta);
        ray.scale_differentials(1.0 / Float::sqrt(sampler.get_samples_per_pixel() as Float));
        let pdf_dir = if let Some((_pdf_pos, pdf_dir)) = camera.pdf_we(&ray.ray) {
            pdf_dir
        } else {
            0.0
        };
        let new_vertex = Arc::new(RwLock::new(Vertex::create_camera_from_ray(
            &camera, &ray.ray, &beta,
        )));
        path[0] = Some(new_vertex);
        if max_depth > 1 {
            random_walk(
                scene,
                &ray,
                sampler,
                arena,
                &beta,
                pdf_dir,
                TransportMode::Radiance,
                1,
                &mut path,
            );
        }
    }
    return path;
}

// GenerateLightSubpath
fn generate_light_subpath_core(
    scene: &Scene,
    sampler: &mut dyn Sampler,
    arena: &mut MemoryArena,
    max_depth: usize,
    time: Float,
    light_distr: &Distribution1D,
    light_to_index: &LightIndexMap,
) -> Vec<Option<Arc<RwLock<Vertex>>>> {
    assert!(max_depth > 0);
    let mut path = vec![None; max_depth];
    let (light_num, light_pdf, _remapped) = light_distr.sample_discrete(sampler.get_1d());
    assert!(light_num < scene.lights.len());
    assert!(light_pdf > 0.0);
    let light = scene.lights[light_num].clone();

    if let Some((le, ray, n_light, pdf_pos, pdf_dir)) =
        light.sample_le(&sampler.get_2d(), &sampler.get_2d(), time)
    {
        if pdf_pos <= 0.0 || pdf_dir <= 0.0 || le.is_black() {
            return path;
        }

        let vertex = Arc::new(RwLock::new(Vertex::create_light_from_ray(
            &light,
            &ray,
            &n_light,
            &le,
            pdf_pos * light_pdf,
        )));
        path[0] = Some(vertex.clone());
        if max_depth == 1 {
            return path;
        }

        let beta = le * n_light.abs_dot(&ray.d) * (1.0 / (light_pdf * pdf_pos * pdf_dir));

        let ray = RayDifferential::from(&ray);
        let n_vertices = random_walk(
            scene,
            &ray,
            sampler,
            arena,
            &beta,
            pdf_dir,
            TransportMode::Importance,
            1,
            &mut path,
        );

        let vertex = vertex.read().unwrap();
        if vertex.is_infinite_light() {
            // Set spatial density of _path[1]_ for infinite area light
            if n_vertices > 0 {
                let next = path[1].as_ref().unwrap();
                let next = next.read().unwrap();
                let mut pdf_fwd = pdf_pos;
                if next.is_on_surface() {
                    pdf_fwd *= Vector3f::abs_dot(&ray.ray.d, &next.get_ng());
                }
                next.pdf_fwd.set(pdf_fwd);
            }
            // Set spatial density of _path[0]_ for infinite area light
            {
                let pdf_fwd =
                    infinite_light_density(scene, light_distr, light_to_index, &ray.ray.d);
                vertex.pdf_fwd.set(pdf_fwd);
            }
        }
    }

    return path;
}

// PealSubpath
fn peal_subpath(
    opath: &mut Vec<Arc<RwLock<Vertex>>>,
    ipath: &Vec<Option<Arc<RwLock<Vertex>>>>,
) -> usize {
    let mut actual_length = 0;
    for i in 0..ipath.len() {
        if let Some(vertex) = ipath[i].as_ref() {
            opath.push(vertex.clone());
            actual_length += 1;
        } else {
            break;
        }
    }
    return actual_length;
}

pub fn generate_camera_subpath(
    scene: &Scene,
    sampler: &mut dyn Sampler,
    arena: &mut MemoryArena,
    max_depth: usize,
    camera: &Arc<dyn Camera>,
    p_film: &Point2f,
    path: &mut Vec<Arc<RwLock<Vertex>>>,
) -> usize {
    if max_depth == 0 {
        return 0;
    }
    let ipath = generate_camera_subpath_core(scene, sampler, arena, max_depth, camera, p_film);
    return peal_subpath(path, &ipath);
}

pub fn generate_light_subpath(
    scene: &Scene,
    sampler: &mut dyn Sampler,
    arena: &mut MemoryArena,
    max_depth: usize,
    time: Float,
    light_distr: &Distribution1D,
    light_to_index: &LightIndexMap,
    path: &mut Vec<Arc<RwLock<Vertex>>>,
) -> usize {
    if max_depth == 0 {
        return 0;
    }
    let ipath = generate_light_subpath_core(
        scene,
        sampler,
        arena,
        max_depth,
        time,
        light_distr,
        light_to_index,
    );
    return peal_subpath(path, &ipath);
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

/*
#[inline]
fn cushion_pdf(pdf: Float) -> Float {
    if !pdf.is_finite() {
        return 1.0;
    } else {
        return Float::max(0.0, pdf);
    }
}
*/

fn mis_weight(
    scene: &Scene,
    light_vertices: &[Arc<RwLock<Vertex>>],
    camera_vertices: &[Arc<RwLock<Vertex>>],
    sampled: &Arc<RwLock<Vertex>>,
    s: i32,
    t: i32,
    light_pdf: &Distribution1D,
    light_to_index: &LightIndexMap,
) -> Float {
    if (s + t) == 2 {
        return 1.0;
    }

    // Temporarily update vertex properties for current strategy

    // Look up connection vertices and their predecessors
    let qs = if s > 0 {
        Some(light_vertices[(s - 1) as usize].clone())
    } else {
        None
    };
    let pt = if t > 0 {
        Some(camera_vertices[(t - 1) as usize].clone())
    } else {
        None
    };

    let qs_minus = if s > 1 {
        Some(light_vertices[(s - 2) as usize].clone())
    } else {
        None
    };

    let pt_minus = if t > 1 {
        Some(camera_vertices[(t - 2) as usize].clone())
    } else {
        None
    };

    // Update sampled vertex for $s=1$ or $t=1$ strategy
    let _a1 = if s == 1 {
        assert!(qs.is_some());
        let qs = qs.as_ref().unwrap();
        let qs = qs.read().unwrap();
        let qs = qs.as_tuple();
        let sampled = sampled.read().unwrap();
        let sampled = sampled.as_tuple();
        let v0 = ScopedAssignment::new(&qs.0, &sampled.0.read().unwrap());
        let v1 = ScopedAssignment::new(&qs.1, &sampled.1.read().unwrap());
        let v2 = ScopedAssignment::new(&qs.2, &sampled.2.read().unwrap());
        let v3 = ScopedAssignment::new(&qs.3, &sampled.3.read().unwrap());
        let v4 = ScopedAssignment::new(&qs.4, &sampled.4.read().unwrap());
        Some((v0, v1, v2, v3, v4))
    } else if t == 1 {
        assert!(pt.is_some());
        let pt = pt.as_ref().unwrap();
        let pt = pt.read().unwrap();
        let pt = pt.as_tuple();
        let sampled = sampled.read().unwrap();
        let sampled = sampled.as_tuple();
        let v0 = ScopedAssignment::new(&pt.0, &sampled.0.read().unwrap());
        let v1 = ScopedAssignment::new(&pt.1, &sampled.1.read().unwrap());
        let v2 = ScopedAssignment::new(&pt.2, &sampled.2.read().unwrap());
        let v3 = ScopedAssignment::new(&pt.3, &sampled.3.read().unwrap());
        let v4 = ScopedAssignment::new(&pt.4, &sampled.4.read().unwrap());
        Some((v0, v1, v2, v3, v4))
    } else {
        None
    };

    // Mark connection vertices as non-degenerate
    let _a2 = if pt.is_some() {
        let pt = pt.as_ref().unwrap();
        let pt = pt.read().unwrap();
        Some(ScopedAssignment::new(&pt.delta.value, &false))
    } else {
        None
    };

    let _a3 = if qs.is_some() {
        let qs = qs.as_ref().unwrap();
        let qs = qs.read().unwrap();
        Some(ScopedAssignment::new(&qs.delta.value, &false))
    } else {
        None
    };

    let _a4 = if pt.is_some() {
        assert!(t > 0);
        let pt = pt.as_ref().unwrap();
        let pt = pt.read().unwrap();
        let pdf = if s > 0 {
            assert!(qs.is_some());
            let qs = qs.as_ref().unwrap();
            let qs = qs.read().unwrap();
            let pdf = qs.pdf(scene, &qs_minus, pt.deref());
            assert!(pdf >= 0.0);
            pdf
        } else {
            assert!(pt_minus.is_some());
            let pt_minus = pt_minus.as_ref().unwrap();
            let pt_minus = pt_minus.read().unwrap();
            let pdf = pt.pdf_light_origin(scene, pt_minus.deref(), light_pdf, light_to_index);
            assert!(pdf >= 0.0);
            pdf
        };
        assert!(pdf >= 0.0);
        Some(ScopedAssignment::new(&pt.pdf_rev.value, &pdf))
    } else {
        None
    };

    let _a5 = if pt_minus.is_some() {
        assert!(pt.is_some());
        let pt_minus = pt_minus.as_ref().unwrap();
        let pt_minus = pt_minus.read().unwrap();
        let pdf = if s > 0 {
            assert!(qs.is_some());
            let pt = pt.as_ref().unwrap();
            let pt = pt.read().unwrap();
            pt.pdf(scene, &qs, pt_minus.deref())
        } else {
            let pt = pt.as_ref().unwrap();
            let pt = pt.read().unwrap();
            pt.pdf_light(scene, pt_minus.deref())
        };
        assert!(pdf >= 0.0);
        Some(ScopedAssignment::new(&pt_minus.pdf_rev.value, &pdf))
    } else {
        None
    };

    let _a6 = if qs.is_some() && pt.is_some() {
        let qs = qs.as_ref().unwrap();
        let qs = qs.read().unwrap();
        let pt = pt.as_ref().unwrap();
        let pt = pt.read().unwrap();
        let pdf = pt.pdf(scene, &pt_minus, qs.deref());
        assert!(pdf >= 0.0);
        Some(ScopedAssignment::new(&qs.pdf_rev.value, &pdf))
    } else {
        None
    };

    let _a7 = if qs_minus.is_some() {
        assert!(qs.is_some());
        let qs_minus = qs_minus.unwrap();
        let qs_minus = qs_minus.read().unwrap();
        let qs = qs.unwrap();
        let qs = qs.read().unwrap();
        let pdf = qs.pdf(scene, &pt, qs_minus.deref());
        assert!(pdf >= 0.0);
        Some(ScopedAssignment::new(&qs_minus.pdf_rev.value, &pdf))
    } else {
        None
    };

    //--------------

    let mut sum_ri = 0.0;
    // Consider hypothetical connection strategies along the camera subpath
    {
        let mut ri = 1.0;
        let mut i = t - 1; //0..
        while i > 0 {
            let vert = &camera_vertices[i as usize];
            let prev = &camera_vertices[i as usize - 1];
            let vert = vert.read().unwrap();
            let prev = prev.read().unwrap();

            let pdf_rev = vert.pdf_rev.get();
            let pdf_fwd = vert.pdf_fwd.get();
            assert!(pdf_fwd >= 0.0);
            assert!(pdf_rev >= 0.0);

            let pdf_delta = remap0(pdf_rev) / remap0(pdf_fwd);
            assert!(pdf_delta.is_finite());
            assert!(pdf_delta >= 0.0);
            // pbrt-r3
            // let pdf_delta = cushion_pdf(pdf_delta);
            // pbrt-r3
            ri *= pdf_delta;

            let cur_delta = vert.delta.get();
            let prv_delta = prev.delta.get();
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
            let vert = &light_vertices[i as usize];
            //let prev = &light_vertices[(i - 1) as usize];
            let vert = vert.read().unwrap();

            let pdf_rev = vert.pdf_rev.get();
            let pdf_fwd = vert.pdf_fwd.get();
            assert!(pdf_fwd >= 0.0);
            assert!(pdf_rev >= 0.0);

            let pdf_delta = remap0(pdf_rev) / remap0(pdf_fwd);
            assert!(pdf_delta.is_finite());
            assert!(pdf_delta >= 0.0);
            // pbrt-r3
            // let pdf_delta = cushion_pdf(pdf_delta);
            // pbrt-r3
            ri *= pdf_delta;

            let cur_delta = vert.delta.get();

            let delta_light_vertex = if i > 0 {
                let prev = &light_vertices[(i - 1) as usize];
                let prev = prev.read().unwrap();
                let prev_delta = prev.delta.get();
                prev_delta
            } else {
                let v0 = &light_vertices[0];
                let v0 = v0.read().unwrap();
                let is_delta_light = v0.is_delta_light();
                is_delta_light
            };
            if !cur_delta && !delta_light_vertex {
                sum_ri += ri;
            }
            i -= 1;
        }
    }

    //println!("sum_ri: {}", sum_ri);

    return 1.0 / (1.0 + sum_ri);
}

// ConnectBDPT
pub fn connect_bdpt(
    scene: &Scene,
    light_vertices: &[Arc<RwLock<Vertex>>],
    camera_vertices: &[Arc<RwLock<Vertex>>],
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
        let t = camera_vertices[(t - 1) as usize].read().unwrap().get_type();
        if t == VertexType::Light {
            return None;
        }
    }

    let mut l = Spectrum::zero();
    let mut sampled = Arc::new(RwLock::new(Vertex::default())); //TODO;  = nullptr;
                                                                // Perform connection and write contribution to _L_
    if s == 0 {
        assert!(t >= 1);
        // Interpret the camera subpath as a complete path
        let pt = camera_vertices[(t - 1) as usize].read().unwrap();
        if pt.is_light() {
            assert!(t >= 2);
            let target = camera_vertices[(t - 2) as usize].read().unwrap();
            let pt_beta = pt.beta.get();
            l = pt.le(scene, &target) * pt_beta;
        }
    } else if t == 1 {
        assert!(s >= 1);
        // Sample a point on the camera and connect it to the light subpath
        let qs = light_vertices[(s - 1) as usize].read().unwrap();
        if qs.is_connectible() {
            let intersection = &qs.get_interaction();
            if let Some((spec, wi, pdf, pr, vis)) =
                camera.as_ref().sample_wi(intersection, &sampler.get_2d())
            {
                if pdf > 0.0 && !spec.is_black() {
                    p_raster = pr;
                    // Initialize dynamically sampled vertex and _L_ for $t=1$ case
                    let sampled_beta = spec * (1.0 / pdf);
                    let sampled_v =
                        Vertex::create_camera_from_interaction(camera, &vis.p1, &sampled_beta);
                    let qs_beta = qs.beta.get();
                    l = qs_beta * qs.f(&sampled_v, TransportMode::Importance) * sampled_beta;
                    if qs.is_on_surface() {
                        l *= Vector3f::abs_dot(&wi, &qs.get_ns());
                    }
                    assert!(l.is_valid());
                    // Only check visibility after we know that the path would
                    // make a non-zero contribution.
                    if !l.is_black() {
                        l *= vis.tr(scene, sampler);
                    }

                    sampled = Arc::new(RwLock::new(sampled_v));
                }
            }
        }
    } else if s == 1 {
        assert!(t >= 1);
        // Sample a point on a light and connect it to the camera subpath
        let pt = camera_vertices[(t - 1) as usize].read().unwrap();
        if pt.is_connectible() {
            let (light_num, light_pdf, _remapped) = light_distr.sample_discrete(sampler.get_1d());
            let light = scene.lights[light_num].clone();
            let inter = pt.get_interaction();
            if let Some((light_weight, wi, pdf, vis)) = light.sample_li(&inter, &sampler.get_2d()) {
                if pdf > 0.0 && !light_weight.is_black() {
                    let ei = EndpointInteraction::from_light_interaction(&light, &vis.p1);
                    let sampled_beta = light_weight * (1.0 / (pdf * light_pdf));
                    let sampled_v = Vertex::create_light_from_endpoint(&ei, &sampled_beta, 0.0);
                    {
                        let pdf_fwd =
                            sampled_v.pdf_light_origin(scene, &pt, light_distr, light_to_index);
                        assert!(pdf_fwd >= 0.0);
                        sampled_v.pdf_fwd.set(pdf_fwd);
                        let pt_beta = pt.beta.get();
                        l = pt_beta * pt.f(&sampled_v, TransportMode::Radiance) * sampled_beta;
                    }
                    if pt.is_on_surface() {
                        l *= Vector3f::abs_dot(&wi, &pt.get_ns());
                    }
                    if !l.is_black() {
                        l *= vis.tr(scene, sampler);
                    }

                    sampled = Arc::new(RwLock::new(sampled_v));
                }
            }
        }
    } else {
        assert!(s >= 1);
        assert!(t >= 1);
        // Handle all other bidirectional connection cases
        let qs = light_vertices[(s - 1) as usize].read().unwrap();
        let pt = camera_vertices[(t - 1) as usize].read().unwrap();
        if qs.is_connectible() && pt.is_connectible() {
            let qs_beta = qs.beta.get();
            let pt_beta = pt.beta.get();
            l = qs_beta
                * qs.f(&pt, TransportMode::Importance)
                * pt_beta
                * pt.f(&qs, TransportMode::Radiance);
            if !l.is_black() {
                l *= g(scene, sampler, &qs, &pt);
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
    if l.is_black() || !l.is_valid() {
        return None;
    } else {
        assert!(l.is_valid() && !l.is_black());
        let mis_weight = mis_weight(
            scene,
            light_vertices,
            camera_vertices,
            &sampled,
            s,
            t,
            light_distr,
            light_to_index,
        );
        return if mis_weight <= 0.0 || !mis_weight.is_finite() {
            None
        } else {
            assert!(mis_weight.is_finite() && mis_weight > 0.0);
            let l = l * mis_weight;
            if l.is_black() || !l.is_valid() {
                None
            } else {
                Some((l, mis_weight, p_raster))
            }
        };
    }
}
