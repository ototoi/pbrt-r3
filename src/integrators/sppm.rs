use crate::core::prelude::*;
use crate::samplers::*;

use std::cell::UnsafeCell;
use std::collections::HashMap;
use std::ops::DerefMut;
use std::sync::Arc;
use std::sync::RwLock;
use std::time::Instant;

use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;

thread_local!(static PHOTON_PATHS: StatCounter = StatCounter::new("Stochastic Progressive Photon Mapping/Photon paths followed"));
thread_local!(static GRID_CELLS_PER_VISIBLE_POINT: StatIntDistribution = StatIntDistribution::new("Stochastic Progressive Photon Mapping/Grid cells per visible point"));
thread_local!(static PIXEL_MEMORY_BYTES: StatMemoryCounter = StatMemoryCounter::new("Memory/SPPM Pixels"));

#[derive(Clone)]
struct SPPMVisiblePoint {
    p: Point3f,
    wo: Vector3f,
    bsdf: Option<Arc<BSDF>>,
    beta: Spectrum,
}

struct SPPMPixel {
    // SPPMPixel Public Methods
    radius: Float,
    ld: Spectrum,
    vp: SPPMVisiblePoint, // VisiblePoint vp;
    n: Float,
    tau: Spectrum,
}

#[derive(Default)]
struct SPPMPixels {
    pixels: Vec<UnsafeCell<SPPMPixel>>,
}
unsafe impl Sync for SPPMPixels {}
unsafe impl Send for SPPMPixels {}
impl SPPMPixels {
    #[inline]
    unsafe fn get(&self, i: usize) -> &SPPMPixel {
        &*self.pixels[i].get()
    }
    #[inline]
    unsafe fn get_mut(&self, i: usize) -> &mut SPPMPixel {
        &mut *self.pixels[i].get()
    }
}

#[derive(Clone)]
struct SPPMHashGrid {
    nodes: Vec<Vec<usize>>,
}
unsafe impl Sync for SPPMHashGrid {}
unsafe impl Send for SPPMHashGrid {}

pub struct SPPMIntegrator {
    pub camera: Arc<dyn Camera>,
    pub initial_search_radius: Float,
    pub n_iterations: u32,
    pub max_depth: i32,
    pub photons_per_iteration: i32,
    pub write_frequency: i32,
}

impl SPPMIntegrator {
    pub fn new(
        camera: Arc<dyn Camera>,
        initial_search_radius: Float,
        n_iterations: i32,
        max_depth: i32,
        photons_per_iteration: i32,
        write_frequency: i32,
    ) -> Self {
        let photons_per_iteration = if photons_per_iteration > 0 {
            photons_per_iteration
        } else {
            let camera = camera.as_ref();
            let film = camera.get_film();
            let film = film.read().unwrap();
            film.cropped_pixel_bounds.area() as i32
        };
        let n_iterations = n_iterations as u32;
        Self {
            camera,
            initial_search_radius,
            n_iterations,
            max_depth,
            photons_per_iteration,
            write_frequency,
        }
    }

    fn get_film(&self) -> Arc<RwLock<Film>> {
        let camera = self.camera.as_ref();
        camera.get_film()
    }
}

#[inline]
fn to_grid(p: &Point3f, bounds: &Bounds3f, grid_res: &[i32; 3]) -> (bool, Point3i) {
    let mut in_bounds = true;
    let pg = bounds.offset(p);
    let mut pi = Point3i::default();
    for i in 0..3 {
        if pg[i] < 0.0 || pg[i] > 1.0 {
            in_bounds = false;
        }
        pi[i] = (grid_res[i] as Float * pg[i]) as i32;
        pi[i] = i32::clamp(pi[i], 0, grid_res[i] - 1);
    }
    return (in_bounds, pi);
}

#[inline]
fn hash(p: &Point3i, hash_size: usize) -> usize {
    let x = p.x as usize * 73856093;
    let y = p.y as usize * 19349663;
    let z = p.z as usize * 83492791;
    let h = (x ^ y ^ z) % hash_size;
    return h;
}

impl Integrator for SPPMIntegrator {
    fn render(&mut self, scene: &Scene) {
        let pixel_bounds = {
            let film = self.get_film();
            let film = film.read().unwrap();
            film.cropped_pixel_bounds
        };
        let timing_enabled = std::env::var_os("PBRT_R3_SPPM_TIMING").is_some();
        let render_t0 = Instant::now();
        let mut t_camera = 0.0f64;
        let mut t_grid = 0.0f64;
        let mut t_photon = 0.0f64;
        let mut t_update = 0.0f64;
        let mut t_image = 0.0f64;
        let mut t_grid_bounds = 0.0f64;
        let mut t_progress = 0.0f64;
        let mut t_finalize = 0.0f64;
        let mut t_setup_film = 0.0f64;
        let mut t_setup_pixels = 0.0f64;
        let mut t_setup_light_distr = 0.0f64;
        let mut t_setup_tiles = 0.0f64;
        let mut t_drop_grid = 0.0f64;

        let setup_film_t0 = Instant::now();
        {
            let film = self.get_film();
            let mut film = film.write().unwrap();
            film.render_start();
        }
        t_setup_film += setup_film_t0.elapsed().as_secs_f64();

        let n_pixels = pixel_bounds.area() as usize;
        let initial_search_radius = self.initial_search_radius;
        let setup_pixels_t0 = Instant::now();
        let mut pixels_vec = Vec::with_capacity(n_pixels);

        PIXEL_MEMORY_BYTES.with(|c| c.add(n_pixels * std::mem::size_of::<SPPMPixel>()));

        for _ in 0..n_pixels {
            let pixel = UnsafeCell::new(SPPMPixel {
                radius: initial_search_radius,
                ld: Spectrum::default(),
                vp: SPPMVisiblePoint {
                    p: Point3f::default(),
                    wo: Vector3f::default(),
                    bsdf: None,
                    beta: Spectrum::zero(),
                },
                n: 0.0,
                tau: Spectrum::zero(),
            });
            pixels_vec.push(pixel);
        }
        let pixels = Arc::new(SPPMPixels { pixels: pixels_vec });
        t_setup_pixels += setup_pixels_t0.elapsed().as_secs_f64();

        // Compute _lightDistr_ for sampling lights proportional to power
        let setup_light_t0 = Instant::now();
        let light_distr = compute_light_power_distribution(scene);
        t_setup_light_distr += setup_light_t0.elapsed().as_secs_f64();

        let n_iterations = self.n_iterations;
        let inv_sqrt_spp = 1.0 / Float::sqrt(n_iterations as Float);
        // Perform _nIterations_ of SPPM integration
        let sampler = HaltonSampler::new(n_iterations, &pixel_bounds, false);

        // Compute number of tiles to use for SPPM camera pass
        let pixel_extent = pixel_bounds.diagonal();
        let tile_size = 16;
        let n_tiles = Point2i::new(
            (pixel_extent.x + tile_size - 1) / tile_size,
            (pixel_extent.y + tile_size - 1) / tile_size,
        );

        let setup_tiles_t0 = Instant::now();
        let mut tiles = Vec::with_capacity(n_tiles.x as usize * n_tiles.y as usize);
        //let mut samplers = Vec::with_capacity(n_tiles.x as usize * n_tiles.y as usize);
        for y in 0..n_tiles.y {
            for x in 0..n_tiles.x {
                let x0 = pixel_bounds.min.x + x * tile_size;
                let x1 = std::cmp::min(x0 + tile_size, pixel_bounds.max.x);
                let y0 = pixel_bounds.min.y + y * tile_size;
                let y1 = std::cmp::min(y0 + tile_size, pixel_bounds.max.y);
                let tile_bounds = Bounds2i::new(&Point2i::new(x0, y0), &Point2i::new(x1, y1));
                let tile_index = y * n_tiles.x + x;
                let sampler = sampler.clone_with_seed(tile_index as u32);
                //
                tiles.push((x, y, tile_bounds, sampler));
                //samplers.push(sampler);
            }
        }
        t_setup_tiles += setup_tiles_t0.elapsed().as_secs_f64();

        let progress = Arc::new(RwLock::new(ProgressReporter::new(
            2 * n_iterations as usize,
            "Rendering",
        )));

        let camera = self.get_camera();
        //let camera = camera.deref();
        let max_depth = self.max_depth;
        for iter in 0..n_iterations {
            // Generate SPPM visible points
            let phase_t0 = Instant::now();
            {
                tiles.par_iter().for_each(|tile| {
                    let mut arena = MemoryArena::new();
                    //let pixel_bounds = tile.5;
                    let pixels = pixels.clone();
                    let camera = camera.as_ref();
                    //let tile = Point2i::new(tile.0, tile.1);
                    let tile_bounds = tile.2;
                    //let tile_index = tile.1 * n_tiles.x + tile.0;
                    let tile_sampler = tile.3.clone();
                    let mut tile_sampler = tile_sampler.write().unwrap();
                    // Follow camera paths for _tile_ in image for SPPM
                    for y in tile_bounds.min.y..tile_bounds.max.y {
                        for x in tile_bounds.min.x..tile_bounds.max.x {
                            let p_pixel = Point2i::new(x, y);
                            tile_sampler.start_pixel(&p_pixel);
                            tile_sampler.set_sample_number(iter);

                            // Generate camera ray for pixel for SPPM
                            let camera_sample = tile_sampler.get_camera_sample(&p_pixel);
                            if let Some((beta, mut ray)) =
                                camera.generate_ray_differential(&camera_sample)
                            {
                                if beta <= 0.0 {
                                    continue;
                                }
                                ray.scale_differentials(inv_sqrt_spp);

                                let mut beta = Spectrum::from(beta);

                                // Follow camera ray path until a visible point is created

                                // Get _SPPMPixel_ for _pPixel_
                                let p_pixel_o = p_pixel - pixel_bounds.min;
                                let pixel_offset = p_pixel_o.x
                                    + p_pixel_o.y * (pixel_bounds.max.x - pixel_bounds.min.x);
                                let pixel_offset = pixel_offset as usize;
                                let pixel = unsafe { pixels.get_mut(pixel_offset) };
                                // Reset visible point at the start of each iteration for this pixel.
                                pixel.vp.beta = Spectrum::zero();
                                pixel.vp.bsdf = None;
                                let mut specular_bounce = false;
                                let mut pixel_ld = Spectrum::zero();
                                let mut visible_point: Option<SPPMVisiblePoint> = None;
                                let mut depth = 0;
                                while depth < max_depth {
                                    if let Some(mut isect) = scene.intersect(&ray.ray) {
                                        // Process SPPM camera ray intersection

                                        // Compute BSDF at SPPM camera ray intersection
                                        {
                                            isect.compute_scattering_functions(
                                                &ray,
                                                &mut arena,
                                                TransportMode::Radiance,
                                                true,
                                            );
                                        }
                                        //isect.ComputeScatteringFunctions(ray, arena, true);
                                        if let Some(bsdf) = isect.bsdf.as_ref() {
                                            // Accumulate direct illumination at SPPM camera ray
                                            // intersection
                                            let wo = -ray.ray.d;
                                            let mut ld = beta
                                                * uniform_sample_one_light_surface(
                                                    &isect,
                                                    scene,
                                                    &mut arena,
                                                    tile_sampler.deref_mut() as &mut dyn Sampler,
                                                    false,
                                                    None,
                                                );
                                            if depth == 0 || specular_bounce {
                                                ld += beta * isect.le(&wo);
                                            }
                                            pixel_ld += ld;
                                            // Possibly create visible point and end camera path

                                            let (is_diffuse, is_glossy) = {
                                                let bsdf = bsdf.as_ref();
                                                let b1 = bsdf.num_components(
                                                    BSDF_DIFFUSE
                                                        | BSDF_REFLECTION
                                                        | BSDF_TRANSMISSION,
                                                ) > 0;
                                                let b2 = bsdf.num_components(
                                                    BSDF_GLOSSY
                                                        | BSDF_REFLECTION
                                                        | BSDF_TRANSMISSION,
                                                ) > 0;
                                                (b1, b2)
                                            };
                                            if is_diffuse || (is_glossy && depth == max_depth - 1)
                                            {
                                                visible_point = Some(SPPMVisiblePoint {
                                                    p: isect.p,
                                                    wo: wo,
                                                    bsdf: Some(bsdf.clone()),
                                                    beta: beta.clone(),
                                                });
                                                break;
                                            }

                                            // Spawn ray from SPPM camera path vertex
                                            if depth < max_depth - 1 {
                                                let bsdf = bsdf.as_ref();
                                                if let Some((f, wi, pdf, t)) = bsdf.sample_f(
                                                    &wo,
                                                    &tile_sampler.get_2d(),
                                                    BSDF_ALL,
                                                ) {
                                                    if pdf <= 0.0 || f.is_black() {
                                                        break;
                                                    }
                                                    specular_bounce = (t & BSDF_SPECULAR) != 0;
                                                    beta *= f
                                                        * (Vector3f::abs_dot(
                                                            &wi,
                                                            &isect.shading.n,
                                                        ) / pdf);

                                                    let beta_y = beta.y();
                                                    if beta_y < 0.25 {
                                                        let continue_prob = Float::min(1.0, beta_y);
                                                        if tile_sampler.get_1d() > continue_prob {
                                                            break;
                                                        }
                                                        beta /= continue_prob;
                                                    }
                                                    ray = isect.spawn_ray(&wi).into();
                                                } else {
                                                    break;
                                                }
                                            }
                                        } else {
                                            ray = isect.spawn_ray(&ray.ray.d).into();
                                            continue;
                                        }
                                    } else {
                                        for light in scene.lights.iter() {
                                            let light = light.as_ref();
                                            pixel_ld += beta * light.le(&ray);
                                        }
                                        break;
                                    }
                                    depth += 1;
                                }
                                pixel.ld += pixel_ld;
                                if let Some(vp) = visible_point {
                                    pixel.vp = vp;
                                }
                            }
                        }
                    }
                });
            }
            t_camera += phase_t0.elapsed().as_secs_f64();

            let progress_t0 = Instant::now();
            {
                let mut progress = progress.write().unwrap();
                progress.update(1);
            }
            t_progress += progress_t0.elapsed().as_secs_f64();

            // Create grid of all SPPM visible points

            // Compute grid bounds for SPPM visible points
            let grid_t0 = Instant::now();
            let grid_bounds_t0 = Instant::now();
            let (grid_bounds, max_radius) = {
                let p0 = unsafe { pixels.get(0) }.vp.p;
                let mut max_radius = unsafe { pixels.get(0) }.radius;
                let mut min = p0;
                let mut max = p0;
                for i in 0..n_pixels {
                    let pixel = unsafe { pixels.get(i) };
                    if pixel.vp.beta.is_black() {
                        continue;
                    }
                    let radius = pixel.radius;
                    let p = pixel.vp.p;
                    for i in 0..3 {
                        min[i] = Float::min(min[i], p[i] - radius);
                        max[i] = Float::max(max[i], p[i] + radius);
                    }
                    max_radius = Float::max(max_radius, radius);
                }
                (Bounds3f::new(&min, &max), max_radius)
            };
            t_grid_bounds += grid_bounds_t0.elapsed().as_secs_f64();
            // Compute resolution of SPPM grid in each dimension
            let diag = grid_bounds.diagonal();
            let max_diag = max_component(&diag);
            let base_grid_res = Float::floor(max_diag / max_radius);
            let grid_res = [
                Float::max(1.0, base_grid_res * diag[0] / max_diag) as i32,
                Float::max(1.0, base_grid_res * diag[1] / max_diag) as i32,
                Float::max(1.0, base_grid_res * diag[2] / max_diag) as i32,
            ];

            // Allocate grid for SPPM visible points
            let hash_size = n_pixels;
            let grid = Arc::new(RwLock::new(SPPMHashGrid {
                nodes: vec![Vec::new(); hash_size],
            }));
            {
                let pixel_indices: Vec<usize> = (0..n_pixels).collect();
                pixel_indices.par_iter().for_each(|pixel_index| {
                    let pixel = unsafe { pixels.get(*pixel_index) };
                    if !pixel.vp.beta.is_black() {
                        // Add pixel's visible point to applicable grid cells
                        let p = pixel.vp.p;
                        let radius = pixel.radius;
                        let (_, pmin) = to_grid(
                            &(p - Vector3f::new(radius, radius, radius)),
                            &grid_bounds,
                            &grid_res,
                        );
                        let (_, pmax) = to_grid(
                            &(p + Vector3f::new(radius, radius, radius)),
                            &grid_bounds,
                            &grid_res,
                        );
                        {
                            let mut grid = grid.write().unwrap();
                            let grid = &mut grid.nodes;
                            for z in pmin.z..=pmax.z {
                                for y in pmin.y..=pmax.y {
                                    for x in pmin.x..=pmax.x {
                                        // Add visible point to grid cell $(x, y, z)$
                                        let h = hash(&Point3i::new(x, y, z), hash_size);
                                        grid[h].push(*pixel_index);
                                    }
                                }
                            }

                            // Update statistics on grid cell sizes
                            {
                                let grid_cells_per_visible_point = (1 + pmax.x - pmin.x)
                                    * (1 + pmax.y - pmin.y)
                                    * (1 + pmax.z - pmin.z);
                                GRID_CELLS_PER_VISIBLE_POINT
                                    .with(|c| c.add(grid_cells_per_visible_point as u64));
                            }
                        }
                    }
                });
            }
            t_grid += grid_t0.elapsed().as_secs_f64();

            // Trace photons and accumulate contributions
            let phase_t0 = Instant::now();
            let photon_accum = {
                let photons_per_iteration = self.photons_per_iteration as usize;
                let photon_accum = (0..photons_per_iteration)
                    .into_par_iter()
                    .fold(
                        || {
                            (
                                MemoryArena::new(),
                                HashMap::<usize, (Spectrum, i64)>::new(),
                            )
                        },
                        |(mut arena, mut local), photon_index| {
                        let light_distr = light_distr.clone();
                        // Follow photon path for _photonIndex_
                        let halton_index =
                            iter as u64 * photons_per_iteration as u64 + photon_index as u64;

                        let mut halton_dim = 0;

                        // Choose light to shoot photon from
                        let light_sample = radical_inverse(halton_dim, halton_index);
                        let (light_num, light_pdf, _remapped) =
                            light_distr.sample_discrete(light_sample);
                        let light = scene.lights[light_num].clone();

                        // Compute sample values for photon from light
                        let u_light0 = Point2f::new(
                            radical_inverse(halton_dim, halton_index),
                            radical_inverse(halton_dim + 1, halton_index),
                        );
                        let u_light1 = Point2f::new(
                            radical_inverse(halton_dim + 2, halton_index),
                            radical_inverse(halton_dim + 3, halton_index),
                        );
                        let camera = camera.as_ref();
                        let shutter = camera.get_shutter();
                        let u_light_time = lerp(
                            radical_inverse(halton_dim + 4, halton_index),
                            shutter.0,
                            shutter.1,
                        );
                        halton_dim += 5;

                        // Generate _photonRay_ from light source and initialize _beta_
                        let light = light.as_ref();
                        if let Some((le, photon_ray, n_light, pdf_pos, pdf_dir)) =
                            light.sample_le(&u_light0, &u_light1, u_light_time)
                        {
                            if le.is_black() || pdf_pos == 0.0 || pdf_dir == 0.0 {
                                return (arena, local);
                            }
                            let mut beta = (le * n_light.abs_dot(&photon_ray.d))
                                / (light_pdf * pdf_pos * pdf_dir);
                            if beta.is_black() {
                                return (arena, local);
                            }

                            let mut photon_ray: RayDifferential = photon_ray.into();

                            // Follow photon path through scene and record intersections
                            let mut depth = 0;
                            while depth < max_depth {
                                if let Some(mut isect) = scene.intersect(&photon_ray.ray) {
                                    if depth > 0 {
                                        // Add photon contribution to nearby visible points
                                        let (in_bounds, photon_grid_index) =
                                            to_grid(&isect.p, &grid_bounds, &grid_res);
                                        if in_bounds {
                                            let h = hash(&photon_grid_index, hash_size);
                                            // Add photon contribution to visible points in
                                            // _grid[h]_
                                            let grid = grid.read().unwrap();
                                            for pixel_index in grid.nodes[h].iter() {
                                                let pixel = unsafe { pixels.get(*pixel_index) };
                                                let wi = -photon_ray.ray.d;
                                                let radius = pixel.radius;
                                                if Vector3f::distance_squared(&pixel.vp.p, &isect.p)
                                                    <= (radius * radius)
                                                {
                                                    if let Some(bsdf) = pixel.vp.bsdf.as_ref() {
                                                        let phi =
                                                            beta * bsdf.f(&pixel.vp.wo, &wi, BSDF_ALL);
                                                        let e = local
                                                            .entry(*pixel_index)
                                                            .or_insert((Spectrum::zero(), 0));
                                                        e.0 += phi;
                                                        e.1 += 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    // Sample new photon ray direction

                                    // Compute BSDF at photon intersection point
                                    isect.compute_scattering_functions(
                                        &photon_ray,
                                        &mut arena,
                                        TransportMode::Importance,
                                        true,
                                    );
                                    if isect.bsdf.is_none() {
                                        // Skip over medium boundaries for photon tracing
                                        photon_ray = isect.spawn_ray(&photon_ray.ray.d).into();
                                        continue;
                                    }

                                    if let Some(photon_bsdf) = isect.bsdf.as_ref() {
                                        // Sample BSDF for photon scattering
                                        let wo = -photon_ray.ray.d;
                                        // Generate _bsdfSample_ for outgoing photon sample
                                        let bsdf_sample = Point2f::new(
                                            radical_inverse(halton_dim, halton_index),
                                            radical_inverse(halton_dim + 1, halton_index),
                                        );
                                        halton_dim += 2;

                                        let photon_bsdf = photon_bsdf.as_ref();
                                        if let Some((fr, wi, pdf, _flags)) =
                                            photon_bsdf.sample_f(&wo, &bsdf_sample, BSDF_ALL)
                                        {
                                            if fr.is_black() || pdf == 0.0 {
                                                break;
                                            }
                                            //
                                            let bnew = beta
                                                * fr
                                                * (Vector3f::abs_dot(&wi, &isect.shading.n) / pdf);

                                            // Possibly terminate photon path with Russian roulette
                                            let q = Float::max(0.0, 1.0 - bnew.y() / beta.y());
                                            let t = radical_inverse(halton_dim, halton_index);
                                            halton_dim += 1;
                                            if t < q {
                                                break;
                                            }
                                            beta = bnew / (1.0 - q);
                                            photon_ray = isect.spawn_ray(&wi).into();
                                        } else {
                                            break;
                                        }
                                    } else {
                                        break;
                                    }
                                } else {
                                    break;
                                }
                                depth += 1;
                            }
                        } else {
                            return (arena, local);
                        }
                        arena.reset();
                        (arena, local)
                    },
                    )
                    .map(|(_, local)| local)
                    .reduce(HashMap::new, |mut a, b| {
                        for (pixel_index, (phi, m)) in b {
                            let e = a.entry(pixel_index).or_insert((Spectrum::zero(), 0));
                            e.0 += phi;
                            e.1 += m;
                        }
                        a
                    });
                PHOTON_PATHS.with(|c| c.add(photons_per_iteration as u64));
                photon_accum
            };
            t_photon += phase_t0.elapsed().as_secs_f64();
            let drop_grid_t0 = Instant::now();
            drop(grid);
            t_drop_grid += drop_grid_t0.elapsed().as_secs_f64();

            // Update pixel values from this pass's photons
            let phase_t0 = Instant::now();
            {
                for (pixel_index, (phi, m)) in photon_accum {
                    // Update pixel photon count, search radius, and tau only for
                    // pixels that actually received photons in this iteration.
                    let p = unsafe { pixels.get_mut(pixel_index) };
                    let gamma = 2.0 / 3.0;
                    let nnew = p.n + gamma * m as Float;
                    let rnew = p.radius * Float::sqrt(nnew / (p.n + m as Float));
                    p.tau = (p.tau + p.vp.beta * phi) * ((rnew * rnew) / (p.radius * p.radius));
                    p.n = nnew;
                    p.radius = rnew;
                }
            }
            t_update += phase_t0.elapsed().as_secs_f64();

            // Periodically store SPPM image in film and write image
            let phase_t0 = Instant::now();
            {
                let photons_per_iteration = self.photons_per_iteration as u32;
                let write_frequency = self.write_frequency as u32;
                let write_image = /*(iter + 1) == n_iterations || */((iter + 1) % write_frequency) == 0;
                {
                    let x0 = pixel_bounds.min.x;
                    let x1 = pixel_bounds.max.x;
                    let y0 = pixel_bounds.min.y;
                    let y1 = pixel_bounds.max.y;
                    let np = (iter + 1) * photons_per_iteration;

                    let area = pixel_bounds.area() as usize;
                    assert_eq!(area, (x1 - x0) as usize * (y1 - y0) as usize);
                    let mut image = vec![Spectrum::default(); area];

                    let width = (x1 - x0) as usize;
                    let mut offset = 0;
                    for y in y0..y1 {
                        for x in x0..x1 {
                            // Compute radiance _L_ for SPPM pixel _pixel_
                            let yy = (y - pixel_bounds.min.y) as usize;
                            let pixel_index = yy * width + (x - x0) as usize;
                            let pixel = unsafe { pixels.get(pixel_index) };
                            let mut l = pixel.ld / ((iter + 1) as Float);
                            l += pixel.tau / (np as Float * PI * pixel.radius * pixel.radius);
                            image[offset] = l;
                            offset += 1;
                        }
                    }
                    {
                        let film = self.get_film();
                        let mut film = film.write().unwrap();
                        film.set_image(&image);
                        if write_image {
                            film.write_image();
                        }
                        film.update_display(&pixel_bounds);
                    }
                }
            }
            t_image += phase_t0.elapsed().as_secs_f64();

            let progress_t0 = Instant::now();
            {
                let mut progress = progress.write().unwrap();
                progress.update(1);
            }
            t_progress += progress_t0.elapsed().as_secs_f64();
        }

        let progress_t0 = Instant::now();
        {
            let mut progress = progress.write().unwrap();
            progress.done();
        }
        t_progress += progress_t0.elapsed().as_secs_f64();

        let finalize_t0 = Instant::now();
        {
            let film = self.get_film();
            let mut film = film.write().unwrap();
            film.render_end();
            film.write_image();
        }
        t_finalize += finalize_t0.elapsed().as_secs_f64();
        if timing_enabled {
            let total = render_t0.elapsed().as_secs_f64();
            let known = t_camera
                + t_grid
                + t_photon
                + t_update
                + t_image
                + t_grid_bounds
                + t_progress
                + t_finalize
                + t_setup_film
                + t_setup_pixels
                + t_setup_light_distr
                + t_setup_tiles
                + t_drop_grid;
            let other = (total - known).max(0.0);
            let denom = if total > 0.0 { total } else { 1.0 };
            eprintln!(
                "[SPPM timing] total={:.3}s setup_film={:.3}s ({:.1}%) setup_pixels={:.3}s ({:.1}%) setup_light_distr={:.3}s ({:.1}%) setup_tiles={:.3}s ({:.1}%) camera={:.3}s ({:.1}%) grid={:.3}s ({:.1}%) grid_bounds={:.3}s ({:.1}%) photon={:.3}s ({:.1}%) drop_grid={:.3}s ({:.1}%) update={:.3}s ({:.1}%) image={:.3}s ({:.1}%) progress={:.3}s ({:.1}%) finalize={:.3}s ({:.1}%) other={:.3}s ({:.1}%)",
                total,
                t_setup_film,
                100.0 * t_setup_film / denom,
                t_setup_pixels,
                100.0 * t_setup_pixels / denom,
                t_setup_light_distr,
                100.0 * t_setup_light_distr / denom,
                t_setup_tiles,
                100.0 * t_setup_tiles / denom,
                t_camera,
                100.0 * t_camera / denom,
                t_grid,
                100.0 * t_grid / denom,
                t_grid_bounds,
                100.0 * t_grid_bounds / denom,
                t_photon,
                100.0 * t_photon / denom,
                t_drop_grid,
                100.0 * t_drop_grid / denom,
                t_update,
                100.0 * t_update / denom,
                t_image,
                100.0 * t_image / denom,
                t_progress,
                100.0 * t_progress / denom,
                t_finalize,
                100.0 * t_finalize / denom,
                other,
                100.0 * other / denom
            );
        }
    }

    fn get_camera(&self) -> Arc<dyn Camera> {
        Arc::clone(&self.camera)
    }
}

pub fn create_sppm_integrator(
    params: &ParamSet,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    let mut n_iterations =
        params.find_one_int("iterations", params.find_one_int("numiterations", 64));
    let max_depth = params.find_one_int("maxdepth", 5);
    let photons_per_iteration = params.find_one_int("photonsperiteration", -1);
    let write_frequency = params.find_one_int("imagewritefrequency", 1 << 31);
    let initial_search_radius = params.find_one_float("radius", 1.0);
    {
        let options = PbrtOptions::get();
        if options.quick_render {
            n_iterations = (n_iterations / 16).max(1);
        }
    }
    //if (PbrtOptions.quickRender) nIterations = std::max(1, nIterations / 16);
    return Ok(Arc::new(RwLock::new(SPPMIntegrator::new(
        Arc::clone(camera),
        initial_search_radius,
        n_iterations,
        max_depth,
        photons_per_iteration,
        write_frequency,
    ))));
}
