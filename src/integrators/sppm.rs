use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;

use crate::core::pbrt::*;
use crate::samplers::HaltonSampler;
use std::ops::DerefMut;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::RwLock;

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
    phi: Spectrum,
    m: i64,
    n: Float,
    tau: Spectrum,
}

#[derive(Clone, Default)]
struct SPPMTile {
    pixels: Vec<Arc<RwLock<SPPMPixel>>>,
}
impl SPPMTile {
    pub fn push(&mut self, pixel: Arc<RwLock<SPPMPixel>>) {
        self.pixels.push(pixel);
    }
}
unsafe impl Sync for SPPMTile {}
unsafe impl Send for SPPMTile {}

#[derive(Clone)]
struct SPPMPixelListNode {
    pixel: Arc<RwLock<SPPMPixel>>,
    next: Option<Arc<SPPMPixelListNode>>,
}

#[derive(Clone)]
struct SPPMHashGrid {
    nodes: Vec<Option<Arc<SPPMPixelListNode>>>,
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

        {
            let film = self.get_film();
            let mut film = film.write().unwrap();
            film.render_start();
        }

        let n_pixels = pixel_bounds.area() as usize;
        let initial_search_radius = self.initial_search_radius;
        let pixels = Arc::new(RwLock::new(SPPMTile::default()));

        PIXEL_MEMORY_BYTES.with(|c| c.add(n_pixels * std::mem::size_of::<SPPMPixel>()));

        {
            let mut pixels = pixels.write().unwrap();
            for _ in 0..n_pixels {
                let pixel = Arc::new(RwLock::new(SPPMPixel {
                    radius: initial_search_radius,
                    ld: Spectrum::default(),
                    vp: SPPMVisiblePoint {
                        p: Point3f::default(),
                        wo: Vector3f::default(),
                        bsdf: None,
                        beta: Spectrum::zero(),
                    },
                    phi: Spectrum::zero(),
                    m: 0,
                    n: 0.0,
                    tau: Spectrum::zero(),
                }));
                pixels.push(pixel);
            }
        }

        // Compute _lightDistr_ for sampling lights proportional to power
        let light_distr = compute_light_power_distribution(scene);

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
                let sampler = Arc::new(Mutex::new(ProxySampler::new(&sampler)));
                //
                tiles.push((x, y, tile_bounds, sampler));
                //samplers.push(sampler);
            }
        }

        let progress = Arc::new(RwLock::new(ProgressReporter::new(
            2 * n_iterations as usize,
            "Rendering",
        )));

        let camera = self.get_camera();
        //let camera = camera.deref();
        let max_depth = self.max_depth;
        for iter in 0..n_iterations {
            // Generate SPPM visible points
            {
                tiles.par_iter().for_each(|tile| {
                    let mut arena = MemoryArena::new();
                    //let pixel_bounds = tile.5;
                    let pixels = pixels.clone();
                    let pixels = pixels.read().unwrap();
                    let camera = camera.as_ref();
                    //let tile = Point2i::new(tile.0, tile.1);
                    let tile_bounds = tile.2;
                    //let tile_index = tile.1 * n_tiles.x + tile.0;
                    let tile_sampler = tile.3.clone();
                    let mut tile_sampler = tile_sampler.lock().unwrap();
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
                                let pixel = pixels.pixels[pixel_offset as usize].clone();
                                let mut specular_bounce = false;
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
                                            let tisect = Interaction::from(&isect);
                                            // Accumulate direct illumination at SPPM camera ray
                                            // intersection
                                            let wo = -ray.ray.d;
                                            if depth == 0 || specular_bounce {
                                                let mut pixel = pixel.write().unwrap();
                                                pixel.ld += beta * isect.le(&wo);
                                            }
                                            {
                                                let mut pixel = pixel.write().unwrap();
                                                pixel.ld += beta
                                                    * uniform_sample_one_light(
                                                        &tisect,
                                                        scene,
                                                        &mut arena,
                                                        tile_sampler.deref_mut()
                                                            as &mut dyn Sampler,
                                                        false,
                                                        None,
                                                    );
                                            }
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
                                            if is_diffuse || (is_glossy && depth == max_depth - 1) {
                                                let mut pixel = pixel.write().unwrap();
                                                pixel.vp = SPPMVisiblePoint {
                                                    p: isect.p,
                                                    wo: wo,
                                                    bsdf: Some(bsdf.clone()),
                                                    beta: beta.clone(),
                                                };
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
                                                        beta *= 1.0 / continue_prob;
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
                                        let mut pixel = pixel.write().unwrap();
                                        for light in scene.lights.iter() {
                                            let light = light.as_ref();
                                            pixel.ld += beta * light.le(&ray);
                                        }
                                        break;
                                    }
                                    depth += 1;
                                }
                            }
                        }
                    }
                });
            }

            {
                let mut progress = progress.write().unwrap();
                progress.update(1);
            }

            // Create grid of all SPPM visible points

            // Compute grid bounds for SPPM visible points
            let (grid_bounds, max_radius) = {
                let pixels = pixels.read().unwrap();
                let p0 = pixels.pixels[0].read().unwrap().vp.p;
                let mut max_radius = pixels.pixels[0].read().unwrap().radius;
                let mut min = p0;
                let mut max = p0;
                for i in 0..n_pixels {
                    let pixel = &pixels.pixels[i];
                    let pixel = pixel.read().unwrap();
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
                nodes: vec![None; hash_size],
            }));
            {
                let pixel_indices: Vec<usize> = (0..n_pixels).collect();
                pixel_indices.par_iter().for_each(|pixel_index| {
                    let pixels = pixels.read().unwrap();
                    let pixel_ref = &pixels.pixels[*pixel_index];
                    let pixel = pixel_ref.read().unwrap();
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
                                        {
                                            let node = Arc::new(SPPMPixelListNode {
                                                pixel: pixel_ref.clone(),
                                                next: grid[h].clone(),
                                            });
                                            grid[h] = Some(node);
                                        }
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

            // Trace photons and accumulate contributions
            {
                let arena = Arc::new(RwLock::new(MemoryArena::new()));
                let photons_per_iteration = self.photons_per_iteration as usize;
                (0..photons_per_iteration)
                    .into_par_iter()
                    .for_each(|photon_index| {
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
                                return;
                            }
                            let mut beta = (le * n_light.abs_dot(&photon_ray.d))
                                / (light_pdf * pdf_pos * pdf_dir);
                            if beta.is_black() {
                                return;
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
                                            let mut node = {
                                                let grid = grid.read().unwrap();
                                                grid.nodes[h].clone()
                                            };
                                            while let Some(n) = node {
                                                let pixel = n.pixel.as_ref();
                                                {
                                                    let pixel = pixel.read().unwrap();
                                                    let radius = pixel.radius;
                                                    if Vector3f::distance_squared(
                                                        &pixel.vp.p,
                                                        &isect.p,
                                                    ) > (radius * radius)
                                                    {
                                                        node = n.next.clone();
                                                        continue;
                                                    }
                                                }
                                                {
                                                    let pixel = pixel.write().unwrap();
                                                    assert!(pixel.vp.bsdf.is_some());
                                                }
                                                // Update _pixel_ $\Phi$ and $M$ for nearby
                                                // photon
                                                let phi = {
                                                    let pixel = pixel.read().unwrap();
                                                    assert!(pixel.vp.bsdf.is_some());
                                                    let bsdf = pixel.vp.bsdf.as_ref().unwrap();
                                                    let wi = -photon_ray.ray.d;
                                                    beta * bsdf.f(&pixel.vp.wo, &wi, BSDF_ALL)
                                                };
                                                {
                                                    let mut pixel = pixel.write().unwrap();
                                                    pixel.phi += phi;
                                                    pixel.m += 1;
                                                }
                                                node = n.next.clone();
                                            }
                                        }
                                    }
                                    // Sample new photon ray direction

                                    // Compute BSDF at photon intersection point
                                    {
                                        let arena = arena.clone();
                                        let mut arena = arena.write().unwrap();
                                        isect.compute_scattering_functions(
                                            &photon_ray,
                                            &mut arena,
                                            TransportMode::Importance,
                                            true,
                                        );
                                    }
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
                            return;
                        }
                    });
                PHOTON_PATHS.with(|c| c.add(photons_per_iteration as u64));
            }

            // Update pixel values from this pass's photons
            {
                (0..n_pixels).into_par_iter().for_each(|i| {
                    let pixels = pixels.clone();
                    let pixels = pixels.read().unwrap();
                    let pixel = &pixels.pixels[i];
                    let mut p = pixel.write().unwrap();
                    if p.m > 0 {
                        // Update pixel photon count, search radius, and $\tau$ from
                        // photons
                        let gamma = 2.0 / 3.0;
                        let nnew = p.n as Float + gamma * p.m as Float;
                        let rnew = p.radius * Float::sqrt(nnew / (p.n as Float + p.m as Float));
                        let phi = p.phi;

                        p.tau = (p.tau + p.vp.beta * phi) * ((rnew * rnew) / (p.radius * p.radius));
                        p.n = nnew;
                        p.radius = rnew;
                        p.m = 0;
                        p.phi = Spectrum::zero();
                    }
                    // Reset _VisiblePoint_ in _pixel_
                    p.vp.beta = Spectrum::zero();
                    p.vp.bsdf = None;
                });
            }

            // Periodically store SPPM image in film and write image
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

                    let pixels = pixels.clone();
                    let pixels = pixels.read().unwrap();
                    let width = (x1 - x0) as usize;
                    let mut offset = 0;
                    for y in y0..y1 {
                        for x in x0..x1 {
                            // Compute radiance _L_ for SPPM pixel _pixel_
                            let yy = (y - pixel_bounds.min.y) as usize;
                            let pixel_index = yy * width + (x - x0) as usize;
                            let pixel = &pixels.pixels[pixel_index];
                            let pixel = pixel.read().unwrap();
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

            {
                let mut progress = progress.write().unwrap();
                progress.update(1);
            }
        }

        {
            let mut progress = progress.write().unwrap();
            progress.done();
        }

        {
            let film = self.get_film();
            let mut film = film.write().unwrap();
            film.render_end();
            film.write_image();
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
