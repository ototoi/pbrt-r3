use super::subpath::*;
use super::vertex::*;
use crate::core::pbrt::*;
use log::*;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use std::ops::DerefMut;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::RwLock;

pub struct BDPTIntegrator {
    sampler: Arc<RwLock<dyn Sampler>>,
    camera: Arc<dyn Camera>,
    max_depth: i32,
    visualize_strategies: bool,
    visualize_weights: bool,
    pixel_bounds: Bounds2i,
    light_sample_strategy: String,
}

impl BDPTIntegrator {
    pub fn new(
        sampler: &Arc<RwLock<dyn Sampler>>,
        camera: &Arc<dyn Camera>,
        max_depth: i32,
        visualize_strategies: bool,
        visualize_weights: bool,
        pixel_bounds: &Bounds2i,
        light_sample_sterategy: &str,
    ) -> Self {
        BDPTIntegrator {
            sampler: sampler.clone(),
            camera: camera.clone(),
            max_depth,
            visualize_strategies: visualize_strategies,
            visualize_weights: visualize_weights,
            pixel_bounds: *pixel_bounds,
            light_sample_strategy: light_sample_sterategy.to_string(),
        }
    }
}

impl Integrator for BDPTIntegrator {
    fn render(&mut self, scene: &Scene) {
        let light_distribution =
            create_light_sample_distribution(&self.light_sample_strategy, scene);
        if let Err(e) = light_distribution {
            error!("Unable to create light sample distribution: {}", e);
            return;
        }
        let light_distribution = light_distribution.unwrap();
        // Compute a reverse mapping from light pointers to offsets into the
        // scene lights vector (and, equivalently, offsets into
        // lightDistr). Added after book text was finalized; this is critical
        // to reasonable performance with 100s+ of light sources.
        let mut light_to_index = LightIndexMap::default();
        for (i, light) in scene.lights.iter().enumerate() {
            let light = light.as_ref();
            let light_ptr = light as *const dyn Light;
            let key = LightKeyType::from(light_ptr);
            light_to_index.insert(key, i);
        }

        // Partition the image into tiles
        let mut tile_indices = Vec::new();
        {
            //let film = self.camera.as_ref().get_film();
            //let film = film.read().unwrap();
            let sampler = self.sampler.clone();
            //let sample_bounds = film.get_sample_bounds();
            let sample_bounds = self.pixel_bounds;
            let sample_extent = sample_bounds.diagonal();
            const TILE_SIZE: i32 = 16;
            let n_tiles = Point2i::from((
                (sample_extent.x + TILE_SIZE - 1) / TILE_SIZE,
                (sample_extent.y + TILE_SIZE - 1) / TILE_SIZE,
            ));
            tile_indices.reserve((n_tiles.x * n_tiles.y) as usize);
            for y in 0..n_tiles.y {
                for x in 0..n_tiles.x {
                    let x0 = sample_bounds.min.x + x * TILE_SIZE;
                    let x1 = i32::min(x0 + TILE_SIZE, sample_bounds.max.x);
                    let y0 = sample_bounds.min.y + y * TILE_SIZE;
                    let y1 = i32::min(y0 + TILE_SIZE, sample_bounds.max.y);
                    let tile_bounds = Bounds2i::from(((x0, y0), (x1, y1)));

                    let seed = (y * n_tiles.x + x) as u32;
                    let s = sampler.read().unwrap().clone_with_seed(seed);
                    let proxy = Arc::new(Mutex::new(ProxySampler::new(&s)));
                    //let camera = self.camera.clone();
                    //let camera = Arc::downgrade(&camera);
                    tile_indices.push((tile_bounds, proxy));
                }
            }
        }

        {
            let camera = self.camera.as_ref();
            let film = camera.get_film();
            let mut film = film.write().unwrap();
            film.render_start();
        }

        let total_tiles = tile_indices.len();
        let reporter = Arc::new(RwLock::new(ProgressReporter::new(total_tiles, "Rendering")));

        let max_depth = self.max_depth as usize;
        // Render and write the output image to disk
        if !scene.lights.is_empty() {
            let camera = self.camera.clone();
            let film: Arc<RwLock<Film>> = camera.as_ref().get_film();
            let proxy_film = Arc::new(Mutex::new(ProxyFilm::new(&film)));
            //let filename = proxy_film.as_ref().lock().unwrap().get_filename();

            //let total = tile_indices.len();

            let visualize_strategies = self.visualize_strategies;
            let visualize_weights = self.visualize_weights;
            let visualize_info = visualize_strategies || visualize_weights;

            let samples_per_pixel = self.sampler.read().unwrap().get_samples_per_pixel();
            {
                tile_indices.par_iter().for_each(|(tile_bounds, sampler)| {
                    //let camera = self.camera.clone();
                    let mut arena = MemoryArena::new();

                    let x0 = tile_bounds.min.x;
                    let x1 = tile_bounds.max.x;
                    let y0 = tile_bounds.min.y;
                    let y1 = tile_bounds.max.y;

                    let mut tile_sampler = sampler.lock().unwrap();

                    let mut film_tile = proxy_film
                        .as_ref()
                        .lock()
                        .unwrap()
                        .get_film_tile(tile_bounds);

                    for yy in y0..y1 {
                        for xx in x0..x1 {
                            let p_pixel = Point2i::from((xx, yy));
                            tile_sampler.start_pixel(&p_pixel);

                            loop {
                                //println!("p_pixel: {:?}", p_pixel);
                                // Generate a single sample using BDPT
                                let p_film = Point2f::new(p_pixel.x as Float, p_pixel.y as Float)
                                    + tile_sampler.get_2d();

                                // Trace the camera subpath
                                let mut camera_vertices = Vec::with_capacity(max_depth + 2);
                                let mut light_vertices = Vec::with_capacity(max_depth + 1);

                                let n_camera = generate_camera_subpath(
                                    scene,
                                    tile_sampler.deref_mut(),
                                    &mut arena,
                                    max_depth + 2,
                                    &camera,
                                    &p_film,
                                    &mut camera_vertices,
                                );
                                // Get a distribution for sampling the light at the
                                // start of the light subpath. Because the light path
                                // follows multiple bounces, basing the sampling
                                // distribution on any of the vertices of the camera
                                // path is unlikely to be a good strategy. We use the
                                // PowerLightDistribution by default here, which
                                // doesn't use the point passed to it.
                                if n_camera <= 0 {
                                    break;
                                }

                                let camera_vertex = camera_vertices[0].as_ref().read().unwrap();

                                let p = camera_vertex.get_p();
                                let light_distribution = light_distribution.as_ref();
                                let light_distr = light_distribution.lookup(&p);
                                // Now trace the light subpath

                                let time = camera_vertex.get_time();
                                let n_light = generate_light_subpath(
                                    scene,
                                    tile_sampler.deref_mut(),
                                    &mut arena,
                                    max_depth + 1,
                                    time,
                                    &light_distr,
                                    &light_to_index,
                                    &mut light_vertices,
                                );

                                let mut l = Spectrum::zero();
                                let n_camera = n_camera as i32;
                                let n_light = n_light as i32;
                                for t in 1..=n_camera {
                                    for s in 0..=n_light {
                                        let depth = (s + t) - 2;

                                        if (s == 1 && t == 1)
                                            || depth < 0
                                            || depth > max_depth as i32
                                        {
                                            continue;
                                        }
                                        // Execute the $(s, t)$ connection strategy and
                                        // update _L_
                                        let (l_path, mis_weight, p_film_new) = connect_bdpt(
                                            scene,
                                            &light_vertices[0..],
                                            &camera_vertices[0..],
                                            s,
                                            t,
                                            &light_distr,
                                            &light_to_index,
                                            &camera,
                                            tile_sampler.deref_mut(),
                                            &p_film,
                                        );
                                        if visualize_info {
                                            let value = if visualize_strategies {
                                                if mis_weight <= 0.0 {
                                                    Spectrum::zero()
                                                } else {
                                                    l_path * (1.0 / mis_weight)
                                                }
                                            } else {
                                                l_path
                                            };
                                            let mut film = proxy_film.as_ref().lock().unwrap();
                                            film.add_splat(&p_film_new, &value);
                                        }
                                        if t != 1 {
                                            l += l_path;
                                        } else {
                                            if !l_path.is_black() && !visualize_info {
                                                //println!("l_path: {:?}", l_path);
                                                let mut film = proxy_film.as_ref().lock().unwrap();
                                                film.add_splat(&p_film_new, &l_path);
                                            }
                                        }
                                    }
                                }

                                if !visualize_info {
                                    film_tile.add_sample(&p_film, &l, 1.0);
                                }
                                arena.reset();

                                if !tile_sampler.start_next_sample() {
                                    break;
                                }
                            }
                        }
                    }

                    {
                        let mut film = proxy_film.as_ref().lock().unwrap();
                        film.merge_splats(1.0 / samples_per_pixel as Float);
                        film.merge_film_tile(&film_tile);
                    }
                    {
                        let mut reporter = reporter.write().unwrap();
                        reporter.update(1);
                    }
                });
            }
        }

        {
            let camera = self.camera.as_ref();
            let film = camera.get_film();
            let mut film = film.as_ref().write().unwrap();
            film.render_end();
            film.write_image();
        }

        {
            let mut reporter = reporter.write().unwrap();
            reporter.done();
        }
    }

    fn get_camera(&self) -> Arc<dyn Camera> {
        return self.camera.clone();
    }
}

unsafe impl Sync for BDPTIntegrator {}

pub fn create_bdpt_integrator(
    params: &ParamSet,
    sampler: &Arc<RwLock<dyn Sampler>>,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    let max_depth = params.find_one_int("maxdepth", 5);
    let visualize_strategies = params.find_one_bool("visualizestrategies", false);
    let visualize_weights = params.find_one_bool("visualizeweights", false);
    //let p =
    let pixel_bounds = camera.get_film().read().unwrap().get_sample_bounds();
    //todo
    let light_sample_strategy = params.find_one_string("lightsamplestrategy", "power");
    return Ok(Arc::new(RwLock::new(BDPTIntegrator::new(
        sampler,
        camera,
        max_depth,
        visualize_strategies,
        visualize_weights,
        &pixel_bounds,
        &light_sample_strategy,
    ))));
}
