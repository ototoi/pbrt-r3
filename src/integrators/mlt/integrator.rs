use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;

use crate::core::pbrt::*;
use crate::integrators::bdpt::subpath::*;
use crate::integrators::bdpt::vertex::*;

use std::sync::atomic::*;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::RwLock;
use std::time::Instant;

use super::sampler::MLTSampler;

// MLTSampler Constants
const CAMERA_STREAM_INDEX: u64 = 0;
const LIGHT_STREAM_INDEX: u64 = 1;
const CONNECTION_STREAM_INDEX: u64 = 2;
const N_SAMPLE_STREAMS: u64 = 3;
const PROGRESS_FREQUENCY: u32 = 32768;
const UPDATE_DISPLAY_INTERVAL: u128 = 1000;

pub struct MLTIntegrator {
    // MLTIntegrator Private Data
    camera: Arc<dyn Camera>,
    max_depth: u32,
    n_bootstrap: u32,
    n_chains: u32,
    mutations_per_pixel: u32,
    sigma: Float,
    large_step_probability: Float,
    sample_bounds: Bounds2i,
}

impl MLTIntegrator {
    pub fn new(
        camera: &Arc<dyn Camera>,
        max_depth: u32,
        n_bootstrap: u32,
        n_chains: u32,
        mutations_per_pixel: u32,
        sigma: Float,
        large_step_probability: Float,
    ) -> Self {
        let film = camera.get_film();
        let sample_bounds = film.read().unwrap().get_sample_bounds();
        MLTIntegrator {
            camera: camera.clone(),
            max_depth,
            n_bootstrap,
            n_chains,
            mutations_per_pixel,
            sigma,
            large_step_probability,
            sample_bounds,
        }
    }

    fn l(
        &self,
        scene: &Scene,
        arena: &mut MemoryArena,
        light_distr: &Distribution1D,
        light_to_index: &LightIndexMap,
        sampler: &mut MLTSampler,
        depth: u32,
    ) -> (Spectrum, Point2f) {
        let camera = self.camera.clone();
        sampler.start_stream(CAMERA_STREAM_INDEX);
        // Determine the number of available strategies for connecting the
        let (s, t, n_strategies) = if depth == 0 {
            (0, 2, 1)
        } else {
            let n_strategies = depth + 2;
            let s = u32::min(
                (sampler.get_1d() * n_strategies as Float) as u32,
                n_strategies as u32 - 1,
            );
            let t = n_strategies as u32 - s;
            (s, t, n_strategies)
        };

        let sample_bounds = self.sample_bounds;
        let sample_bounds = Bounds2f::new(
            &Point2f::new(sample_bounds.min.x as Float, sample_bounds.min.y as Float),
            &Point2f::new(sample_bounds.max.x as Float, sample_bounds.max.y as Float),
        );
        let p_raster = sample_bounds.lerp(&sampler.get_2d());
        // Generate a camera subpath with exactly _t_ vertices
        let mut camera_vertices = Vec::with_capacity(t as usize);
        if generate_camera_subpath(
            scene,
            sampler,
            arena,
            t as usize,
            &camera,
            &p_raster,
            &mut camera_vertices,
        ) != t as usize
        {
            return (Spectrum::zero(), p_raster);
        };

        // Generate a light subpath with exactly _s_ vertices
        sampler.start_stream(LIGHT_STREAM_INDEX);
        let time = camera_vertices[0].read().unwrap().get_time();
        let mut light_vertices = Vec::with_capacity(s as usize);
        if generate_light_subpath(
            scene,
            sampler,
            arena,
            s as usize,
            time,
            light_distr,
            light_to_index,
            &mut light_vertices,
        ) != s as usize
        {
            return (Spectrum::zero(), p_raster);
        };
        //let t = camera_vertices.len() as u32;
        //let s = light_vertices.len() as u32;

        // Execute connection strategy and return the radiance estimate
        sampler.start_stream(CONNECTION_STREAM_INDEX);
        if let Some((spec, _, p_raster_new)) = connect_bdpt(
            scene,
            &light_vertices,
            &camera_vertices,
            s as i32,
            t as i32,
            light_distr,
            light_to_index,
            &camera,
            sampler,
            &p_raster,
        ) {
            return (spec * n_strategies as Float, p_raster_new);
        } else {
            return (Spectrum::zero(), p_raster);
        }
    }
}

impl Integrator for MLTIntegrator {
    fn render(&mut self, scene: &Scene) {
        {
            let camera = self.camera.clone();
            let film = camera.get_film();
            let mut film = film.write().unwrap();
            film.render_start();
        }

        let light_distr = compute_light_power_distribution(scene);

        // Compute a reverse mapping from light pointers to offsets into the
        // scene lights vector (and, equivalently, offsets into
        // lightDistr). Added after book text was finalized; this is critical
        // to reasonable performance with 100s+ of light sources.
        let mut light_to_index = LightIndexMap::new();
        for (i, light) in scene.lights.iter().enumerate() {
            let light = light.as_ref();
            let light_ptr = light as *const dyn Light;
            let key = LightKeyType::from(light_ptr);
            light_to_index.insert(key, i);
        }

        // Generate bootstrap samples and compute normalization constant $b$
        let mutations_per_pixel = self.mutations_per_pixel;
        let n_bootstrap = self.n_bootstrap as u32;
        let max_depth = self.max_depth as u32;
        let n_bootstrap_samples = n_bootstrap * (max_depth + 1);

        let bootstrap_weights = Mutex::new(vec![0.0; n_bootstrap_samples as usize]);
        if scene.lights.len() > 0 {
            let reporter = Arc::new(RwLock::new(ProgressReporter::new(
                n_bootstrap as usize,
                "Generating bootstrap paths",
            )));

            let sigma = self.sigma;
            let large_step_probability = self.large_step_probability;
            let n_sample_streams = N_SAMPLE_STREAMS;
            (0..n_bootstrap).into_par_iter().for_each(|i| {
                // Generate _i_th bootstrap sample
                let mut arena = MemoryArena::new();
                for depth in 0..=max_depth {
                    let rng_index = i * (max_depth + 1) + depth;
                    let mut sampler = MLTSampler::new(
                        mutations_per_pixel,
                        rng_index as u64,
                        sigma,
                        large_step_probability,
                        n_sample_streams,
                    );
                    let (l, _p_raster) = self.l(
                        scene,
                        &mut arena,
                        &light_distr,
                        &light_to_index,
                        &mut sampler,
                        depth,
                    );
                    let y = l.y();
                    {
                        let mut bootstrap_weights = bootstrap_weights.lock().unwrap();
                        bootstrap_weights[rng_index as usize] = y;
                    }
                }
                {
                    let mut reporter = reporter.write().unwrap();
                    reporter.update(1);
                }
            });
            {
                let mut reporter = reporter.write().unwrap();
                reporter.done();
            }
        }

        let bootstrap = {
            let bootstrap_weights = bootstrap_weights.lock().unwrap();
            Distribution1D::new(&bootstrap_weights)
        };
        let b = bootstrap.func_int * (max_depth + 1) as Float;
        let splat_scale = b as Float / mutations_per_pixel as Float;
        //println!("b: {}, {}", bootstrap.func_int, max_depth);

        // Run _nChains_ Markov chains in parallel
        if scene.lights.len() > 0 {
            let camera = self.camera.clone();
            let film = camera.get_film();
            let pixel_bounds = film.read().unwrap().cropped_pixel_bounds;
            let n_total_mutations =
                self.mutations_per_pixel * film.read().unwrap().get_sample_bounds().area() as u32;
            let progress_frequency = PROGRESS_FREQUENCY as usize;
            let n_chains = self.n_chains;
            let mutations_per_pixel = self.mutations_per_pixel;
            let sigma = self.sigma;
            let large_step_probability = self.large_step_probability;
            let n_sample_streams = N_SAMPLE_STREAMS;

            //let n_total_work = (n_total_mutations / progress_frequency) as usize;
            let n_chain_mutations_list: Vec<usize> = (0..n_chains)
                .map(|i| {
                    i64::max(
                        1,
                        i64::min(
                            ((i + 1) * n_total_mutations / n_chains) as i64,
                            n_total_mutations as i64,
                        ) - (i * n_total_mutations / n_chains) as i64,
                    ) as usize
                })
                .collect();
            let n_total_iterations = n_chain_mutations_list.iter().sum::<usize>();
            //println!(
            //    "total mutations: {}, n_total_iterations: {}",
            //    n_total_mutations,
            //    n_total_iterations
            //);
            assert!(n_total_iterations <= n_total_mutations as usize);
            let reporter = Arc::new(RwLock::new(ProgressReporter::new(
                n_total_iterations,
                "Rendering",
            )));
            let count_iteration = AtomicUsize::new(0);
            let prev_time = Arc::new(Mutex::new(Instant::now()));
            //let count_iteration = Arc::new(RwLock::new(0 as usize));
            (0..n_chains).into_par_iter().for_each(|i| {
                let reporter = reporter.clone();
                let n_chain_mutations = n_chain_mutations_list[i as usize];

                // Follow {i}th Markov chain for _nChainMutations_
                let mut arena = MemoryArena::new();

                // Select initial state from the set of bootstrap samples
                let mut rng = RNG::new_sequence(i as u64);
                let (bootstrap_index, _, _) = bootstrap.sample_discrete(rng.uniform_float());
                let depth = bootstrap_index as u32 % (max_depth + 1);

                // Initialize local variables for selected state
                let mut sampler = MLTSampler::new(
                    mutations_per_pixel,
                    bootstrap_index as u64,
                    sigma,
                    large_step_probability,
                    n_sample_streams,
                );

                let (mut l_current, mut p_current) = self.l(
                    scene,
                    &mut arena,
                    &light_distr,
                    &light_to_index,
                    &mut sampler,
                    depth,
                );

                // Run the Markov chain for _nChainMutations_ steps
                for _ in 0..n_chain_mutations {
                    sampler.start_iteration();
                    let (l_proposed, p_proposed) = self.l(
                        scene,
                        &mut arena,
                        &light_distr,
                        &light_to_index,
                        &mut sampler,
                        depth,
                    );
                    // Compute acceptance probability for proposed sample
                    let accept = Float::min(1.0, l_proposed.y() / l_current.y());
                    //assert!(accept >= 0.0);

                    // Splat both current and proposed samples to _film_
                    {
                        let mut film = film.write().unwrap();
                        if accept > 0.0 {
                            let spec = l_proposed * (accept / l_proposed.y());
                            film.add_splat(&p_proposed, &spec);
                        }
                        let spec = l_current * ((1.0 - accept) / l_current.y());
                        film.add_splat(&p_current, &spec);
                    }

                    // Accept or reject the proposal
                    if rng.uniform_float() < accept {
                        l_current = l_proposed;
                        p_current = p_proposed;
                        sampler.accept();
                        //acceptedMutations++;
                    } else {
                        sampler.reject();
                    }
                    //totalMutations++;

                    {
                        let mut reporter = reporter.write().unwrap();
                        reporter.update(1);
                    }

                    let iteration = count_iteration.fetch_add(1, Ordering::Acquire);
                    {
                        if (iteration % progress_frequency) == 0 {
                            let mut prev_time = prev_time.lock().unwrap();
                            let elapsed = prev_time.elapsed();
                            if elapsed.as_millis() > UPDATE_DISPLAY_INTERVAL {
                                let mut film = film.write().unwrap();
                                film.merge_splats(splat_scale);
                                let display_scale =
                                    n_total_iterations as Float / iteration as Float;
                                film.update_display_scale(&pixel_bounds, display_scale);
                                *prev_time = Instant::now();
                            }
                        }
                    }
                }
                arena.reset();
            });

            // Store final image computed with MLT
            {
                let camera = self.camera.clone();
                let film = camera.get_film();
                let mut film = film.write().unwrap();
                film.merge_splats(splat_scale);
                film.update_display(&pixel_bounds);
                film.render_end();
                film.write_image();
            }

            {
                let mut reporter = reporter.write().unwrap();
                reporter.done();
            }
        }
    }

    fn get_camera(&self) -> Arc<dyn Camera> {
        return self.camera.clone();
    }
}

unsafe impl Sync for MLTIntegrator {}

pub fn create_mlt_integrator(
    params: &ParamSet,
    camera: &Arc<dyn Camera>,
) -> Result<Arc<RwLock<dyn Integrator>>, PbrtError> {
    let max_depth = params.find_one_int("maxdepth", 5);
    let mut n_bootstrap = params.find_one_int("bootstrapsamples", 100000);
    let n_chains = params.find_one_int("chains", 1000);
    let mut mutations_per_pixel = params.find_one_int("mutationsperpixel", 100);
    let large_step_probability = params.find_one_float("largestepprobability", 0.3);
    let sigma = params.find_one_float("sigma", 0.01);
    {
        let options = PbrtOptions::get();
        if options.quick_render {
            mutations_per_pixel = (mutations_per_pixel / 16).max(1);
            n_bootstrap = (n_bootstrap / 16).max(1);
        }
    }

    return Ok(Arc::new(RwLock::new(MLTIntegrator::new(
        &camera,
        max_depth as u32,
        n_bootstrap as u32,
        n_chains as u32,
        mutations_per_pixel as u32,
        sigma,
        large_step_probability,
    ))));
}
