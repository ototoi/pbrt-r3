use super::subpath::*;
use crate::core::camera::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::interaction::*;
use crate::core::light::*;
use crate::core::lightdistrib::*;
use crate::core::material::*;
use crate::core::memory::*;
use crate::core::misc::*;
use crate::core::options::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::profile::*;
use crate::core::sampler::*;
use crate::core::sampling::*;
use crate::core::scene::*;
use crate::core::spectrum::*;

use std::sync::atomic::*;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::RwLock;
use std::time::Instant;

use rayon::iter::*;

// MLTSampler Constants
const CAMERA_STREAM_INDEX: u64 = 0;
const LIGHT_STREAM_INDEX: u64 = 1;
const CONNECTION_STREAM_INDEX: u64 = 2;
const N_SAMPLE_STREAMS: u64 = 3;

// MLT Constants
const PROGRESS_FREQUENCY: u32 = 32768;
const UPDATE_DISPLAY_INTERVAL: u128 = 1000;

#[derive(Debug, PartialEq, Default, Clone)]
pub struct PrimarySample {
    // PrimarySample Public Data
    pub value: Float,
    pub last_modification_iteration: i64,
    pub value_backup: Float,
    pub modify_backup: i64,
}

impl PrimarySample {
    // PrimarySample Public Data
    pub fn backup(&mut self) {
        self.value_backup = self.value;
        self.modify_backup = self.last_modification_iteration;
    }
    pub fn restore(&mut self) {
        self.value = self.value_backup;
        self.last_modification_iteration = self.modify_backup;
    }
}

#[derive(Debug, PartialEq, Default, Clone)]
pub struct MLTSampler {
    base: BaseSampler,
    rng: RNG,
    sigma: Float,
    large_step_probability: Float,
    stream_count: u64,

    x: Vec<PrimarySample>,
    current_iteration: i64,
    large_step: bool,
    last_large_step_iteration: i64,
    stream_index: u64,
    sample_index: u64,
}

impl MLTSampler {
    // MLTSampler Public Methods
    //
    pub fn new(
        mutation_per_pixel: u32,
        rng_sequence_index: u64,
        sigma: Float,
        large_step_probability: Float,
        stream_count: u64,
    ) -> Self {
        let x = Vec::new();
        let rng = RNG::new_sequence(rng_sequence_index);
        MLTSampler {
            base: BaseSampler::new(mutation_per_pixel),
            rng,
            sigma,
            large_step_probability,
            stream_count: stream_count,
            x,
            current_iteration: 0,
            large_step: true,
            last_large_step_iteration: 0,
            stream_index: 0,
            sample_index: 0,
        }
    }

    pub fn start_iteration(&mut self) {
        self.current_iteration += 1;
        self.large_step = self.rng.uniform_float() < self.large_step_probability;
    }

    pub fn accept(&mut self) {
        if self.large_step {
            self.last_large_step_iteration = self.current_iteration;
        }
    }

    pub fn reject(&mut self) {
        for xi in &mut self.x {
            if xi.last_modification_iteration == self.current_iteration {
                xi.restore();
            }
        }
        self.current_iteration -= 1;
    }

    pub fn start_stream(&mut self, index: u64) {
        assert!(index < self.stream_count);
        self.stream_index = index;
        self.sample_index = 0;
    }

    pub fn get_next_index(&mut self) -> u64 {
        let index = self.stream_index + self.stream_count * self.sample_index;
        self.sample_index += 1;
        return index;
    }

    fn ensure_ready(&mut self, index: usize) {
        // Enlarge _MLTSampler::X_ if necessary and get current $\VEC{X}_i$
        if index >= self.x.len() {
            self.x.resize(index + 1, PrimarySample::default());
        }
        let xi = &mut self.x[index];

        // Reset $\VEC{X}_i$ if a large step took place in the meantime
        if xi.last_modification_iteration < self.last_large_step_iteration {
            xi.value = self.rng.uniform_float();
            xi.last_modification_iteration = self.last_large_step_iteration;
        }

        // Apply remaining sequence of mutations to _xi_
        xi.backup();
        if self.large_step {
            xi.value = self.rng.uniform_float();
        } else {
            let n_small = self.current_iteration - xi.last_modification_iteration;
            // Apply _n_small_ small step mutations

            // Sample the standard normal distribution $N(0, 1)$
            let normal_sample = SQRT_2 * erf_inv(2.0 * self.rng.uniform_float() - 1.0);

            // Compute the effective standard deviation and apply perturbation to $\VEC{X}_i$
            // $\VEC{X}_i$
            let eff_sigma = self.sigma * Float::sqrt(n_small as Float);
            xi.value += normal_sample * eff_sigma;
            xi.value -= Float::floor(xi.value);
        }
        xi.last_modification_iteration = self.current_iteration;
    }
}

impl Sampler for MLTSampler {
    // Sampler Public Methods
    fn start_pixel(&mut self, p: &Point2i) {
        self.base.start_pixel(p);
    }

    fn get_1d(&mut self) -> Float {
        let _p = ProfilePhase::new(Prof::GetSample);
        let index = self.get_next_index() as usize;
        self.ensure_ready(index);
        return self.x[index].value;
    }

    fn get_2d(&mut self) -> Point2f {
        let x = self.get_1d();
        let y = self.get_1d();
        return Point2f { x, y };
    }

    fn request_1d_array(&mut self, n: u32) {
        self.base.request_1d_array(n);
    }

    fn request_2d_array(&mut self, n: u32) {
        self.base.request_2d_array(n);
    }

    fn get_1d_array(&mut self, n: u32) -> Option<Vec<Float>> {
        return self.base.get_1d_array(n);
    }

    fn get_2d_array(&mut self, n: u32) -> Option<Vec<Vector2f>> {
        return self.base.get_2d_array(n);
    }

    fn start_next_sample(&mut self) -> bool {
        return self.base.start_next_sample();
    }

    fn clone_with_seed(&self, seed: u32) -> Arc<RwLock<dyn Sampler>> {
        let sampler = MLTSampler {
            base: self.base.clone(),
            rng: RNG::new_sequence(seed as u64),
            sigma: self.sigma,
            large_step_probability: self.large_step_probability,
            stream_count: self.stream_count,
            x: self.x.clone(),
            current_iteration: self.current_iteration,
            large_step: self.large_step,
            last_large_step_iteration: self.last_large_step_iteration,
            stream_index: self.stream_index,
            sample_index: self.sample_index,
        };
        return Arc::new(RwLock::new(sampler));
    }

    fn get_samples_per_pixel(&self) -> u32 {
        return self.base.get_samples_per_pixel();
    }
}

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
        camera_vertices: &mut Vec<Arc<Vertex>>,
        light_vertices: &mut Vec<Arc<Vertex>>,
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
        //let mut camera_vertices = Vec::with_capacity(t as usize);
        camera_vertices.clear();
        if generate_camera_subpath(
            scene,
            sampler,
            arena,
            t as usize,
            &camera,
            &p_raster,
            camera_vertices,
        ) != t as usize
        {
            return (Spectrum::zero(), p_raster);
        };

        // Generate a light subpath with exactly _s_ vertices
        sampler.start_stream(LIGHT_STREAM_INDEX);
        let time = camera_vertices[0].get_time();
        //let mut light_vertices = Vec::with_capacity(s as usize);
        light_vertices.clear();
        if generate_light_subpath(
            scene,
            sampler,
            arena,
            s as usize,
            time,
            light_distr,
            light_to_index,
            light_vertices,
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

#[inline]
fn safe_div(x: Float, y: Float) -> Float {
    if x == 0.0 && y == 0.0 {
        return 0.0;
    } else if y == 0.0 {
        return Float::INFINITY;
    } else {
        return x / y;
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
                let mut camera_vertices = Vec::with_capacity(max_depth as usize + 2);
                let mut light_vertices = Vec::with_capacity(max_depth as usize + 1);

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
                        &mut camera_vertices,
                        &mut light_vertices,
                        depth,
                    );
                    let y = l.y();
                    if y > 0.0 {
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
            let pixel_bounds: Bounds2<i32> = film.read().unwrap().cropped_pixel_bounds;
            let n_total_mutations = (self.mutations_per_pixel as i64)
                * (film.read().unwrap().get_sample_bounds().area() as i64);
            let progress_frequency = PROGRESS_FREQUENCY as usize;
            let n_chains = self.n_chains as i64;
            let mutations_per_pixel = self.mutations_per_pixel;
            let sigma = self.sigma;
            let large_step_probability = self.large_step_probability;
            let n_sample_streams = N_SAMPLE_STREAMS;

            //let n_total_work = (n_total_mutations / progress_frequency) as usize;
            let permutations_per_chain = n_total_mutations / n_chains;
            //println!("n_total_mutations: {}", n_total_mutations);
            //println!("n_chains: {}", n_chains);
            //println!("permutations_per_chain: {}", permutations_per_chain);

            let mut n_chain_mutations_list: Vec<usize> = (0..n_chains)
                .map(|i| {
                    (i64::min((i + 1) * permutations_per_chain, n_total_mutations)
                        - (i * permutations_per_chain)) as usize
                })
                .collect();
            let n_total_iterations = n_chain_mutations_list.iter().sum::<usize>();
            //let l = n_chain_mutations_list.len();
            //let nn = &n_chain_mutations_list[l-20..];
            //println!("n_chain_mutations_list: {:?}", nn);
            assert!(n_total_iterations <= n_total_mutations as usize);
            let last = n_chain_mutations_list.len() - 1;
            n_chain_mutations_list[last] +=
                usize::max(0, n_total_mutations as usize - n_total_iterations); //pbrt-r3
            let n_total_iterations = n_chain_mutations_list.iter().sum::<usize>();
            assert!(n_total_iterations == n_total_mutations as usize);
            //println!(
            //    "total mutations: {}, n_total_iterations: {}",
            //    n_total_mutations,
            //    n_total_iterations
            //);

            //let n_total_mutations = n_total_iterations;
            let reporter = Arc::new(RwLock::new(ProgressReporter::new(
                n_total_iterations,
                "Rendering",
            )));
            let count_iteration = AtomicUsize::new(0);
            let prev_time = Arc::new(Mutex::new(Instant::now()));
            //let count_iteration = Arc::new(RwLock::new(0 as usize));
            let n_chains = self.n_chains as usize;
            (0..n_chains).into_par_iter().for_each(|i| {
                let reporter = reporter.clone();
                let n_chain_mutations = n_chain_mutations_list[i];

                // Follow {i}th Markov chain for _nChainMutations_
                let mut arena = MemoryArena::new();
                let mut camera_vertices = Vec::with_capacity(max_depth as usize + 2);
                let mut light_vertices = Vec::with_capacity(max_depth as usize + 1);

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
                    &mut camera_vertices,
                    &mut light_vertices,
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
                        &mut camera_vertices,
                        &mut light_vertices,
                        depth,
                    );
                    // Compute acceptance probability for proposed sample
                    let accept = Float::min(1.0, safe_div(l_proposed.y(), l_current.y()));
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
        camera,
        max_depth as u32,
        n_bootstrap as u32,
        n_chains as u32,
        mutations_per_pixel as u32,
        sigma,
        large_step_probability,
    ))));
}
