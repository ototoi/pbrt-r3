use std::cell::Cell;
use std::fmt::{Display, Formatter};
use std::sync::atomic::{AtomicI32, Ordering};
use std::sync::{Arc, RwLock};

use super::stats_accumlator::*;
use super::stat_reporter::*;

thread_local!(static PROFILER_STATE: Cell<u64> = Cell::new(0));

#[derive(Debug, Default, Clone, Copy)]
pub struct Prof(u32);

#[allow(non_upper_case_globals)]
impl Prof {
    pub const SceneConstruction: Prof = Prof(0);
    pub const AccelConstruction: Prof = Prof(1);
    pub const TextureLoading: Prof = Prof(2);
    pub const MIPMapCreation: Prof = Prof(3);

    pub const IntegratorRender: Prof = Prof(4);
    pub const SamplerIntegratorLi: Prof = Prof(5);
    pub const SPPMCameraPass: Prof = Prof(6);
    pub const SPPMGridConstruction: Prof = Prof(7);
    pub const SPPMPhotonPass: Prof = Prof(8);
    pub const SPPMStatsUpdate: Prof = Prof(9);
    pub const BDPTGenerateSubpath: Prof = Prof(10);
    pub const BDPTConnectSubpaths: Prof = Prof(11);
    pub const LightDistribLookup: Prof = Prof(12);
    pub const LightDistribSpinWait: Prof = Prof(13);
    pub const LightDistribCreation: Prof = Prof(14);
    pub const DirectLighting: Prof = Prof(15);
    pub const BSDFEvaluation: Prof = Prof(16);
    pub const BSDFSampling: Prof = Prof(17);
    pub const BSDFPdf: Prof = Prof(18);
    pub const BSSRDFEvaluation: Prof = Prof(19);
    pub const BSSRDFSampling: Prof = Prof(20);
    pub const PhaseFuncEvaluation: Prof = Prof(21);
    pub const PhaseFuncSampling: Prof = Prof(22);
    pub const AccelIntersect: Prof = Prof(23);
    pub const AccelIntersectP: Prof = Prof(24);
    pub const LightSample: Prof = Prof(25);
    pub const LightPdf: Prof = Prof(26);
    pub const MediumSample: Prof = Prof(27);
    pub const MediumTr: Prof = Prof(28);
    pub const TriIntersect: Prof = Prof(29);
    pub const TriIntersectP: Prof = Prof(30);
    pub const CurveIntersect: Prof = Prof(31);
    pub const CurveIntersectP: Prof = Prof(32);
    pub const ShapeIntersect: Prof = Prof(33);
    pub const ShapeIntersectP: Prof = Prof(34);
    pub const ComputeScatteringFuncs: Prof = Prof(35);
    pub const GenerateCameraRay: Prof = Prof(36);
    pub const MergeFilmTile: Prof = Prof(37);
    pub const SplatFilm: Prof = Prof(38);
    pub const AddFilmSample: Prof = Prof(39);
    pub const StartPixel: Prof = Prof(40);
    pub const GetSample: Prof = Prof(41);
    pub const TexFiltTrilerp: Prof = Prof(42);
    pub const TexFiltEWA: Prof = Prof(43);
    pub const TexFiltPtex: Prof = Prof(44);
    pub const NumProfCategories: Prof = Prof(45);
}
const PROF_NAMES: [&str; 46] = [
    "Scene parsing and creation",
    "Acceleration structure creation",
    "Texture loading",
    "MIP map generation",
    "Integrator::Render()",
    "SamplerIntegrator::Li()",
    "SPPM camera pass",
    "SPPM grid construction",
    "SPPM photon pass",
    "SPPM photon statistics update",
    "BDPT subpath generation",
    "BDPT subpath connections",
    "SpatialLightDistribution lookup",
    "SpatialLightDistribution spin wait",
    "SpatialLightDistribution creation",
    "Direct lighting",
    "BSDF::f()",
    "BSDF::Sample_f()",
    "BSDF::PDF()",
    "BSSRDF::f()",
    "BSSRDF::Sample_f()",
    "PhaseFunction::p()",
    "PhaseFunction::Sample_p()",
    "Accelerator::Intersect()",
    "Accelerator::IntersectP()",
    "Light::Sample_*()",
    "Light::Pdf()",
    "Medium::Sample()",
    "Medium::Tr()",
    "Triangle::Intersect()",
    "Triangle::IntersectP()",
    "Curve::Intersect()",
    "Curve::IntersectP()",
    "Other Shape::Intersect()",
    "Other Shape::IntersectP()",
    "Material::ComputeScatteringFunctions()",
    "Camera::GenerateRay[Differential]()",
    "Film::MergeTile()",
    "Film::AddSplat()",
    "Film::AddSample()",
    "Sampler::StartPixelSample()",
    "Sampler::GetSample[12]D()",
    "MIPMap::Lookup() (trilinear)",
    "MIPMap::Lookup() (EWA)",
    "Ptex lookup",
    "Num prof categories",
];

impl Display for Prof {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let index = self.0 as usize;
        if index < PROF_NAMES.len() {
            return write!(f, "{}", PROF_NAMES[index]);
        } else {
            return write!(f, "Unknown");
        }
    }
}

fn prof_to_bits(p: Prof) -> u64 {
    return 1u64 << p.0;
}

#[derive(Debug)]
pub struct ProfilePhase {
    pub reset: bool,
    pub category_bit: u64,
}

impl ProfilePhase {
    pub fn new(category: Prof) -> Self {
        let category_bit = prof_to_bits(category);
        let reset = (PROFILER_STATE.get() & category_bit) == 0;
        PROFILER_STATE.with(|state| {
            state.set(state.get() | category_bit);
        });
        ProfilePhase {
            reset,
            category_bit,
        }
    }
}

impl Drop for ProfilePhase {
    fn drop(&mut self) {
        if self.reset {
            PROFILER_STATE.with(|state| {
                state.set(state.get() & !self.category_bit);
            });
        }
    }
}

// For a given profiler state (i.e., a set of "on" bits corresponding to
// profiling categories that are active), ProfileSample stores a count of
// the number of times that state has been active when the timer interrupt
// to record a profiling sample has fired.
struct ProfileSample {
    pub profile_state: u64,
    pub count: u64,
}

struct ProfleReporter {
    pub samples: Vec<ProfileSample>,
}

impl ProfleReporter {
    pub fn new() -> Self {
        ProfleReporter { samples: Vec::new() }
    }

    pub fn sample(&mut self, state: u64) {
        let mut found = false;
        for sample in self.samples.iter_mut() {
            if sample.profile_state == state {
                sample.count += 1;
                found = true;
                break;
            }
        }
        if !found {
            self.samples.push(ProfileSample {
                profile_state: state,
                count: 1,
            });
        }
    }
}

impl StatReporter for ProfleReporter {
    fn report(&self, accum: &mut StatsAccumulator) {
        //
    }

    fn clear(&mut self) {
        //
    }
}

struct ProfileSampler {
    reporter: Arc<RwLock<ProfleReporter>>,
}

impl ProfileSampler {
    pub fn new() -> Self {
        let reporter = Arc::new(RwLock::new(ProfleReporter::new()));
        register_stat_reporter(reporter.clone());
        ProfileSampler { reporter }
    }

    pub fn sample(&mut self) {
        let mut reporter = self.reporter.write().unwrap();
        let state = PROFILER_STATE.with(|state| state.get());
        reporter.sample(state);
    }
}

struct Profiler {
    pub suspend_count: AtomicI32,
}

impl Profiler {
    pub fn new() -> Self {
        PROFILER_STATE.with(|state| {
            state.set(0);
        });
        
        Profiler {
            suspend_count: AtomicI32::new(0),
        }
    }

    pub fn suspend(&mut self) {
        self.suspend_count.fetch_add(1, Ordering::Relaxed);
    }

    pub fn resume(&mut self) {
        self.suspend_count.fetch_sub(1, Ordering::Relaxed);
    }

    pub fn is_suspended(&self) -> bool {
        return self.suspend_count.load(std::sync::atomic::Ordering::Relaxed) > 0;
    }

    pub fn clear(&mut self) {
        PROFILER_STATE.with(|state| {
            state.set(0);
        });
    }
}

pub fn init_profiler() {
    //
}

pub fn suspend_profiler() {
    //
}
