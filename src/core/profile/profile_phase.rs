use std::cell::Cell;
use std::fmt::{Display, Formatter};

use crate::core::options::PbrtOptions;

thread_local!(static PROFILER_STATE: Cell<u64> = Cell::new(0));

#[derive(Debug, Clone, Copy)]
#[repr(u64)]
pub enum Prof {
    SceneConstruction = 0,
    AccelConstruction = 1,
    TextureLoading = 2,
    MIPMapCreation = 3,
    IntegratorRender = 4,
    SamplerIntegratorLi = 5,
    SPPMCameraPass = 6,
    SPPMGridConstruction = 7,
    SPPMPhotonPass = 8,
    SPPMStatsUpdate = 9,
    BDPTGenerateSubpath = 10,
    BDPTConnectSubpaths = 11,
    LightDistribLookup = 12,
    LightDistribSpinWait = 13,
    LightDistribCreation = 14,
    DirectLighting = 15,
    BSDFEvaluation = 16,
    BSDFSampling = 17,
    BSDFPdf = 18,
    BSSRDFEvaluation = 19,
    BSSRDFSampling = 20,
    PhaseFuncEvaluation = 21,
    PhaseFuncSampling = 22,
    AccelIntersect = 23,
    AccelIntersectP = 24,
    LightSample = 25,
    LightPdf = 26,
    MediumSample = 27,
    MediumTr = 28,
    TriIntersect = 29,
    TriIntersectP = 30,
    CurveIntersect = 31,
    CurveIntersectP = 32,
    ShapeIntersect = 33,
    ShapeIntersectP = 34,
    ComputeScatteringFuncs = 35,
    GenerateCameraRay = 36,
    MergeFilmTile = 37,
    SplatFilm = 38,
    AddFilmSample = 39,
    StartPixel = 40,
    GetSample = 41,
    TexFiltTrilerp = 42,
    TexFiltEWA = 43,
    TexFiltPtex = 44,
    NumProfCategories = 45,
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
        let index = *self as usize;
        if index < PROF_NAMES.len() {
            return write!(f, "{}", PROF_NAMES[index]);
        } else {
            return write!(f, "Unknown");
        }
    }
}

fn prof_to_bits(p: Prof) -> u64 {
    return 1u64 << p as u64;
}

#[derive(Debug)]
pub struct ProfilePhase {
    pub reset: bool,
    pub category_bit: u64,
}

impl ProfilePhase {
    pub fn new(category: Prof) -> Self {
        let options = PbrtOptions::get();
        if options.no_profile {
            ProfilePhase {
                reset: false,
                category_bit: 0,
            }
        } else {
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
}

impl Drop for ProfilePhase {
    fn drop(&mut self) {
        let options = PbrtOptions::get();
        if !options.no_profile {
            if self.reset {
                PROFILER_STATE.with(|state| {
                    state.set(state.get() & !self.category_bit);
                });
            }
        }
    }
}
