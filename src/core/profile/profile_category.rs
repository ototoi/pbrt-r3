use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, Copy)]
#[repr(u64)]
pub enum ProfileCategory {
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

impl Display for ProfileCategory {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let index = *self as usize;
        if index < PROF_NAMES.len() {
            return write!(f, "{}", PROF_NAMES[index]);
        } else {
            return write!(f, "Unknown");
        }
    }
}

fn prof_to_bits(p: ProfileCategory) -> u64 {
    return 1u64 << p as u64;
}

impl ProfileCategory {
    pub fn to_bits(&self) -> u64 {
        return prof_to_bits(*self);
    }
}

pub type Prof = ProfileCategory;
