use std::fmt::{Display, Formatter};

/*
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
*/

#[derive(Debug, Default, Clone, Copy)]
pub struct ProfileCategory(pub u32);

#[allow(non_upper_case_globals)]
impl ProfileCategory {
    pub const SceneConstruction: Self = Self(0);
    pub const AccelConstruction: Self = Self(1);
    pub const TextureLoading: Self = Self(2);
    pub const MIPMapCreation: Self = Self(3);

    pub const IntegratorRender: Self = Self(4);
    pub const SamplerIntegratorLi: Self = Self(5);
    pub const SPPMCameraPass: Self = Self(6);
    pub const SPPMGridConstruction: Self = Self(7);
    pub const SPPMPhotonPass: Self = Self(8);
    pub const SPPMStatsUpdate: Self = Self(9);
    pub const BDPTGenerateSubpath: Self = Self(10);
    pub const BDPTConnectSubpaths: Self = Self(11);
    pub const LightDistribLookup: Self = Self(12);
    pub const LightDistribSpinWait: Self = Self(13);
    pub const LightDistribCreation: Self = Self(14);
    pub const DirectLighting: Self = Self(15);
    pub const BSDFEvaluation: Self = Self(16);
    pub const BSDFSampling: Self = Self(17);
    pub const BSDFPdf: Self = Self(18);
    pub const BSSRDFEvaluation: Self = Self(19);
    pub const BSSRDFSampling: Self = Self(20);
    pub const PhaseFuncEvaluation: Self = Self(21);
    pub const PhaseFuncSampling: Self = Self(22);
    pub const AccelIntersect: Self = Self(23);
    pub const AccelIntersectP: Self = Self(24);
    pub const LightSample: Self = Self(25);
    pub const LightPdf: Self = Self(26);
    pub const MediumSample: Self = Self(27);
    pub const MediumTr: Self = Self(28);
    pub const TriIntersect: Self = Self(29);
    pub const TriIntersectP: Self = Self(30);
    pub const CurveIntersect: Self = Self(31);
    pub const CurveIntersectP: Self = Self(32);
    pub const ShapeIntersect: Self = Self(33);
    pub const ShapeIntersectP: Self = Self(34);
    pub const ComputeScatteringFuncs: Self = Self(35);
    pub const GenerateCameraRay: Self = Self(36);
    pub const MergeFilmTile: Self = Self(37);
    pub const SplatFilm: Self = Self(38);
    pub const AddFilmSample: Self = Self(39);
    pub const StartPixel: Self = Self(40);
    pub const GetSample: Self = Self(41);
    pub const TexFiltTrilerp: Self = Self(42);
    pub const TexFiltEWA: Self = Self(43);
    pub const TexFiltPtex: Self = Self(44);
    pub const NumProfCategories: Self = Self(45);
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
        let index = self.0 as usize;
        if index < PROF_NAMES.len() {
            return write!(f, "{}", PROF_NAMES[index]);
        } else {
            return write!(f, "Unknown");
        }
    }
}

fn prof_to_bits(p: ProfileCategory) -> u64 {
    return 1u64 << p.0 as u64;
}

impl ProfileCategory {
    pub fn to_bits(&self) -> u64 {
        return prof_to_bits(*self);
    }
}

pub type Prof = ProfileCategory;
