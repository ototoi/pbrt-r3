use crate::core::base::*;
use crate::core::reflection::*;
use crate::core::spectrum::*;
use std::sync::Arc;

pub enum BxDFEnum {
    Dyn(Arc<dyn BxDF>),
    LambertianReflection(LambertianReflection),
    LambertianTransmission(LambertianTransmission),
    OrenNayar(OrenNayar),
    MicrofacetReflection(MicrofacetReflection),
    MicrofacetTransmission(MicrofacetTransmission),
    FresnelBlend(FresnelBlend),
    ScaledBxDF(ScaledBxDF),
    SpecularReflection(SpecularReflection),
    SpecularTransmission(SpecularTransmission),
    FresnelSpecular(FresnelSpecular),
}

impl BxDFEnum {
    #[inline]
    pub fn matches_flags(&self, t: BxDFType) -> bool {
        match self {
            BxDFEnum::Dyn(v) => v.matches_flags(t),
            BxDFEnum::LambertianReflection(v) => v.matches_flags(t),
            BxDFEnum::LambertianTransmission(v) => v.matches_flags(t),
            BxDFEnum::OrenNayar(v) => v.matches_flags(t),
            BxDFEnum::MicrofacetReflection(v) => v.matches_flags(t),
            BxDFEnum::MicrofacetTransmission(v) => v.matches_flags(t),
            BxDFEnum::FresnelBlend(v) => v.matches_flags(t),
            BxDFEnum::ScaledBxDF(v) => v.matches_flags(t),
            BxDFEnum::SpecularReflection(v) => v.matches_flags(t),
            BxDFEnum::SpecularTransmission(v) => v.matches_flags(t),
            BxDFEnum::FresnelSpecular(v) => v.matches_flags(t),
        }
    }

    #[inline]
    pub fn get_type(&self) -> BxDFType {
        match self {
            BxDFEnum::Dyn(v) => v.get_type(),
            BxDFEnum::LambertianReflection(v) => v.get_type(),
            BxDFEnum::LambertianTransmission(v) => v.get_type(),
            BxDFEnum::OrenNayar(v) => v.get_type(),
            BxDFEnum::MicrofacetReflection(v) => v.get_type(),
            BxDFEnum::MicrofacetTransmission(v) => v.get_type(),
            BxDFEnum::FresnelBlend(v) => v.get_type(),
            BxDFEnum::ScaledBxDF(v) => v.get_type(),
            BxDFEnum::SpecularReflection(v) => v.get_type(),
            BxDFEnum::SpecularTransmission(v) => v.get_type(),
            BxDFEnum::FresnelSpecular(v) => v.get_type(),
        }
    }

    #[inline]
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        match self {
            BxDFEnum::Dyn(v) => v.f(wo, wi),
            BxDFEnum::LambertianReflection(v) => v.f(wo, wi),
            BxDFEnum::LambertianTransmission(v) => v.f(wo, wi),
            BxDFEnum::OrenNayar(v) => v.f(wo, wi),
            BxDFEnum::MicrofacetReflection(v) => v.f(wo, wi),
            BxDFEnum::MicrofacetTransmission(v) => v.f(wo, wi),
            BxDFEnum::FresnelBlend(v) => v.f(wo, wi),
            BxDFEnum::ScaledBxDF(v) => v.f(wo, wi),
            BxDFEnum::SpecularReflection(v) => v.f(wo, wi),
            BxDFEnum::SpecularTransmission(v) => v.f(wo, wi),
            BxDFEnum::FresnelSpecular(v) => v.f(wo, wi),
        }
    }

    #[inline]
    pub fn rho(&self, wo: &Vector3f, samples: &[Point2f]) -> Spectrum {
        match self {
            BxDFEnum::Dyn(v) => v.rho(wo, samples),
            BxDFEnum::LambertianReflection(v) => v.rho(wo, samples),
            BxDFEnum::LambertianTransmission(v) => v.rho(wo, samples),
            BxDFEnum::OrenNayar(v) => v.rho(wo, samples),
            BxDFEnum::MicrofacetReflection(v) => v.rho(wo, samples),
            BxDFEnum::MicrofacetTransmission(v) => v.rho(wo, samples),
            BxDFEnum::FresnelBlend(v) => v.rho(wo, samples),
            BxDFEnum::ScaledBxDF(v) => v.rho(wo, samples),
            BxDFEnum::SpecularReflection(v) => v.rho(wo, samples),
            BxDFEnum::SpecularTransmission(v) => v.rho(wo, samples),
            BxDFEnum::FresnelSpecular(v) => v.rho(wo, samples),
        }
    }

    #[inline]
    pub fn rho2(&self, samples: &[(Point2f, Point2f)]) -> Spectrum {
        match self {
            BxDFEnum::Dyn(v) => v.rho2(samples),
            BxDFEnum::LambertianReflection(v) => v.rho2(samples),
            BxDFEnum::LambertianTransmission(v) => v.rho2(samples),
            BxDFEnum::OrenNayar(v) => v.rho2(samples),
            BxDFEnum::MicrofacetReflection(v) => v.rho2(samples),
            BxDFEnum::MicrofacetTransmission(v) => v.rho2(samples),
            BxDFEnum::FresnelBlend(v) => v.rho2(samples),
            BxDFEnum::ScaledBxDF(v) => v.rho2(samples),
            BxDFEnum::SpecularReflection(v) => v.rho2(samples),
            BxDFEnum::SpecularTransmission(v) => v.rho2(samples),
            BxDFEnum::FresnelSpecular(v) => v.rho2(samples),
        }
    }

    #[inline]
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        sample: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        match self {
            BxDFEnum::Dyn(v) => v.sample_f(wo, sample),
            BxDFEnum::LambertianReflection(v) => v.sample_f(wo, sample),
            BxDFEnum::LambertianTransmission(v) => v.sample_f(wo, sample),
            BxDFEnum::OrenNayar(v) => v.sample_f(wo, sample),
            BxDFEnum::MicrofacetReflection(v) => v.sample_f(wo, sample),
            BxDFEnum::MicrofacetTransmission(v) => v.sample_f(wo, sample),
            BxDFEnum::FresnelBlend(v) => v.sample_f(wo, sample),
            BxDFEnum::ScaledBxDF(v) => v.sample_f(wo, sample),
            BxDFEnum::SpecularReflection(v) => v.sample_f(wo, sample),
            BxDFEnum::SpecularTransmission(v) => v.sample_f(wo, sample),
            BxDFEnum::FresnelSpecular(v) => v.sample_f(wo, sample),
        }
    }

    #[inline]
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        match self {
            BxDFEnum::Dyn(v) => v.pdf(wo, wi),
            BxDFEnum::LambertianReflection(v) => v.pdf(wo, wi),
            BxDFEnum::LambertianTransmission(v) => v.pdf(wo, wi),
            BxDFEnum::OrenNayar(v) => v.pdf(wo, wi),
            BxDFEnum::MicrofacetReflection(v) => v.pdf(wo, wi),
            BxDFEnum::MicrofacetTransmission(v) => v.pdf(wo, wi),
            BxDFEnum::FresnelBlend(v) => v.pdf(wo, wi),
            BxDFEnum::ScaledBxDF(v) => v.pdf(wo, wi),
            BxDFEnum::SpecularReflection(v) => v.pdf(wo, wi),
            BxDFEnum::SpecularTransmission(v) => v.pdf(wo, wi),
            BxDFEnum::FresnelSpecular(v) => v.pdf(wo, wi),
        }
    }
}

impl From<Arc<dyn BxDF>> for BxDFEnum {
    fn from(value: Arc<dyn BxDF>) -> Self {
        BxDFEnum::Dyn(value)
    }
}

impl From<LambertianReflection> for BxDFEnum {
    fn from(value: LambertianReflection) -> Self {
        BxDFEnum::LambertianReflection(value)
    }
}

impl From<LambertianTransmission> for BxDFEnum {
    fn from(value: LambertianTransmission) -> Self {
        BxDFEnum::LambertianTransmission(value)
    }
}

impl From<OrenNayar> for BxDFEnum {
    fn from(value: OrenNayar) -> Self {
        BxDFEnum::OrenNayar(value)
    }
}

impl From<MicrofacetReflection> for BxDFEnum {
    fn from(value: MicrofacetReflection) -> Self {
        BxDFEnum::MicrofacetReflection(value)
    }
}

impl From<MicrofacetTransmission> for BxDFEnum {
    fn from(value: MicrofacetTransmission) -> Self {
        BxDFEnum::MicrofacetTransmission(value)
    }
}

impl From<FresnelBlend> for BxDFEnum {
    fn from(value: FresnelBlend) -> Self {
        BxDFEnum::FresnelBlend(value)
    }
}

impl From<ScaledBxDF> for BxDFEnum {
    fn from(value: ScaledBxDF) -> Self {
        BxDFEnum::ScaledBxDF(value)
    }
}

impl From<SpecularReflection> for BxDFEnum {
    fn from(value: SpecularReflection) -> Self {
        BxDFEnum::SpecularReflection(value)
    }
}

impl From<SpecularTransmission> for BxDFEnum {
    fn from(value: SpecularTransmission) -> Self {
        BxDFEnum::SpecularTransmission(value)
    }
}

impl From<FresnelSpecular> for BxDFEnum {
    fn from(value: FresnelSpecular) -> Self {
        BxDFEnum::FresnelSpecular(value)
    }
}
