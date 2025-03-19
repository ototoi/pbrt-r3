use crate::core::distribution::*;
use crate::core::pbrt::*;
use crate::core::refrection::*;
use crate::core::sampling::*;
use crate::core::spectrum::*;

pub struct FresnelBlend {
    rd: Spectrum,
    rs: Spectrum,
    distribution: Box<dyn MicrofacetDistribution>,
}

fn pow5(v: Float) -> Float {
    return v.powf(5.0);
}

impl FresnelBlend {
    pub fn new(
        rd: &Spectrum,
        rs: &Spectrum,
        distribution: Box<dyn MicrofacetDistribution>,
    ) -> Self {
        FresnelBlend {
            rd: *rd,
            rs: *rs,
            distribution: distribution,
        }
    }

    pub fn schlick_fresnel(&self, cos_theta: Float) -> Spectrum {
        return self.rs + (Spectrum::one() - self.rs) * pow5(1.0 - cos_theta);
    }
}

impl BxDF for FresnelBlend {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let diffuse = self.rd
            * (Spectrum::one() - self.rs)
            * (1.0 - pow5(1.0 - 0.5 * abs_cos_theta(wi)))
            * (1.0 - pow5(1.0 - 0.5 * abs_cos_theta(wo)))
            * (28.0 / (23.0 * PI));
        let wh = *wi + *wo;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::zero();
        }
        let wh = wh.normalize();
        let distribution = self.distribution.as_ref();
        let specular = self.schlick_fresnel(Vector3f::dot(wi, &wh))
            * (distribution.d(&wh)
                / (4.0
                    * Vector3f::abs_dot(wi, &wh)
                    * Float::max(abs_cos_theta(wi), abs_cos_theta(wo))));

        return diffuse + specular;
    }

    fn sample_f(
        &self,
        wo: &Vector3f,
        u_orig: &Vector2f,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        let (mut ux, uy) = (u_orig.x, u_orig.y);
        let wi = if ux < 0.5 {
            ux = Float::min(2.0 * ux, ONE_MINUS_EPSILON);
            let mut wi = cosine_sample_hemisphere(&Point2f::new(ux, uy));
            if wo.z < 0.0 {
                wi.z *= -1.0;
            }
            wi
        } else {
            ux = Float::min(2.0 * (ux - 0.5), ONE_MINUS_EPSILON);
            // Sample microfacet orientation $\wh$ and reflected direction $\wi$
            let distribution = self.distribution.as_ref();
            let wh = distribution.sample_wh(wo, &Point2f::new(ux, uy));
            let wi = reflect(wo, &wh);
            if !same_hemisphere(wo, &wi) {
                return None;
            }
            wi
        };

        let pdf = self.pdf(wo, &wi);
        if pdf > 0.0 {
            let f = self.f(wo, &wi);
            return Some((f, wi, pdf, 0));
        } else {
            return None;
        }
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !same_hemisphere(wo, wi) {
            return 0.0;
        }
        let wh = (*wo + *wi).normalize();
        // pbrt-r3:
        if Vector3f::dot(wo, &wh) < 0.0 {
            return 0.0; // Should be rare
        }
        // pbrt-r3:
        let distribution = self.distribution.as_ref();
        let pdf_wh = distribution.pdf(wo, &wh);
        return 0.5 * (abs_cos_theta(wi) * INV_PI + pdf_wh / (4.0 * Vector3f::dot(wo, &wh)));
    }

    fn get_type(&self) -> BxDFType {
        return BSDF_REFLECTION | BSDF_GLOSSY;
    }

    fn to_string(&self) -> String {
        return format!("FresnelBlend {:?}", self.get_type());
    }
}
