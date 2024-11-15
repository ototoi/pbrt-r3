use crate::core::pbrt::*;

use std::fmt::Formatter;
use std::fmt::Result;
use std::{fmt::Debug, sync::Arc};

#[derive(Clone)]
pub struct BSDF {
    pub eta: Float,
    pub ns: Normal3f,
    pub ng: Normal3f,
    pub ss: Vector3f,
    pub ts: Vector3f,
    pub bxdfs: Vec<Arc<dyn BxDF>>,
}

impl Debug for BSDF {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.debug_struct("BSDF")
            .field("eta", &self.eta)
            .field("ns", &self.ns)
            .field("ng", &self.ng)
            .field("ss", &self.ss)
            .field("ts", &self.ts)
            .field("bxdfs", &self.bxdfs.len())
            .finish()
    }
}

impl BSDF {
    pub fn new(si: &SurfaceInteraction, eta: Float) -> Self {
        let ns = si.shading.n;
        let ng = si.n;
        let ss = si.shading.dpdu.normalize();
        let ts = Vector3f::cross(&ns, &ss).normalize();
        BSDF {
            eta,
            ns,
            ng,
            ss,
            ts,
            bxdfs: Vec::new(),
        }
    }

    pub fn num_components(&self, t: BxDFType) -> i32 {
        let mut num = 0;
        for it in self.bxdfs.iter() {
            let bxdf = it.as_ref();
            if bxdf.matches_flags(t) {
                num += 1;
            }
        }
        return num;
    }

    pub fn world_to_local(&self, v: &Vector3f) -> Vector3f {
        let x = Vector3f::dot(v, &self.ss);
        let y = Vector3f::dot(v, &self.ts);
        let z = Vector3f::dot(v, &self.ns);
        return Vector3f::new(x, y, z);
    }

    pub fn local_to_world(&self, v: &Vector3f) -> Vector3f {
        let x = self.ss.x * v.x + self.ts.x * v.y + self.ns.x * v.z;
        let y = self.ss.y * v.x + self.ts.y * v.y + self.ns.y * v.z;
        let z = self.ss.z * v.x + self.ts.z * v.y + self.ns.z * v.z;
        return Vector3f::new(x, y, z);
    }

    pub fn sample_f(
        &self,
        wo_w: &Vector3f,
        u: &Point2f,
        flags: BxDFType,
    ) -> Option<(Spectrum, Vector3f, Float, BxDFType)> {
        let matching_comps = self.num_components(flags);
        if matching_comps == 0 {
            return None;
        }
        let comp = i32::min(
            Float::floor(u[0] * matching_comps as Float) as i32,
            matching_comps - 1,
        );
        //if matching_comps >= 2 {
        //    println!("comp/matching_comps:{}/{}", comp+1, matching_comps);
        //}
        let n_bxdfs = self.bxdfs.len();
        let mut target_index: i32 = -1;
        {
            //2/4/5 == comp/matching_comps/len(bxdfs)
            //a _ a b a
            //2 1 1 0 0
            let mut count = comp;
            for i in 0..n_bxdfs {
                let bxdf = self.bxdfs[i].as_ref();
                if bxdf.matches_flags(flags) && count == 0 {
                    target_index = i as i32;
                    break;
                }
                count -= 1;
            }
        }
        //assert!(target_index >= 0);
        assert!(target_index >= 0);
        let index = target_index as usize;
        let found_bxdf = self.bxdfs[index].as_ref();
        // Remap _BxDF_ sample _u_ to $[0,1)^2$
        let remapped = Vector2f::from((
            Float::min(
                (u[0] * matching_comps as Float) - comp as Float,
                ONE_MINUS_EPSILON,
            ),
            u[1],
        ));

        let wo = self.world_to_local(wo_w);
        if wo.z == 0.0 {
            return None;
        }
        let mut sampled_type = found_bxdf.get_type();
        if let Some((f, wi, pdf, t)) = found_bxdf.sample_f(&wo, &remapped) {
            assert!(Float::is_finite(f.y()));
            if pdf == 0.0 {
                //t = 0;
                return None;
            }
            let mut f = f;
            let mut pdf = pdf;

            if t != 0 {
                sampled_type = t;
            }

            let bxdf_type = found_bxdf.get_type();

            let wi_world = self.local_to_world(&wi);

            // Compute overall PDF with all matching _BxDF_s
            if (bxdf_type & BSDF_SPECULAR) == 0 && matching_comps > 1 {
                for i in 0..n_bxdfs {
                    if i != index {
                        let bxdf = self.bxdfs[i].as_ref();
                        if bxdf.matches_flags(flags) {
                            pdf += bxdf.pdf(&wo, &wi);
                        }
                    }
                }
            }

            if matching_comps > 1 {
                pdf /= matching_comps as Float;
            }

            // Compute value of BSDF for sampled direction
            if (bxdf_type & BSDF_SPECULAR) == 0 {
                let ng = self.ng;
                let reflect = (Vector3f::dot(&wi_world, &ng) * Vector3f::dot(wo_w, &ng)) > 0.0;
                f = self
                    .bxdfs
                    .iter()
                    .filter(|bxdf| -> bool {
                        let tp = bxdf.get_type();
                        let b1 = bxdf.matches_flags(flags);
                        let b2 = reflect && ((tp & BSDF_REFLECTION) != 0);
                        let b3 = !reflect && ((tp & BSDF_TRANSMISSION) != 0);
                        return b1 && (b2 || b3);
                    })
                    .map(|bxdf| -> Spectrum {
                        return bxdf.f(&wo, &wi);
                    })
                    .fold(Spectrum::zero(), |a, b| a + b);
            }
            assert!(pdf >= 0.0);
            return Some((f, wi_world, pdf, sampled_type));
        } else {
            return None;
            //let f = Spectrum::zero();
            //let wi_w = Vector3f::zero();
            //return Some((f, wi_w, 0.0, sampled_type));
        }
    }

    pub fn f(&self, wo_w: &Vector3f, wi_w: &Vector3f, flags: BxDFType) -> Spectrum {
        let wi = self.world_to_local(wi_w);
        let wo = self.world_to_local(wo_w);

        if wo.z == 0.0 {
            return Spectrum::zero();
        }

        let ng = self.ng;
        let reflect = (Vector3f::dot(wi_w, &ng) * Vector3f::dot(wo_w, &ng)) > 0.0;

        let f = self
            .bxdfs
            .iter()
            .filter(|bxdf| -> bool {
                let tp = bxdf.get_type();
                let b1 = bxdf.matches_flags(flags);
                let b2 = reflect && ((tp & BSDF_REFLECTION) != 0);
                let b3 = !reflect && ((tp & BSDF_TRANSMISSION) != 0);
                return b1 && (b2 || b3);
            })
            .map(|bxdf| -> Spectrum {
                return bxdf.f(&wo, &wi);
            })
            .fold(Spectrum::zero(), |a, b| a + b);
        return f;
    }

    pub fn pdf(&self, wo_w: &Vector3f, wi_w: &Vector3f, flags: BxDFType) -> Float {
        let wi = self.world_to_local(wi_w);
        let wo = self.world_to_local(wo_w);

        if wo.z == 0.0 {
            return 0.0;
        }

        let targets: Vec<_> = self
            .bxdfs
            .iter()
            .filter(|b| -> bool {
                let bxdf = b.as_ref();
                return bxdf.matches_flags(flags);
            })
            .collect();
        let count = targets.len();
        if count > 0 {
            let pdf = targets
                .iter()
                .map(|b| -> Float {
                    let bxdf = b.as_ref();
                    return bxdf.pdf(&wo, &wi);
                })
                .fold(0.0 as Float, |a, b| a + b)
                / (count as Float);
            return pdf;
        } else {
            return 0.0;
        }
    }

    pub fn add(&mut self, bxdf: &Arc<dyn BxDF>) {
        self.bxdfs.push(Arc::clone(bxdf));
    }
}
