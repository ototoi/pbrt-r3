use crate::core::base::*;

#[derive(Debug, Default, Clone)]
pub struct Distribution1D {
    pub func: Vec<Float>,
    pub cdf: Vec<Float>,
    pub func_int: Float,
    inv_count: Float,
}

#[inline(always)]
fn find_interval_cdf(cdf: &[Float], u: Float) -> usize {
    let mut first = 0usize;
    let mut len = cdf.len();
    while len > 0 {
        let half = len >> 1;
        let middle = first + half;
        if cdf[middle] <= u {
            first = middle + 1;
            len -= half + 1;
        } else {
            len = half;
        }
    }
    let idx = first.saturating_sub(1);
    if idx > cdf.len() - 2 {
        cdf.len() - 2
    } else {
        idx
    }
}

impl Distribution1D {
    pub fn new(f: &[Float]) -> Self {
        let n = f.len();
        let func = Vec::from(f);
        let mut cdf = vec![0.0; n + 1];
        cdf[0] = 0.0;
        for i in 1..(n + 1) {
            cdf[i] = cdf[i - 1] + func[i - 1] / (n as Float);
        }
        let func_int = cdf[n];
        if func_int == 0.0 {
            for i in 1..(n + 1) {
                cdf[i] = (i as Float) / (n as Float);
            }
        } else {
            for i in 1..(n + 1) {
                cdf[i] /= func_int;
            }
        }
        Distribution1D {
            func,
            cdf,
            func_int,
            inv_count: 1.0 / (n as Float),
        }
    }

    pub fn count(&self) -> usize {
        return self.func.len();
    }

    //value, pdf, offset
    #[inline(always)]
    pub fn sample_continuous(&self, u: Float) -> (Float, Float, usize) {
        let offset = find_interval_cdf(&self.cdf, u);

        let cdf0 = self.cdf[offset];
        let cdf1 = self.cdf[offset + 1];
        let mut du = u - cdf0;
        let cdf_span = cdf1 - cdf0;
        if cdf_span > 0.0 {
            du /= cdf_span;
        }
        // Compute PDF for sampled offset
        let pdf = if self.func_int > 0.0 {
            self.func[offset] / self.func_int
        } else {
            0.0
        };
        let r = ((offset as Float) + du) * self.inv_count;
        return (r, pdf, offset);
    }

    //offset, pdf, remapped
    #[inline(always)]
    pub fn sample_discrete(&self, u: Float) -> (usize, Float, Float) {
        let offset = find_interval_cdf(&self.cdf, u);
        let cdf0 = self.cdf[offset];
        let cdf1 = self.cdf[offset + 1];
        assert!(cdf0 <= u);
        assert!(u <= cdf1);
        let pdf = if self.func_int > 0.0 {
            self.func[offset] * self.inv_count / self.func_int
        } else {
            0.0
        };
        let remapped = (u - cdf0) / (cdf1 - cdf0);
        assert!(remapped >= 0.0 && remapped <= 1.0);
        return (offset, pdf, remapped);
    }

    pub fn discrete_pdf(&self, index: usize) -> Float {
        return self.func[index] / (self.func_int * self.count() as Float);
    }
}

pub struct Distribution2D {
    pub conditional_v: Vec<Distribution1D>,
    pub marginal: Distribution1D,
}

impl Distribution2D {
    pub fn new(data: &[Float], nu: usize, nv: usize) -> Self {
        let mut conditional_v = Vec::with_capacity(nv);
        for v in 0..nv {
            let a = v * nu;
            let b = a + nu;
            let s = &data[a..b];
            conditional_v.push(Distribution1D::new(s));
        }
        let mut marginal_func = Vec::with_capacity(nv);
        for v in 0..nv {
            marginal_func.push(conditional_v[v].func_int);
        }
        Distribution2D {
            conditional_v,
            marginal: Distribution1D::new(&marginal_func),
        }
    }

    pub fn sample_continuous(&self, u: &Point2f) -> (Point2f, Float) {
        let (d1, pdf1, v) = self.marginal.sample_continuous(u[1]);
        let (d0, pdf0, _) = self.conditional_v[v].sample_continuous(u[0]);
        return (Point2f::new(d0, d1), pdf0 * pdf1);
    }
    pub fn pdf(&self, p: &Point2f) -> Float {
        let ucount = self.conditional_v[0].count();
        let vcount = self.marginal.count();
        let iu = usize::clamp((p[0] * ucount as Float) as usize, 0, ucount - 1);
        let iv = usize::clamp((p[1] * vcount as Float) as usize, 0, vcount - 1);
        return self.conditional_v[iv].func[iu] / self.marginal.func_int;
    }
}
