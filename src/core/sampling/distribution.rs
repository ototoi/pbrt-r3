use crate::core::base::*;

#[derive(Debug, Default, Clone)]
pub struct Distribution1D {
    pub func: Vec<Float>,
    pub cdf: Vec<Float>,
    pub func_int: Float,
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
        }
    }

    pub fn count(&self) -> usize {
        return self.func.len();
    }

    //value, pdf, offset
    pub fn sample_continuous(&self, u: Float) -> (Float, Float, usize) {
        let offset = find_interval(&self.cdf, &|vv, i| -> bool {
            return vv[i] <= u;
        });

        let mut du = u - self.cdf[offset];
        if (self.cdf[offset + 1] - self.cdf[offset]) > 0.0 {
            du /= self.cdf[offset + 1] - self.cdf[offset];
        }
        // Compute PDF for sampled offset
        let pdf = if self.func_int > 0.0 {
            self.func[offset] / self.func_int
        } else {
            0.0
        };
        let r = ((offset as Float) + du) / (self.count() as Float);
        return (r, pdf, offset);
    }

    //offset, pdf, remapped
    pub fn sample_discrete(&self, u: Float) -> (usize, Float, Float) {
        let offset = find_interval(&self.cdf, &|vv, i| -> bool {
            return vv[i] <= u;
        });
        assert!(self.cdf[offset] <= u);
        assert!(u <= self.cdf[offset + 1]);
        let pdf = if self.func_int > 0.0 {
            self.func[offset] / (self.func_int * self.count() as Float)
        } else {
            0.0
        };
        let remapped = (u - self.cdf[offset]) / (self.cdf[offset + 1] - self.cdf[offset]);
        assert!(remapped >= 0.0 && remapped <= 1.0);
        return (offset, pdf, remapped);
    }

    pub fn discrete_pdf(&self, index: usize) -> Float {
        return self.func[index] / (self.func_int * self.count() as Float);
    }
}

pub struct Distribution2D {
    pub conditional_v: Vec<Box<Distribution1D>>,
    pub marginal: Box<Distribution1D>,
}

impl Distribution2D {
    pub fn new(data: &[Float], nu: usize, nv: usize) -> Self {
        let mut conditional_v = Vec::with_capacity(nv);
        for v in 0..nv {
            let a = v * nu;
            let b = a + nu;
            let s = &data[a..b];
            conditional_v.push(Box::new(Distribution1D::new(s)));
        }
        let mut marginal_func = Vec::with_capacity(nv);
        for v in 0..nv {
            marginal_func.push(conditional_v[v].func_int);
        }
        Distribution2D {
            conditional_v,
            marginal: Box::new(Distribution1D::new(&marginal_func)),
        }
    }

    pub fn sample_continuous(&self, u: &Point2f) -> (Point2f, Float) {
        let (d1, pdf1, v) = self.marginal.sample_continuous(u[1]);
        let (d0, pdf0, _) = self.conditional_v[v].as_ref().sample_continuous(u[0]);
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
