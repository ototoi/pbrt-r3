use crate::core::medium::*;
use crate::core::pbrt::*;
use std::cell::Cell;
use std::sync::Arc;

#[derive(Default, Clone)]
pub struct Ray {
    pub o: Point3f,
    pub d: Vector3f,

    pub t_max: Cell<Float>,

    pub time: Float,
    pub medium: Option<Arc<dyn Medium>>,
}

impl Ray {
    pub fn new(o: &Point3f, d: &Vector3f, t_max: Float, time: Float) -> Self {
        //let d_len = d.length();
        //let t_max = Float::min(t_max * d_len, Float::INFINITY);
        //let d = if d_len != 0.0 { d.normalize() } else { *d };
        let o = *o;
        let d = *d;
        let t_max = Cell::<Float>::new(t_max);
        Ray {
            o,
            d,
            t_max,
            time,
            medium: None,
        }
    }

    pub fn zero() -> Self {
        Ray {
            o: Vector3f::new(0.0, 0.0, 0.0),
            d: Vector3f::new(0.0, 0.0, 0.0),
            t_max: Cell::<Float>::new(10000.0),
            time: 0.0,
            medium: None,
        }
    }

    pub fn position(&self, t: Float) -> Point3f {
        return self.o + self.d * t;
    }
}

impl From<(&Point3f, &Vector3f, Float, Float, &Option<Arc<dyn Medium>>)> for Ray {
    fn from(value: (&Point3f, &Vector3f, Float, Float, &Option<Arc<dyn Medium>>)) -> Self {
        Ray {
            o: *value.0,
            d: *value.1,
            t_max: Cell::<Float>::new(value.2),
            time: value.3,
            medium: value.4.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let r1 = Ray::new(
            &Point3f::new(1.0, 2.0, 3.0),
            &Vector3f::new(1.0, 0.0, 0.0),
            1000.0,
            0.0,
        );
        let p1 = r1.position(4.0);
        let p2 = Point3f::new(5.0, 2.0, 3.0);
        assert_eq!(p1, p2);
    }
}
