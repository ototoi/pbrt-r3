use super::intersect::*;
use crate::core::pbrt::*;

#[derive(Debug, PartialEq, Default, Copy, Clone)]
pub struct Bounds3<T> {
    pub min: Vector3<T>,
    pub max: Vector3<T>,
}

pub type Bounds3f = Bounds3<f32>;
pub type Bounds3d = Bounds3<f64>;
pub type Bounds3i = Bounds3<i32>;

impl<T: Copy> Bounds3<T> {
    pub fn new(min: &Vector3<T>, max: &Vector3<T>) -> Self {
        Bounds3::<T> {
            min: *min,
            max: *max,
        }
    }

    pub fn corner(&self, i: usize) -> Vector3<T> {
        assert!(i < 8);
        return match i {
            0 => Vector3::<T>::new(self.min.x, self.min.y, self.min.z),
            1 => Vector3::<T>::new(self.max.x, self.min.y, self.min.z),
            2 => Vector3::<T>::new(self.min.x, self.max.y, self.min.z),
            3 => Vector3::<T>::new(self.max.x, self.max.y, self.min.z),
            4 => Vector3::<T>::new(self.min.x, self.min.y, self.max.z),
            5 => Vector3::<T>::new(self.max.x, self.min.y, self.max.z),
            6 => Vector3::<T>::new(self.min.x, self.max.y, self.max.z),
            7 => Vector3::<T>::new(self.max.x, self.max.y, self.max.z),
            _ => Vector3::<T>::new(self.min.x, self.min.y, self.min.z),
        };
    }
}

fn min_<T: Copy + PartialOrd>(a: T, b: T) -> T {
    return if a <= b { a } else { b };
}

fn max_<T: Copy + PartialOrd>(a: T, b: T) -> T {
    return if a >= b { a } else { b };
}

impl<
        T: Copy
            + PartialOrd
            + std::ops::Add<Output = T>
            + std::ops::Sub<Output = T>
            + std::ops::Mul<Output = T>
            + std::ops::Div<Output = T>,
    > Bounds3<T>
{
    pub fn area(&self) -> T {
        return (self.max.x - self.min.x) * (self.max.y - self.min.y);
    }
    pub fn diagonal(&self) -> Vector3<T> {
        return self.max - self.min;
    }
    pub fn maximum_extent(&self) -> usize {
        let d = self.diagonal();
        if d.x > d.y && d.x > d.z {
            return 0;
        } else if d.y > d.z {
            return 1;
        } else {
            return 2;
        }
    }

    pub fn offset(&self, p: &Vector3<T>) -> Vector3<T> {
        let mut o = *p - self.min;
        if self.max.x > self.min.x {
            o.x = o.x / (self.max.x - self.min.x);
        }
        if self.max.y > self.min.y {
            o.y = o.y / (self.max.y - self.min.y);
        }
        if self.max.z > self.min.z {
            o.z = o.z / (self.max.z - self.min.z);
        }
        return o;
    }

    pub fn expand(&self, delta: T) -> Self {
        let delta = Vector3::<T>::new(delta, delta, delta);
        return Bounds3 {
            min: self.min - delta,
            max: self.max + delta,
        };
    }

    pub fn union(&self, other: &Self) -> Self {
        let min = Vector3::<T>::new(
            min_(self.min.x, other.min.x),
            min_(self.min.y, other.min.y),
            min_(self.min.z, other.min.z),
        );
        let max = Vector3::<T>::new(
            max_(self.max.x, other.max.x),
            max_(self.max.y, other.max.y),
            max_(self.max.z, other.max.z),
        );
        return Bounds3 { min, max };
    }

    pub fn intersect(&self, other: &Self) -> Self {
        let min = Vector3::<T>::new(
            max_(self.min.x, other.min.x),
            max_(self.min.y, other.min.y),
            max_(self.min.z, other.min.z),
        );
        let max = Vector3::<T>::new(
            min_(self.max.x, other.max.x),
            min_(self.max.y, other.max.y),
            min_(self.max.z, other.max.z),
        );
        return Bounds3 { min, max };
    }

    pub fn union_p(&self, other: &Vector3<T>) -> Self {
        let min = Vector3::<T>::new(
            min_(self.min.x, other.x),
            min_(self.min.y, other.y),
            min_(self.min.z, other.z),
        );
        let max = Vector3::<T>::new(
            max_(self.max.x, other.x),
            max_(self.max.y, other.y),
            max_(self.max.z, other.z),
        );
        return Bounds3 { min, max };
    }

    pub fn inside_exclusive(&self, p: &Vector3<T>) -> bool {
        return p.x >= self.min.x
            && p.x < self.max.x
            && p.y >= self.min.y
            && p.y < self.max.y
            && p.z >= self.min.z
            && p.z < self.max.z;
    }
}

impl Bounds3f {
    /*
    pub fn intersect_p(&self, ray: &Ray) -> (bool, Float, Float) {
        let t0 = 0.0;
        let t1 = ray.t_max.get();
        return intersect_box(&self.min, &self.max, &ray.o, &ray.d, t0, t1);
    }
    */
    pub fn intersect_p(&self, ray: &Ray) -> Option<(Float, Float)> {
        let t0 = 0.0;
        let t1 = ray.t_max.get();
        let (b, tmin, tmax) = intersect_box(&self.min, &self.max, &ray.o, &ray.d, t0, t1);
        if b {
            return Some((tmin, tmax));
        } else {
            return None;
        }
    }

    pub fn lerp(&self, t: &Vector3f) -> Vector3f {
        return Vector3f::new(
            lerp(t.x, self.min.x, self.max.x),
            lerp(t.y, self.min.y, self.max.y),
            lerp(t.z, self.min.z, self.max.z),
        );
    }

    pub fn bounding_sphere(&self) -> (Vector3f, Float) {
        let center = (self.min + self.max) * 0.5;
        let radius = (self.max - self.min).length() * 0.5;
        return (center, radius);
    }

    pub fn surface_area(&self) -> Float {
        let d = self.diagonal();
        return 2.0 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
}

impl<T: Copy> From<((T, T, T), (T, T, T))> for Bounds3<T> {
    fn from(value: ((T, T, T), (T, T, T))) -> Self {
        Bounds3::<T> {
            min: Vector3::<T>::from(value.0),
            max: Vector3::<T>::from(value.1),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let v1 = Vector3f::new(1.0, 2.0, 3.0);
        let v2 = Vector3f::new(4.0, 5.0, 6.0);
        let b1 = Bounds3f::new(&v1, &v2);
        assert_eq!(b1.min, v1);
        assert_eq!(b1.max, v2);
    }
}
