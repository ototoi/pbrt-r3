use crate::core::pbrt::*;

#[derive(Debug, PartialEq, Default, Copy, Clone)]
pub struct Bounds2<T> {
    pub min: Vector2<T>,
    pub max: Vector2<T>,
}

pub type Bounds2f = Bounds2<f32>;
pub type Bounds2d = Bounds2<f64>;
pub type Bounds2i = Bounds2<i32>;

impl<T: Copy> Bounds2<T> {
    pub fn new(min: &Vector2<T>, max: &Vector2<T>) -> Self {
        Bounds2::<T> {
            min: *min,
            max: *max,
        }
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
    > Bounds2<T>
{
    pub fn area(&self) -> T {
        return (self.max.x - self.min.x) * (self.max.y - self.min.y);
    }
    pub fn diagonal(&self) -> Vector2<T> {
        return self.max - self.min;
    }
    pub fn offset(&self, p: &Vector2<T>) -> Vector2<T> {
        let mut o = *p - self.min;
        if self.max.x > self.min.x {
            o.x = o.x / (self.max.x - self.min.x);
        }
        if self.max.y > self.min.y {
            o.y = o.y / (self.max.y - self.min.y);
        }
        return o;
    }

    pub fn expand(&self, delta: T) -> Self {
        let delta = Vector2::<T>::new(delta, delta);
        return Bounds2 {
            min: self.min - delta,
            max: self.max + delta,
        };
    }

    pub fn union(&self, other: &Self) -> Self {
        let min = Vector2::<T>::new(min_(self.min.x, other.min.x), min_(self.min.y, other.min.y));
        let max = Vector2::<T>::new(max_(self.max.x, other.max.x), max_(self.max.y, other.max.y));
        return Bounds2 { min, max };
    }

    pub fn intersect(&self, other: &Self) -> Self {
        let min = Vector2::<T>::new(max_(self.min.x, other.min.x), max_(self.min.y, other.min.y));
        let max = Vector2::<T>::new(min_(self.max.x, other.max.x), min_(self.max.y, other.max.y));
        return Bounds2 { min, max };
    }

    pub fn union_p(&self, other: &Vector2<T>) -> Self {
        let min = Vector2::<T>::new(min_(self.min.x, other.x), min_(self.min.y, other.y));
        let max = Vector2::<T>::new(max_(self.max.x, other.x), max_(self.max.y, other.y));
        return Bounds2 { min, max };
    }

    pub fn inside(&self, p: &Vector2<T>) -> bool {
        return p.x >= self.min.x && p.x <= self.max.x && p.y >= self.min.y && p.y <= self.max.y;
    }

    pub fn inside_exclusive(&self, p: &Vector2<T>) -> bool {
        return p.x >= self.min.x && p.x < self.max.x && p.y >= self.min.y && p.y < self.max.y;
    }
}

impl Bounds2f {
    pub fn lerp(&self, t: &Vector2f) -> Vector2f {
        return Vector2f::new(
            lerp(t.x, self.min.x, self.max.x),
            lerp(t.y, self.min.y, self.max.y),
        );
    }
}

impl<T: Copy> From<((T, T), (T, T))> for Bounds2<T> {
    fn from(value: ((T, T), (T, T))) -> Self {
        Bounds2::<T> {
            min: Vector2::<T>::from(value.0),
            max: Vector2::<T>::from(value.1),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let v1 = Vector2f::new(1.0, 2.0);
        let v2 = Vector2f::new(4.0, 5.0);
        let b1 = Bounds2f::new(&v1, &v2);
        assert_eq!(b1.min, v1);
        assert_eq!(b1.max, v2);
    }
}
