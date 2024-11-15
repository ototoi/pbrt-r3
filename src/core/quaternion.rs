use crate::core::pbrt::*;
use std::ops;

#[derive(Debug, PartialEq, Default, Copy, Clone)]
pub struct Quaternion {
    pub x: Float,
    pub y: Float,
    pub z: Float,
    pub w: Float,
}

impl Quaternion {
    pub fn new(x: Float, y: Float, z: Float, w: Float) -> Self {
        Quaternion { x, y, z, w }
    }

    pub fn identity() -> Self {
        Quaternion {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 1.0,
        }
    }

    pub fn normalize(&self) -> Self {
        let l = Float::sqrt(Quaternion::dot(self, self));
        return Quaternion::new(self.x / l, self.y / l, self.z / l, self.w / l);
    }

    pub fn dot(q1: &Quaternion, q2: &Quaternion) -> Float {
        return (q1.x * q2.x) + (q1.y * q2.y) + (q1.z * q2.z) + (q1.w * q2.w);
    }

    pub fn slerp(t: Float, q1: &Quaternion, q2: &Quaternion) -> Quaternion {
        const T: Float = 1.0 - Float::EPSILON;

        let c = Self::dot(q1, q2);
        if c >= T {
            return (*q1) * (1.0 - t) + (*q2) * t;
        } else {
            let theta = Float::acos(c);
            let s = Float::recip(Float::sin(theta));
            return ((*q1) * Float::sin((1.0 - t) * theta) + (*q2) * Float::sin(t * theta)) * s;
        }
    }

    pub fn to_matrix(&self) -> Matrix4x4 {
        let x = self.x;
        let y = self.y;
        let z = self.z;
        let w = self.w;

        let xx = x * x;
        let yy = y * y;
        let zz = z * z;
        let xy = x * y;
        let xz = x * z;
        let yz = y * z;
        let wx = x * w;
        let wy = y * w;
        let wz = z * w;

        let mut m = Matrix4x4::identity();
        m.m[4 * 0 + 0] = 1.0 - 2.0 * (yy + zz);
        m.m[4 * 0 + 1] = 2.0 * (xy + wz);
        m.m[4 * 0 + 2] = 2.0 * (xz - wy);
        m.m[4 * 1 + 0] = 2.0 * (xy - wz);
        m.m[4 * 1 + 1] = 1.0 - 2.0 * (xx + zz);
        m.m[4 * 1 + 2] = 2.0 * (yz + wx);
        m.m[4 * 2 + 0] = 2.0 * (xz + wy);
        m.m[4 * 2 + 1] = 2.0 * (yz - wx);
        m.m[4 * 2 + 2] = 1.0 - 2.0 * (xx + yy);
        return m.transpose();
    }
}

impl ops::Add<Quaternion> for Quaternion {
    type Output = Quaternion;
    fn add(self, rhs: Quaternion) -> Quaternion {
        Quaternion {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            w: self.w + rhs.w,
        }
    }
}

impl ops::Sub<Quaternion> for Quaternion {
    type Output = Quaternion;
    fn sub(self, rhs: Quaternion) -> Quaternion {
        Quaternion {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            w: self.w - rhs.w,
        }
    }
}

impl ops::Mul<Quaternion> for Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: Quaternion) -> Quaternion {
        Quaternion {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
            w: self.w * rhs.w,
        }
    }
}

impl ops::Mul<Float> for Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: Float) -> Quaternion {
        Quaternion {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
            w: self.w * rhs,
        }
    }
}

impl ops::Neg for Quaternion {
    type Output = Quaternion;
    fn neg(self) -> Quaternion {
        return Quaternion {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        };
    }
}

impl From<Matrix4x4> for Quaternion {
    fn from(m: Matrix4x4) -> Self {
        let trace = m.m[0] + m.m[5] + m.m[10];
        if trace > 0.0 {
            // Compute w from matrix trace, then xyz
            // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
            let s = Float::sqrt(trace + 1.0);
            let w = s / 2.0;
            let s2 = 0.5 / s;
            let x = (m.m[4 * 2 + 1] - m.m[4 * 1 + 2]) * s2; //21 12
            let y = (m.m[4 * 0 + 2] - m.m[4 * 2 + 0]) * s2; //02 20
            let z = (m.m[4 * 1 + 0] - m.m[4 * 0 + 1]) * s2; //10 01
            return Quaternion::new(x, y, z, w);
        } else {
            // Compute largest of $x$, $y$, or $z$, then remaining components
            let nxt = [1, 2, 0];
            let mut q = [0.0; 3];
            let mut i = 0;
            if m.m[4 * 1 + 1] > m.m[4 * 0 + 0] {
                i = 1;
            }
            if m.m[4 * 2 + 2] > m.m[4 * i + i] {
                i = 2;
            }

            let j = nxt[i];
            let k = nxt[j];
            let mut s = Float::sqrt((m.m[4 * i + i] - (m.m[4 * j + j] + m.m[4 * k + k])) + 1.0);
            q[i] = s * 0.5;
            if s != 0.0 {
                s = 0.5 / s;
            }
            let w = (m.m[4 * k + j] - m.m[4 * j + k]) * s;
            q[j] = (m.m[4 * j + i] + m.m[4 * i + j]) * s;
            q[k] = (m.m[4 * k + i] + m.m[4 + i + k]) * s;
            let x = q[0];
            let y = q[1];
            let z = q[2];
            return Quaternion::new(x, y, z, w);
        }
    }
}
