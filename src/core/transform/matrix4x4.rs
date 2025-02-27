use crate::core::pbrt::*;
use std::ops;

#[inline]
fn radians(x: Float) -> Float {
    return x * (PI / 180.0);
}

pub fn solve_linear_system_2x2(a: &[[Float; 2]; 2], b: &[Float; 2]) -> Option<(Float, Float)> {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if Float::abs(det) < 1e-10 {
        return None;
    }
    let x0 = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
    let x1 = (a[0][0] * b[1] - a[1][0] * b[0]) / det;
    if x0.is_nan() || x1.is_nan() {
        return None;
    }
    return Some((x0, x1));
}

#[derive(Debug, PartialEq, Default, Copy, Clone)]
pub struct Matrix4x4 {
    pub m: [Float; 16],
}

impl Matrix4x4 {
    pub fn new(
        e0: Float,
        e1: Float,
        e2: Float,
        e3: Float,
        e4: Float,
        e5: Float,
        e6: Float,
        e7: Float,
        e8: Float,
        e9: Float,
        e10: Float,
        e11: Float,
        e12: Float,
        e13: Float,
        e14: Float,
        e15: Float,
    ) -> Self {
        Matrix4x4 {
            m: [
                e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15,
            ],
        }
    }

    pub fn identity() -> Self {
        Matrix4x4 {
            m: [
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        }
    }

    pub fn translate(x: Float, y: Float, z: Float) -> Self {
        return Matrix4x4::new(
            1.0, 0.0, 0.0, x, 0.0, 1.0, 0.0, y, 0.0, 0.0, 1.0, z, 0.0, 0.0, 0.0, 1.0,
        );
    }

    pub fn rotate_x(theta: Float) -> Self {
        let s = Float::sin(radians(theta));
        let c = Float::cos(radians(theta));
        Matrix4x4 {
            m: [
                1.0, 0.0, 0.0, 0.0, 0.0, c, -s, 0.0, 0.0, s, c, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        }
    }

    pub fn rotate_y(theta: Float) -> Self {
        let s = Float::sin(radians(theta));
        let c = Float::cos(radians(theta));
        Matrix4x4 {
            m: [
                c, 0.0, s, 0.0, 0.0, 1.0, 0.0, 0.0, -s, 0.0, c, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        }
    }

    pub fn rotate_z(theta: Float) -> Self {
        let s = Float::sin(radians(theta));
        let c = Float::cos(radians(theta));
        Matrix4x4 {
            m: [
                c, -s, 0.0, 0.0, s, c, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        }
    }

    pub fn rotate(theta: Float, x: Float, y: Float, z: Float) -> Self {
        let a = Vector3f::new(x, y, z).normalize();
        let s = Float::sin(radians(theta));
        let c = Float::cos(radians(theta));

        let _m00 = a.x * a.x + (1.0 - a.x * a.x) * c;
        let _m01 = a.x * a.y + (1.0 - a.x * a.x) * c;
        /*
        m.m[0][2] = a.x * a.y * (1 - c) + a.y * s;
        m.m[0][3] = 0;
        m.m[0][1] = x;
        m.m[0][2] = a.x * a.z * (1 - c) + a.y * s;
        m.m[0][3] = 0;
        */
        Matrix4x4 {
            m: [
                c, -s, 0.0, 0.0, s, c, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        }
    }

    pub fn scale(x: Float, y: Float, z: Float) -> Self {
        Matrix4x4 {
            m: [
                x, 0.0, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, 0.0, z, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        }
    }

    pub fn look_at(
        ex: Float,
        ey: Float,
        ez: Float,
        lx: Float,
        ly: Float,
        lz: Float,
        ux: Float,
        uy: Float,
        uz: Float,
    ) -> Self {
        return Self::camera_to_world(ex, ey, ez, lx, ly, lz, ux, uy, uz)
            .inverse()
            .unwrap();
    }

    pub fn camera_to_world(
        ex: Float,
        ey: Float,
        ez: Float,
        lx: Float,
        ly: Float,
        lz: Float,
        ux: Float,
        uy: Float,
        uz: Float,
    ) -> Self {
        let pos = Vector3f::new(ex, ey, ez);
        let look = Vector3f::new(lx, ly, lz);
        let up = Vector3f::new(ux, uy, uz).normalize();

        /*
        camera_to_world[3] = pos.x;
        camera_to_world[7] = pos.y;
        camera_to_world[11] = pos.z;
        camera_to_world[15] = 1.0;
        */

        let dir = (look - pos).normalize(); //z
                                            //println!("z:{:?}", dir);

        let mut right = Vector3f::cross(&up, &dir); //x
                                                    //println!("x:{:?}", right);
        assert!(right.length() != 0.0);
        right = right.normalize();
        let new_up = Vector3f::cross(&dir, &right).normalize();
        /*
        camera_to_world[0] = right.x;
        camera_to_world[4] = right.y;
        camera_to_world[8] = right.z;
        camera_to_world[12] = 0.;
        camera_to_world[1] = new_up.x;
        camera_to_world[5] = new_up.y;
        camera_to_world[9] = new_up.z;
        camera_to_world[13] = 0.0;
        camera_to_world[2] = dir.x;
        camera_to_world[6] = dir.y;
        camera_to_world[10] = dir.z;
        camera_to_world[14] = 0.;
        */
        let m: [Float; 16] = [
            right.x, new_up.x, dir.x, pos.x, //
            right.y, new_up.y, dir.y, pos.y, //
            right.z, new_up.z, dir.z, pos.z, //
            0.0, 0.0, 0.0, 1.0,
        ];
        Matrix4x4 { m }
    }

    pub fn transpose(&self) -> Self {
        Matrix4x4 {
            m: [
                self.m[0], self.m[4], self.m[8], self.m[12], self.m[1], self.m[5], self.m[9],
                self.m[13], self.m[2], self.m[6], self.m[10], self.m[14], self.m[3], self.m[7],
                self.m[11], self.m[15],
            ],
        }
    }

    pub fn inverse(&self) -> Option<Self> {
        let mut indxc = [0; 4];
        let mut indxr = [0; 4];
        let mut ipiv = [0; 4];
        let mut minv: [Float; 16] = self.m;
        for i in 0..4 {
            let mut irow = 0;
            let mut icol = 0;
            let mut big: Float = 0.0;
            // Choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for k in 0..4 {
                        if ipiv[k] == 0 {
                            if Float::abs(minv[4 * j + k]) >= big {
                                big = Float::abs(minv[4 * j + k]);
                                irow = j;
                                icol = k;
                            }
                        } else if ipiv[k] > 1 {
                            //Error("Singular matrix in MatrixInvert");
                            return None;
                        }
                    }
                }
            }
            ipiv[icol] += 1;
            if irow != icol {
                for k in 0..4 {
                    //let tmp = minv[4 * irow + k];
                    //minv[4 * irow + k] = minv[4 * icol + k];
                    //minv[4 * icol + k] = tmp;
                    minv.swap(4 * irow + k, 4 * icol + k);
                }
            }
            indxr[i] = irow;
            indxc[i] = icol;
            if minv[4 * icol + icol] == 0.0 {
                //Error("Singular matrix in MatrixInvert");
                return None;
            }

            // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
            let pivinv = 1.0 / minv[4 * icol + icol];
            minv[4 * icol + icol] = 1.0;
            for j in 0..4 {
                minv[4 * icol + j] *= pivinv;
            }

            // Subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save = minv[4 * j + icol];
                    minv[4 * j + icol] = 0.0;
                    for k in 0..4 {
                        minv[4 * j + k] -= minv[4 * icol + k] * save;
                    }
                }
            }
        }

        // Swap columns to reflect permutation
        for j in [3, 2, 1, 0] {
            if indxr[j] != indxc[j] {
                for k in 0..4 {
                    let src = 4 * k + indxr[j];
                    let dst = 4 * k + indxc[j];
                    minv.swap(src, dst);
                }
            }
        }

        return Some(Matrix4x4 { m: minv });
    }

    pub fn transform_point(&self, p: &Point3f) -> Point3f {
        let x = p.x;
        let y = p.y;
        let z = p.z;
        let xp = self.m[0] * x + self.m[1] * y + self.m[2] * z + self.m[3];
        let yp = self.m[4] * x + self.m[5] * y + self.m[6] * z + self.m[7];
        let zp = self.m[8] * x + self.m[9] * y + self.m[10] * z + self.m[11];
        let wp = self.m[12] * x + self.m[13] * y + self.m[14] * z + self.m[15];
        if wp == 1.0 {
            return Point3f::new(xp, yp, zp);
        } else {
            return Point3f::new(xp / wp, yp / wp, zp / wp);
        }
    }

    pub fn transform_vector(&self, p: &Vector3f) -> Vector3f {
        let x = p.x;
        let y = p.y;
        let z = p.z;
        let xp = self.m[0] * x + self.m[1] * y + self.m[2] * z;
        let yp = self.m[4] * x + self.m[5] * y + self.m[6] * z;
        let zp = self.m[8] * x + self.m[9] * y + self.m[10] * z;
        return Vector3f::new(xp, yp, zp);
    }

    pub fn transform_normal(&self, p: &Normal3f) -> Normal3f {
        let x = p.x;
        let y = p.y;
        let z = p.z;
        let xp = self.m[0] * x + self.m[4] * y + self.m[8] * z;
        let yp = self.m[1] * x + self.m[5] * y + self.m[9] * z;
        let zp = self.m[2] * x + self.m[6] * y + self.m[10] * z;
        return Normal3f::new(xp, yp, zp);
    }
}

fn mul4x4(a: &[Float], b: &[Float]) -> Float {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

impl ops::Mul<Matrix4x4> for Matrix4x4 {
    type Output = Matrix4x4;
    fn mul(self, rhs: Matrix4x4) -> Matrix4x4 {
        Matrix4x4::new(
            mul4x4(&self.m[0..4], &[rhs.m[0], rhs.m[4], rhs.m[8], rhs.m[12]]),
            mul4x4(&self.m[0..4], &[rhs.m[1], rhs.m[5], rhs.m[9], rhs.m[13]]),
            mul4x4(&self.m[0..4], &[rhs.m[2], rhs.m[6], rhs.m[10], rhs.m[14]]),
            mul4x4(&self.m[0..4], &[rhs.m[3], rhs.m[7], rhs.m[11], rhs.m[15]]),
            mul4x4(&self.m[4..8], &[rhs.m[0], rhs.m[4], rhs.m[8], rhs.m[12]]),
            mul4x4(&self.m[4..8], &[rhs.m[1], rhs.m[5], rhs.m[9], rhs.m[13]]),
            mul4x4(&self.m[4..8], &[rhs.m[2], rhs.m[6], rhs.m[10], rhs.m[14]]),
            mul4x4(&self.m[4..8], &[rhs.m[3], rhs.m[7], rhs.m[11], rhs.m[15]]),
            mul4x4(&self.m[8..12], &[rhs.m[0], rhs.m[4], rhs.m[8], rhs.m[12]]),
            mul4x4(&self.m[8..12], &[rhs.m[1], rhs.m[5], rhs.m[9], rhs.m[13]]),
            mul4x4(&self.m[8..12], &[rhs.m[2], rhs.m[6], rhs.m[10], rhs.m[14]]),
            mul4x4(&self.m[8..12], &[rhs.m[3], rhs.m[7], rhs.m[11], rhs.m[15]]),
            mul4x4(&self.m[12..16], &[rhs.m[0], rhs.m[4], rhs.m[8], rhs.m[12]]),
            mul4x4(&self.m[12..16], &[rhs.m[1], rhs.m[5], rhs.m[9], rhs.m[13]]),
            mul4x4(&self.m[12..16], &[rhs.m[2], rhs.m[6], rhs.m[10], rhs.m[14]]),
            mul4x4(&self.m[12..16], &[rhs.m[3], rhs.m[7], rhs.m[11], rhs.m[15]]),
        )
    }
}

impl From<[Float; 16]> for Matrix4x4 {
    fn from(v: [Float; 16]) -> Self {
        Matrix4x4 { m: v }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let m1 = Matrix4x4::scale(4.0, 4.0, 4.0);
        let m2 = m1.inverse().unwrap();
        let m3 = Matrix4x4::scale(0.25, 0.25, 0.25);
        assert_eq!(m2, m3);
    }

    #[test]
    fn test_002() {
        let m1 = Matrix4x4::translate(4.0, 4.0, 4.0);
        let m2 = m1.inverse().unwrap();
        let m3 = Matrix4x4::translate(-4.0, -4.0, -4.0);
        assert_eq!(m2, m3);
    }

    #[test]
    fn test_003() {
        let m1 = Matrix4x4::translate(4.0, 4.0, 4.0);
        let m2 = Matrix4x4::inverse(&m1).unwrap();
        let m3 = Matrix4x4::translate(-4.0, -4.0, -4.0);
        assert_eq!(m2, m3);
    }

    #[test]
    fn test_004() {
        let m1 = Matrix4x4::rotate_x(90.0);
        let m2 = Matrix4x4::inverse(&m1).unwrap();
        let m3 = Matrix4x4::rotate_x(-90.0);
        assert_eq!(m2, m3);
    }

    #[test]
    fn test_005() {
        let m1 = Matrix4x4::rotate_y(90.0);
        let m2 = Matrix4x4::inverse(&m1).unwrap();
        let m3 = Matrix4x4::rotate_y(-90.0);
        assert_eq!(m2, m3);
    }

    #[test]
    fn test_006() {
        let m1 = Matrix4x4::rotate_z(90.0);
        let m2 = Matrix4x4::inverse(&m1).unwrap();
        let m3 = Matrix4x4::rotate_z(-90.0);
        assert_eq!(m2, m3);
    }

    #[test]
    fn test_007() {
        let m1 = Matrix4x4::rotate_z(90.0);
        let m2 = Matrix4x4::transpose(&m1);
        let m3 = Matrix4x4::rotate_z(-90.0);
        assert_eq!(m2, m3);
    }

    #[test]
    fn test_009() {
        let v1 = Vector3f::new(0.0, 0.0, 0.0);
        let m1 = Matrix4x4::look_at(0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0); //world2camera
        let m2 = m1.inverse().unwrap(); //camera2world
                                        /*
                                        -1.0, 0.0, 0.0, 0.0,
                                         0.0, 0.7071068, -0.70710677, -4.214685e-8,
                                        -0.0, -0.70710677, -0.7071068, 1.4142135,
                                         0.0, 0.0, 0.0, 1.0
                                        */
        let v2 = m2.transform_point(&v1);

        assert!(Vector3f::distance_squared(&v2, &Vector3f::new(0.0, 1.0, 1.0)) < 0.01);
    }
}
