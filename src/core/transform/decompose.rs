use super::matrix4x4::*;
use crate::core::base::*;
use crate::core::quaternion::*;

fn invertible(m: &Matrix4x4) -> bool {
    m.inverse().is_some()
}

fn length(x: Float, y: Float, z: Float) -> Float {
    (x * x + y * y + z * z).sqrt()
}

fn suppress_for_scale(m: Matrix4x4) -> Matrix4x4 {
    let mut mm = Matrix4x4::identity();
    for i in 0..3 {
        mm.m[4 * i + i] = m.m[4 * i + i];
    }
    return mm;
}

// Decompose a 4x4 matrix into translation, rotation, and scale.
pub fn decompose(
    m: &Matrix4x4,
    epsilon: Float,
    max_count: i32,
) -> Option<(Vector3f, Quaternion, Vector3f)> {
    // Extract translation _T_ from transformation matrix
    let t = Vector3f::new(m.m[4 * 0 + 3], m.m[4 * 1 + 3], m.m[4 * 2 + 3]);

    // Compute new transformation matrix _M_ without translation
    let mut mm = *m;
    for i in 0..3 {
        mm.m[4 * i + 3] = 0.0;
    }
    mm.m[15] = 1.0;

    // Extract rotation _R_ from transformation matrix
    let mut r = mm;
    // pbrt-r3
    let sx = length(r.m[4 * 0 + 0], r.m[4 * 1 + 0], r.m[4 * 2 + 0]);
    let sy = length(r.m[4 * 0 + 1], r.m[4 * 1 + 1], r.m[4 * 2 + 1]);
    let sz = length(r.m[4 * 0 + 2], r.m[4 * 1 + 2], r.m[4 * 2 + 2]);
    if sx != 0.0 {
        r.m[4 * 0 + 0] /= sx;
        r.m[4 * 0 + 1] /= sx;
        r.m[4 * 0 + 2] /= sx;
    }
    if sy != 0.0 {
        r.m[4 * 1 + 0] /= sy;
        r.m[4 * 1 + 1] /= sy;
        r.m[4 * 1 + 2] /= sy;
    }
    if sz != 0.0 {
        r.m[4 * 2 + 0] /= sz;
        r.m[4 * 2 + 1] /= sz;
        r.m[4 * 2 + 2] /= sz;
    }
    // pbrt-r3
    let mut count = 0;
    let mut norm: Float = 0.0;
    loop {
        // Compute inverse of _R_ and check for singularity
        assert!(invertible(&r));
        if let Some(r_it) = r.transpose().inverse() {
            // Compute next matrix _Rnext_ in series
            let mut r_next = r;
            for i in 0..4 {
                for j in 0..4 {
                    r_next.m[4 * i + j] = lerp(0.5, r.m[4 * i + j], r_it.m[4 * i + j]);
                }
            }

            // pbrt-r3
            assert!(invertible(&r_next));
            if let Some(ir_next) = r_next.inverse() {
                let s = suppress_for_scale(ir_next * mm);
                if let Some(is) = s.inverse() {
                    r_next = mm * is;
                }
            }
            let q = Quaternion::from(r_next);
            let q = q.normalize();
            r_next = q.to_matrix();
            // pbrt-r3

            // Compute norm of difference between _R_ and _Rnext_
            for i in 0..3 {
                let n = Float::abs(r.m[4 * i + 0] - r_next.m[4 * i + 0])
                    + Float::abs(r.m[4 * i + 1] - r_next.m[4 * i + 1])
                    + Float::abs(r.m[4 * i + 2] - r_next.m[4 * i + 2]);
                norm = Float::max(norm, n);
            }
            r = r_next;
        } else {
            return None;
        }

        count += 1;
        if !(count < max_count && norm > epsilon) {
            break;
        }
    }

    if let Some(ir) = r.inverse() {
        let s = suppress_for_scale(ir * mm);
        if let Some(is) = s.inverse() {
            r = mm * is;
        }
        let q = Quaternion::from(r);
        let q = q.normalize(); //pbrt-r3
        let s = Vector3f::new(s.m[0], s.m[5], s.m[10]);
        return Some((t, q, s));
    }
    return None;
}
