use crate::core::pbrt::*;

pub fn decompose(
    m: &Matrix4x4,
    epsilon: Float,
    max_count: i32,
) -> Option<(Vector3f, Quaternion, Matrix4x4)> {
    let t = Vector3f::new(m.m[3], m.m[7], m.m[11]);
    let mut mm = *m;
    for i in 0..3 {
        mm.m[4 * i + 3] = 0.0;
    }
    mm.m[15] = 1.0;
    let mut r = mm;
    let mut count = 0;
    let mut norm: Float = 0.0;
    loop {
        let mut r_next = Matrix4x4::identity();
        if let Some(r_it) = r.transpose().inverse() {
            for i in 0..4 {
                for j in 0..4 {
                    r_next.m[4 * i + j] = 0.5 * (r.m[4 * i + j] + r_it.m[4 * i + j]);
                }
            }
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

    if let Some(r_it) = r.inverse() {
        let q = Quaternion::from(r);
        let s = r_it * mm;
        return Some((t, q, s));
    }
    return None;
}
