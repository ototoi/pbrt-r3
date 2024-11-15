use super::sobolmatrices::*;
use crate::core::pbrt::*;

#[inline]
pub fn sobol_interval_to_index(m: u32, frame: u64, p: &Point2i) -> u64 {
    let mut frame = frame;
    if m == 0 {
        return 0;
    }
    let m2 = m.wrapping_shl(1);
    let mut index = frame.wrapping_shl(m2);
    let mut delta: u64 = 0;
    let mut c = 0;
    while frame != 0 {
        if (frame & 1) != 0 {
            delta ^= VDC_SOBOL_MATRICES[(m - 1) as usize][c];
        }
        c += 1;
        frame = frame.wrapping_shr(1);
    }
    let mut b = ((((p.x as u32) as u64).wrapping_shl(m)) | (p.y as u64)) ^ delta;
    let mut c = 0;
    while b != 0 {
        if (b & 1) != 0 {
            index ^= VDC_SOBOL_MATRICES_INV[(m - 1) as usize][c];
        }
        c += 1;
        b = b.wrapping_shr(1);
    }
    return index;
}

#[inline]
pub fn sobol_sample(a: i64, dimension: u32, scramble: u32) -> Float {
    let mut a = a;
    let mut v = scramble;
    //let dimension = dimension % ((NUM_SOBOL_DIMENSIONS) as u32);
    let mut i = usize::min(
        dimension as usize * SOBOL_MATRIX_SIZE,
        SOBOL_MATRICES_32.len() - 1,
    );
    while a != 0 {
        if (a & 1) != 0 {
            v ^= SOBOL_MATRICES_32[i] as u32;
        }
        a >>= 1;
        i += 1;
        i %= SOBOL_MATRICES_32.len();
    }
    let fv = ((v as f64) * 2.3283064365386963e-10) as Float;
    return Float::min(fv, ONE_MINUS_EPSILON);
}
