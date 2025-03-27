use crate::core::base::*;

// To convert matrices from the "(t, m, s)-Nets and Maximized Minimum
// Distance" paper:
// - Recall that each uint32_t represents a column of the matrix, so
//   copy columns of bits.
// - Reverse the bits in each column (so we don't need to reverse the
//   result after the matrix multiply.)
pub const CMAXMIN_DIST: [[u32; 32]; 17] = [
    [
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x80000000,
    ],
    [
        0xc0000000, 0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0xa0000000, 0x40000000, 0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0xc0000000, 0x50000000, 0x20000000, 0x30000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x88000000, 0x58000000, 0x20000000, 0x40000000, 0x80000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0xe0000000, 0x60000000, 0x28000000, 0x10000000, 0x18000000, 0x04000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x82000000, 0x44000000, 0x2c000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x70000000, 0x30000000, 0x14000000, 0x08000000, 0x0c000000, 0x02000000,
        0x01000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0xc0000000, 0x41000000, 0x22000000, 0x16000000, 0x08000000, 0x10000000, 0x20000000,
        0x40800000, 0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x18000000, 0x08000000, 0x1c000000, 0x1e000000,
        0x03000000, 0x00800000, 0x00400000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x60400000, 0x20800000, 0x11000000, 0x0b000000, 0x04000000, 0x08000000,
        0x10000000, 0x20000000, 0x40000000, 0x00200000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x1c000000, 0x0c000000, 0x05000000, 0x02000000,
        0x03000000, 0x00800000, 0x00400000, 0x00200000, 0x00100000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x30200000, 0x10400000, 0x08800000, 0x05800000, 0x02000000,
        0x04000000, 0x08000000, 0x10000000, 0x20000000, 0x00100000, 0x00080000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x0e000000, 0x06000000, 0x02800000,
        0x01000000, 0x01800000, 0x00400000, 0x00200000, 0x00100000, 0x00080000, 0x00040000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x18100000, 0x08200000, 0x04400000, 0x02c00000,
        0x01000000, 0x02000000, 0x04000000, 0x08000000, 0x10000000, 0x00080000, 0x00040000,
        0x00020000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x08000000, 0x07000000, 0x03000000,
        0x01400000, 0x00800000, 0x00c00000, 0x00200000, 0x00100000, 0x00080000, 0x00040000,
        0x00020000, 0x00010000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
];

#[inline]
pub fn multiply_generator(c: &[u32], a: u32) -> u32 {
    let mut v = 0;
    let mut a = a;
    let mut i = 0;
    while a != 0 {
        if (a & 1) != 0 {
            v ^= c[i];
        }
        i += 1;
        a >>= 1;
    }
    return v;
}

#[inline]
pub fn sample_generator_matrix(c: &[u32], a: u32, scramble: u32) -> Float {
    return Float::min(
        (multiply_generator(c, a) ^ scramble) as Float * 2.3283064365386963e-10,
        ONE_MINUS_EPSILON,
    );
}
