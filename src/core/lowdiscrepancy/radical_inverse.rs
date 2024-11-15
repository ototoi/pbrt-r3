use super::primes::PRIMES;
use crate::core::pbrt::*;

/*
pub fn reverse_bits32(mut n: u32) -> u32 {
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    return n;
}

pub fn reverse_bits64(n: u64) -> u64 {
    let n0 = reverse_bits32(n as u32) as u64;
    let n1 = reverse_bits32((n >> 32) as u32) as u64;
    return (n0 << 32) | n1;
}
*/

fn reverse_bits32(mut n: u32) -> u32 {
    n = (n.wrapping_shl(16)) | (n.wrapping_shr(16));
    n = ((n & 0x00ff00ff).wrapping_shl(8)) | ((n & 0xff00ff00).wrapping_shr(8));
    n = ((n & 0x0f0f0f0f).wrapping_shl(4)) | ((n & 0xf0f0f0f0).wrapping_shr(4));
    n = ((n & 0x33333333).wrapping_shl(2)) | ((n & 0xcccccccc).wrapping_shr(2));
    n = ((n & 0x55555555).wrapping_shl(1)) | ((n & 0xaaaaaaaa).wrapping_shr(1));
    return n;
}

fn reverse_bits64(n: u64) -> u64 {
    let n0 = reverse_bits32(n as u32) as u64;
    let n1 = reverse_bits32((n.wrapping_shr(32)) as u32) as u64;
    return (n0.wrapping_shl(32)) | n1;
}

fn radical_inverse_specialized(base: u64, mut a: u64) -> Float {
    let inv_base = 1.0 / base as Float;
    let mut reversed_digits = 0;
    let mut inv_base_n = 1.0;
    while a != 0 {
        let next = a / base;
        let digit = a - next * base;
        reversed_digits = reversed_digits * base + digit;
        inv_base_n *= inv_base;
        a = next;
    }
    return Float::min(reversed_digits as Float * inv_base_n, ONE_MINUS_EPSILON);
}

pub fn radical_inverse(base_index: u32, a: u64) -> Float {
    assert!(base_index < 1024);
    return match base_index {
        0 => reverse_bits64(a) as Float * 5.4210108624275222e-20,
        _ => radical_inverse_specialized(PRIMES[base_index as usize], a),
    };
}

pub fn inverse_radical_inverse(base: u64, inverse: u64, ndigits: usize) -> u64 {
    let mut inverse = inverse;
    let mut index = 0;
    for _ in 0..ndigits {
        let digit = inverse % base;
        inverse /= base;
        index = index * base + digit;
    }
    return index;
}

fn scrambled_radical_inverse_specialized(base: u64, perm: &[u16], a: u64) -> Float {
    let mut a = a;
    let inv_base = Float::recip(base as Float);
    let mut reverse_digits = 0u64;
    let mut inv_base_n = 1.0;
    while a != 0 {
        let next = a / base;
        let digit = a - next * base;
        reverse_digits = reverse_digits * base + perm[digit as usize] as u64;
        inv_base_n *= inv_base;
        a = next;
    }
    return Float::min(
        inv_base_n * (reverse_digits as Float + inv_base * perm[0] as Float / (1.0 - inv_base)),
        ONE_MINUS_EPSILON,
    );
}

pub fn scrambled_radical_inverse(base_index: u32, a: u64, perm: &[u16]) -> Float {
    return scrambled_radical_inverse_specialized(PRIMES[base_index as usize], perm, a);
}

pub fn compute_radical_inverse_permutations(perms: &mut Vec<u16>, rng: &mut RNG) {
    // Allocate space in _perms_ for radical inverse permutations
    let mut perm_array_size = 0;
    for i in 0..PRIMES.len() {
        perm_array_size += PRIMES[i] as usize;
    }
    perms.resize(perm_array_size, 0);
    let mut offset = 0;
    for i in 0..PRIMES.len() {
        let len = PRIMES[i] as usize;
        for j in 0..len {
            perms[offset + j] = j as u16;
        }
        shuffle_array(&mut perms[offset..], len, 1, rng);
        offset += len;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let u1 = 1;
        let u2 = reverse_bits64(u1);
        let u3 = reverse_bits64(u2);
        assert_eq!(u1, u3);
    }
}
