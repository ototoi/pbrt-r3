// Imported from sampling.cpp

use pbrt_r3::core::lowdiscrepancy::maxmin::*;
use pbrt_r3::core::lowdiscrepancy::primes::*;
use pbrt_r3::core::lowdiscrepancy::radical_inverse::*;
use pbrt_r3::core::lowdiscrepancy::*;
use pbrt_r3::core::prelude::*;
use pbrt_r3::samplers::*;

fn near_equal(a: Float, b: Float, e: Float) -> bool {
    (a - b).abs() < e
}

#[test]
fn lowdiscrepancy_radical_inverse() {
    for a in 0..1024 {
        assert_eq!(
            reverse_bits32(a) as Float * 2.3283064365386963e-10,
            radical_inverse(0, a as u64)
        );
    }
}

#[test]
fn lowdiscrepancy_scrambled_radical_inverse() {
    for dim in 0..128 {
        let mut rng = RNG::new_sequence(dim as u64);
        // Random permutation table
        let base = PRIMES[dim as usize] as u16;
        let mut perm: Vec<u16> = (0..base).map(|i| base - 1 - i).collect();
        shuffle_array(&mut perm, base as usize, 1, &mut rng);
        for index in [0, 1, 2, 1151, 32351, 4363211, 681122] {
            // First, compare to the pbrt-v2 implementation.
            {
                let mut val = 0.0;
                let inv_base = 1.0 / base as Float;
                let mut inv_bi = inv_base as Float;
                let mut n = index as u64;
                while n > 0 {
                    let d_i = perm[(n % base as u64) as usize];
                    val += d_i as Float * inv_bi;
                    n = (n as Float * inv_base) as u64;
                    inv_bi = inv_bi * inv_base;
                }
                // For the case where the permutation table permutes the digit 0
                // to
                // another digit, account for the infinite sequence of that
                // digit
                // trailing at the end of the radical inverse value.
                val += perm[0] as Float * (base as Float) / (base as Float - 1.0) * inv_bi;

                let scrampled_val = scrambled_radical_inverse(dim, index, &perm);
                assert!(
                    near_equal(val, scrampled_val, 1e-5),
                    "val: {}, scrampled_val:{}, dim: {}, index: {}",
                    val,
                    scrampled_val,
                    dim,
                    index
                );
            }
        }
    }
}

#[test]
fn lowdiscrepancy_generator_matrix() {
    // Identity matrix, column-wise
    let c = (0..32).map(|i| 1 << i).collect::<Vec<u32>>();
    let c_rev = (0..32).map(|i| reverse_bits32(c[i])).collect::<Vec<u32>>();

    for a in 0..128 {
        // Make sure identity generator matrix matches van der Corput
        assert_eq!(a, multiply_generator(&c, a));
        assert_eq!(
            radical_inverse(0, a as u64),
            reverse_bits32(multiply_generator(&c, a)) as Float * 2.3283064365386963e-10
        );
        assert_eq!(
            radical_inverse(0, a as u64),
            sample_generator_matrix(&c_rev, a, 0)
        );
    }

    // Random / goofball generator matrix
    let mut rng = RNG::new();
    let c = (0..32).map(|_| rng.uniform_uint32()).collect::<Vec<u32>>();
    let c_rev = (0..32).map(|i| reverse_bits32(c[i])).collect::<Vec<u32>>();
    for a in 0..1024 {
        assert_eq!(
            reverse_bits32(multiply_generator(&c, a)),
            multiply_generator(&c_rev, a)
        );
    }
}

#[test]
fn lowdiscrepancy_gray_code_sample() {
    // Identity matrix, column-wise
    let c = (0..32).map(|i| 1 << i).collect::<Vec<u32>>();

    let mut v = [0.0; 64];
    gray_code_sample_1d(&c, 64, 0, &mut v);
    for a in 0..64 {
        let u = multiply_generator(&c, a) as Float * 2.3283064365386963e-10;
        assert!(v.contains(&u));
    }
}

#[test]
fn lowdiscrepancy_sobol() {
    // Check that float and double variants match (as float values).
    for i in 0..256 {
        for dim in 0..100 {
            let val_f = sobol_sample_float(i, dim, 0) as Float;
            let val_d = sobol_sample_double(i, dim, 0) as Float;
            assert!(
                near_equal(val_f, val_d, 1e-5),
                "val_f: {}, val_d: {}, i: {}, dim: {}",
                val_f,
                val_d,
                i,
                dim
            );
        }
    }

    // Make sure first dimension is the regular base 2 radical inverse
    for i in 0..8192 {
        assert_eq!(
            sobol_sample_float(i as i64, 0, 0),
            reverse_bits32(i as u32) as f32 * 2.3283064365386963e-10
        );
    }
}

fn check_sampler(name: &str, sampler: &mut dyn Sampler, log_samples: usize) {
    // Get all of the samples for a pixel.
    sampler.start_pixel(&Point2i::new(0, 0));
    let mut samples = Vec::new();
    loop {
        samples.push(sampler.get_2d());
        if !sampler.start_next_sample() {
            break;
        }
    }

    for i in 0..log_samples {
        // Check one set of elementary intervals: number of intervals
        // in each dimension.
        let nx = 1 << i;
        let ny = 1 << (log_samples - i);

        let mut count = vec![0; (1 << log_samples) as usize];
        for s in samples.iter() {
            // Map the sample to an interval
            let x = nx as Float * s.x;
            let y = ny as Float * s.y;
            assert!(x >= 0.0);
            assert!(x < nx as Float);
            assert!(y >= 0.0);
            assert!(y < ny as Float);
            let xx = x.floor() as i32;
            let yy = y.floor() as i32;
            let index = yy * nx + xx;
            assert!(index >= 0);
            assert!(index < count.len() as i32);
            let index = index as usize;

            // This should be the first time a sample has landed in its
            // interval.
            assert_eq!(
                count[index], 0,
                "Sampler {}: count[{}] = {}",
                name, index, count[index]
            );
            count[index] += 1;
        }
    }
}
// Make sure samplers that are supposed to generate a single sample in
// each of the elementary intervals actually do so.
// TODO: check Halton (where the elementary intervals are (2^i, 3^j)).
#[test]
fn lowdiscrepancy_elementary_intervals() {
    for log_samples in 2..=10 {
        {
            let mut sampler = MaxMinDistSampler::new(1 << log_samples, 2);
            check_sampler("MaxMinDistSampler", &mut sampler, log_samples);
        }
        {
            let mut sampler = ZeroTwoSequenceSampler::new(1 << log_samples, 2);
            check_sampler("ZeroTwoSequenceSampler", &mut sampler, log_samples);
        }
        {
            let bounds = Bounds2i::from(((0, 0), (10, 10)));
            let mut sampler = SobolSampler::new(1 << log_samples, &bounds);
            check_sampler("SobolSampler", &mut sampler, log_samples);
        }
    }
}

#[test]
fn max_min_dist_min_dist() {
    // Expected minimum distances from Gruenschloss et al.'s paper.
    const EXPECTED_MIN_DIST: [Float; 17] = [
        0.0, /* not checked */
        0.0, /* not checked */
        0.35355, 0.35355, 0.22534, 0.16829, 0.11267, 0.07812, 0.05644, 0.03906, 0.02816, 0.01953,
        0.01408, 0.00975, 0.00704, 0.00486, 0.00352,
    ];

    // We use a silly O(n^2) distance check below, so don't go all the way up
    // to 2^16 samples.
    for log_samples in 2..=10 {
        // Store a pixel's worth of samples in the vector s.
        let mut mm = MaxMinDistSampler::new(1 << log_samples, 2);
        mm.start_pixel(&Point2i::new(0, 0));
        let mut s = Vec::new();
        loop {
            s.push(mm.get_2d());
            if !mm.start_next_sample() {
                break;
            }
        }

        // Distance with toroidal topology
        let dist = |p0: &Point2f, p1: &Point2f| -> Float {
            let mut d = (*p1 - *p0).abs();
            if d.x > 0.5 {
                d.x = 1.0 - d.x;
            }
            if d.y > 0.5 {
                d.y = 1.0 - d.y;
            }
            return d.length();
        };

        let mut min_dist = Float::INFINITY;
        for i in 0..s.len() {
            for j in 0..s.len() {
                if i == j {
                    continue;
                }
                min_dist = min_dist.min(dist(&s[i], &s[j]));
            }
        }

        // Increase the tolerance by a small slop factor.
        assert!(min_dist > 0.99 * EXPECTED_MIN_DIST[log_samples as usize]);
    }
}

#[test]
fn distribution_1d_discrete() {
    // Carefully chosen distribution so that transitions line up with
    // (inverse) powers of 2.
    let func = [0.0, 1.0, 0.0, 3.0];
    let dist = Distribution1D::new(&func);
    assert_eq!(4, dist.func.len());

    assert_eq!(0.0, dist.discrete_pdf(0));
    assert_eq!(0.25, dist.discrete_pdf(1));
    assert_eq!(0.0, dist.discrete_pdf(2));
    assert_eq!(0.75, dist.discrete_pdf(3));

    let (offset, pdf, _u_remapped) = dist.sample_discrete(0.0);
    assert_eq!(1, offset);
    assert_eq!(0.25, pdf);

    let (offset, pdf, u_remapped) = dist.sample_discrete(0.125);
    assert_eq!(1, offset);
    assert_eq!(0.25, pdf);
    assert_eq!(0.5, u_remapped);

    let (offset, pdf, _u_remapped) = dist.sample_discrete(0.24999);
    assert_eq!(1, offset);
    assert_eq!(0.25, pdf);

    let (offset, pdf, _u_remapped) = dist.sample_discrete(0.250001);
    assert_eq!(3, offset);
    assert_eq!(0.75, pdf);

    let (offset, pdf, u_remapped) = dist.sample_discrete(0.625);
    assert_eq!(3, offset);
    assert_eq!(0.75, pdf);
    assert_eq!(0.5, u_remapped);

    let (offset, pdf, _u_remapped) = dist.sample_discrete(ONE_MINUS_EPSILON);
    assert_eq!(3, offset);
    assert_eq!(0.75, pdf);

    let (offset, pdf, _u_remapped) = dist.sample_discrete(1.0);
    assert_eq!(3, offset);
    assert_eq!(0.75, pdf);

    // Compute the interval to test over.
    let mut u = 0.25;
    let mut u_max = 0.25;
    for _ in 0..20 {
        u = next_float_down(u);
        u_max = next_float_up(u_max);
    }
    // We should get a stream of hits in the first interval, up until the
    // cross-over point at 0.25 (plus/minus fp slop).
    while u < u_max {
        let (interval, _pdf, _u_remapped) = dist.sample_discrete(u);
        if interval == 3 {
            break;
        }
        assert_eq!(1, interval);
        u = next_float_up(u);
    }
    assert!(u < u_max);
    // And then all the rest should be in the third interval
    while u <= u_max {
        let (interval, _pdf, _u_remapped) = dist.sample_discrete(u);
        assert_eq!(3, interval);
        u = next_float_up(u);
    }
}

#[test]
fn distribution_1d_continuous() {
    let func = [1.0, 1.0, 2.0, 4.0, 8.0];
    let dist = Distribution1D::new(&func);
    assert_eq!(5, dist.count());

    let (value, pdf, offset) = dist.sample_continuous(0.0);
    assert_eq!(0.0, value);
    assert_eq!(dist.count() as Float * 1.0 / 16.0, pdf);
    assert_eq!(0, offset);

    // Right at the bounary between the 4 and the 8 segments.
    let (value, _pdf, _offset) = dist.sample_continuous(0.5);
    assert_eq!(0.8, value);

    // Middle of the 8 segment
    let (value, pdf, offset) = dist.sample_continuous(0.75);
    assert_eq!(0.9, value);
    assert_eq!(dist.count() as Float * 8.0 / 16.0, pdf);
    assert_eq!(4, offset);

    let (value, _pdf, _offset) = dist.sample_continuous(0.0);
    assert_eq!(0.0, value);

    let (value, _pdf, _offset) = dist.sample_continuous(1.0);
    assert_eq!(1.0, value);
}
