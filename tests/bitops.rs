// Imported from bitops.cpp

use pbrt_r3::core::pbrt::*;

#[test]
fn log2_basics() {
    for i in 0..32 {
        let ui = 1 << i;
        assert_eq!(i, log2int(ui));
        assert_eq!(i, log2_int64(ui as u64));
    }

    for i in 1..31 {
        let ui = 1 << i;
        assert_eq!(i, log2int(ui + 1));
        assert_eq!(i, log2_int64((ui + 1) as u64));
    }

    for i in 0..64 {
        let ui = 1 << i;
        assert_eq!(i, log2_int64(ui));
    }
}

#[test]
fn pow2_basics() {
    for i in 0..32 {
        let ui = 1 << i;
        assert_eq!(true, is_power_of_2(ui));
        if ui > 1 {
            assert_eq!(false, is_power_of_2(ui + 1));
        }
        if ui > 2 {
            assert_eq!(false, is_power_of_2(ui - 1));
        }
    }
}

#[test]
fn count_trailing_zeros_basics() {
    for i in 0..32 {
        let ui = 1 << i;
        assert_eq!(i, count_trailing_zeros(ui));
    }
}

#[test]
fn round_up_pow2_basics() {
    for i in 0..32 {
        assert_eq!(8, round_up_pow2(7));
    }
}

