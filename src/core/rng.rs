use crate::core::pbrt::*;

const PCG32_DEFAULT_STATE: u64 = 0x853c49e6748fea9b;
const PCG32_DEFAULT_STREAM: u64 = 0xda3e39cb94b95bdb;
const PCG32_MULT: u64 = 0x5851f42d4c957f2d;

#[derive(Debug, PartialEq, Clone)]
pub struct RNG {
    pub state: u64,
    pub inc: u64,
}

impl RNG {
    pub fn new() -> Self {
        RNG {
            state: PCG32_DEFAULT_STATE,
            inc: PCG32_DEFAULT_STREAM,
        }
    }

    pub fn new_sequence(initseq: u64) -> Self {
        let mut r = Self::new();
        r.set_sequence(initseq);
        return r;
    }

    pub fn set_sequence(&mut self, initseq: u64) {
        self.state = 0;
        self.inc = (initseq << 1) | 1;
        self.uniform_uint32();
        self.state += PCG32_DEFAULT_STATE;
        self.uniform_uint32();
    }

    #[inline]
    pub fn uniform_uint32(&mut self) -> u32 {
        let oldstate: u64 = self.state;
        //self.state = oldstate * PCG32_MULT + self.inc;
        self.state = oldstate.wrapping_mul(PCG32_MULT).wrapping_add(self.inc);
        //let xorshifted: u32 = (((oldstate >> 18) ^ oldstate) >> 27) as u32;
        let xorshifted: u32 = ((oldstate.wrapping_shr(18) ^ oldstate).wrapping_shr(27)) as u32;
        //let rot: u32 = (oldstate >> 59) as u32;
        let rot: u32 = (oldstate.wrapping_shr(59)) as u32;
        return (xorshifted.wrapping_shr(rot))
            | (xorshifted.wrapping_shl(((!rot).wrapping_add(1)) & 31));
    }

    pub fn uniform_uint32_threshold(&mut self, b: u32) -> u32 {
        let threshold = (!b + 1) % b;
        loop {
            let r = self.uniform_uint32();
            if r >= threshold {
                return r % b;
            }
        }
    }

    #[inline]
    pub fn uniform_float(&mut self) -> Float {
        return self.uniform_float32() as Float;
    }

    pub fn uniform_float32(&mut self) -> f32 {
        let f: f32 = self.uniform_uint32() as f32 * 2.3283064365386963e-10;
        return FLOAT_ONE_MINUS_EPSILON.min(f);
    }
}

impl Default for RNG {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let a: u32 = 0b10100000;
        let b: u32 = 0b01011111;
        assert_eq!((!a) & 0b11111111, b);
    }

    #[test]
    fn test_002() {
        let mut rng = RNG::new();
        let a: f32 = rng.uniform_float32();
        let astate = rng.state;
        let b: f32 = rng.uniform_float32();
        let bstate = rng.state;
        assert_ne!(a, b);
        assert_ne!(astate, bstate);
    }
}
