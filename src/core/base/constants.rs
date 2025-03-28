use super::types::Float;
/*
static PBRT_CONSTEXPR Float ShadowEpsilon = 0.0001f;
static PBRT_CONSTEXPR Float Pi = 3.14159265358979323846;
static PBRT_CONSTEXPR Float InvPi = 0.31830988618379067154;
static PBRT_CONSTEXPR Float Inv2Pi = 0.15915494309189533577;
static PBRT_CONSTEXPR Float Inv4Pi = 0.07957747154594766788;
static PBRT_CONSTEXPR Float PiOver2 = 1.57079632679489661923;
static PBRT_CONSTEXPR Float PiOver4 = 0.78539816339744830961;
static PBRT_CONSTEXPR Float Sqrt2 = 1.41421356237309504880;
*/
pub const DOUBLE_ONE_MINUS_EPSILON: f64 = 0.99999999999999989;
pub const FLOAT_ONE_MINUS_EPSILON: f32 = 0.99999994;

#[cfg(not(feature = "float-as-double"))]
mod detail {
    use super::*;

    pub const SHADOW_EPSILON: Float = 0.0001;
    pub const PI: Float = std::f32::consts::PI; //3.14159265358979323846;
    pub const INV_PI: Float = std::f32::consts::FRAC_1_PI; //0.31830988618379067154;
    pub const INV_2_PI: Float = INV_PI * 0.5;
    pub const INV_4_PI: Float = INV_PI * 0.25;

    pub const INV_SQRT_PI: Float = 0.5 * std::f32::consts::FRAC_2_SQRT_PI; //1 / sqrt(pi)

    pub const PI_OVER_2: Float = PI / 2.0; //1.57079632679489661923
    pub const PI_OVER_4: Float = PI / 4.0; //0.78539816339744830961

    pub const SQRT_2: Float = std::f32::consts::SQRT_2;

    pub const ONE_MINUS_EPSILON: f32 = FLOAT_ONE_MINUS_EPSILON;
}

#[cfg(feature = "float-as-double")]
mod detail {
    use super::*;

    pub const SHADOW_EPSILON: Float = 0.0001;
    pub const PI: Float = std::f64::consts::PI; //3.14159265358979323846;
    pub const INV_PI: Float = std::f64::consts::FRAC_1_PI; //0.31830988618379067154;
    pub const INV_2_PI: Float = INV_PI * 0.5;
    pub const INV_4_PI: Float = INV_PI * 0.25;

    pub const INV_SQRT_PI: Float = 0.5 * std::f64::consts::FRAC_2_SQRT_PI; //1 / sqrt(pi)

    pub const PI_OVER_2: Float = PI / 2.0; //1.57079632679489661923
    pub const PI_OVER_4: Float = PI / 4.0; //0.78539816339744830961

    pub const SQRT_2: Float = std::f64::consts::SQRT_2;

    pub const ONE_MINUS_EPSILON: f64 = DOUBLE_ONE_MINUS_EPSILON;
}

pub use detail::*;
