use crate::core::base::Float;

pub fn xyz_to_rgb(xyz: &[Float]) -> [Float; 3] {
    let mut rgb: [Float; 3] = [0.0; 3];
    rgb[0] = 3.240479 * xyz[0] - 1.537150 * xyz[1] - 0.498535 * xyz[2];
    rgb[1] = -0.969256 * xyz[0] + 1.875991 * xyz[1] + 0.041556 * xyz[2];
    rgb[2] = 0.055648 * xyz[0] - 0.204043 * xyz[1] + 1.057311 * xyz[2];
    return rgb;
}

pub fn rgb_to_xyz(rgb: &[Float]) -> [Float; 3] {
    let mut xyz: [Float; 3] = [0.0; 3];
    xyz[0] = 0.412453 * rgb[0] + 0.357580 * rgb[1] + 0.180423 * rgb[2];
    xyz[1] = 0.212671 * rgb[0] + 0.715160 * rgb[1] + 0.072169 * rgb[2];
    xyz[2] = 0.019334 * rgb[0] + 0.119193 * rgb[1] + 0.950227 * rgb[2];
    return xyz;
}
