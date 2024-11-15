use super::beam_diffusion::BSSRDFTable;
use crate::core::{interpolation::invert_catmull_rom, pbrt::*};

// BSSRDF Utility Functions
pub fn fresnel_moment1(eta: Float) -> Float {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let eta5 = eta4 * eta;
    if eta < 1.0 {
        return 0.45966 - 1.73965 * eta + 3.37668 * eta2 - 3.904945 * eta3 + 2.49277 * eta4
            - 0.68441 * eta5;
    } else {
        return -4.61686 + 11.1136 * eta - 10.4646 * eta2 + 5.11455 * eta3 - 1.27198 * eta4
            + 0.12746 * eta5;
    }
}

pub fn fresnel_moment2(eta: Float) -> Float {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let eta5 = eta4 * eta;
    if eta < 1.0 {
        return 0.27614 - 0.87350 * eta + 1.12077 * eta2 - 0.65095 * eta3
            + 0.07883 * eta4
            + 0.04860 * eta5;
    } else {
        let r_eta = 1.0 / eta;
        let r_eta2 = r_eta * r_eta;
        let r_eta3 = r_eta2 * r_eta;

        return -547.033 + 45.3087 * r_eta3 - 218.725 * r_eta2 + 458.843 * r_eta + 404.557 * eta
            - 189.519 * eta2
            + 54.9327 * eta3
            - 9.00603 * eta4
            + 0.63942 * eta5;
    }
}

pub fn subsurface_from_diffuse(
    t: &BSSRDFTable,
    rho_eff: &Spectrum,
    mfp: &Spectrum,
) -> (Spectrum, Spectrum) {
    let rho_eff = rho_eff.to_vec();
    let mfp = mfp.to_vec();
    let mut sigma_a = vec![0.0; mfp.len()];
    let mut sigma_s = vec![0.0; mfp.len()];
    for c in 0..mfp.len() {
        let rho = invert_catmull_rom(&t.rho_samples, &t.rho_eff, rho_eff[c]);
        sigma_s[c] = rho / mfp[c];
        sigma_a[c] = (1.0 - rho) / mfp[c];
    }
    //println!("{:?}, {:?}", sigma_a, sigma_s);
    return (Spectrum::from(sigma_a), Spectrum::from(sigma_s));
}
