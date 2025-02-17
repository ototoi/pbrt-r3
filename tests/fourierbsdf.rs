// Imported from fourierbsdf.cpp

use std::sync::Arc;

use pbrt_r3::core::pbrt::*;

// Serialized version of bsdfs/roughgold_alpha_0.2.bsdf.
const FOURIER_DATA: &'static [u8] = include_bytes!("bsdfs/roughgold_alpha_0.2.bsdf");

#[test]
fn bsdfs_fourier() {
    let table = FourierBSDFTable::read_from_bytes(&FOURIER_DATA);
    assert!(table.is_ok());
    let table = table.unwrap();
    let table = Arc::new(table);
    let bsdf = FourierBSDF::new(&table, TransportMode::Radiance);

    // Relative error
    let err = |a: Float, b: Float| (a - b).abs() / b;

    // Check evaluation.
    let wo = Vector3f::new(-0.5, -0.5, 0.8).normalize();
    // Sort-of close to the specular reflection direction.
    let wi = Vector3f::new(0.4, 0.52, 0.7).normalize();
    assert!(err(bsdf.f(&wo, &wi).y(), 2.679294) < 0.001);

    // Pdf
    assert!(err(bsdf.pdf(&wo, &wi), 2.438230) < 0.001);
    assert!(err(bsdf.pdf(&wi, &wo), 2.503326) < 0.001);

    // Sampling
    if let Some((fr, wi, pdf, _)) = bsdf.sample_f(&wo, &Point2f::new(0.1, 0.8)) {
        assert!(err(fr.y(), 2.596391) < 0.001);
        assert!(err(pdf, 1.855472) < 0.001);
        assert!(err(wi.x, 0.539052) < 0.001);
        assert!(err(wi.y, 0.617347) < 0.001);
        assert!(err(wi.z, 0.572980) < 0.001);
    } else {
        assert!(false);
    }
}
