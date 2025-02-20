// Imported from bsdfs.cpp

use pbrt_r3::core::pbrt::*;
use pbrt_r3::shapes::*;

use std::sync::Arc;

static CHI2_THETA_RES: usize = 10;
static CHI2_PHI_RES: usize = 2 * CHI2_THETA_RES;
static CHI2_SAMPLECOUNT: usize = 1000000;
static CHI2_RUNS: usize = 5;


fn frequency_table(bsdf: & BSDF, wo: &Vector3f, rng: &mut RNG, sample_count: usize, theta_res: usize, phi_res: usize) -> Vec<Float> {
    let mut table = vec![0.0; (theta_res * phi_res) as usize];
    for bxdf in bsdf.bxdfs.iter() {
        if bxdf.matches_flags(BSDF_REFLECTION | BSDF_TRANSMISSION) {
            continue;
        }
        let flags = bxdf.get_type();
        let index = match flags {
            BSDF_DIFFUSE => 0,
            BSDF_GLOSSY => 1,
            BSDF_SPECULAR => 2,
            _ => 3,
        };
        let f = bxdf.f(wo, &Vector3f::default());
        if f.y() > 0.0 {
            table[index] += 1.0;
        }
    }
    table
}

fn chi2_test(frequencies: &[Float], exp_frequencies: &[Float]) -> Result<(), String> {
    let mut sum = 0.0;
    for i in 0..frequencies.len() {
        let diff = frequencies[i] - exp_frequencies[i];
        sum += diff * diff / exp_frequencies[i];
    }
    Ok(())
}

fn test_bsdf(create_bsdf: &(dyn Fn(&mut BSDF, &mut MemoryArena)), description: &str) {
    let mut arena = MemoryArena::new();

    let theta_res = CHI2_THETA_RES;
    let phi_res = CHI2_PHI_RES;
    let sample_count = CHI2_SAMPLECOUNT;

    let mut rng = RNG::new();

    let t = Transform::rotate_x(-90.0);
    let t_inv = t.inverse();

    let mut bsdf = None;
    {
        let reverse_orientation = false;
        let disk = Disk::new(
            &t,
            &t_inv,
            reverse_orientation,
            0.0,
            1.0,
            0.0,
            360.0,
        );
        let origin = Point3f::new(0.1, 1.0, 0.0);
        let direction = Vector3f::new(0.0, -1.0, 0.0);
        let r = Ray::new(&origin, &direction, Float::INFINITY, 0.0);
        if let Some((_t_hit, isect)) = disk.intersect(&r) {
            let mut b = BSDF::new(&isect, 1.0);
            create_bsdf(&mut b, &mut arena);
            bsdf = Some(b);
        }
    }
    if let Some(bsdf) = bsdf.as_ref() {
        for _ in 0..CHI2_RUNS {
            /* Randomly pick an outgoing direction on the hemisphere */
            let sample = Point2f::new(rng.uniform_float(), rng.uniform_float());
            let wo_l = cosine_sample_hemisphere(&sample);
            let wo = bsdf.world_to_local(&wo_l);
            let frequencies = frequency_table(bsdf, &wo, &mut rng, sample_count, theta_res, phi_res);
            //
            
        }
    }
}
    
fn create_lambertian(bsdf: &mut BSDF, _arena: &mut MemoryArena) {
    let kd = Spectrum::from(1.0);
    let r: Arc<dyn BxDF> = Arc::new(LambertianReflection::new(&kd));
    bsdf.add(&r);
}


#[test]
fn bsdf_sampling_lambertian() {
    test_bsdf(&create_lambertian, "Lambertian");
}