use crate::core::pbrt::*;
use std::sync::Arc;

thread_local!(pub static DENSITY_BYTES: StatMemoryCounter = StatMemoryCounter::new("Memory/Volume density grid"));

#[derive(Debug, Clone)]
pub struct GridDensityMedium {
    #[allow(dead_code)]
    sigma_a: Spectrum,
    sigma_s: Spectrum,
    g: Float,
    nx: u32,
    ny: u32,
    nz: u32,
    world_to_medium: Transform,
    density: Vec<f32>,
    sigma_t: Float,
    inv_max_density: Float,
}

impl GridDensityMedium {
    pub fn new(
        sigma_a: &Spectrum,
        sigma_s: &Spectrum,
        g: Float,
        nx: u32,
        ny: u32,
        nz: u32,
        medium_to_world: &Transform,
        data: &[f32],
    ) -> Self {
        let sigma_a = *sigma_a;
        let sigma_s = *sigma_s;
        let sigma_t = (sigma_a + sigma_s)[0];
        //
        let mut max_density = 0.0;
        for i in 0..data.len() {
            max_density = Float::max(max_density, data[i]);
        }
        let inv_max_density = 1.0 / max_density;
        let density = data.to_vec();
        /*
        println!(
            "sigma_a: {:?}, sigma_s: {:?}, n: [{}, {}, {}], max density: {}, inv_max_density: {}",
            sigma_a, sigma_s,
            nx, ny, nz,
            max_density, inv_max_density
        );
        */

        DENSITY_BYTES.with(|s| {
            s.add(
                std::mem::size_of::<GridDensityMedium>() + std::mem::size_of::<f32>() * data.len(),
            );
        });

        GridDensityMedium {
            sigma_a,
            sigma_s,
            g,
            nx,
            ny,
            nz,
            world_to_medium: Transform::inverse(medium_to_world),
            density,
            sigma_t,
            inv_max_density,
        }
    }

    fn d(&self, p: &Point3i) -> f32 {
        let nx = self.nx as i32;
        let ny = self.ny as i32;
        let nz = self.nz as i32;
        let sample_bounds = Bounds3i::new(&Point3i::new(0, 0, 0), &Point3i::new(nx, ny, nz));
        if !sample_bounds.inside_exclusive(p) {
            return 0.0;
        }
        let index = (p.z * ny + p.y) * nx + p.x;
        self.density[index as usize]
    }

    fn density(&self, p: &Point3f) -> f32 {
        // Compute voxel coordinates and offsets for _p_
        let nx = self.nx as i32;
        let ny = self.ny as i32;
        let nz = self.nz as i32;
        let p_samples = Point3f::new(
            p.x * nx as Float - 0.5,
            p.y * ny as Float - 0.5,
            p.z * nz as Float - 0.5,
        );
        let ix = Float::floor(p_samples.x) as i32;
        let iy = Float::floor(p_samples.y) as i32;
        let iz = Float::floor(p_samples.z) as i32;
        let pi = Point3i::new(ix, iy, iz);
        let dx = p_samples.x - ix as Float;
        let dy = p_samples.y - iy as Float;
        let dz = p_samples.z - iz as Float;

        // Trilinearly interpolate density values to compute local density
        let d00 = lerp(
            dx,
            self.d(&(pi + Vector3i::new(0, 0, 0))),
            self.d(&(pi + Vector3i::new(1, 0, 0))),
        );
        let d10 = lerp(
            dx,
            self.d(&(pi + Vector3i::new(0, 1, 0))),
            self.d(&(pi + Vector3i::new(1, 1, 0))),
        );
        let d01 = lerp(
            dx,
            self.d(&(pi + Vector3i::new(0, 0, 1))),
            self.d(&(pi + Vector3i::new(1, 0, 1))),
        );
        let d11 = lerp(
            dx,
            self.d(&(pi + Vector3i::new(0, 1, 1))),
            self.d(&(pi + Vector3i::new(1, 1, 1))),
        );
        let d0 = lerp(dy, d00, d10);
        let d1 = lerp(dy, d01, d11);
        return lerp(dz, d0, d1);
    }
}

impl Medium for GridDensityMedium {
    fn sample(
        &self,
        r_world: &Ray,
        sampler: &mut dyn Sampler,
        _arena: &mut MemoryArena,
    ) -> (Spectrum, Option<MediumInteraction>) {
        let _p = ProfilePhase::new(Prof::MediumSample);
        let t_max = r_world.t_max.get();
        let ray = Ray::new(
            &r_world.o,
            &r_world.d.normalize(),
            t_max * r_world.d.length(),
            r_world.time,
        );
        let (ray, _, _) = self.world_to_medium.transform_ray(&ray);
        // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
        let b = Bounds3f::new(&Point3f::new(0.0, 0.0, 0.0), &Point3f::new(1.0, 1.0, 1.0));
        if let Some((tmin, tmax)) = b.intersect_p(&ray) {
            // Run delta-tracking iterations to sample a medium interaction
            let g = self.g;
            let sigma_t = self.sigma_t;
            let inv_max_density = self.inv_max_density;
            let inv_max_density_sigma_t = inv_max_density / sigma_t;
            let sigma_s = self.sigma_s;
            let sigma_s_t = sigma_s * (1.0 / sigma_t);
            let mut t = tmin;
            loop {
                t -= Float::ln(1.0 - sampler.get_1d()) * inv_max_density_sigma_t;
                if t >= tmax {
                    break;
                }
                let p = ray.position(t);
                let density = self.density(&p);
                if density * inv_max_density > sampler.get_1d() {
                    // Populate _mi_ with medium interaction information and return
                    let phase = HenyeyGreenstein::new(g);
                    let p_world = r_world.position(t);

                    let mi = MediumInteraction::new(
                        &p_world,
                        &(-r_world.d),
                        r_world.time,
                        &Some(Arc::new(phase)),
                    );
                    return (sigma_s_t, Some(mi));
                }
            }
        }
        return (Spectrum::one(), None);
    }

    fn tr(&self, r_world: &Ray, sampler: &mut dyn Sampler) -> Spectrum {
        let _p = ProfilePhase::new(Prof::MediumTr);
        let t_max = r_world.t_max.get();
        let ray = Ray::new(
            &r_world.o,
            &r_world.d.normalize(),
            t_max * r_world.d.length(),
            r_world.time,
        );
        let (ray, _, _) = self.world_to_medium.transform_ray(&ray);
        // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
        let b = Bounds3f::new(&Point3f::new(0.0, 0.0, 0.0), &Point3f::new(1.0, 1.0, 1.0));
        if let Some((tmin, tmax)) = b.intersect_p(&ray) {
            // Run delta-tracking iterations to sample a medium interaction
            //let g = self.g;
            let sigma_t = self.sigma_t;
            let inv_max_density = self.inv_max_density;
            let inv_max_density_sigma_t = inv_max_density / sigma_t;
            //let sigma_s = self.sigma_s;
            //let sigma_s_t = sigma_s * (1.0 / sigma_t);

            let mut tr = 1.0;
            let mut t = tmin;
            loop {
                t -= Float::ln(1.0 - sampler.get_1d()) * inv_max_density_sigma_t;
                if t >= tmax {
                    break;
                }
                let p = ray.position(t);
                let density = self.density(&p);
                tr *= 1.0 - Float::max(0.0, density * inv_max_density);
                // Added after book publication: when transmittance gets low,
                // start applying Russian roulette to terminate sampling.
                const RR_THRESHOLD: Float = 0.1;
                if tr < RR_THRESHOLD {
                    let q = Float::max(0.05, 1.0 - tr);
                    if sampler.get_1d() < q {
                        return Spectrum::zero();
                    }
                    tr /= 1.0 - q;
                }
            }
            return Spectrum::from(tr);
        }
        return Spectrum::one();
    }
}
