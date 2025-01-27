use crate::core::lowdiscrepancy::maxmin::*;
use crate::core::pbrt::*;
use log::*;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, PartialEq, Default, Clone)]
pub struct MaxMinDistSampler {
    base: BasePixelSampler,
    cpixel: [u32; 32],
}

impl MaxMinDistSampler {
    pub fn new(samples_per_pixel: u32, n_sampled_dimensions: u32) -> Self {
        let mut spp = samples_per_pixel;
        let cindex = log2int(spp) as usize;
        let maxminlen = CMAXMIN_DIST.len();
        assert_eq!(maxminlen, 17);
        if cindex >= maxminlen {
            //spp = (1 << maxminlen) - 1; //original code
            spp = 1 << (maxminlen - 1); //pbrt-r3 patched

            //Warning
            warn!("No more than {} samples per pixel are supported with MaxMinDistSampler. Rounding down.", spp);
        }
        if !is_power_of_2(spp as i32) {
            spp = round_up_pow2(spp as i32) as u32;

            //Warning
            warn!(
                "Non power-of-two sample count rounded up to {} for MaxMinDistSampler.",
                spp
            );
        }
        let cindex = log2int(spp) as usize;
        assert!(cindex < maxminlen);
        let cpixel = CMAXMIN_DIST[cindex];
        MaxMinDistSampler {
            base: BasePixelSampler::new(spp, n_sampled_dimensions),
            cpixel,
        }
    }
}

impl Sampler for MaxMinDistSampler {
    fn start_pixel(&mut self, p: &Point2i) {
        let _p = ProfilePhase::new(Prof::StartPixel);

        let samples_per_pixel = self.base.base.samples_per_pixel as usize;
        let inv_spp = 1.0 / samples_per_pixel as Float;
        for i in 0..samples_per_pixel {
            self.base.samples2d[0][i] = Point2f::new(
                (i as Float) * inv_spp,
                sample_generator_matrix(&self.cpixel, i as u32, 0),
            );
        }
        shuffle_array(
            &mut self.base.samples2d[0],
            samples_per_pixel,
            1,
            &mut self.base.rng,
        );
        // Generate remaining samples for _MaxMinDistSampler_
        for i in 0..self.base.samples1d.len() {
            van_der_corput(
                1,
                samples_per_pixel as u32,
                &mut self.base.samples1d[i],
                &mut self.base.rng,
            );
        }

        for i in 1..self.base.samples2d.len() {
            sobol_2d(
                1,
                samples_per_pixel as u32,
                &mut self.base.samples2d[i],
                &mut self.base.rng,
            );
        }

        for i in 0..self.base.base.samples1d_array_sizes.len() {
            let count = self.base.base.samples1d_array_sizes[i];
            van_der_corput(
                count,
                samples_per_pixel as u32,
                &mut self.base.base.sample_array1d[i],
                &mut self.base.rng,
            );
        }

        for i in 0..self.base.base.samples2d_array_sizes.len() {
            let count = self.base.base.samples2d_array_sizes[i];
            sobol_2d(
                count,
                samples_per_pixel as u32,
                &mut self.base.base.sample_array2d[i],
                &mut self.base.rng,
            );
        }
        self.base.start_pixel(p);
    }

    fn start_next_sample(&mut self) -> bool {
        return self.base.start_next_sample();
    }

    fn set_sample_number(&mut self, sample_num: u32) -> bool {
        return self.base.set_sample_number(sample_num);
    }

    fn get_1d(&mut self) -> Float {
        return self.base.get_1d();
    }

    fn get_2d(&mut self) -> Vector2f {
        return self.base.get_2d();
    }

    fn request_1d_array(&mut self, n: u32) {
        self.base.base.request_1d_array(n);
    }

    fn request_2d_array(&mut self, n: u32) {
        self.base.base.request_2d_array(n);
    }

    fn get_1d_array(&mut self, n: u32) -> Option<Vec<Float>> {
        return self.base.base.get_1d_array(n);
    }

    fn get_2d_array(&mut self, n: u32) -> Option<Vec<Vector2f>> {
        return self.base.base.get_2d_array(n);
    }

    fn round_count(&self, n: u32) -> u32 {
        return round_up_pow2(n as i32) as u32;
    }

    fn clone_with_seed(&self, seed: u32) -> Arc<RwLock<dyn Sampler>> {
        let mut ss = self.clone();
        ss.base.rng.set_sequence(seed as u64);
        return Arc::new(RwLock::new(ss));
    }

    fn get_samples_per_pixel(&self) -> u32 {
        return self.base.base.samples_per_pixel;
    }
}

impl PixelSampler for MaxMinDistSampler {}

pub fn create_maxmindist_sampler(params: &ParamSet) -> Result<Arc<RwLock<dyn Sampler>>, PbrtError> {
    let mut nsamp = params.find_one_int("pixelsamples", 16) as u32;
    let sd = params.find_one_int("dimensions", 4) as u32;
    {
        let options = PbrtOptions::get();
        if options.quick_render {
            nsamp = 1;
        }
    }
    return Ok(Arc::new(RwLock::new(MaxMinDistSampler::new(nsamp, sd))));
}
