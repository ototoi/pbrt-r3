include!(concat!(env!("OUT_DIR"), "/spectrum_utils.rs"));
use crate::core::misc::find_interval;

pub fn spectrum_samples_sorted(lambda: &[f32], _vals: &[f32]) -> bool {
    for i in 0..(lambda.len() - 1) {
        if lambda[i] > lambda[i + 1] {
            return false;
        }
    }
    return true;
}

pub fn sort_spectrum_samples(lambda: &mut [f32], vals: &mut [f32]) {
    let n = lambda.len();
    let mut sort_vec: Vec<(f32, f32)> = Vec::with_capacity(n);
    for i in 0..n {
        sort_vec.push((lambda[i], vals[i]));
    }
    sort_vec.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for i in 0..n {
        lambda[i] = sort_vec[i].0;
        vals[i] = sort_vec[i].1;
    }
}

pub fn interpolate_spectrum_samples(lambda: &[f32], vals: &[f32], l: f32) -> f32 {
    //for (int i = 0; i < n - 1; ++i) CHECK_GT(lambda[i + 1], lambda[i]);
    let n = lambda.len();
    if l <= lambda[0] {
        return vals[0];
    }
    if l >= lambda[n - 1] {
        return vals[n - 1];
    }
    let offset = find_interval(lambda, &|v, index| -> bool {
        return v[index] <= l;
    });
    //CHECK(l >= lambda[offset] && l <= lambda[offset + 1]);
    let t = (l - lambda[offset]) / (lambda[offset + 1] - lambda[offset]);
    return lerp(t, vals[offset], vals[offset + 1]);
}

// Given a piecewise-linear SPD with values in vIn[] at corresponding
// wavelengths lambdaIn[], where lambdaIn is assumed to be sorted but may
// be irregularly spaced, resample the spectrum over the range of
// wavelengths [lambdaMin, lambdaMax], with a total of nOut wavelength
// samples between lambdaMin and lamdbaMax (including those at
// endpoints). The resampled spectrum values are written to vOut.
//
// In general, this is a standard sampling and reconstruction problem, with
// the complication that for any given invocation, some of the
// reconstruction points may involve upsampling the input distribution and
// others may involve downsampling. For upsampling, we just point-sample,
// and for downsampling, we apply a box filter centered around the
// destination wavelength with total width equal to the sample spacing.
pub fn resample_linear_spectrum(
    lambda_in: &[f32],
    v_in: &[f32],
    lambda_min: f32,
    lambda_max: f32,
    v_out: &mut [f32],
) {
    let n_in = lambda_in.len();
    let n_out = v_out.len();
    assert!(n_out > 2);
    for i in 0..(n_in - 1) {
        assert!(lambda_in[i] <= lambda_in[i + 1]);
    }
    assert!(lambda_min < lambda_max);

    // Spacing between samples in the output distribution.
    let delta = (lambda_max - lambda_min) / (n_out - 1) as f32;

    // We assume that the SPD is constant outside of the specified
    // wavelength range, taking on the respectively first/last SPD value
    // for out-of-range wavelengths.
    //
    // To make this convention fit cleanly into the code below, we create
    // virtual samples in the input distribution with index -1 for the
    // sample before the first valid sample and index nIn for the sample
    // after the last valid sample. In turn, can place those virtual
    // samples beyond the endpoints of the target range so that we can
    // always assume that the source range is broader than the target
    // range, which in turn lets us not worry about various boundary cases
    // below.

    // The wavelengths of the virtual samples at the endpoints are set so
    // that they are one destination sample spacing beyond the destination
    // endpoints.  (Note that this potentially means that if we swept along
    // indices from -1 to nIn, we wouldn't always see a monotonically
    // increasing set of wavelength values. However, this isn't a problem
    // since we only access these virtual samples if the destination range
    // is wider than the source range.)
    let lambda_in_clamped = |index: i32| -> f32 {
        assert!(index >= -1 && index <= n_in as i32);
        if index == -1 {
            assert!(lambda_min - delta < lambda_in[0]);
            return lambda_min - delta;
        } else if index == n_in as i32 {
            assert!(lambda_max + delta > lambda_in[n_in - 1]);
            return lambda_max + delta;
        } else {
            return lambda_in[index as usize];
        }
    };

    // Due to the piecewise-constant assumption, the SPD values outside the
    // specified range are given by the valid endpoints.
    let v_in_clamped = |index: i32| -> f32 {
        assert!(index >= -1 && index <= n_in as i32);
        return v_in[index.clamp(0, n_in as i32 - 1) as usize];
    };

    // Helper that upsamples ors downsample the given SPD at the given
    // wavelength lambda.
    let resample = |lambda: f32| -> f32 {
        // Handle the edge cases first so that we don't need to worry about
        // them in the following.
        //
        // First, if the entire filtering range for the destination is
        // outside of the range of valid samples, we can just return the
        // endpoint value.

        if lambda + delta / 2.0 <= lambda_in[0] {
            return v_in[0];
        }
        if lambda - delta / 2.0 >= lambda_in[n_in - 1] {
            return v_in[n_in - 1];
        }
        // Second, if there's only one sample, then the SPD has the same
        // value everywhere, and we're done.
        if n_in == 1 {
            return v_in[0];
        }

        // Otherwise, find indices into the input SPD that bracket the
        // wavelength range [lambda-delta, lambda+delta]. Note that this is
        // a 2x wider range than we will actually filter over in the end.
        let start;
        let mut end;
        if lambda - delta < lambda_in[0] {
            // Virtual sample at the start, as described above.
            start = -1;
        } else {
            start = find_interval(lambda_in, &|v, index| -> bool {
                return v[index] <= lambda - delta;
            }) as i32;
            assert!(start >= 0 && start < n_in as i32);
        }

        if lambda + delta > lambda_in[n_in - 1] {
            // Virtual sample at the end, as described above.
            end = n_in as i32;
        } else {
            // Linear search from the starting point. (Presumably more
            // efficient than a binary search from scratch, or doesn't
            // matter either way.)
            end = if start > 0 { start } else { 0 };
            while end < n_in as i32 && lambda + delta > lambda_in[end as usize] {
                end += 1;
            }
        }

        if end - start == 2
            && lambda_in_clamped(start) <= lambda - delta
            && lambda_in[(start + 1) as usize] == lambda
            && lambda_in_clamped(end) >= lambda + delta
        {
            // Downsampling: special case where the input and output
            // wavelengths line up perfectly, so just return the
            // corresponding point sample at lambda.
            return v_in[(start + 1) as usize];
        } else if end - start == 1 {
            // Downsampling: evaluate the piecewise-linear function at
            // lambda.
            let t = (lambda - lambda_in_clamped(start))
                / (lambda_in_clamped(end) - lambda_in_clamped(start));
            assert!(t >= 0.0 && t <= 1.0);
            return lerp(t, v_in_clamped(start), v_in_clamped(end));
        } else {
            // Upsampling: use a box filter and average all values in the
            // input spectrum from lambda +/- delta / 2.
            return average_spectrum_samples(
                lambda_in,
                v_in,
                lambda - delta / 2.0,
                lambda + delta / 2.0,
            );
        }
    };

    // For each destination sample, compute the wavelength lambda for the
    // sample and then resample the source SPD distribution at that point.
    for out_offset in 0..n_out {
        // TODO: Currently, resample() does a binary search each time,
        // even though we could do a single sweep across the input array,
        // since we're resampling it at a regular and increasing set of
        // lambdas. It would be nice to polish that up.
        let lambda = lerp(
            out_offset as f32 / (n_out as f32 - 1.0),
            lambda_min,
            lambda_max,
        );
        v_out[out_offset] = resample(lambda);
    }
}
