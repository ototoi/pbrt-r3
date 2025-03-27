use crate::core::prelude::*;

use std::sync::Arc;

struct MarbleTexture {
    mapping: Box<dyn TextureMapping3D>,
    octaves: u32,
    omega: Float,
    scale: Float,
    variation: Float,
}

fn lerps(c0: &Spectrum, c1: &Spectrum, t: Float) -> Spectrum {
    return *c0 * (1.0 - t) + *c1 * t;
}

impl MarbleTexture {
    pub fn new(
        mapping: Box<dyn TextureMapping3D>,
        octaves: u32,
        omega: Float,
        scale: Float,
        variation: Float,
    ) -> Self {
        Self {
            mapping,
            octaves,
            omega,
            scale,
            variation,
        }
    }

    pub fn evaluate_colors(&self, colors: &[[Float; 3]], si: &SurfaceInteraction) -> Spectrum {
        let variation = self.variation;
        let scale = self.scale;
        let omega = self.omega;
        let octaves = self.octaves;
        let (p, dpdx, dpdy) = self.mapping.map(si);
        let p = scale * p;
        let marble = p.y + variation * fbm(&p, &(scale * dpdx), &(scale * dpdy), omega, octaves);
        let t = 0.5 + 0.5 * Float::sin(marble);
        // Evaluate marble spline at _t_
        let nc = colors.len();
        let nseg = nc - 3;
        let first = usize::min(1, Float::floor(t * nseg as Float) as usize);
        let t = t * nseg as Float - first as Float;
        let c0 = Spectrum::from_rgb(&colors[first], SpectrumType::Reflectance);
        let c1 = Spectrum::from_rgb(&colors[first + 1], SpectrumType::Reflectance);
        let c2 = Spectrum::from_rgb(&colors[first + 2], SpectrumType::Reflectance);
        let c3 = Spectrum::from_rgb(&colors[first + 3], SpectrumType::Reflectance);
        // Bezier spline evaluated with de Castilejau's algorithm
        let s0 = lerps(&c0, &c1, t);
        let s1 = lerps(&c1, &c2, t);
        let s2 = lerps(&c2, &c3, t);
        let s0 = lerps(&s0, &s1, t);
        let s1 = lerps(&s1, &s2, t);
        // Extra scale of 1.5 to increase variation among colors
        return lerps(&s0, &s1, t) * 1.5;
    }
}

const C: [[Float; 3]; 9] = [
    [0.58, 0.58, 0.6],
    [0.58, 0.58, 0.6],
    [0.58, 0.58, 0.6],
    [0.5, 0.5, 0.5],
    [0.6, 0.59, 0.58],
    [0.58, 0.58, 0.6],
    [0.58, 0.58, 0.6],
    [0.2, 0.2, 0.33],
    [0.58, 0.58, 0.6],
];

impl Texture<Spectrum> for MarbleTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        return self.evaluate_colors(&C, si);
    }
}

//

pub fn create_marble_spectrum_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let map: Box<dyn TextureMapping3D> = Box::new(IdentityMapping3D::new(tex2world));
    let octaves = tp.find_int("octaves", 8) as u32;
    let roughness = tp.find_float("roughness", 0.5);
    let scale = tp.find_float("scale", 1.0);
    let variation = tp.find_float("variation", 0.2);
    return Ok(Arc::new(MarbleTexture::new(
        map, octaves, roughness, scale, variation,
    )));
}
