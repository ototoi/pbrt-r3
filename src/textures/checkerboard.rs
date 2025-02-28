use crate::core::pbrt::*;

use std::ops::*;
use std::sync::Arc;

// AAMethod Declaration
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum AAMethod {
    None,
    ClosedForm,
}

pub struct Checkerboard2DTexture<T> {
    mapping: Box<dyn TextureMapping2D>,
    tex1: Arc<dyn Texture<T>>,
    tex2: Arc<dyn Texture<T>>,
    aa_method: AAMethod,
}

impl<T> Checkerboard2DTexture<T> {
    pub fn new(
        mapping: Box<dyn TextureMapping2D>,
        tex1: &Arc<dyn Texture<T>>,
        tex2: &Arc<dyn Texture<T>>,
        aa_method: AAMethod,
    ) -> Self {
        return Checkerboard2DTexture::<T> {
            mapping,
            tex1: Arc::clone(tex1),
            tex2: Arc::clone(tex2),
            aa_method,
        };
    }
}

#[inline]
fn bump_int(x: Float) -> Float {
    return Float::floor(x / 2.0) + 2.0 * Float::max(x / 2.0 - Float::floor(x / 2.0) - 0.5, 0.0);
}

impl<T: Copy + Add<T, Output = T> + Mul<Float, Output = T>> Texture<T>
    for Checkerboard2DTexture<T>
{
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let (st, dstdx, dstdy) = self.mapping.map(si);
        if self.aa_method == AAMethod::None {
            // Point sample _Checkerboard2DTexture_
            if ((Float::floor(st[0]) as i32) + (Float::floor(st[1]) as i32)) % 2 == 0 {
                return self.tex1.as_ref().evaluate(si);
            } else {
                return self.tex2.as_ref().evaluate(si);
            }
        } else {
            // Compute closed-form box-filtered _Checkerboard2DTexture_ value

            // Evaluate single check if filter is entirely inside one of them
            let ds = Float::max(Float::abs(dstdx[0]), Float::abs(dstdy[0]));
            let dt = Float::max(Float::abs(dstdx[1]), Float::abs(dstdy[1]));
            let s0 = st[0] - ds;
            let s1 = st[0] + ds;
            let t0 = st[1] - dt;
            let t1 = st[1] + dt;
            if Float::floor(s0) == Float::floor(s1) && Float::floor(t0) == Float::floor(t1) {
                // Point sample _Checkerboard2DTexture_
                if ((Float::floor(st[0]) as i32) + (Float::floor(st[1]) as i32)) % 2 == 0 {
                    return self.tex1.as_ref().evaluate(si);
                } else {
                    return self.tex2.as_ref().evaluate(si);
                }
            }

            let sint = (bump_int(s1) - bump_int(s0)) / (2.0 * ds);
            let tint = (bump_int(t1) - bump_int(t0)) / (2.0 * dt);
            let mut area2 = sint + tint - 2.0 * sint * tint;
            if ds > 1.0 || dt > 1.0 {
                area2 = 0.5;
            }
            return self.tex1.as_ref().evaluate(si) * (1.0 - area2)
                + self.tex2.as_ref().evaluate(si) * area2;
        }
    }
}

pub struct Checkerboard3DTexture<T> {
    mapping: Box<dyn TextureMapping3D>,
    tex1: Arc<dyn Texture<T>>,
    tex2: Arc<dyn Texture<T>>,
}

impl<T> Checkerboard3DTexture<T> {
    pub fn new(
        mapping: Box<dyn TextureMapping3D>,
        tex1: &Arc<dyn Texture<T>>,
        tex2: &Arc<dyn Texture<T>>,
    ) -> Self {
        return Checkerboard3DTexture::<T> {
            mapping,
            tex1: Arc::clone(tex1),
            tex2: Arc::clone(tex2),
        };
    }
}

impl<T: Copy> Texture<T> for Checkerboard3DTexture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let (st, _dstdx, _dstdy) = self.mapping.map(si);
        if ((Float::floor(st[0]) as i32)
            + (Float::floor(st[1]) as i32)
            + (Float::floor(st[2]) as i32))
            % 2
            == 0
        {
            return self.tex1.as_ref().evaluate(si);
        } else {
            return self.tex2.as_ref().evaluate(si);
        }
    }
}

fn create_aa_method(aa: &str) -> Result<AAMethod, PbrtError> {
    match aa {
        "none" => Ok(AAMethod::None),
        "closedform" => Ok(AAMethod::ClosedForm),
        _ => {
            let msg = format!("Antialiasing mode \"{}\" not understood by \"Checkerboard2DTexture\"; using \"closedform\"", aa);
            return Err(PbrtError::error(&msg));
        }
    }
}

pub fn create_checkerboard_float_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let dim = tp.find_int("dimension", 2);
    if dim != 2 && dim != 3 {
        let msg = format!("{} dimensional checkerboard texture not supported", dim);
        return Err(PbrtError::error(&msg));
    }
    let tex1 = tp.get_float_texture("tex1", 1.0);
    let tex2 = tp.get_float_texture("tex2", 0.0);
    if dim == 2 {
        let map = create_texture_mapping2d(tex2world, tp)?;
        // Compute _aaMethod_ for _CheckerboardTexture_
        let aa = tp.find_string("aamode", "closedform");
        let aa_method = create_aa_method(&aa)?;
        return Ok(Arc::new(Checkerboard2DTexture::<Float>::new(
            map, &tex1, &tex2, aa_method,
        )));
    } else {
        // Initialize 3D texture mapping _map_ from _tp_
        let map = Box::new(IdentityMapping3D::new(tex2world));
        return Ok(Arc::new(Checkerboard3DTexture::<Float>::new(
            map, &tex1, &tex2,
        )));
    }
}

pub fn create_checkerboard_spectrum_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let dim = tp.find_int("dimension", 2);
    if dim != 2 && dim != 3 {
        let msg = format!("{} dimensional checkerboard texture not supported", dim);
        return Err(PbrtError::error(&msg));
    }
    let tex1 = tp.get_spectrum_texture("tex1", &Spectrum::from(1.0));
    let tex2 = tp.get_spectrum_texture("tex2", &Spectrum::from(0.0));
    if dim == 2 {
        let map = create_texture_mapping2d(tex2world, tp)?;
        // Compute _aaMethod_ for _CheckerboardTexture_
        let aa = tp.find_string("aamode", "closedform");
        let aa_method = create_aa_method(&aa)?;
        return Ok(Arc::new(Checkerboard2DTexture::<Spectrum>::new(
            map, &tex1, &tex2, aa_method,
        )));
    } else {
        // Initialize 3D texture mapping _map_ from _tp_
        let map = Box::new(IdentityMapping3D::new(tex2world));
        return Ok(Arc::new(Checkerboard3DTexture::<Spectrum>::new(
            map, &tex1, &tex2,
        )));
    }
}
