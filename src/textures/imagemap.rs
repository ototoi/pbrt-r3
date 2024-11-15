use crate::core::imageio::*;
use crate::core::pbrt::*;
use crate::core::texture::TexInfo;

use std::marker::PhantomData;
use std::path::Path;
use std::sync::Arc;

pub struct ImageTexture<Tmemory, Treturn> {
    mapping: Box<dyn TextureMapping2D>,
    mipmap: MIPMap<Tmemory>,
    phantom: PhantomData<Treturn>, //non allocate
}

impl ImageTexture<Float, Float> {
    pub fn new(mapping: Box<dyn TextureMapping2D>, mipmap: MIPMap<Float>) -> Self {
        Self {
            mapping,
            mipmap,
            phantom: PhantomData,
        }
    }

    pub fn convert_in(from: &RGBSpectrum, scale: Float, gamma: bool) -> Float {
        return scale
            * (if gamma {
                inverse_gamma_correct(from.y())
            } else {
                from.y()
            });
    }
}

impl ImageTexture<RGBSpectrum, Spectrum> {
    pub fn new(mapping: Box<dyn TextureMapping2D>, mipmap: MIPMap<RGBSpectrum>) -> Self {
        Self {
            mapping,
            mipmap,
            phantom: PhantomData,
        }
    }

    fn igc(from: &RGBSpectrum) -> RGBSpectrum {
        let rgb = from.to_rgb();
        let rgb: Vec<Float> = rgb.iter().map(|x| inverse_gamma_correct(*x)).collect();
        return RGBSpectrum::from(rgb);
    }

    pub fn convert_in(from: &RGBSpectrum, scale: Float, gamma: bool) -> RGBSpectrum {
        return (if gamma { Self::igc(from) } else { *from }) * scale;
    }

    pub fn convert_out(from: &RGBSpectrum) -> Spectrum {
        return Spectrum::from_rgb(&from.to_rgb(), SpectrumType::Reflectance);
    }
}

impl Texture<Float> for ImageTexture<Float, Float> {
    fn evaluate(&self, si: &SurfaceInteraction) -> Float {
        let (st, dstdx, dstdy) = self.mapping.as_ref().map(si);
        return self.mipmap.lookup_delta(&st, &dstdx, &dstdy);
    }
}

impl Texture<Spectrum> for ImageTexture<RGBSpectrum, Spectrum> {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        let (st, dstdx, dstdy) = self.mapping.as_ref().map(si);
        return Self::convert_out(&self.mipmap.lookup_delta(&st, &dstdx, &dstdy));
    }
}

fn has_extension(path: &str, ext: &str) -> bool {
    if let Some(e) = Path::new(path).extension() {
        if let Some(s) = e.to_str() {
            if s == ext {
                return true;
            }
        }
    }
    return false;
}

fn has_annotation(path: &str, annotations: &[&str]) -> bool {
    for anno in annotations {
        if path.contains(*anno) {
            return true;
        }
    }
    return false;
}

fn get_wrap_mode(mode: &str) -> ImageWrap {
    return match mode {
        "repeat" => ImageWrap::Repeat,
        "black" => ImageWrap::Black,
        "clamp" => ImageWrap::Clamp,
        _ => ImageWrap::Repeat,
    };
}

pub fn create_texinfo(tp: &TextureParams) -> Option<TexInfo> {
    // Initialize _ImageTexture_ parameters
    let filename = tp.find_filename("filename", "");
    if filename.is_empty() {
        return None;
    }

    let max_aniso = tp.find_float("maxanisotropy", 8.0);
    let trilinear = tp.find_bool("trilinear", false);
    let wrap = tp.find_string("wrap", "repeat");
    let swrap = tp.find_string("swrap", &wrap);
    let twrap = tp.find_string("twrap", &wrap);

    let swrap_mode = get_wrap_mode(&swrap);
    let twrap_mode = get_wrap_mode(&twrap);

    let scale = tp.find_float("scale", 1.0);
    let mut gamma = !has_extension(&filename, "exr");
    // Fixed in pbrt-r3
    if has_annotation(&filename, &["_alpha", "_bump"]) {
        //println!("has_annotation:{}", filename);
        gamma = false;
    }
    gamma = tp.find_bool("gamma", gamma);

    return Some(TexInfo {
        filename,
        trilinear,
        max_aniso,
        swrap_mode,
        twrap_mode,
        scale,
        gamma,
        flip_y: true,
    });
}

fn flip_y<T: Copy>(data: &mut [T], resolution: &Vector2i) {
    // Flip image in y; texture coordinate space has (0,0) at the lower
    // left corner.
    let w = resolution.x;
    let h = resolution.y;
    for y in 0..(h / 2) {
        for x in 0..w {
            let o1 = (y * w + x) as usize;
            let o2 = ((h - 1 - y) * w + x) as usize;
            data.swap(o1, o2);
        }
    }
}

pub type FloatImageTexture = ImageTexture<Float, Float>;
pub type SpectrumImageTexture = ImageTexture<RGBSpectrum, Spectrum>;

pub fn create_image_float_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Float>>, PbrtError> {
    let mapping = create_texture_mapping2d(tex2world, tp)?;
    let texinfo = create_texinfo(tp).unwrap();
    let scale = texinfo.scale;
    let gamma = texinfo.gamma;

    let cache_path = get_mipmap_cache_dir(&texinfo);
    let file_hash = get_mipmap_file_hash(&texinfo.filename)?;
    if let Ok(mipmap) = load_float_mipmap_cache(&cache_path, &file_hash) {
        //println!("load from cache:{}", cache_path);
        return Ok(Arc::new(FloatImageTexture::new(mapping, mipmap)));
    } else {
        //println!("make from scratch");
        let (data, resolution) = read_image_gamma_correct(&texinfo.filename, false)?;
        let mut data: Vec<Float> = data
            .iter()
            .map(|spec| -> Float { FloatImageTexture::convert_in(spec, scale, gamma) })
            .collect();
        if texinfo.flip_y {
            flip_y(&mut data, &resolution);
        }
        let mipmap = MIPMap::<Float>::new(
            &resolution,
            &data,
            texinfo.trilinear,
            texinfo.max_aniso,
            texinfo.swrap_mode,
            texinfo.twrap_mode,
        );
        let _ = save_float_mipmap_cache(&cache_path, &texinfo, &file_hash, &mipmap);
        return Ok(Arc::new(FloatImageTexture::new(mapping, mipmap)));
    }
}

pub fn create_image_spectrum_texture(
    tex2world: &Transform,
    tp: &TextureParams,
) -> Result<Arc<dyn Texture<Spectrum>>, PbrtError> {
    let mapping = create_texture_mapping2d(tex2world, tp)?;
    let texinfo = create_texinfo(tp).unwrap();
    let scale = texinfo.scale;
    let gamma = texinfo.gamma;

    let cache_path = get_mipmap_cache_dir(&texinfo);
    let file_hash = get_mipmap_file_hash(&texinfo.filename)?;
    if let Ok(mipmap) = load_spectrum_mipmap_cache(&cache_path, &file_hash) {
        //println!("load from cache:{}", cache_path);
        return Ok(Arc::new(SpectrumImageTexture::new(mapping, mipmap)));
    } else {
        //println!("make from scratch");
        let (data, resolution) = read_image_gamma_correct(&texinfo.filename, false)?;
        let mut data: Vec<RGBSpectrum> = data
            .iter()
            .map(|spec| -> RGBSpectrum { SpectrumImageTexture::convert_in(spec, scale, gamma) })
            .collect();
        if texinfo.flip_y {
            flip_y(&mut data, &resolution);
        }
        let mipmap = MIPMap::<RGBSpectrum>::new(
            &resolution,
            &data,
            texinfo.trilinear,
            texinfo.max_aniso,
            texinfo.swrap_mode,
            texinfo.twrap_mode,
        );
        let _ = save_spectrum_mipmap_cache(&cache_path, &texinfo, &file_hash, &mipmap);
        return Ok(Arc::new(SpectrumImageTexture::new(mapping, mipmap)));
    }
}
