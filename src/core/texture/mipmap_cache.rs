use std::env;
use std::error::Error;
use std::fs::*;
use std::io::BufReader;
use std::io::Read;
use std::path::Path;

use super::mipmap::*;
use super::texture::*;
use crate::core::base::*;
use crate::core::error::*;
use crate::core::imageio::*;
use crate::core::spectrum::*;

use crypto::digest::Digest;
use crypto::md5::Md5;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct TexCacheInfo {
    pub texinfo: TexInfo,
    pub image_file_hash: String,
    pub mip_levels: u32,
}

fn hash(s: &str) -> String {
    let mut md5 = Md5::new();
    md5.input(s.as_bytes());
    return md5.result_str();
}

pub fn get_mipmap_file_hash(path: &str) -> Result<String, PbrtError> {
    let mut file = File::open(path)?;
    let mut buf = Vec::new();
    let _ = file.read_to_end(&mut buf)?;
    let mut md5 = Md5::new();
    md5.input(&buf);
    return Ok(md5.result_str());
}

pub fn get_mipmap_cache_dir(texinfo: &TexInfo) -> String {
    //temp_dir/pbrt/textures/filename/hash(fullpath)/hash(jsonify(texinfo))/
    let temp_dir = env::temp_dir();
    let fullpath = texinfo.filename.clone();
    let filename = Path::new(&fullpath).file_stem().unwrap().to_str().unwrap();
    let path = temp_dir
        .join("pbrt")
        .join("textures")
        .join(filename)
        .join(hash(&fullpath))
        .join(hash(&texinfo.to_string()));
    return String::from(path.to_str().unwrap());
}

fn load_cacheinfo_from_json(path: &str) -> Result<TexCacheInfo, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let cache_info = serde_json::from_reader(reader)?;
    return Ok(cache_info);
}

fn save_cacheinfo_to_json(path: &str, cache_info: &TexCacheInfo) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, cache_info)?;
    return Ok(());
}

fn convert_float_image(img: Vec<f32>, resolution: &Point2i) -> Vec<f32> {
    let channels = img.len() / (resolution.x * resolution.y) as usize;
    if channels == 1 {
        return img;
    } else {
        let mut v = vec![0.0; (resolution.x * resolution.y) as usize];
        for i in 0..v.len() {
            v[i] = img[channels * i];
        }
        return v;
    }
}

fn convert_spectrum_image(img: Vec<f32>, _resolution: &Point2i) -> Vec<f32> {
    return img;
}

pub fn load_float_mipmap_cache(dir: &str, file_hash: &str) -> Result<MIPMap<Float>, PbrtError> {
    let dir = Path::new(&dir);
    let json_path = dir.join("cache_info.json");
    let cache_info = load_cacheinfo_from_json(json_path.to_str().unwrap())?;
    if cache_info.image_file_hash != file_hash {
        return Err(PbrtError::error("File hash of image is not same."));
    }
    let mut pyramid = Vec::with_capacity(cache_info.mip_levels as usize);
    for i in 0..cache_info.mip_levels {
        let mip_image_path = dir.join(format!("{}.exr", i));
        let (image, resolution) = read_cache_image(mip_image_path.to_str().unwrap())?;
        let image = convert_float_image(image, &resolution);
        let mip_image = F32MIPMapImage::new(
            image,
            (resolution.x as usize, resolution.y as usize),
        );
        pyramid.push(mip_image);
    }
    let mipmap = MIPMap::<Float>::make_from_pyramid(
        pyramid,
        cache_info.texinfo.trilinear,
        cache_info.texinfo.max_aniso,
        cache_info.texinfo.swrap_mode,
        cache_info.texinfo.twrap_mode,
    );
    return Ok(mipmap);
}

pub fn save_float_mipmap_cache(
    dir: &str,
    texinfo: &TexInfo,
    file_hash: &str,
    mipmap: &MIPMap<Float>,
) -> Result<(), PbrtError> {
    if !Path::new(dir).exists() {
        create_dir_all(dir)?;
    }
    let dir = Path::new(&dir);
    let json_path = dir.join("cache_info.json");
    let cache_info = TexCacheInfo {
        texinfo: texinfo.clone(),
        image_file_hash: String::from(file_hash),
        mip_levels: mipmap.levels() as u32,
    };
    save_cacheinfo_to_json(json_path.to_str().unwrap(), &cache_info)?;
    for i in 0..cache_info.mip_levels {
        let mip_image_path = dir.join(format!("{}.exr", i));
        let image = &mipmap.pyramid[i as usize];
        write_cache_image(
            mip_image_path.to_str().unwrap(),
            &image.data,
            &Point2i::new(image.resolution.0 as i32, image.resolution.1 as i32),
        )?;
    }
    return Ok(());
}

pub fn load_spectrum_mipmap_cache(
    dir: &str,
    file_hash: &str,
) -> Result<MIPMap<RGBSpectrum>, PbrtError> {
    let dir = Path::new(&dir);
    let json_path = dir.join("cache_info.json");
    let cache_info = load_cacheinfo_from_json(json_path.to_str().unwrap())?;
    if cache_info.image_file_hash != file_hash {
        return Err(PbrtError::error("File hash of image is not same."));
    }
    let mut pyramid = Vec::with_capacity(cache_info.mip_levels as usize);
    for i in 0..cache_info.mip_levels {
        let mip_image_path = dir.join(format!("{}.exr", i));
        let (image, resolution) = read_cache_image(mip_image_path.to_str().unwrap())?;
        let image = convert_spectrum_image(image, &resolution);
        let mip_image = F32MIPMapImage::new(
            image,
            (resolution.x as usize, resolution.y as usize),
        );
        pyramid.push(mip_image);
    }
    let mipmap = MIPMap::<RGBSpectrum>::make_from_pyramid(
        pyramid,
        cache_info.texinfo.trilinear,
        cache_info.texinfo.max_aniso,
        cache_info.texinfo.swrap_mode,
        cache_info.texinfo.twrap_mode,
    );
    return Ok(mipmap);
}

pub fn save_spectrum_mipmap_cache(
    dir: &str,
    texinfo: &TexInfo,
    file_hash: &str,
    mipmap: &MIPMap<RGBSpectrum>,
) -> Result<(), PbrtError> {
    if !Path::new(dir).exists() {
        create_dir_all(dir)?;
    }
    let dir = Path::new(&dir);
    let json_path = dir.join("cache_info.json");
    let cache_info = TexCacheInfo {
        texinfo: texinfo.clone(),
        image_file_hash: String::from(file_hash),
        mip_levels: mipmap.levels() as u32,
    };
    save_cacheinfo_to_json(json_path.to_str().unwrap(), &cache_info)?;
    for i in 0..cache_info.mip_levels {
        let mip_image_path = dir.join(format!("{}.exr", i));
        let image = &mipmap.pyramid[i as usize];
        write_cache_image(
            mip_image_path.to_str().unwrap(),
            &image.data,
            &Point2i::new(image.resolution.0 as i32, image.resolution.1 as i32),
        )?;
    }
    return Ok(());
}
