use crate::core::pbrt::*;
use nom::character::complete::{alphanumeric1, digit1, multispace0};
use nom::error::*;
use nom::number::complete::*;
use nom::IResult;
use nom::{character, sequence};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/*
PF
512 256
-1
*/

fn read_word(s: &[u8]) -> IResult<&[u8], &[u8]> {
    let (s, word) = sequence::delimited(multispace0, alphanumeric1, multispace0)(s)?;
    return Ok((s, word));
}

fn read_value(s: &[u8]) -> IResult<&[u8], &[u8]> {
    let (s, word) = sequence::delimited(multispace0, recognize_float, multispace0)(s)?;
    return Ok((s, word));
}

fn read_cc(input: &[u8]) -> IResult<&[u8], &[u8]> {
    let (input, cc) = read_word(input)?;
    return Ok((input, cc));
}

fn read_header(input: &[u8]) -> IResult<&[u8], (u32, u32, f32)> {
    let (input, width) = read_value(input)?;
    let (input, height) = read_value(input)?;
    let (input, scale) = read_value(input)?;
    let width = std::str::from_utf8(width).unwrap().parse::<u32>().unwrap();
    let height = std::str::from_utf8(height).unwrap().parse::<u32>().unwrap();
    let scale = std::str::from_utf8(scale).unwrap().parse::<f32>().unwrap();
    return Ok((input, (width, height, scale)));
}

fn read_image_pfm_core(input: &[u8]) -> IResult<&[u8], (Vec<RGBSpectrum>, Point2i)> {
    let (input, cc) = read_cc(input)?;
    let cc = std::str::from_utf8(cc).unwrap();
    let (input, (width, height, scale)) = read_header(input)?;

    let n_channels;
    if cc == "Pf" {
        n_channels = 1;
    } else if cc == "PF" {
        n_channels = 3;
    } else {
        return Err(nom::Err::Error(Error::new(input, ErrorKind::Tag)));
    }

    let width = width as usize;
    let height = height as usize;
    let n_channels = n_channels as usize;

    //println!("width: {}, height: {}, scale: {}, n_channels: {}", width, height, scale, n_channels);

    // apply endian conversian and scale if appropriate
    let file_little_endian = scale < 0.0;
    // read the data
    let n_floats = n_channels * width * height;
    let mut data = vec![0.0; n_floats];
    {
        let mut input = input;
        for y in 0..height {
            let yy = height - y - 1;
            for x in 0..width {
                let index = (yy * width + x) as usize;
                for c in 0..n_channels {
                    let (inp, f) = if file_little_endian {
                        le_f32(input)?
                    } else {
                        be_f32(input)?
                    };
                    data[n_channels * index + c] = f;
                    input = inp;
                }
            }
        }
    }
    let scale = scale.abs();
    if scale != 1.0 {
        for f in &mut data {
            *f *= scale;
        }
    }

    let mut spcs = Vec::with_capacity(width as usize * height as usize);
    if n_channels == 1 {
        for f in data {
            spcs.push(RGBSpectrum::rgb_from_rgb(&[f, f, f]));
        }
    } else {
        for i in 0..(width * height) as usize {
            let r = data[i * 3];
            let g = data[i * 3 + 1];
            let b = data[i * 3 + 2];
            spcs.push(RGBSpectrum::rgb_from_rgb(&[r, g, b]));
        }
    }
    assert!(spcs.len() == (width * height) as usize);
    return Ok((input, (spcs, Point2i::from((width as i32, height as i32)))));
}

pub fn read_image_pfm(name: &str) -> Result<(Vec<RGBSpectrum>, Point2i), PbrtError> {
    let path = Path::new(name);
    let mut reader = BufReader::new(File::open(path)?);
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes)?;

    match read_image_pfm_core(&bytes) {
        Ok((_, v)) => {
            return Ok(v);
        }
        Err(e) => {
            let msg = format!("Error reading PFM file \"{}\": {}", path.display(), e);
            return Err(PbrtError::error(&msg));
        }
    }
}
