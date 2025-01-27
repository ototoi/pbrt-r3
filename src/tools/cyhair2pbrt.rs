use clap::*;

use std::ops::{Add, Mul};
use std::path::Path;
use std::path::PathBuf;
use std::process;

use nom::number::complete::*;
use nom::IResult;

pub enum CyHairError {
    Message(String),
    Error(Box<dyn std::error::Error>),
}

impl std::fmt::Display for CyHairError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            CyHairError::Message(msg) => write!(f, "{}", msg),
            CyHairError::Error(err) => write!(f, "{}", err),
        }
    }
}

impl From<std::io::Error> for CyHairError {
    fn from(error: std::io::Error) -> Self {
        return CyHairError::Error(Box::new(error));
    }
}

impl From<nom::error::Error<&[u8]>> for CyHairError {
    fn from(error: nom::error::Error<&[u8]>) -> Self {
        return CyHairError::Message(format!("Parse error: {:?}", error));
    }
}

impl From<Box<dyn std::error::Error>> for CyHairError {
    fn from(error: Box<dyn std::error::Error>) -> Self {
        return CyHairError::Error(error);
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct CyHairHeader {
    magic: [u8; 4],
    num_strands: u32,
    total_points: u32,
    flags: u32,
    default_segments: u32,
    default_thickness: f32,
    default_transparency: f32,
    default_color: [f32; 3],
    infomation: [u8; 88],
}

#[derive(Debug)]
struct CyHair {
    segments: Vec<u16>,
    points: Vec<f32>,
    thicknesses: Vec<f32>,
    transparencies: Vec<f32>,
    colors: Vec<f32>,

    header: CyHairHeader,

    // Processed CyHair values
    strand_offsets: Vec<usize>,
}

impl CyHair {
    pub fn load(path: &Path) -> Result<CyHair, CyHairError> {
        let file = std::fs::File::open(path)?;
        let mut reader = std::io::BufReader::new(file);
        let header = CyHair::load_header(&mut reader)?;
        //println!("header: {:?}", header);
        let magic = header.magic.iter().map(|&c| c as char).collect::<String>();
        if magic != "HAIR" {
            return Err(CyHairError::Message("Not a CyHair file".to_string()));
        }
        return CyHair::load_strands(&mut reader, &header);
    }

    fn load_header(reader: &mut dyn std::io::Read) -> Result<CyHairHeader, CyHairError> {
        let mut buffer: [u8; 128] = [0; 128];
        reader.read_exact(&mut buffer)?;

        let input = &buffer[..];
        match CyHair::parse_header(input) {
            Ok((_, header)) => return Ok(header),
            Err(err) => return Err(CyHairError::Message(format!("Parse error: {:?}", err))),
        }
    }

    fn parse_header(input: &[u8]) -> IResult<&[u8], CyHairHeader> {
        //let original_input = input;
        let infomation: [u8; 88] = input[128 - 88..].try_into().unwrap();
        let (input, magic0) = le_u8(input)?;
        let (input, magic1) = le_u8(input)?;
        let (input, magic2) = le_u8(input)?;
        let (input, magic3) = le_u8(input)?;

        let (input, num_strands) = le_u32(input)?;
        let (input, total_points) = le_u32(input)?;
        let (input, flags) = le_u32(input)?;
        let (input, default_segments) = le_u32(input)?;
        let (input, default_thickness) = le_f32(input)?;
        let (input, default_transparency) = le_f32(input)?;
        let (input, default_color0) = le_f32(input)?;
        let (input, default_color1) = le_f32(input)?;
        let (input, default_color2) = le_f32(input)?;

        let has_segments = (flags & 1) != 0;
        let has_points = (flags & 2) != 0;
        let has_thickness = (flags & 4) != 0;
        let has_transparency = (flags & 8) != 0;
        let has_color = (flags & 16) != 0;

        //eprintln!("magic = [{}, {}, {}, {}]", magic0, magic1, magic2, magic3);
        eprintln!(
            "magic = {}",
            String::from_utf8_lossy(&[magic0, magic1, magic2, magic3])
        );
        eprintln!("flags = {}", flags);
        eprintln!("  has_segments = {}", has_segments);
        eprintln!("  has_points = {}", has_points);
        eprintln!("  has_thickness = {}", has_thickness);
        eprintln!("  has_transparency = {}", has_transparency);
        eprintln!("  has_color = {}", has_color);
        eprintln!("num_strands = {}", num_strands);
        eprintln!("total_points = {}", total_points);

        eprintln!("default_segments = {}", default_segments);
        eprintln!("default_thickness = {}", default_thickness);
        eprintln!("default_transparency = {}", default_transparency);
        eprintln!(
            "default_color = [{}, {}, {}]",
            default_color0, default_color1, default_color2
        );

        return Ok((
            input,
            CyHairHeader {
                magic: [magic0, magic1, magic2, magic3],
                num_strands: num_strands,
                total_points: total_points,
                flags: flags,
                default_segments: default_segments,
                default_thickness: default_thickness,
                default_transparency: default_transparency,
                default_color: [default_color0, default_color1, default_color2],
                infomation: infomation,
            },
        ));
    }

    fn load_strands(
        reader: &mut dyn std::io::Read,
        header: &CyHairHeader,
    ) -> Result<CyHair, CyHairError> {
        let default_segments = header.default_segments;

        let flags = header.flags;
        let has_segments = (flags & 1) != 0;
        let has_points = (flags & 2) != 0;
        if !has_points {
            return Err(CyHairError::Message("No point data in CyHair.".to_string()));
        }

        if default_segments < 1 && !has_segments {
            return Err(CyHairError::Message(
                "No valid segment information in CyHair.".to_string(),
            ));
        }

        let mut buffer = Vec::new();
        reader.read_to_end(&mut buffer)?;
        match CyHair::parse_strands(&buffer, header) {
            Ok((_, mut cyhair)) => {
                let num_strands = cyhair.header.num_strands as usize;
                let mut strand_offsets = vec![0; num_strands];
                for i in 1..num_strands {
                    let num_segments = if cyhair.segments.is_empty() {
                        cyhair.header.default_segments as usize
                    } else {
                        cyhair.segments[i - 1] as usize
                    };
                    strand_offsets[i] = strand_offsets[i - 1] as usize + (num_segments + 1);
                }
                cyhair.strand_offsets = strand_offsets;
                return Ok(cyhair);
            }
            Err(err) => return Err(CyHairError::Message(format!("Parse error: {:?}", err))),
        }
    }

    fn parse_strands<'a>(input: &'a [u8], header: &'a CyHairHeader) -> IResult<&'a [u8], CyHair> {
        let flags = header.flags;
        //println!("flags = {}", flags);
        let has_segments = (flags & 1) != 0;
        let has_points = (flags & 2) != 0;
        let has_thickness = (flags & 4) != 0;
        let has_transparency = (flags & 8) != 0;
        let has_color = (flags & 16) != 0;

        let num_strands = header.num_strands as usize;
        let total_points = header.total_points as usize;

        let mut input = input;
        // First read all strand data from a file.
        let mut segments = Vec::new();
        if has_segments {
            let (_input, _segments) = nom::multi::count(le_u16, num_strands)(input)?;
            input = _input;
            segments = _segments;
        }
        let mut points = Vec::new();
        if has_points {
            let (_input, _points) = nom::multi::count(le_f32, total_points * 3)(input)?;
            input = _input;
            points = _points;
        }
        let mut thicknesses = Vec::new();
        if has_thickness {
            let (_input, _thicknesses) = nom::multi::count(le_f32, total_points)(input)?;
            input = _input;
            thicknesses = _thicknesses;
        }
        let mut transparencies = Vec::new();
        if has_transparency {
            let (_input, _transparencies) = nom::multi::count(le_f32, total_points)(input)?;
            input = _input;
            transparencies = _transparencies;
        }
        let mut colors = Vec::new();
        if has_color {
            let (_input, _colors) = nom::multi::count(le_f32, total_points * 3)(input)?;
            input = _input;
            colors = _colors;
        }

        let cyhair = CyHair {
            segments: segments,
            points: points,
            thicknesses: thicknesses,
            transparencies: transparencies,
            colors: colors,
            header: *header,
            strand_offsets: Vec::new(),
        };
        return Ok((input, cyhair));
    }
}

//type Real3 = (f32, f32, f32);
#[derive(Debug, Default, Copy, Clone, PartialEq)]
struct Real3(pub (f32, f32, f32));
impl From<(f32, f32, f32)> for Real3 {
    fn from(t: (f32, f32, f32)) -> Self {
        Real3(t)
    }
}

impl Add for Real3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Real3((
            self.0 .0 + other.0 .0,
            self.0 .1 + other.0 .1,
            self.0 .2 + other.0 .2,
        ))
    }
}

impl Mul<Real3> for f32 {
    type Output = Real3;
    fn mul(self, other: Real3) -> Real3 {
        Real3((self * other.0 .0, self * other.0 .1, self * other.0 .2))
    }
}

impl Mul<f32> for Real3 {
    type Output = Self;
    fn mul(self, other: f32) -> Self {
        Real3((self.0 .0 * other, self.0 .1 * other, self.0 .2 * other))
    }
}

impl Mul<Real3> for Real3 {
    type Output = Self;
    fn mul(self, other: Real3) -> Self {
        Real3((
            self.0 .0 * other.0 .0,
            self.0 .1 * other.0 .1,
            self.0 .2 * other.0 .2,
        ))
    }
}

const TO_C2B: [[f32; 4]; 4] = [
    [0.0, 6.0 / 6.0, 0.0, 0.0],
    [-1.0 / 6.0, 6.0 / 6.0, 1.0 / 6.0, 0.0],
    [0.0, 1.0 / 6.0, 6.0 / 6.0, -1.0 / 6.0],
    [0.0, 0.0, 6.0 / 6.0, 0.0],
];

const TO_C2B0: [[f32; 4]; 4] = [
    [0.0, 6.0 / 6.0, 0.0, 0.0],
    [0.0, 3.0 / 6.0, 4.0 / 6.0, -1.0 / 6.0],
    [0.0, 1.0 / 6.0, 6.0 / 6.0, -1.0 / 6.0],
    [0.0, 0.0, 6.0 / 6.0, 0.0],
];

const TO_C2B1: [[f32; 4]; 4] = [
    [0.0, 6.0 / 6.0, 0.0, 0.0],
    [-1.0 / 6.0, 6.0 / 6.0, 1.0 / 6.0, 0.0],
    [-1.0 / 6.0, 4.0 / 6.0, 3.0 / 6.0, 0.0],
    [0.0, 0.0, 6.0 / 6.0, 0.0],
];

fn mul_matrix(mat: &[[f32; 4]; 4], pt: &[Real3; 4]) -> [Real3; 4] {
    let mut out = [Real3::default(); 4];
    for i in 0..4 {
        out[i] = mat[i][0] * pt[0] + mat[i][1] * pt[1] + mat[i][2] * pt[2] + mat[i][3] * pt[3];
    }
    return out;
}

fn camull_rom_to_cubic_bezier(cps: &[Real3], seg_idx: usize) -> [Real3; 4] {
    let cps_size = cps.len();
    if cps_size == 2 {
        let i0 = seg_idx;
        let i1 = seg_idx + 1;
        let p0 = cps[i0];
        let p1 = cps[i0] * (2.0 / 3.0) + cps[i1] * (1.0 / 3.0);
        let p2 = cps[i0] * (1.0 / 3.0) + cps[i1] * (2.0 / 3.0);
        let p3 = cps[i1];
        return [p0, p1, p2, p3];
    } else {
        if seg_idx == 0 {
            let i0 = seg_idx;
            let i1 = seg_idx + 1;
            let i2 = seg_idx + 2;
            let p = [Real3::default(), cps[i0], cps[i1], cps[i2]];
            return mul_matrix(&TO_C2B0, &p);
        } else if seg_idx == cps_size - 2 {
            let i0 = seg_idx - 1;
            let i1 = seg_idx;
            let i2 = seg_idx + 1;
            let p = [cps[i0], cps[i1], cps[i2], Real3::default()];
            return mul_matrix(&TO_C2B1, &p);
        } else {
            let i0 = seg_idx - 1;
            let i1 = seg_idx;
            let i2 = seg_idx + 1;
            let i3 = seg_idx + 2;
            let p = [cps[i0], cps[i1], cps[i2], cps[i3]];
            return mul_matrix(&TO_C2B, &p);
        }
    }
}

fn to_cubic_bezier_curves(
    cyhair: &CyHair,
    vertex_scale: Real3,
    vertex_translate: Real3,
    max_strands: Option<usize>,
    user_thickness: Option<f32>,
) -> Result<(Vec<f32>, Vec<f32>), CyHairError> {
    if cyhair.points.is_empty() || cyhair.strand_offsets.is_empty() {
        return Err(CyHairError::Message("No valid CyHair data.".to_string()));
    }
    let mut vertices = Vec::new();
    let mut radiuss = Vec::new();

    let mut num_strands = cyhair.header.num_strands as usize;
    if let Some(max_strands) = max_strands {
        num_strands = usize::min(num_strands, max_strands);
    }

    for i in 0..num_strands {
        let num_segments = if cyhair.segments.is_empty() {
            cyhair.header.default_segments as usize
        } else {
            cyhair.segments[i] as usize
        };

        if num_segments < 2 {
            continue;
        }

        let mut segment_points: Vec<Real3> = Vec::new();
        for k in 0..num_segments {
            // Zup -> Yup
            let x = cyhair.points[3 * (cyhair.strand_offsets[i] + k) + 0]; //0
            let y = cyhair.points[3 * (cyhair.strand_offsets[i] + k) + 2]; //2
            let z = cyhair.points[3 * (cyhair.strand_offsets[i] + k) + 1]; //1
            let p = (x, y, z).into();
            segment_points.push(p);
        }

        // Skip both endpoints
        for s in 1..num_segments - 1 {
            let seg_idx = s - 1;
            let q = camull_rom_to_cubic_bezier(&segment_points, seg_idx);
            for j in 0..4 {
                let qq = q[j];
                let qqt = qq * vertex_scale + vertex_translate;
                let qqt: (f32, f32, f32) = qqt.0;
                vertices.push(qqt.0);
                vertices.push(qqt.1);
                vertices.push(qqt.2);
            }

            if let Some(user_thickness) = user_thickness {
                for _ in 0..4 {
                    radiuss.push(user_thickness);
                }
            } else {
                let thickness = cyhair.header.default_thickness; //todo
                for _ in 0..4 {
                    radiuss.push(thickness);
                }
            }
        }
    }
    return Ok((vertices, radiuss));
}

fn write_pbrt(
    writer: &mut dyn std::io::Write,
    vertices: &Vec<f32>,
    radiuss: &Vec<f32>,
    input_file: &str,
    user_thickness: f32,
) -> Result<(), CyHairError> {
    let mut bounds: [Real3; 2] = [
        Real3::from((1e30, 1e30, 1e30)),
        Real3::from((-1e30, -1e30, -1e30)),
    ];
    for i in 0..vertices.len() / 3 {
        let thickness = radiuss[i];
        let x = vertices[3 * i + 0];
        let y = vertices[3 * i + 1];
        let z = vertices[3 * i + 2];

        bounds[0].0 .0 = f32::min(bounds[0].0 .0, x - thickness);
        bounds[0].0 .1 = f32::min(bounds[0].0 .1, y - thickness);
        bounds[0].0 .2 = f32::min(bounds[0].0 .2, z - thickness);
        bounds[1].0 .0 = f32::max(bounds[1].0 .0, x + thickness);
        bounds[1].0 .1 = f32::max(bounds[1].0 .1, y + thickness);
        bounds[1].0 .2 = f32::max(bounds[1].0 .2, z + thickness);
    }
    let num_strands = radiuss.len() / 4;
    writer.write_all(format!("# Converted from \"{}\" by cyhair2pbrt\n", input_file).as_bytes())?;
    writer.write_all(
        format!(
            "# The number of strands = {}. user_thickness = {}\n",
            num_strands, user_thickness
        )
        .as_bytes(),
    )?;
    let x0 = bounds[0].0 .0;
    let y0 = bounds[0].0 .1;
    let z0 = bounds[0].0 .2;
    let x1 = bounds[1].0 .0;
    let y1 = bounds[1].0 .1;
    let z1 = bounds[1].0 .2;
    writer.write_all(
        format!(
            "# Scene bounds: ({}, {}, {}) - ({}, {}, {})\n\n\n",
            x0, y0, z0, x1, y1, z1
        )
        .as_bytes(),
    )?;

    let num_curves = radiuss.len() / 4;
    for i in 0..num_curves {
        writer.write_all(b"Shape \"curve\" \"string type\" [ \"cylinder\" ] \"point P\" [ ")?;
        for j in 0..12 {
            let idx = i * 12 + j;
            writer.write_all(format!("{} ", vertices[idx]).as_bytes())?;
        }
        writer.write_all(
            format!(
                " ] \"float width0\" [ {} ] \"float width1\" [ {} ]\n",
                radiuss[4 * i + 0],
                radiuss[4 * i + 3]
            )
            .as_bytes(),
        )?;
    }
    writer.flush()?;
    return Ok(());
}

#[derive(Debug, Parser)]
#[clap(author, about, version, disable_help_flag = true)]
struct CommandOptions {
    /// Input .cyhair file.
    #[arg(short, long, value_name = "filename")]
    pub infile: Option<PathBuf>,

    /// Write the final image to the given filename.
    #[arg(short, long, value_name = "filename")]
    pub outfile: Option<PathBuf>,

    /// Print this help text.
    #[arg(short, long, action = clap::ArgAction::HelpLong)]
    pub help: Option<bool>,

    /// Number of strands to generate.
    #[arg(long, value_name = "count")]
    pub max_strands: Option<usize>,

    ///Tickness of the hair strands.
    #[arg(long, value_name = "thickness")]
    pub user_thickness: Option<f32>,

    #[arg(value_name = "options")]
    pub options: Option<Vec<String>>,
}

pub fn main() {
    let opts = CommandOptions::parse();

    let mut other_options = Vec::new();
    if let Some(args) = opts.options.as_ref() {
        for arg in args {
            other_options.push(arg.clone());
        }
    }

    let mut infile: Option<PathBuf> = None;
    if let Some(path) = opts.infile.as_ref() {
        infile = Some(path.clone());
    } else if other_options.len() > 0 {
        infile = Some(PathBuf::from(other_options.remove(0)));
    }
    let mut outfile: Option<PathBuf> = None;
    if let Some(path) = opts.outfile.as_ref() {
        outfile = Some(path.clone());
    } else if other_options.len() > 0 {
        outfile = Some(PathBuf::from(other_options.remove(0)));
    }

    let max_strands = opts.max_strands;
    let user_thickness = opts.user_thickness;

    if infile.is_none() {
        eprintln!(
            "usage: cyhair2pbrt [CyHair filename] [pbrt output filename] (max strands) (thickness)"
        );
        process::exit(1);
    }

    let infile = infile.as_ref().unwrap();

    let cyhair = match CyHair::load(&infile) {
        Ok(cyhair) => cyhair,
        Err(err) => {
            eprintln!("Failed to load CyHair file [ {} ]", err);
            process::exit(1);
        }
    };

    let vertex_scale = Real3((1.0, 1.0, 1.0));
    let vertex_translate = Real3((0.0, 0.0, 0.0));

    let (vertices, radiuss) = match to_cubic_bezier_curves(
        &cyhair,
        vertex_scale,
        vertex_translate,
        max_strands,
        user_thickness,
    ) {
        Ok((vertices, radiuss)) => (vertices, radiuss),
        Err(_err) => {
            eprintln!("Failed to convert CyHair data");
            process::exit(1);
        }
    };

    let mut writer: Box<dyn std::io::Write> = match outfile {
        Some(path) => {
            let file = std::fs::File::create(path).unwrap();
            let writer = std::io::BufWriter::new(file);
            Box::new(writer)
        }
        None => Box::new(std::io::stdout()),
    };

    let basename = infile.file_name().unwrap().to_str().unwrap();
    let user_thickness = user_thickness.unwrap_or(-1.0);
    match write_pbrt(&mut *writer, &vertices, &radiuss, basename, user_thickness) {
        Ok(_) => {}
        Err(err) => {
            eprintln!("Failed to write pbrt file [ {} ]", err);
            process::exit(1);
        }
    }

    eprintln!("Converted {} strands.", radiuss.len() / 4);

    process::exit(0);
}
