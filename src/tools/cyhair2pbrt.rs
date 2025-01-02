use clap::*;

use core::num;
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
        let has_segments = (flags & 1) != 0;
        //let has_points = (flags & 2) != 0;
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
        {
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

fn to_cubic_bezier_curves(cyhair: &CyHair, max_strands: Option<usize>, user_thickness: Option<f32>) -> Result<(Vec<f32>, Vec<f32>), CyHairError> {
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

        let mut segment_points = Vec::new();
        for k in 0..num_segments {
            // Zup -> Yup
            let x = cyhair.points[3 * (cyhair.strand_offsets[i] + k) + 0];//0
            let y = cyhair.points[3 * (cyhair.strand_offsets[i] + k) + 2];//1
            let z = cyhair.points[3 * (cyhair.strand_offsets[i] + k) + 1];//
            segment_points.push((x, y, z));
        }
    }

    // Skip both endpoints
    for s in 1..num_strands-1 {
        
    }


    return Ok((vertices, radiuss));
}

fn write_pbrt(
    writer: &mut dyn std::io::Write,
    vertices: &Vec<f32>,
    radiuss: &Vec<f32>,
) -> Result<(), CyHairError> {
    let num_curves = radiuss.len() / 4;
    for i in 0..num_curves {
        writer.write_all(b"Shape \"curve\" \"string type\" \"cylinder\" \"point P\" [")?;
        for j in 0..12 {
            let idx = i * 12 + j;
            writer.write_all(format!("{} ", vertices[idx]).as_bytes())?;
        }
    }

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
    let thickness = opts.user_thickness;

    if infile.is_none() {
        eprintln!(
            "usage: cyhair2pbrt [CyHair filename] [pbrt output filename] (max strands) (thickness)"
        );
        process::exit(1);
    }

    let cyhair = match CyHair::load(infile.as_ref().unwrap()) {
        Ok(cyhair) => cyhair,
        Err(err) => {
            eprintln!("Failed to load CyHair file [ {} ]", err);
            process::exit(1);
        }
    };

    let (vertices, radiuss) = match to_cubic_bezier_curves(&cyhair, max_strands, thickness) {
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

    match write_pbrt(&mut *writer, &vertices, &radiuss, ) {
        Ok(_) => {}
        Err(err) => {
            eprintln!("Failed to write pbrt file [ {} ]", err);
            process::exit(1);
        }
    }
    //

    process::exit(0);
}
