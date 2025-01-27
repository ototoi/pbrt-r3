use super::common::*;
use super::read_file::{read_file_with_include, read_file_without_include};
use super::remove_comment::remove_comment;
use crate::core::api::ParseContext;
use crate::core::param_set::ParamSet;
use crate::core::pbrt::PbrtError;

use nom::bytes;
use nom::character;
use nom::multi;
use nom::number;
use nom::sequence;
use nom::IResult;

pub fn pbrt_parse_file(filename: &str, context: &mut dyn ParseContext) -> Result<(), PbrtError> {
    let s = read_file_with_include(filename)?;
    return pbrt_parse_string_core(&s, context);
}

pub fn pbrt_parse_string(s: &str, context: &mut dyn ParseContext) -> Result<(), PbrtError> {
    let ops = parse_opnodes(s)?;
    return evaluate_opnodes(&ops, context);
}

pub fn pbrt_parse_file_without_include(
    filename: &str,
    context: &mut dyn ParseContext,
) -> Result<(), PbrtError> {
    let s = read_file_without_include(filename)?;
    return pbrt_parse_string_core(&s, context);
}
//-----------------------------------

pub fn pbrt_parse_string_core(s: &str, context: &mut dyn ParseContext) -> Result<(), PbrtError> {
    let ops = parse_opnodes_core(s)?;
    return evaluate_opnodes(&ops, context);
}

fn parse_opnodes(s: &str) -> Result<Vec<OPNode>, PbrtError> {
    let r = remove_comment(s);
    match r {
        Ok((_, s)) => {
            return parse_opnodes_core(&s);
        }
        Err(e) => {
            return Err(PbrtError::from(e.to_string()));
        }
    }
}

fn parse_opnodes_core(s: &str) -> Result<Vec<OPNode>, PbrtError> {
    let result = nom::combinator::all_consuming(multi::many0(sequence::delimited(
        space0,
        parse_operation,
        space0,
    )))(s);
    match result {
        Ok((_, nodes)) => {
            return Ok(nodes);
        }
        Err(e) => {
            return Err(PbrtError::from(e.to_string()));
        }
    }
}

struct OPNode {
    pub name: String,
    pub args: Option<ParamSet>,
    pub params: Option<ParamSet>,
}

impl OPNode {
    pub fn new(name: &str, args: Option<ParamSet>, params: Option<ParamSet>) -> Self {
        OPNode {
            name: String::from(name),
            args,
            params,
        }
    }
}

fn evaluate_opnodes(ops: &[OPNode], context: &mut dyn ParseContext) -> Result<(), PbrtError> {
    for op in ops {
        //println!("{}", op.name);
        let opname: &str = &op.name;
        match opname {
            "Identity" => {
                //fn pbrt_identity(&mut self);
                context.pbrt_identity();
            }
            "Translate" => {
                //fn pbrt_translate(&mut self, dx: Float, dy: Float, dz: Float);
                let args = op.args.as_ref().unwrap();
                let vec = args.get_floats("args");
                if vec.len() != 3 {
                    let msg = format!("{} required {} arguments", opname, 3);
                    return Err(PbrtError::error(&msg));
                }
                context.pbrt_translate(vec[0], vec[1], vec[2]);
            }
            "Rotate" => {
                //fn pbrt_rotate(&mut self, angle: Float, ax: Float, ay: Float, az: Float);
                let args = op.args.as_ref().unwrap();
                let vec = args.get_floats("args");
                if vec.len() != 4 {
                    let msg = format!("{} required {} arguments", opname, 4);
                    return Err(PbrtError::error(&msg));
                }
                context.pbrt_rotate(vec[0], vec[1], vec[2], vec[3]);
            }
            "Scale" => {
                //fn pbrt_scale(&mut self, sx: Float, sy: Float, sz: Float);
                let args = op.args.as_ref().unwrap();
                let vec = args.get_floats("args");
                if vec.len() != 3 {
                    let msg = format!("{} required {} arguments", opname, 3);
                    return Err(PbrtError::error(&msg));
                }
                context.pbrt_scale(vec[0], vec[1], vec[2]);
            }
            "LookAt" => {
                /*
                fn pbrt_look_at(
                    &mut self,
                    ex: Float,
                    ey: Float,
                    ez: Float,
                    lx: Float,
                    ly: Float,
                    lz: Float,
                    ux: Float,
                    uy: Float,
                    uz: Float,
                );
                */
                let args = op.args.as_ref().unwrap();
                let vec = args.get_floats("args");
                if vec.len() != 9 {
                    let msg = format!("{} required {} arguments", opname, 9);
                    return Err(PbrtError::error(&msg));
                }
                context.pbrt_look_at(
                    vec[0], vec[1], vec[2], vec[3], vec[4], vec[5], vec[6], vec[7], vec[8],
                );
            }
            "ConcatTransform" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_floats("arg1");
                context.pbrt_concat_transform(&vec);
            }
            "Transform" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_floats("arg1");
                context.pbrt_transform(&vec);
            }
            "CoordinateSystem" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                context.pbrt_coordinate_system(name);
            }
            "CoordSysTransform" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                context.pbrt_coord_sys_transform(name);
            }
            "ActiveTransformAll" => {
                context.pbrt_active_transform_all();
            }
            "ActiveTransformEndTime" => {
                context.pbrt_active_transform_end_time();
            }
            "ActiveTransformStartTime" => {
                context.pbrt_active_transform_start_time();
            }
            /*
            fn pbrt_transform_times(&mut self, start: Float, end: Float);
            */
            "PixelFilter" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_pixel_filter(name, params);
            }
            "Film" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_film(name, params);
            }
            "Sampler" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_sampler(name, params);
            }
            "Accelerator" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_accelerator(name, params);
            }
            "Integrator" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_integrator(name, params);
            }
            "Camera" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_camera(name, params);
            }
            "MakeNamedMedium" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_make_named_medium(name, params);
            }
            "MediumInterface" => {
                let args = op.args.as_ref().unwrap();
                let vec1 = args.get_strings("arg1");
                let vec2 = args.get_strings("arg2");
                let inside_name = vec1.first().unwrap();
                let outside_name = vec2.first().unwrap();
                context.pbrt_medium_interface(inside_name, outside_name);
            }
            "WorldBegin" => {
                context.pbrt_world_begin();
            }
            "AttributeBegin" => {
                context.pbrt_attribute_begin();
            }
            "AttributeEnd" => {
                context.pbrt_attribute_end();
            }
            "TransformBegin" => {
                context.pbrt_transform_begin();
            }
            "TransformEnd" => {
                context.pbrt_transform_end();
            }
            "Texture" => {
                let args = op.args.as_ref().unwrap();
                let name = String::from(args.get_strings("arg1").first().unwrap());
                let tp = String::from(args.get_strings("arg2").first().unwrap());
                let tex_name = String::from(args.get_strings("arg3").first().unwrap());
                let params = op.params.as_ref().unwrap();
                context.pbrt_texture(&name, &tp, &tex_name, params);
            }
            "Material" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_material(name, params);
            }
            "MakeNamedMaterial" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_make_named_material(name, params);
            }
            "NamedMaterial" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                context.pbrt_named_material(name);
            }
            "LightSource" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_light_source(name, params);
            }
            "AreaLightSource" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_area_light_source(name, params);
            }
            "Shape" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_shape(name, params);
            }
            "ReverseOrientation" => {
                context.pbrt_reverse_orientation();
            }
            "ObjectBegin" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                context.pbrt_object_begin(name);
            }
            "ObjectEnd" => {
                context.pbrt_object_end();
            }
            "ObjectInstance" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let name: &str = vec.first().unwrap();
                context.pbrt_object_instance(name);
            }
            "WorldEnd" => {
                context.pbrt_world_end();
            }
            "WorkDirBegin" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let path: &str = vec.first().unwrap();
                context.pbrt_work_dir_begin(path);
            }
            "WorkDirEnd" => {
                context.pbrt_work_dir_end();
            }
            "Include" => {
                let args = op.args.as_ref().unwrap();
                let vec = args.get_strings("arg1");
                let filename: &str = vec.first().unwrap();
                let params = op.params.as_ref().unwrap();
                context.pbrt_include(filename, params);
            }
            _ => {
                let msg = format!("Unexpected token: {}", opname);
                return Err(PbrtError::error(&msg));
            }
        }
    }
    return Ok(());
}

/*
#[warn(dead_code)]
fn parse_opnodes_debug(_: &str) -> Result<Vec<OPNode>, String> {
    let mut ops = Vec::<OPNode>::new();
    ops.push(OPNode::new("Identity", None, None));
    ops.push(OPNode::new("Identity", None, None));
    ops.push(OPNode::new("WorldBegin", None, None));
    ops.push(OPNode::new("WorldEnd", None, None));
    {
        let mut args = ParamSet::new();
        args.add_string("name", "a_object");
        ops.push(OPNode::new("ObjectBegin", Some(args), None));
    }
    ops.push(OPNode::new("ObjectEnd", None, None));
    {
        let mut args = ParamSet::new();
        args.add_string("name", "a_texture");
        args.add_string("type", "a_type");
        args.add_string("tex_name", "a_tex_name");
        let params = ParamSet::new();
        ops.push(OPNode::new("Texture", Some(args), Some(params)));
    }
    return Ok(ops);
}
*/

//-----------------------------------

fn parse_operation(s: &str) -> IResult<&str, OPNode> {
    return nom::branch::alt((
        nom::branch::alt((
            parse_identity,
            parse_translate,
            parse_rotate,
            parse_scale,
            parse_look_at,
            parse_concat_transform,
            parse_transform,
            parse_coordinate_system,
            parse_coord_sys_transform,
            parse_active_transform,
            parse_transform_times,
        )),
        nom::branch::alt((
            parse_pixel_filter,
            parse_film,
            parse_sampler,
            parse_accelerator,
            parse_integrator,
            parse_camera,
            parse_make_named_medium,
            parse_medium_interface,
        )),
        nom::branch::alt((
            parse_world_begin,
            parse_attribute_begin,
            parse_attribute_end,
            parse_transform_begin,
            parse_transform_end,
            parse_texture,
            parse_material,
            parse_make_named_material,
            parse_named_material,
            parse_light_source,
            parse_area_light_source,
            parse_shape,
            parse_reverse_orientation,
            parse_object_begin,
            parse_object_end,
            parse_object_instance,
            parse_world_end,
            parse_include,
        )),
        nom::branch::alt((parse_work_dir_begin, parse_work_dir_end)),
    ))(s);
}

fn parse_op_void<'a>(s: &'a str, name: &str) -> IResult<&'a str, OPNode> {
    let (s, _) = sequence::terminated(bytes::complete::tag(name), space0)(s)?;
    return Ok((s, OPNode::new(name, None, None)));
}

fn parse_op_float_n<'a>(s: &'a str, opname: &str, n: usize) -> IResult<&'a str, OPNode> {
    let (s, (op, a)) = nom::branch::permutation((
        sequence::terminated(bytes::complete::tag(opname), space1),
        multi::count(
            sequence::terminated(number::complete::recognize_float, space0),
            n,
        ),
    ))(s)?;
    let mut args = ParamSet::new();
    let v: Vec<f32> = a.iter().map(|x| (*x).parse::<f32>().unwrap()).collect();
    args.add_floats("args", &v);
    return Ok((s, OPNode::new(op, Some(args), None)));
}

fn parse_op_floats<'a>(s: &'a str, opname: &str) -> IResult<&'a str, OPNode> {
    let (s, (op, a)) = nom::branch::permutation((
        sequence::terminated(bytes::complete::tag(opname), space1),
        sequence::delimited(
            character::complete::char('['),
            sequence::delimited(
                character::complete::multispace0,
                multi::separated_list1(
                    character::complete::multispace1,
                    number::complete::recognize_float,
                ),
                character::complete::multispace0,
            ),
            character::complete::char(']'),
        ),
    ))(s)?;
    let mut args = ParamSet::new();
    let v: Vec<f32> = a.iter().map(|x| (*x).parse::<f32>().unwrap()).collect();
    args.add_floats("arg1", &v);
    return Ok((s, OPNode::new(op, Some(args), None)));
}

fn parse_op_string<'a>(s: &'a str, opname: &'a str) -> IResult<&'a str, OPNode> {
    let (s, (op, name)) = nom::branch::permutation((
        sequence::terminated(bytes::complete::tag(opname), space1),
        sequence::terminated(string_literal, space0),
    ))(s)?;
    let mut args = ParamSet::new();
    args.add_string("arg1", name);
    return Ok((s, OPNode::new(op, Some(args), None)));
}

fn parse_op_string_string<'a>(s: &'a str, opname: &str) -> IResult<&'a str, OPNode> {
    let (s, (op, b, c)) = nom::branch::permutation((
        sequence::terminated(bytes::complete::tag(opname), space1),
        sequence::terminated(string_literal, space1),
        sequence::terminated(string_literal, space0),
    ))(s)?;
    let mut args = ParamSet::new();
    args.add_string("arg1", b);
    args.add_string("arg2", c);
    return Ok((s, OPNode::new(op, Some(args), None)));
}

fn parse_op_string_params<'a>(s: &'a str, opname: &str) -> IResult<&'a str, OPNode> {
    let (s, (op, a, params)) = nom::branch::permutation((
        sequence::terminated(bytes::complete::tag(opname), space1),
        sequence::terminated(string_literal, space0),
        parse_params,
    ))(s)?;
    let mut args = ParamSet::new();
    args.add_string("arg1", a);
    return Ok((s, OPNode::new(op, Some(args), Some(params))));
}

fn parse_op_string_string_string_params<'a>(s: &'a str, opname: &str) -> IResult<&'a str, OPNode> {
    let (s, (op, a, params)) = nom::branch::permutation((
        sequence::terminated(bytes::complete::tag(opname), space1),
        multi::count(sequence::terminated(string_literal, space0), 3),
        parse_params,
    ))(s)?;
    let mut args = ParamSet::new();
    args.add_string("arg1", a[0]);
    args.add_string("arg2", a[1]);
    args.add_string("arg3", a[2]);
    return Ok((s, OPNode::new(op, Some(args), Some(params))));
}

fn parse_identity(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "Identity");
}
//fn pbrt_translate(&mut self, dx: Float, dy: Float, dz: Float);
fn parse_translate(s: &str) -> IResult<&str, OPNode> {
    return parse_op_float_n(s, "Translate", 3);
}

//fn pbrt_rotate(&mut self, angle: Float, ax: Float, ay: Float, az: Float);
fn parse_rotate(s: &str) -> IResult<&str, OPNode> {
    return parse_op_float_n(s, "Rotate", 4);
}

//fn pbrt_scale(&mut self, sx: Float, sy: Float, sz: Float);
fn parse_scale(s: &str) -> IResult<&str, OPNode> {
    return parse_op_float_n(s, "Scale", 3);
}

fn parse_look_at(s: &str) -> IResult<&str, OPNode> {
    return parse_op_float_n(s, "LookAt", 9);
}

fn parse_concat_transform(s: &str) -> IResult<&str, OPNode> {
    return parse_op_floats(s, "ConcatTransform");
}

fn parse_transform(s: &str) -> IResult<&str, OPNode> {
    return parse_op_floats(s, "Transform");
}

//fn pbrt_coordinate_system(&mut self, name: &str);
fn parse_coordinate_system(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string(s, "CoordinateSystem");
}

//fn pbrt_coord_sys_transform(&mut self, name: &str);
fn parse_coord_sys_transform(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string(s, "CoordSysTransform");
}

fn parse_active_transform(s: &str) -> IResult<&str, OPNode> {
    let (s, (op, t)) = nom::branch::permutation((
        sequence::terminated(bytes::complete::tag("ActiveTransform"), space1),
        sequence::terminated(
            nom::branch::alt((
                bytes::complete::tag("All"),
                bytes::complete::tag("EndTime"),
                bytes::complete::tag("StartTime"),
            )),
            space0,
        ),
    ))(s)?;

    //fn pbrt_active_transform_all(&mut self);
    //fn pbrt_active_transform_end_time(&mut self);
    //fn pbrt_active_transform_start_time(&mut self);
    let name = String::from(op) + t;
    return Ok((s, OPNode::new(&name, None, None)));
}

//fn pbrt_transform_times(&mut self, start: Float, end: Float);
fn parse_transform_times(s: &str) -> IResult<&str, OPNode> {
    return parse_op_float_n(s, "TransformTimes", 2);
}

fn parse_pixel_filter(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "PixelFilter");
}

fn parse_film(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "Film");
}

fn parse_sampler(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "Sampler");
}

fn parse_accelerator(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "Accelerator");
}

fn parse_integrator(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "Integrator");
}

fn parse_camera(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "Camera");
}

fn parse_make_named_medium(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "MakeNamedMedium");
}

fn parse_medium_interface(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_string(s, "MediumInterface");
}

fn parse_world_begin(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "WorldBegin");
}

fn parse_attribute_begin(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "AttributeBegin");
}

fn parse_attribute_end(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "AttributeEnd");
}

fn parse_transform_begin(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "TransformBegin");
}

fn parse_transform_end(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "TransformEnd");
}

fn parse_texture(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_string_string_params(s, "Texture");
}

fn parse_material(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "Material");
}

fn parse_make_named_material(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "MakeNamedMaterial");
}

fn parse_named_material(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string(s, "NamedMaterial");
}

fn parse_light_source(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "LightSource");
}

fn parse_area_light_source(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "AreaLightSource");
}

fn parse_shape(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "Shape");
}

fn parse_reverse_orientation(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "ReverseOrientation");
}

fn parse_object_begin(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string(s, "ObjectBegin");
}

fn parse_object_end(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "ObjectEnd");
}

fn parse_object_instance(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string(s, "ObjectInstance");
}

fn parse_world_end(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "WorldEnd");
}

fn parse_work_dir_begin(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string(s, "WorkDirBegin");
}

fn parse_work_dir_end(s: &str) -> IResult<&str, OPNode> {
    return parse_op_void(s, "WorkDirEnd");
}

fn parse_include(s: &str) -> IResult<&str, OPNode> {
    return parse_op_string_params(s, "Include");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_string_literal() {
        let s = "\"abc\"";
        let (_, b) = string_literal(s).unwrap();
        assert_eq!(b, "abc");
    }

    #[test]
    fn test_parse_ops_001() {
        let s = "Integrator \"path\" \n
        WorldBegin";
        let r = parse_opnodes(s);
        match r {
            Ok(nodes) => {
                assert_eq!(nodes[0].name, "Integrator");
            }
            _ => {
                assert_eq!(1, 0);
            }
        }
    }
    #[test]
    fn test_parse_ops_002() {
        let s = "\n
            Translate 0 0 -140\n
            Rotate 0 1 2 3\n
            LookAt 0 1 2 3 4 5 6 7 8\n

            Transform [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]\n
            Texture \"a\" \"b\" \"c\"
        ";
        let r = parse_opnodes(s);
        match r {
            Ok(nodes) => {
                assert_eq!(nodes[0].name, "Translate");
                {
                    let op = &nodes[0];
                    let p = op.args.as_ref().unwrap();
                    let v = p.get_floats("args");
                    assert_eq!(v.len(), 3);
                }
                assert_eq!(nodes[1].name, "Rotate");
                {
                    let op = &nodes[1];
                    let p = op.args.as_ref().unwrap();
                    let v = p.get_floats("args");
                    assert_eq!(v.len(), 4);
                }
                assert_eq!(nodes[2].name, "LookAt");
                {
                    let op = &nodes[2];
                    let p = op.args.as_ref().unwrap();
                    let v = p.get_floats("args");
                    assert_eq!(v.len(), 9);
                }
                assert_eq!(nodes[3].name, "Transform");
                {
                    let op = &nodes[3];
                    let p = op.args.as_ref().unwrap();
                    let v = p.get_floats("arg1");
                    assert_eq!(v.len(), 16);
                }
                assert_eq!(nodes[4].name, "Texture");
                {
                    let op = &nodes[4];
                    let p = op.args.as_ref().unwrap();
                    let v1 = p.get_strings("arg1");
                    assert_eq!(v1.len(), 1);
                }
            }
            _ => {
                assert!(false);
            }
        }
    }

    #[test]
    fn test_parse_ops_003() {
        let s = "\n
        Film \"image\" \"integer xresolution\" [700] \"integer yresolution\" [700]\n
            \"string filename\" \"killeroo-simple.exr\"\n
        ";
        let r = parse_opnodes(s);
        match r {
            Ok(nodes) => {
                assert_eq!(nodes[0].name, "Film");
                let params = nodes[0].params.as_ref().unwrap();
                {
                    let p1 = params.get_ints("xresolution");
                    assert_eq!(p1[0], 700);
                }
                {
                    let p1 = params.get_ints("yresolution");
                    assert_eq!(p1[0], 700);
                }
                {
                    let p1 = params.get_strings("filename");
                    assert_eq!(p1[0], "killeroo-simple.exr");
                }
            }
            _ => {
                assert!(false);
            }
        }
    }

    #[test]
    fn test_parse_ops_004() {
        let s = "\n
        AttributeBegin # A\n
            Material \"matte\" \"color Kd\" [.5 .5 .8]\n
            Translate 0 0 -140\n
            Shape \"trianglemesh\" \"point P\" [ -1000 -1000 0 1000 -1000 0 1000 1000 0 -1000 1000 0 ]\n
                \"float uv\" [ 0 0 5 0 5 5 0 5 ]\n
                \"integer indices\" [ 0 1 2 2 3 0]\n
            Shape \"trianglemesh\" \"point P\" [ -400 -1000 -1000   -400 1000 -1000   -400 1000 1000 -400 -1000 1000 ]\n
                \"float uv\" [ 0 0 5 0 5 5 0 5 ]\n
                \"integer indices\" [ 0 1 2 2 3 0]\n
        AttributeEnd\n
        ";
        let r = parse_opnodes(s);
        match r {
            Ok(nodes) => {
                assert_eq!(nodes[0].name, "AttributeBegin");
                assert_eq!(nodes[1].name, "Material");
                assert_eq!(nodes[2].name, "Translate");
                assert_eq!(nodes[3].name, "Shape");
                assert_eq!(nodes[4].name, "Shape");
                assert_eq!(nodes[5].name, "AttributeEnd");
                {
                    let params = nodes[1].params.as_ref().unwrap();
                    {
                        let p1 = params.get_points("Kd");
                        assert_eq!(p1, vec![0.5, 0.5, 0.8]);
                    }
                }
                {
                    let params = nodes[3].params.as_ref().unwrap();
                    {
                        let p1 = params.get_points("P");
                        assert_eq!(
                            p1,
                            vec![
                                -1000.0, -1000.0, 0., 1000., -1000., 0., 1000., 1000., 0., -1000.,
                                1000., 0.
                            ]
                        );
                    }
                    {
                        let p1 = params.get_floats("uv");
                        assert_eq!(p1, vec![0., 0., 5., 0., 5., 5., 0., 5.]);
                    }
                    {
                        let p1 = params.get_ints("indices");
                        assert_eq!(p1, vec![0, 1, 2, 2, 3, 0]);
                    }
                }
            }
            _ => {
                assert!(false);
            }
        }
    }
}
