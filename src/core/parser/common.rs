use crate::core::base::Float;
use crate::core::param_set::wellknown_params;
use crate::core::param_set::ParamSet;
use nom::bytes;
use nom::character;
use nom::multi;
use nom::number;
use nom::sequence;
use nom::IResult;

pub fn space0(s: &str) -> IResult<&str, &str> {
    return character::complete::multispace0(s);
}

pub fn space1(s: &str) -> IResult<&str, &str> {
    return character::complete::multispace1(s);
}

pub fn bool_literal(s: &str) -> IResult<&str, &str> {
    return nom::branch::alt((
        nom::bytes::complete::tag_no_case("true"),
        nom::bytes::complete::tag_no_case("false"),
    ))(s);
}

pub fn float_literal(s: &str) -> IResult<&str, &str> {
    return number::complete::recognize_float(s);
}

pub fn string_literal(s: &str) -> IResult<&str, &str> {
    return sequence::delimited(
        character::complete::char('"'),
        bytes::complete::take_until("\""),
        character::complete::char('"'),
    )(s);
}

pub fn parse_literal(s: &str) -> IResult<&str, &str> {
    return nom::branch::alt((float_literal, bool_literal, string_literal))(s);
}

pub fn parse_listed_literal(s: &str) -> IResult<&str, Vec<&str>> {
    let (s, r) = parse_literal(s)?;
    return Ok((s, vec![r]));
}

pub fn parse_list(s: &str) -> IResult<&str, Vec<&str>> {
    return sequence::delimited(
        character::complete::char('['),
        sequence::delimited(
            space0,
            multi::separated_list1(space1, parse_literal), //nom::branch::alt((number::complete::recognize_float, token))),
            space0,
        ),
        character::complete::char(']'),
    )(s);
}

pub fn get_param_type(s: &str) -> (&str, &str) {
    let ss: Vec<&str> = s.split_ascii_whitespace().collect();
    if ss.len() == 2 {
        return (ss[0], ss[1]);
    } else if ss.len() == 1 {
        if let Some(t) = wellknown_params::find_type_from_key(ss[0]) {
            return (t, ss[0]);
        } else {
            return ("", ss[0]);
        }
    } else {
        return ("", s);
    }
}

pub fn convert_bool(s: &str) -> Result<bool, std::str::ParseBoolError> {
    let s = String::from(s).to_lowercase();
    let s2: &str = &s;
    match s2 {
        "true" => return Ok(true),
        "false" => return Ok(false),
        "\"true\"" => return Ok(true),
        "\"false\"" => return Ok(false),
        _ => return s2.parse::<bool>(),
    }
}

pub fn parse_params(s: &str) -> IResult<&str, ParamSet> {
    let (s, v) = multi::separated_list0(
        space1,
        nom::branch::permutation((
            sequence::terminated(string_literal, space1),
            nom::branch::alt((parse_list, parse_listed_literal)),
        )),
    )(s)?;
    let mut params = ParamSet::new();
    for vv in v {
        let org_key = vv.0;
        let (t, key) = get_param_type(org_key);
        let new_key = format!("{t} {key}");
        match t {
            "string" => {
                let s_values = vv.1;
                params.add_strings(&new_key, &s_values);
            }
            "texture" => {
                let s_values = vv.1;
                params.add_strings(&new_key, &s_values);
            }
            "spectrum" => {
                let s_values = vv.1;
                params.add_strings(&new_key, &s_values);
            }
            "bool" => {
                let s_values = vv.1;
                let values: Vec<bool> = s_values.iter().map(|s| convert_bool(s).unwrap()).collect();
                params.add_bools(&new_key, &values);
            }
            "integer" => {
                let s_values = vv.1;
                let values: Vec<i32> = s_values.iter().map(|s| s.parse::<i32>().unwrap()).collect();
                params.add_ints(&new_key, &values);
            }
            "color" => {
                let s_values = vv.1;
                let values: Vec<Float> = s_values
                    .iter()
                    .map(|s| s.parse::<Float>().unwrap())
                    .collect();
                params.add_color(&new_key, &values);
            }
            "rgb" => {
                let s_values = vv.1;
                let values: Vec<Float> = s_values
                    .iter()
                    .map(|s| s.parse::<Float>().unwrap())
                    .collect();
                params.add_rgb(&new_key, &values);
            }
            "xyz" => {
                let s_values = vv.1;
                let values: Vec<Float> = s_values
                    .iter()
                    .map(|s| s.parse::<Float>().unwrap())
                    .collect();
                params.add_xyz(&new_key, &values);
            }
            "blackbody" => {
                let s_values = vv.1;
                let values: Vec<Float> = s_values
                    .iter()
                    .map(|s| s.parse::<Float>().unwrap())
                    .collect();
                params.add_blackbody(&new_key, &values);
            }
            "point" => {
                let s_values = vv.1;
                let values: Vec<Float> = s_values
                    .iter()
                    .map(|s| s.parse::<Float>().unwrap())
                    .collect();
                params.add_point(&new_key, &values);
            }
            "vector" => {
                let s_values = vv.1;
                let values: Vec<Float> = s_values
                    .iter()
                    .map(|s| s.parse::<Float>().unwrap())
                    .collect();
                params.add_point(&new_key, &values);
            }
            "normal" => {
                let s_values = vv.1;
                let values: Vec<Float> = s_values
                    .iter()
                    .map(|s| s.parse::<Float>().unwrap())
                    .collect();
                params.add_point(&new_key, &values);
            }
            "float" => {
                let s_values = vv.1;
                let values: Vec<Float> = s_values
                    .iter()
                    .map(|s| s.parse::<Float>().unwrap())
                    .collect();
                params.add_floats(&new_key, &values);
            }
            _ => {
                let s_values = vv.1;
                //println!("{}:{:?}", t, s_values);
                let values: Vec<Float> = s_values
                    .iter()
                    .map(|s| s.parse::<Float>().unwrap())
                    .collect();
                params.add_floats(&new_key, &values);
            }
        }
    }
    return Ok((s, params));
}
