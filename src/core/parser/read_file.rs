use super::common::*;
use super::remove_comment::remove_comment;
use nom::bytes;
use nom::character;
use nom::multi;
use nom::sequence;
use nom::IResult;
use std::fs;
use std::io::Error;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};

pub fn read_file_with_include(path: &str) -> Result<String, Error> {
    let path = Path::new(path);
    if path.exists() {
        let path = path.canonicalize().unwrap();
        let mut dirs = Vec::<PathBuf>::new();
        dirs.push(PathBuf::from(path.parent().unwrap()));
        let mut ss = String::new();
        ss += &print_work_dir_begin(&dirs);
        ss += &read_file_with_include_core(path.as_path(), &mut dirs)?;
        ss += &print_work_dir_end();
        return Ok(ss);
    } else {
        return Err(Error::new(ErrorKind::NotFound, "File is not found."));
    }
}

fn read_file_with_include_core(path: &Path, dirs: &mut Vec<PathBuf>) -> Result<String, Error> {
    let string_result = fs::read_to_string(path);
    match string_result {
        Ok(s) => {
            return evaluate_include(&s, dirs);
        }
        Err(e) => {
            return Err(e);
        }
    }
}

fn get_next_path(filename: &Path, dirs: &[PathBuf]) -> Option<PathBuf> {
    for d in dirs.iter().rev() {
        let dir = d;
        let path = dir.join(Path::new(filename));
        if path.exists() {
            return Some(path);
        }
    }
    return None;
}

fn print_work_dir_begin(dirs: &[PathBuf]) -> String {
    let path = &dirs[dirs.len() - 1];
    let path = path.as_path().to_str().unwrap();
    return format!("WorkDirBegin \"{}\"\n", path);
}

fn print_work_dir_end() -> String {
    return "WorkDirEnd\n".to_string();
}

fn evaluate_include(s: &str, dirs: &mut Vec<PathBuf>) -> Result<String, Error> {
    let r = remove_comment(s);
    match r {
        Ok((_, s)) => {
            let s = nom::combinator::all_consuming(nom::multi::many0(parse_one))(&s);
            match s {
                Ok((_, vs)) => {
                    let mut ss = String::new();
                    //ss += &print_work_dir_begin(dirs);
                    for s in vs {
                        if s.starts_with("Include") {
                            let vv: Vec<&str> = s.split('|').collect();
                            let filename = Path::new(vv[1]);
                            let path_opt = get_next_path(filename, dirs);
                            match path_opt {
                                Some(next_path) => {
                                    dirs.push(PathBuf::from(next_path.parent().unwrap()));
                                    ss += &print_work_dir_begin(dirs);

                                    let rss =
                                        read_file_with_include_core(next_path.as_path(), dirs)?;
                                    ss += &rss;
                                    if vv.len() > 2 {
                                        for i in 2..vv.len() {
                                            ss += &format!(" {}", vv[i]);
                                        }
                                        ss += "\n";
                                    }
                                    dirs.pop();
                                    ss += &print_work_dir_end();
                                }
                                None => {
                                    return Err(Error::new(
                                        ErrorKind::NotFound,
                                        "File is not found.",
                                    ));
                                }
                            }
                        } else {
                            ss += &s;
                        }
                    }
                    //ss += &print_work_dir_end();
                    //print!("ss:{}", ss);
                    return Ok(ss);
                }
                Err(e) => {
                    return Err(Error::new(ErrorKind::Other, e.to_string()));
                }
            }
        }
        Err(e) => {
            return Err(Error::new(ErrorKind::Other, e.to_string()));
        }
    }
}

fn parse_one(s: &str) -> IResult<&str, String> {
    return nom::branch::alt((
        parse_space1,
        parse_string_literal,
        parse_include,
        parse_token,
        parse_float,
        parse_any,
    ))(s);
}

fn parse_token(s: &str) -> IResult<&str, String> {
    let r = nom::branch::permutation((
        character::complete::alpha1,
        bytes::complete::take_while(|c: char| c.is_alphanumeric() || c == '_'),
    ))(s);
    match r {
        Ok((s, (a, b))) => {
            let ss = format!("{}{}", a, b);
            return Ok((s, ss));
        }
        Err(e) => {
            return Err(e);
        }
    }
}

fn parse_float(s: &str) -> IResult<&str, String> {
    let r = nom::number::complete::recognize_float(s);
    match r {
        Ok((s, a)) => {
            let ss = String::from(a);
            return Ok((s, ss));
        }
        Err(e) => {
            return Err(e);
        }
    }
}

fn parse_any(s: &str) -> IResult<&str, String> {
    let r = character::complete::anychar(s);
    match r {
        Ok((s, b)) => {
            return Ok((s, String::from(b)));
        }
        Err(e) => {
            return Err(e);
        }
    }
}

fn parse_space1(s: &str) -> IResult<&str, String> {
    let r = character::complete::multispace1(s);
    match r {
        Ok((s, b)) => {
            return Ok((s, String::from(b)));
        }
        Err(e) => {
            return Err(e);
        }
    }
}

fn parse_string_literal(s: &str) -> IResult<&str, String> {
    let r = string_literal(s);
    match r {
        Ok((s, ss)) => {
            let k = format!("{}{}{}", "\"", ss, "\"");
            return Ok((s, k));
        }
        Err(e) => {
            return Err(e);
        }
    }
}

pub fn parse_params(s: &str) -> IResult<&str, String> {
    let (s, v) = multi::separated_list0(
        space1,
        nom::branch::permutation((
            sequence::terminated(string_literal, space1),
            nom::branch::alt((parse_list, parse_listed_literal)),
        )),
    )(s)?;
    let mut ss = String::new();
    if v.len() > 0 {
        for (key, value) in v {
            let sv = value.join(" ");
            ss += &format!("|\"{}\" [{}]", key, sv);
        }
    }
    return Ok((s, ss));
}

fn parse_include(s: &str) -> IResult<&str, String> {
    let (s, (op, a, params)) = nom::branch::permutation((
        sequence::terminated(bytes::complete::tag("Include"), space1),
        sequence::terminated(string_literal, space0),
        parse_params,
    ))(s)?;
    let ss = format!("{}|{}{}", op, a, params);
    //print!("parse_include:{}", ss);
    return Ok((s, ss));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(dead_code)]
    fn evaluate_include2(s: &str) -> IResult<&str, String> {
        let r = nom::combinator::all_consuming(nom::multi::many0(parse_one))(s);
        match r {
            Ok((s, vs)) => {
                let mut ss = String::new();
                for rs in vs {
                    if rs.starts_with("Include") {
                        ss += &rs;
                    }
                }
                return Ok((s, ss));
            }
            Err(e) => {
                return Err(e);
            }
        }
    }

    #[allow(dead_code)]
    fn make_next_path(path: &Path, s: &str) -> Option<String> {
        let dir = path.parent().unwrap();
        let vv: Vec<&str> = s.split('|').collect();
        let filename = vv[1];
        let path2 = dir.join(Path::new(filename));
        //.to_str().unwrap());
        //let path2 = String::from(path2.to_str().unwrap());
        //let path2 = path::absolute::(path2.);
        //let path2 = String::from(path2.unwrap().to_str().unwrap());
        return Some(String::from(path2.as_path().to_str().unwrap()));
    }

    #[test]
    fn test_001() {
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
    Include \"./ss.txt\"\n
    ";
        let (a, b) = evaluate_include2(s).unwrap();
        assert_eq!(a, "");
        assert_eq!(b, "Include|./ss.txt");
        let next_ = make_next_path(Path::new("~/a.txt"), &b).unwrap();
        assert_eq!(next_, "~/./ss.txt");
    }
}
