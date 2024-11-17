use nom::bytes;
use nom::character;

use nom::sequence;
use nom::IResult;

pub fn remove_comment(s: &str) -> IResult<&str, String> {
    let (rs, vs) = nom::combinator::all_consuming(nom::multi::many0(parse_one))(s)?;
    return Ok((rs, vs.join("")));
}

fn parse_one(s: &str) -> IResult<&str, String> {
    return nom::branch::alt((
        parse_space1,
        parse_string_literal,
        parse_removed_comment,
        parse_token,
        parse_float,
        parse_any,
    ))(s);
}

fn parse_token(s: &str) -> IResult<&str, String> {
    let (s, (a, b)) = nom::branch::permutation((
        character::complete::alpha1,
        bytes::complete::take_while(|c: char| c.is_alphanumeric() || c == '_'),
    ))(s)?;
    return Ok((s, format!("{}{}", a, b)));
}

fn parse_float(s: &str) -> IResult<&str, String> {
    let (s, a) = nom::number::complete::recognize_float(s)?;
    return Ok((s, a.to_string()));
}

fn parse_any(s: &str) -> IResult<&str, String> {
    let (s, a) = character::complete::anychar(s)?;
    return Ok((s, a.to_string()));
}

fn parse_space1(s: &str) -> IResult<&str, String> {
    let (s, a) = character::complete::multispace1(s)?;
    return Ok((s, a.to_string()));
}

fn parse_comment(s: &str) -> IResult<&str, &str> {
    return sequence::preceded(
        character::complete::char('#'),
        bytes::complete::take_till(|c| c == '\n'),
    )(s);
}

fn parse_removed_comment(s: &str) -> IResult<&str, String> {
    let (s, _) = parse_comment(s)?;
    return Ok((s, String::from("")));
}

fn parse_string_literal(s: &str) -> IResult<&str, String> {
    let (s, a) = sequence::delimited(
        character::complete::char('"'),
        bytes::complete::take_until("\""),
        character::complete::char('"'),
    )(s)?;
    return Ok((s, format!("{}{}{}", "\"", a, "\"")));
}

//----------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_001() {
        let s = "Integrator";
        let (a, b) = parse_token(s).unwrap();
        assert_eq!(a, "");
        assert_eq!(b, "Integrator");
    }

    #[test]
    fn test_002() {
        let s = "Integrator ";
        let (a, b) = parse_token(s).unwrap();
        assert_eq!(a, " ");
        assert_eq!(b, "Integrator");
    }

    #[test]
    fn test_003() {
        let s = "#hoge ";
        let (a, b) = parse_removed_comment(s).unwrap();
        assert_eq!(a, "");
        assert_eq!(b, "");
    }

    #[test]
    fn test_004() {
        let s = "aaa llll";
        let (a, b) = parse_token(s).unwrap();
        assert_eq!(a, " llll");
        assert_eq!(b, "aaa");
    }

    #[test]
    fn test_005() {
        let s = "\"aaa\" #1234\n aaa";
        let (a, b) = (nom::multi::many0(parse_one)(s)).unwrap();
        assert_eq!(a, "");
        assert_eq!(b.len(), 5);
        assert_eq!(b[0], "\"aaa\"");
        assert_eq!(b[1], " ");
        assert_eq!(b[2], "");
    }

    #[test]
    fn test_006() {
        let s = "\"aaa\" #1234\n aaa";
        let (_, ss) = remove_comment(s).unwrap();
        assert_eq!(ss, "\"aaa\" \n aaa");
    }
}
