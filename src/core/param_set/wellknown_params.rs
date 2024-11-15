const WELLKNOWN_PARAMS: [(&str, &str, &str); 8] = [
    ("", "integer", "xresolution"),
    ("", "integer", "yresolution"),
    ("", "integer", "maxdepth"),
    ("", "integer", "pixelsamples"),
    ("", "string", "filename"),
    ("", "float", "fov"),
    ("", "float", "radius"),
    ("", "color", "L"),
];

pub fn find_type_from_key(key: &str) -> Option<&str> {
    if let Some((_, t, _k)) = WELLKNOWN_PARAMS.iter().find(|(_a, _tt, kk)| -> bool {
        return *kk == key;
    }) {
        return Some(*t);
    } else {
        return None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let t = find_type_from_key("xresolution");
        assert!(t.is_some());
        assert_eq!("integer", t.unwrap());
    }
}
