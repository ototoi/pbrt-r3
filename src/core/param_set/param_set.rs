use super::wellknown_params;
use crate::core::pbrt::types::*;
use crate::core::spectrum::*;
use std::cell::{Ref, RefCell, RefMut};
use std::collections::HashMap;

pub struct ParamSet {
    pub bools: HashMap<String, RefCell<Vec<bool>>>,
    pub ints: HashMap<String, RefCell<Vec<i32>>>,
    pub floats: HashMap<String, RefCell<Vec<Float>>>,
    pub strings: HashMap<String, RefCell<Vec<String>>>,
    pub spectrums: HashMap<String, RefCell<Vec<Spectrum>>>,
    pub points: HashMap<String, RefCell<Vec<Float>>>,
    pub keys: Vec<String>,
}

fn get_key_type(key: &str) -> String {
    let ss: Vec<&str> = key.split_ascii_whitespace().collect();
    match ss.len() {
        2 => {
            return String::from(ss[0]);
        }
        1 => {
            let name = ss[0];
            if let Some(t) = wellknown_params::find_type_from_key(name) {
                return String::from(t);
            }
        }
        _ => {}
    }
    return String::from("");
}

fn get_key_name(key: &str) -> String {
    let ss: Vec<&str> = key.split_ascii_whitespace().collect();
    match ss.len() {
        2 => String::from(ss[1]),
        _ => String::from(key),
    }
}

fn add_value<T: Clone>(
    k: &mut Vec<String>,
    m: &mut HashMap<String, RefCell<Vec<T>>>,
    key: &str,
    v: T,
) {
    k.push(key.to_string());
    let keyname = get_key_name(key);
    match m.get(&keyname) {
        Some(r) => {
            let mut items = r.borrow_mut();
            items.push(v);
        }
        _ => {
            let r = RefCell::<Vec<T>>::new(Vec::<T>::new());
            {
                let mut items = r.borrow_mut();
                items.push(v);
            }
            m.insert(keyname, r);
        }
    }
}

fn add_values<T: Clone>(
    k: &mut Vec<String>,
    m: &mut HashMap<String, RefCell<Vec<T>>>,
    key: &str,
    v: &[T],
) {
    k.push(key.to_string());
    let keyname = get_key_name(key);
    if let Some(r) = m.get(&keyname) {
        let mut items = r.borrow_mut();
        items.clear(); //
        for x in v.iter() {
            items.push(x.clone());
        }
    } else {
        let r = RefCell::<Vec<T>>::new(Vec::<T>::new());
        {
            let mut items = r.borrow_mut();
            for x in v.iter() {
                items.push(x.clone());
            }
        }
        m.insert(keyname, r);
    }
}

fn add_value_no_key<T: Clone>(m: &mut HashMap<String, RefCell<Vec<T>>>, key: &str, v: T) {
    let keyname = get_key_name(key);
    match m.get(&keyname) {
        Some(r) => {
            let mut items = r.borrow_mut();
            items.push(v);
        }
        _ => {
            let r = RefCell::<Vec<T>>::new(Vec::<T>::new());
            {
                let mut items = r.borrow_mut();
                items.push(v);
            }
            m.insert(keyname, r);
        }
    }
}

fn get_values<T: Clone>(m: &HashMap<String, RefCell<Vec<T>>>, key: &str) -> Vec<T> {
    let keyname = get_key_name(key);
    match m.get(&keyname) {
        Some(r) => r.borrow().clone(),
        _ => Vec::<T>::new(),
    }
}

fn get_values_ref<'a, T: Clone>(
    m: &'a HashMap<String, RefCell<Vec<T>>>,
    key: &str,
) -> Option<Ref<'a, Vec<T>>> {
    let keyname = get_key_name(key);
    let r = m.get(&keyname);
    match r {
        Some(r) => {
            return Some(r.borrow());
        }
        _ => {
            return None;
        }
    }
}

fn get_values_mut<'a, T: Clone>(
    m: &'a HashMap<String, RefCell<Vec<T>>>,
    key: &str,
) -> Option<RefMut<'a, Vec<T>>> {
    let keyname = get_key_name(key);
    let r = m.get(&keyname);
    match r {
        Some(r) => {
            return Some(r.borrow_mut());
        }
        _ => {
            return None;
        }
    }
}

impl ParamSet {
    pub fn new() -> Self {
        ParamSet {
            bools: HashMap::new(),
            ints: HashMap::new(),
            floats: HashMap::new(),
            strings: HashMap::new(),
            spectrums: HashMap::new(),
            points: HashMap::new(),
            keys: Vec::<String>::new(),
        }
    }

    //--------------------

    pub fn add_bool(&mut self, key: &str, v: bool) {
        add_value(&mut self.keys, &mut self.bools, key, v);
    }

    pub fn add_bools(&mut self, key: &str, v: &[bool]) {
        add_values(&mut self.keys, &mut self.bools, key, v);
    }

    pub fn add_int(&mut self, key: &str, v: i32) {
        add_value(&mut self.keys, &mut self.ints, key, v);
    }

    pub fn add_ints(&mut self, key: &str, v: &[i32]) {
        add_values(&mut self.keys, &mut self.ints, key, v);
    }

    pub fn add_float(&mut self, key: &str, v: f32) {
        add_value(&mut self.keys, &mut self.floats, key, v);
    }

    pub fn add_floats(&mut self, key: &str, v: &[Float]) {
        let t = get_key_type(key);
        match &t as &str {
            "point" => add_values(&mut self.keys, &mut self.points, key, v),
            "normal" => add_values(&mut self.keys, &mut self.points, key, v),
            "vector" => add_values(&mut self.keys, &mut self.points, key, v),
            "color" => self.add_color(key, v),
            "rgb" => self.add_color(key, v),
            "blackbody" => add_values(&mut self.keys, &mut self.floats, key, v),
            _ => add_values(&mut self.keys, &mut self.floats, key, v),
        }
    }

    pub fn add_string(&mut self, key: &str, v: &str) {
        add_value(&mut self.keys, &mut self.strings, key, String::from(v));
    }

    pub fn add_strings(&mut self, key: &str, v: &[&str]) {
        let vv: Vec<String> = v.iter().map(|s| String::from(*s)).collect();
        add_values(&mut self.keys, &mut self.strings, key, &vv);
    }

    pub fn add_spectrum(&mut self, key: &str, v: &Spectrum) {
        add_value(&mut self.keys, &mut self.spectrums, key, *v);
    }

    pub fn add_spectrums(&mut self, key: &str, v: &[Spectrum]) {
        add_values(&mut self.keys, &mut self.spectrums, key, v);
    }

    pub fn add_point(&mut self, key: &str, v: &[Float]) {
        add_values(&mut self.keys, &mut self.points, key, v);
    }

    //--------------------

    pub fn add_color(&mut self, key: &str, v: &[Float]) {
        add_values(&mut self.keys, &mut self.points, key, v);
    }

    pub fn add_rgb(&mut self, key: &str, v: &[Float]) {
        add_values(&mut self.keys, &mut self.points, key, v);
    }

    pub fn add_xyz(&mut self, key: &str, v: &[Float]) {
        let xyz = RGBSpectrum::rgb_from_xyz(v);
        let rgb = xyz.to_rgb();
        add_values(&mut self.keys, &mut self.points, key, &rgb);
    }

    pub fn add_blackbody(&mut self, key: &str, v: &[Float]) {
        add_values(&mut self.keys, &mut self.floats, key, v);
    }

    pub fn add_spectrum_no_key(&mut self, key: &str, v: &Spectrum) {
        add_value_no_key(&mut self.spectrums, key, *v);
    }

    //--------------------

    pub fn add_point2f(&mut self, key: &str, v: &Point2f) {
        let tv = [v.x, v.y];
        self.add_floats(key, &tv);
    }

    pub fn add_vector2f(&mut self, key: &str, v: &Vector2f) {
        let tv = [v.x, v.y];
        self.add_floats(key, &tv);
    }

    //--------------------

    pub fn add_point3f(&mut self, key: &str, v: &Point3f) {
        let tv = [v.x, v.y, v.z];
        self.add_point(key, &tv);
    }

    pub fn add_vector3f(&mut self, key: &str, v: &Vector3f) {
        let tv = [v.x, v.y, v.z];
        self.add_point(key, &tv);
    }

    pub fn add_normal3f(&mut self, key: &str, v: &Normal3f) {
        let tv = [v.x, v.y, v.z];
        self.add_point(key, &tv);
    }

    //--------------------

    pub fn get_bools(&self, key: &str) -> Vec<bool> {
        return get_values(&self.bools, key);
    }

    pub fn get_ints(&self, key: &str) -> Vec<i32> {
        return get_values(&self.ints, key);
    }

    pub fn get_floats(&self, key: &str) -> Vec<Float> {
        return get_values(&self.floats, key);
    }

    pub fn get_strings(&self, key: &str) -> Vec<String> {
        return get_values(&self.strings, key);
    }

    pub fn get_spectrums(&self, key: &str) -> Vec<Spectrum> {
        return get_values(&self.spectrums, key);
    }

    pub fn get_points(&self, key: &str) -> Vec<Float> {
        return get_values(&self.points, key);
    }

    //--------------------
    pub fn get_bools_ref(&self, key: &str) -> Option<Ref<Vec<bool>>> {
        return get_values_ref(&self.bools, key);
    }

    pub fn get_ints_ref(&self, key: &str) -> Option<Ref<Vec<i32>>> {
        return get_values_ref(&self.ints, key);
    }

    pub fn get_floats_ref(&self, key: &str) -> Option<Ref<Vec<Float>>> {
        return get_values_ref(&self.floats, key);
    }

    pub fn get_strings_ref(&self, key: &str) -> Option<Ref<Vec<String>>> {
        return get_values_ref(&self.strings, key);
    }

    pub fn get_textures_ref(&self, key: &str) -> Option<Ref<Vec<String>>> {
        return get_values_ref(&self.strings, key);
    }

    pub fn get_spectrums_ref(&self, key: &str) -> Option<Ref<Vec<Spectrum>>> {
        return get_values_ref(&self.spectrums, key);
    }

    pub fn get_points_ref(&self, key: &str) -> Option<Ref<Vec<Float>>> {
        return get_values_ref(&self.points, key);
    }

    //--------------------

    pub fn get_bools_mut(&mut self, key: &str) -> Option<RefMut<Vec<bool>>> {
        return get_values_mut(&self.bools, key);
    }

    pub fn get_ints_mut(&mut self, key: &str) -> Option<RefMut<Vec<i32>>> {
        return get_values_mut(&self.ints, key);
    }

    pub fn get_strings_mut(&mut self, key: &str) -> Option<RefMut<Vec<String>>> {
        return get_values_mut(&self.strings, key);
    }

    //--------------------
    pub fn find_one_bool(&self, key: &str, value: bool) -> bool {
        let r = self.get_bools_ref(key);
        match r {
            Some(v) => v[0],
            None => value,
        }
    }

    pub fn find_one_int(&self, key: &str, value: i32) -> i32 {
        let r = self.get_ints_ref(key);
        match r {
            Some(v) => v[0],
            None => value,
        }
    }

    pub fn find_one_float(&self, key: &str, value: Float) -> Float {
        let r = self.get_floats_ref(key);
        match r {
            Some(v) => v[0],
            None => value,
        }
    }

    pub fn find_one_string(&self, key: &str, value: &str) -> String {
        let r = self.get_strings_ref(key);
        match r {
            Some(v) => v[0].clone(),
            None => String::from(value),
        }
    }

    pub fn find_one_filename(&self, key: &str, value: &str) -> String {
        let r = self.get_strings_ref(key);
        match r {
            Some(v) => v[0].clone(),
            None => String::from(value),
        }
    }

    pub fn find_one_spectrum(&self, key: &str, value: &Spectrum) -> Spectrum {
        if let Some(r) = self.get_points_ref(key) {
            return Spectrum::from_rgb(&r, SpectrumType::Reflectance); //todo
        } else if let Some(r) = self.get_spectrums_ref(key) {
            return r[0];
        } else {
            return *value;
        }
    }

    pub fn find_one_point(&self, key: &str, value: &[Float]) -> Vec<Float> {
        let r = self.get_points_ref(key);
        match r {
            Some(v) => {
                return v.clone();
            }
            None => value.to_vec(),
        }
    }

    pub fn find_one_point3f(&self, key: &str, value: &Point3f) -> Point3f {
        let v = vec![value.x, value.y, value.z];
        let a: &[f32] = &self.find_one_point(key, &v);
        return Point3f::from(a);
    }

    pub fn find_one_vector3f(&self, key: &str, value: &Vector3f) -> Vector3f {
        let v = vec![value.x, value.y, value.z];
        let a: &[f32] = &self.find_one_point(key, &v);
        return Vector3f::from(a);
    }

    pub fn find_one_normal3f(&self, key: &str, value: &Normal3f) -> Normal3f {
        let v = vec![value.x, value.y, value.z];
        let a: &[f32] = &self.find_one_point(key, &v);
        return Normal3f::from(a);
    }

    //--------------------

    pub fn replace_one_bool(&mut self, key: &str, value: bool) {
        let mut replaced = false;
        if let Some(mut r) = self.get_bools_mut(key) {
            if !r.is_empty() {
                r[0] = value;
                replaced = true;
            }
        }
        if !replaced {
            self.add_bool(key, value);
        }
    }

    pub fn replace_one_int(&mut self, key: &str, value: i32) {
        let mut replaced = false;
        if let Some(mut r) = self.get_ints_mut(key) {
            if !r.is_empty() {
                r[0] = value;
                replaced = true;
            }
        }
        if !replaced {
            self.add_int(key, value);
        }
    }

    pub fn replace_one_string(&mut self, key: &str, value: &str) {
        let mut replaced = false;
        if let Some(mut r) = self.get_strings_mut(key) {
            if !r.is_empty() {
                r[0] = String::from(value);
                replaced = true;
            }
        }
        if !replaced {
            self.add_string(key, value);
        }
    }

    //--------------------
    pub fn set(&mut self, other: &ParamSet) {
        self.bools = other.bools.clone();
        self.ints = other.ints.clone();
        self.floats = other.floats.clone();
        self.strings = other.strings.clone();
        self.spectrums = other.spectrums.clone();
        self.points = other.points.clone();
        self.keys = other.keys.clone();
    }
    //--------------------

    pub fn get_keys(&self) -> Vec<String> {
        return self.keys.clone();
    }

    pub fn get_key_name(&self, key: &str) -> String {
        return get_key_name(key);
    }

    pub fn get_key_type(&self, key: &str) -> String {
        return get_key_type(key);
    }
}

impl Clone for ParamSet {
    fn clone(&self) -> Self {
        ParamSet {
            bools: self.bools.clone(),
            ints: self.ints.clone(),
            floats: self.floats.clone(),
            strings: self.strings.clone(),
            spectrums: self.spectrums.clone(),
            points: self.points.clone(),
            keys: self.keys.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let mut params = ParamSet::new();
        params.add_bool("bool is_a", true);
        let val1 = params.find_one_bool("is_a", false);
        let val2 = params.find_one_bool("bool is_a", false);
        let val3 = params.find_one_bool("fuga", false);
        let val4 = params.find_one_bool("bool fuga", false);
        assert_eq!(val1, true);
        assert_eq!(val2, true);
        assert_eq!(val3, false);
        assert_eq!(val4, false);
    }

    #[test]
    fn test_002() {
        let mut params = ParamSet::new();
        params.add_int("integer count", 1234);
        let val1 = params.find_one_int("count", 5678);
        let val2 = params.find_one_int("integer count", 5678);
        let val3 = params.find_one_int("fuga", 5678);
        let val4 = params.find_one_int("integer fuga", 5678);
        assert_eq!(val1, 1234);
        assert_eq!(val2, 1234);
        assert_eq!(val3, 5678);
        assert_eq!(val4, 5678);
    }

    #[test]
    fn test_003() {
        let mut params = ParamSet::new();
        params.add_float("float value1", 1234.0);
        let val1 = params.find_one_float("value1", 5678.0);
        let val2 = params.find_one_float("float value1", 5678.0);
        let val3 = params.find_one_float("fuga", 5678.0);
        let val4 = params.find_one_float("float fuga", 5678.0);
        assert_eq!(val1, 1234.0);
        assert_eq!(val2, 1234.0);
        assert_eq!(val3, 5678.0);
        assert_eq!(val4, 5678.0);
    }

    #[test]
    fn test_004() {
        let mut params = ParamSet::new();
        let s1 = Point3f::from([1.0, 2.0, 3.0]);
        let s2 = Point3f::from([4.0, 5.0, 6.0]);
        params.add_point3f("point P", &s1);
        let val1 = params.find_one_point3f("P", &s2);
        let val2 = params.find_one_point3f("point P", &s2);
        let val3 = params.find_one_point3f("fuga", &s2);
        let val4 = params.find_one_point3f("point fuga", &s2);
        assert_eq!(val1, s1);
        assert_eq!(val2, s1);
        assert_eq!(val3, s2);
        assert_eq!(val4, s2);
    }

    #[test]
    fn test_005() {
        let mut params = ParamSet::new();
        params.add_string("string value1", "hello!");
        let val1 = params.find_one_string("value1", "world!");
        let val2 = params.find_one_string("string value1", "world!");
        let val3 = params.find_one_string("fuga", "world!");
        let val4 = params.find_one_string("string fuga", "world!");
        assert_eq!(val1, "hello!");
        assert_eq!(val2, "hello!");
        assert_eq!(val3, "world!");
        assert_eq!(val4, "world!");
    }

    #[test]
    fn test_006() {
        let mut params = ParamSet::new();
        let s1 = Spectrum::from([1.0, 2.0, 3.0]);
        let s2 = Spectrum::from([4.0, 5.0, 6.0]);
        params.add_spectrum("spectrum value1", &s1);
        let val1 = params.find_one_spectrum("value1", &s2);
        let val2 = params.find_one_spectrum("spectrum value1", &s2);
        let val3 = params.find_one_spectrum("fuga", &s2);
        let val4 = params.find_one_spectrum("spectrum fuga", &s2);
        assert_eq!(val1, s1);
        assert_eq!(val2, s1);
        assert_eq!(val3, s2);
        assert_eq!(val4, s2);
    }

    #[test]
    fn test_007() {
        let mut params1 = ParamSet::new();
        let s1 = 1234;
        let s2 = 4567;
        let s3 = 78910;
        params1.add_int("integer count", s1);
        let mut params2 = params1.clone();
        let val1 = params2.find_one_int("count", s2);
        params2.add_int("integer total", s3);
        let val2 = params1.find_one_int("total", s2);
        let val3 = params2.find_one_int("total", s2);
        params1.set(&params2);
        let val4 = params1.find_one_int("total", s2);

        assert_eq!(val1, s1);
        assert_eq!(val2, s2);
        assert_eq!(val3, s3);
        assert_eq!(val4, s3);
    }
}
