use super::parse_context::*;
use crate::core::base::*;
use crate::core::param_set::*;

use std::cell::{Cell, RefCell};
use std::io::Write;
use std::sync::Arc;

fn get_param_type(s: &str) -> (&str, &str) {
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

fn get_type(s: &str) -> &str {
    let (t, _) = get_param_type(s);
    match t {
        "string" => "s",
        "spectrum" => "s",
        "texture" => "s",
        "bool" => "b",
        "integer" => "i",
        "point" | "point2" | "point3" | "point4" => "p",
        "normal" => "p",
        "vector" | "vector2" | "vector3" | "vector4" => "p",
        "color" => "p",
        "rgb" => "p",
        "blackbody" => "f",
        _ => "f",
    }
}

pub struct PrintContext {
    pub writer: Arc<RefCell<dyn Write>>,
    pub omit_long_values: bool,
    pub indent: Cell<i32>,
}

impl PrintContext {
    pub fn new(w: Arc<RefCell<dyn Write>>) -> Self {
        PrintContext {
            writer: w,
            omit_long_values: false,
            indent: Cell::<i32>::new(0),
        }
    }

    pub fn new_with_params(w: Arc<RefCell<dyn Write>>, omit_long_values: bool) -> Self {
        PrintContext {
            writer: w,
            omit_long_values,
            indent: Cell::<i32>::new(0),
        }
    }

    pub fn new_stdout(omit_long_values: bool) -> Self {
        PrintContext {
            writer: Arc::new(RefCell::new(std::io::stdout())),
            omit_long_values,
            indent: Cell::<i32>::new(0),
        }
    }

    pub fn inc_indent(&mut self) {
        self.indent.set(self.indent.get() + 1);
    }
    pub fn dec_indent(&mut self) {
        self.indent.set(self.indent.get() - 1);
    }

    pub fn get_indent(&self) -> String {
        return self.get_indent_i(self.indent.get());
    }

    pub fn get_indent_i(&self, count: i32) -> String {
        let mut s = String::new();
        for _ in 0..count {
            s += "    ";
        }
        return s;
    }

    fn print(&self, s: &str) {
        _ = self.writer.borrow_mut().write_all(s.as_bytes());
    }

    fn convert_transform(&self, values: &[Float]) -> String {
        let mut s = String::from("[");
        let len = values.len();
        for i in 0..len {
            let v = values[i];
            s += &format!("{v}");
            if i != len - 1 {
                s += " ";
            }
        }
        s += "]";
        return s;
    }

    fn convert_values(&self, params: &ParamSet, key: &str) -> String {
        let mut s = String::from("");
        let t = get_type(key);
        s += "[";
        match t {
            "s" => {
                let values = params.get_strings(key);
                let len = values.len();
                for i in 0..len {
                    let v = &values[i];
                    s += &format!("\"{v}\"");
                    if i != len - 1 {
                        s += " ";
                    }
                }
            }
            "b" => {
                let values = params.get_bools(key);
                let len = values.len();
                for i in 0..len {
                    let v = values[i];
                    if v {
                        s += "\"true\"";
                    } else {
                        s += "\"false\"";
                    }
                    if i != len - 1 {
                        s += " ";
                    }
                }
            }
            "i" => {
                let values = params.get_ints_ref(key).unwrap();
                if !self.omit_long_values || values.len() <= 16 {
                    let len = values.len();
                    for i in 0..len {
                        let v = values[i];
                        s += &format!("{v}");
                        if i != len - 1 {
                            s += " ";
                        }
                    }
                } else {
                    s += "...";
                }
            }
            "p" => {
                let values = params.get_points_ref(key).unwrap();
                if !self.omit_long_values || values.len() <= 16 {
                    let len = values.len();
                    for i in 0..len {
                        let v = values[i];
                        s += &format!("{v}");
                        if i != len - 1 {
                            s += " ";
                        }
                    }
                } else {
                    s += "...";
                }
            }
            _ => {
                let values = params.get_floats_ref(key).unwrap();
                if !self.omit_long_values || values.len() <= 16 {
                    let len = values.len();
                    for i in 0..len {
                        let v = values[i];
                        s += &format!("{v}");
                        if i != len - 1 {
                            s += " ";
                        }
                    }
                } else {
                    s += "...";
                }
            }
        }
        s += "]";
        return s;
    }

    fn convert_params(&self, params: &ParamSet) -> String {
        let mut s = String::from("");
        let keys = &params.keys;
        if !keys.is_empty() {
            let indent = self.get_indent_i(self.indent.get() + 1);
            for key in keys {
                let values = self.convert_values(params, key);
                s += "\n";
                s += &format!("{indent}\"{key}\" {values}");
            }
        }
        return s;
    }

    fn with_params(&self, params: &ParamSet) -> String {
        if !params.keys.is_empty() {
            return format!(" {}", self.convert_params(params));
        } else {
            return String::from("");
        }
    }
}

impl ParseContext for PrintContext {
    fn pbrt_cleanup(&mut self) {
        self.print(&format!("{}Cleanup\n", self.get_indent()));
    }

    fn pbrt_identity(&mut self) {
        self.print(&format!("{}Identity\n", self.get_indent()));
    }

    fn pbrt_translate(&mut self, dx: Float, dy: Float, dz: Float) {
        self.print(&format!("{}Translate {dx} {dy} {dz}\n", self.get_indent()));
    }

    fn pbrt_rotate(&mut self, angle: Float, ax: Float, ay: Float, az: Float) {
        self.print(&format!(
            "{}Rotate {angle} {ax} {ay} {az}\n",
            self.get_indent()
        ));
    }

    fn pbrt_scale(&mut self, sx: Float, sy: Float, sz: Float) {
        self.print(&format!("{}Scale {sx} {sy} {sz}\n", self.get_indent()));
    }

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
    ) {
        self.print(&format!(
            "{}LookAt {ex} {ey} {ez} {lx} {ly} {lz} {ux} {uy} {uz}\n",
            self.get_indent()
        ));
    }

    fn pbrt_concat_transform(&mut self, transform: &[Float]) {
        let s_transform = self.convert_transform(transform);
        self.print(&format!(
            "{}ConcatTransform {s_transform}\n",
            self.get_indent()
        ));
    }

    fn pbrt_transform(&mut self, transform: &[Float]) {
        let s_transform = self.convert_transform(transform);
        self.print(&format!("{}Transform {s_transform}\n", self.get_indent()));
    }

    fn pbrt_coordinate_system(&mut self, name: &str) {
        self.print(&format!(
            "{}CoordinateSystem \"{name}\"\n",
            self.get_indent()
        ));
    }

    fn pbrt_coord_sys_transform(&mut self, name: &str) {
        self.print(&format!(
            "{}CoordSysTransform \"{name}\"\n",
            self.get_indent()
        ));
    }

    fn pbrt_active_transform_all(&mut self) {
        self.print(&format!("{}ActiveTransfrom All\n", self.get_indent()));
    }

    fn pbrt_active_transform_end_time(&mut self) {
        self.print(&format!("{}ActiveTransfrom EndTime\n", self.get_indent()));
    }

    fn pbrt_active_transform_start_time(&mut self) {
        self.print(&format!("{}ActiveTransfrom StartTime\n", self.get_indent()));
    }

    fn pbrt_transform_times(&mut self, start: Float, end: Float) {
        self.print(&format!(
            "{}TransformTimes {start} {end}\n",
            self.get_indent()
        ));
    }

    fn pbrt_pixel_filter(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}PixelFilter \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_film(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!("{}Film \"{name}\"{s_params}\n", self.get_indent()));
    }

    fn pbrt_sampler(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}Sampler \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_accelerator(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}Accelerator \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_integrator(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}Integrator \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_camera(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}Camera \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_make_named_medium(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}NamedMedium \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_medium_interface(&mut self, inside_name: &str, outside_name: &str) {
        self.print(&format!(
            "{}MediumInterface \"{inside_name}\" \"{outside_name}\"\n",
            self.get_indent()
        ));
    }

    fn pbrt_world_begin(&mut self) {
        self.print(&format!("{}WorldBegin\n", self.get_indent()));
        self.inc_indent();
    }

    fn pbrt_attribute_begin(&mut self) {
        self.print(&format!("{}AttributeBegin\n", self.get_indent()));
        self.inc_indent();
    }

    fn pbrt_attribute_end(&mut self) {
        self.dec_indent();
        self.print(&format!("{}AttributeEnd\n", self.get_indent()));
    }

    fn pbrt_transform_begin(&mut self) {
        self.print(&format!("{}TransformBegin\n", self.get_indent()));
        self.inc_indent();
    }

    fn pbrt_transform_end(&mut self) {
        self.dec_indent();
        self.print(&format!("{}TransformEnd\n", self.get_indent()));
    }

    fn pbrt_texture(&mut self, name: &str, t: &str, tex_name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}Texture \"{name}\" \"{t}\" \"{tex_name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_material(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}Material \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_make_named_material(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}MakeNamedMaterial \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_named_material(&mut self, name: &str) {
        self.print(&format!("{}NamedMaterial \"{name}\"\n", self.get_indent()));
    }

    fn pbrt_light_source(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}LightSource \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_area_light_source(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}AreaLightSource \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_shape(&mut self, name: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}Shape \"{name}\"{s_params}\n",
            self.get_indent()
        ));
    }

    fn pbrt_reverse_orientation(&mut self) {
        self.print(&format!("{}ReverseOrientation\n", self.get_indent()));
    }

    fn pbrt_object_begin(&mut self, name: &str) {
        self.print(&format!("{}ObjectBegin \"{name}\"\n", self.get_indent()));
        self.inc_indent();
    }

    fn pbrt_object_end(&mut self) {
        self.dec_indent();
        self.print(&format!("{}ObjectEnd\n", self.get_indent()));
    }

    fn pbrt_object_instance(&mut self, name: &str) {
        self.print(&format!("{}ObjectInstance \"{name}\"\n", self.get_indent()));
    }

    fn pbrt_world_end(&mut self) {
        self.dec_indent();
        self.print(&format!("{}WorldEnd\n", self.get_indent()));
    }

    fn pbrt_parse_file(&mut self, _file_name: &str) {
        /*
        self.print(&format!(
            "{}Include \"{file_name}\"\n",
            self.get_indent()
        ));
        */
    }
    fn pbrt_parse_string(&mut self, _s: &str) {
        /*
        self.print(&format!(
            "{}Include \"{file_name}\"\n",
            self.get_indent()
        ));
        */
    }

    fn pbrt_work_dir_begin(&mut self, _path: &str) {
        //Do not anything!
    }

    fn pbrt_work_dir_end(&mut self) {
        //Do not anything!
    }

    fn pbrt_include(&mut self, filename: &str, params: &ParamSet) {
        let s_params = self.with_params(params);
        self.print(&format!(
            "{}Include \"{filename}\"{s_params}\n",
            self.get_indent()
        ));
    }
}
