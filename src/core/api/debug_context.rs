use super::parse_context::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use std::cell::{Cell, RefCell};

pub struct Operation {
    pub name: String,
    pub args: Vec<String>,
}

impl Operation {
    pub fn new(op: &str, args: &[String]) -> Self {
        Operation {
            name: op.to_string(),
            args: args.to_vec(),
        }
    }
}

#[derive(Default)]
pub struct DebugContext {
    pub operations: RefCell<Vec<Operation>>,
    pub indent: Cell<i32>,
}

impl DebugContext {
    pub fn new() -> Self {
        DebugContext {
            operations: RefCell::<Vec<Operation>>::new(Vec::<Operation>::new()),
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
        let mut s = String::new();
        let count = self.indent.get();
        for _ in 0..count {
            s += "    ";
        }
        return s;
    }
}

impl ParseContext for DebugContext {
    fn pbrt_cleanup(&mut self) {
        println!("{}pbrt_cleanup", self.get_indent());
        let v = vec![];
        self.operations
            .borrow_mut()
            .push(Operation::new("Cleanup", &v));
    }
    fn pbrt_identity(&mut self) {
        println!("{}pbrt_identity", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("Identitiy", &v));
    }
    fn pbrt_translate(&mut self, dx: Float, dy: Float, dz: Float) {
        println!("{}pbrt_translate:[{dx}, {dy}, {dz}]", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("Translate", &v));
    }
    fn pbrt_rotate(&mut self, angle: Float, ax: Float, ay: Float, az: Float) {
        println!(
            "{}pbrt_rotate:[{angle}, {ax}, {ay}, {az}]",
            self.get_indent()
        );
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("Rotate", &v));
    }
    fn pbrt_scale(&mut self, sx: Float, sy: Float, sz: Float) {
        println!("{}pbrt_scale:[{sx}, {sy}, {sz}]", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("Scale", &v));
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
        println!(
            "{}pbrt_look_at:[{ex}, {ey}, {ez}, {lx}, {ly}, {lz}, {ux}, {uy}, {uz}]",
            self.get_indent()
        );
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("LookAt", &v));
    }

    fn pbrt_concat_transform(&mut self, _tansform: &[Float]) {
        println!("{}pbrt_concat_transform", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("ConcatTransform", &v));
    }
    fn pbrt_transform(&mut self, _tansform: &[Float]) {
        println!("{}pbrt_transform", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("Transform", &v));
    }
    fn pbrt_coordinate_system(&mut self, name: &str) {
        println!("{}pbrt_coordinate_system:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("CoordinateSystem", &v));
    }
    fn pbrt_coord_sys_transform(&mut self, name: &str) {
        println!("{}pbrt_coord_sys_transform:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("CoordSysTransform", &v));
    }
    fn pbrt_active_transform_all(&mut self) {
        println!("{}pbrt_active_transform_all", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("ActiveTransformAll", &v));
    }
    fn pbrt_active_transform_end_time(&mut self) {
        println!("{}pbrt_active_transform_end_time", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("ActiveTransformEndTime", &v));
    }
    fn pbrt_active_transform_start_time(&mut self) {
        println!("{}pbrt_active_transform_start_time", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("ActiveTransformStartTime", &v));
    }
    fn pbrt_transform_times(&mut self, _start: Float, _end: Float) {
        println!("{}pbrt_transform_times", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("TransformTimes", &v));
    }
    fn pbrt_pixel_filter(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_pixel_filter:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("PixelFilter", &v));
    }
    fn pbrt_film(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_film:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("Film", &v));
    }
    fn pbrt_sampler(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_sampler:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("Sampler", &v));
    }
    fn pbrt_accelerator(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_accelerator:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("Accelerator", &v));
    }
    fn pbrt_integrator(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_integrator:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("Integrator", &v));
    }
    fn pbrt_camera(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_camera:\"{name}\"", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("Camera", &v));
    }
    fn pbrt_make_named_medium(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_make_named_medium:\"{name}\"", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("NamedMedium", &v));
    }
    fn pbrt_medium_interface(&mut self, inside_name: &str, outside_name: &str) {
        println!(
            "{}pbrt_medium_interface:\"{inside_name}\", \"{outside_name}\"",
            self.get_indent()
        );
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("MediumInterface", &v));
    }
    fn pbrt_world_begin(&mut self) {
        println!("{}pbrt_world_begin", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("WorldBegin", &v));
        self.inc_indent();
    }
    fn pbrt_attribute_begin(&mut self) {
        println!("{}pbrt_attribute_begin", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("AttributeBegin", &v));
        self.inc_indent();
    }
    fn pbrt_attribute_end(&mut self) {
        self.dec_indent();
        println!("{}pbrt_attribute_end", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("AttributeEnd", &v));
    }
    fn pbrt_transform_begin(&mut self) {
        println!("{}pbrt_transform_begin", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("TransformBegin", &v));
        self.inc_indent();
    }
    fn pbrt_transform_end(&mut self) {
        self.dec_indent();
        println!("{}pbrt_transform_end", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("TransformEnd", &v));
    }
    fn pbrt_texture(&mut self, name: &str, t: &str, tex_name: &str, _params: &ParamSet) {
        println!(
            "{}pbrt_texture:\"{}\", \"{}\", \"{}\"",
            self.get_indent(),
            name,
            t,
            tex_name
        );
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("Texture", &v));
    }
    fn pbrt_material(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_material:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("Material", &v));
    }
    fn pbrt_make_named_material(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_make_named_material:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("MakeNamedMaterial", &v));
    }
    fn pbrt_named_material(&mut self, name: &str) {
        println!("{}pbrt_named_material:\"{name}\"", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("NamedMaterial", &v));
    }
    fn pbrt_light_source(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_light_source:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("LightSource", &v));
    }
    fn pbrt_area_light_source(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_area_light_source:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("AreaLightSource", &v));
    }
    fn pbrt_shape(&mut self, name: &str, _params: &ParamSet) {
        println!("{}pbrt_shape:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("Shape", &v));
    }
    fn pbrt_reverse_orientation(&mut self) {
        println!("{}pbrt_reverse_orientation", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("ReverseOrientation", &v));
    }
    fn pbrt_object_begin(&mut self, name: &str) {
        println!("{}pbrt_object_begin:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("ObjectBegin", &v));
        self.inc_indent();
    }
    fn pbrt_object_end(&mut self) {
        self.dec_indent();
        println!("{}pbrt_object_end", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("ObjectEnd", &v));
    }
    fn pbrt_object_instance(&mut self, name: &str) {
        println!("{}pbrt_object_instance:\"{name}\"", self.get_indent());
        let v = vec![String::from(name)];
        self.operations
            .borrow_mut()
            .push(Operation::new("ObjectInstance", &v));
    }
    fn pbrt_world_end(&mut self) {
        self.dec_indent();
        println!("{}pbrt_world_end", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("WorldEnd", &v));
    }
    fn pbrt_parse_file(&mut self, file_name: &str) {
        println!("{}pbrt_parse_file:\"{file_name}\"", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("ParseFile", &v));
    }
    fn pbrt_parse_string(&mut self, _s: &str) {
        println!("{}pbrt_parse_string", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("ParseString", &v));
    }

    fn pbrt_work_dir_begin(&mut self, path: &str) {
        println!("{}pbrt_work_dir_begin:\"{path}\"", self.get_indent());
        let v = vec![String::from(path)];
        self.operations
            .borrow_mut()
            .push(Operation::new("WorkDirBegin", &v));
    }

    fn pbrt_work_dir_end(&mut self) {
        println!("{}pbrt_work_dir_end", self.get_indent());
        let v = vec![String::from("")];
        self.operations
            .borrow_mut()
            .push(Operation::new("WorkDirEnd", &v));
    }

    fn pbrt_include(&mut self, filename: &str, _params: &ParamSet) {
        println!("{}pbrt_include:\"{filename}\"", self.get_indent());
        let v = vec![String::from(filename)];
        self.operations
            .borrow_mut()
            .push(Operation::new("Include", &v));
    }
}
