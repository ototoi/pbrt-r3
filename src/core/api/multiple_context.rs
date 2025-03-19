use super::parse_context::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use std::cell::RefCell;
use std::sync::Arc;

pub struct MutipleContext {
    pub contexts: Vec<Arc<RefCell<dyn ParseContext>>>,
}

impl MutipleContext {
    pub fn new() -> Self {
        MutipleContext {
            contexts: Vec::new(),
        }
    }
    pub fn add(&mut self, context: Arc<RefCell<dyn ParseContext>>) {
        self.contexts.push(context);
    }
}

impl ParseContext for MutipleContext {
    fn pbrt_cleanup(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_cleanup();
        }
    }

    fn pbrt_identity(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_identity();
        }
    }

    fn pbrt_translate(&mut self, dx: Float, dy: Float, dz: Float) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_translate(dx, dy, dz);
        }
    }

    fn pbrt_rotate(&mut self, angle: Float, ax: Float, ay: Float, az: Float) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_rotate(angle, ax, ay, az);
        }
    }

    fn pbrt_scale(&mut self, sx: Float, sy: Float, sz: Float) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_scale(sx, sy, sz);
        }
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
        for r in self.contexts.iter() {
            r.borrow_mut()
                .pbrt_look_at(ex, ey, ez, lx, ly, lz, ux, uy, uz);
        }
    }

    fn pbrt_concat_transform(&mut self, transform: &[Float]) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_concat_transform(transform);
        }
    }

    fn pbrt_transform(&mut self, transform: &[Float]) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_transform(transform);
        }
    }

    fn pbrt_coordinate_system(&mut self, name: &str) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_coordinate_system(name);
        }
    }

    fn pbrt_coord_sys_transform(&mut self, name: &str) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_coord_sys_transform(name);
        }
    }

    fn pbrt_active_transform_all(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_active_transform_all();
        }
    }

    fn pbrt_active_transform_end_time(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_active_transform_end_time();
        }
    }

    fn pbrt_active_transform_start_time(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_active_transform_start_time();
        }
    }

    fn pbrt_transform_times(&mut self, start: Float, end: Float) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_transform_times(start, end);
        }
    }

    fn pbrt_pixel_filter(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_pixel_filter(name, params);
        }
    }

    fn pbrt_film(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_film(name, params);
        }
    }

    fn pbrt_sampler(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_sampler(name, params);
        }
    }

    fn pbrt_accelerator(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_accelerator(name, params);
        }
    }

    fn pbrt_integrator(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_integrator(name, params);
        }
    }
    fn pbrt_camera(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_camera(name, params);
        }
    }

    fn pbrt_make_named_medium(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_make_named_medium(name, params);
        }
    }

    fn pbrt_medium_interface(&mut self, inside_name: &str, outside_name: &str) {
        for r in self.contexts.iter() {
            r.borrow_mut()
                .pbrt_medium_interface(inside_name, outside_name);
        }
    }

    fn pbrt_world_begin(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_world_begin();
        }
    }

    fn pbrt_attribute_begin(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_attribute_begin();
        }
    }

    fn pbrt_attribute_end(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_attribute_end();
        }
    }

    fn pbrt_transform_begin(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_transform_begin();
        }
    }

    fn pbrt_transform_end(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_transform_end();
        }
    }

    fn pbrt_texture(&mut self, name: &str, t: &str, tex_name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_texture(name, t, tex_name, params);
        }
    }

    fn pbrt_material(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_material(name, params);
        }
    }

    fn pbrt_make_named_material(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_make_named_material(name, params);
        }
    }

    fn pbrt_named_material(&mut self, name: &str) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_named_material(name);
        }
    }

    fn pbrt_light_source(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_light_source(name, params);
        }
    }

    fn pbrt_area_light_source(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_area_light_source(name, params);
        }
    }

    fn pbrt_shape(&mut self, name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_shape(name, params);
        }
    }

    fn pbrt_reverse_orientation(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_reverse_orientation();
        }
    }

    fn pbrt_object_begin(&mut self, name: &str) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_object_begin(name);
        }
    }

    fn pbrt_object_end(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_object_end();
        }
    }
    fn pbrt_object_instance(&mut self, name: &str) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_object_instance(name);
        }
    }

    fn pbrt_world_end(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_world_end();
        }
    }

    fn pbrt_parse_file(&mut self, file_name: &str) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_parse_file(file_name);
        }
    }

    fn pbrt_parse_string(&mut self, s: &str) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_parse_string(s);
        }
    }

    fn pbrt_work_dir_begin(&mut self, path: &str) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_work_dir_begin(path);
        }
    }

    fn pbrt_work_dir_end(&mut self) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_work_dir_end();
        }
    }

    fn pbrt_include(&mut self, file_name: &str, params: &ParamSet) {
        for r in self.contexts.iter() {
            r.borrow_mut().pbrt_include(file_name, params);
        }
    }
}
