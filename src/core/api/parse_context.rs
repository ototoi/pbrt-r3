use crate::core::pbrt::*;

//Sized
pub trait ParseContext {
    fn pbrt_cleanup(&mut self);
    fn pbrt_identity(&mut self);
    fn pbrt_translate(&mut self, dx: Float, dy: Float, dz: Float);
    fn pbrt_rotate(&mut self, angle: Float, ax: Float, ay: Float, az: Float);
    fn pbrt_scale(&mut self, sx: Float, sy: Float, sz: Float);
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
    fn pbrt_concat_transform(&mut self, tansform: &[Float]);
    fn pbrt_transform(&mut self, tansform: &[Float]);
    fn pbrt_coordinate_system(&mut self, name: &str);
    fn pbrt_coord_sys_transform(&mut self, name: &str);
    fn pbrt_active_transform_all(&mut self);
    fn pbrt_active_transform_end_time(&mut self);
    fn pbrt_active_transform_start_time(&mut self);
    fn pbrt_transform_times(&mut self, start: Float, end: Float);
    fn pbrt_pixel_filter(&mut self, name: &str, params: &ParamSet);
    fn pbrt_film(&mut self, name: &str, params: &ParamSet);
    fn pbrt_sampler(&mut self, name: &str, params: &ParamSet);
    fn pbrt_accelerator(&mut self, name: &str, params: &ParamSet);
    fn pbrt_integrator(&mut self, name: &str, params: &ParamSet);
    fn pbrt_camera(&mut self, name: &str, params: &ParamSet);
    fn pbrt_make_named_medium(&mut self, name: &str, params: &ParamSet);
    fn pbrt_medium_interface(&mut self, inside_name: &str, outside_name: &str);

    fn pbrt_world_begin(&mut self);
    fn pbrt_attribute_begin(&mut self);
    fn pbrt_attribute_end(&mut self);
    fn pbrt_transform_begin(&mut self);
    fn pbrt_transform_end(&mut self);
    fn pbrt_texture(&mut self, name: &str, _type: &str, tex_name: &str, params: &ParamSet);
    fn pbrt_material(&mut self, name: &str, params: &ParamSet);
    fn pbrt_make_named_material(&mut self, name: &str, params: &ParamSet);
    fn pbrt_named_material(&mut self, name: &str);
    fn pbrt_light_source(&mut self, name: &str, params: &ParamSet);
    fn pbrt_area_light_source(&mut self, name: &str, params: &ParamSet);
    fn pbrt_shape(&mut self, name: &str, params: &ParamSet);
    fn pbrt_reverse_orientation(&mut self);
    fn pbrt_object_begin(&mut self, name: &str);
    fn pbrt_object_end(&mut self);
    fn pbrt_object_instance(&mut self, name: &str);
    fn pbrt_world_end(&mut self);

    fn pbrt_parse_file(&mut self, filename: &str);
    fn pbrt_parse_string(&mut self, s: &str);

    //----------------------------------------
    fn pbrt_work_dir_begin(&mut self, _path: &str) {}
    fn pbrt_work_dir_end(&mut self) {}
    fn pbrt_include(&mut self, _filename: &str, _params: &ParamSet) {}
    //----------------------------------------
}
