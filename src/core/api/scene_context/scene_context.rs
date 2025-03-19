use super::graphics_state::GraphicsState;
use super::graphics_state::MaterialInstance;
use super::render_options::RenderOptions;

use crate::accelerators::*;
use crate::cameras::*;
use crate::core::api::parse_context::*;
use crate::core::camera::*;
use crate::core::error::*;
use crate::core::film::*;
use crate::core::filter::*;
use crate::core::integrator::*;
use crate::core::light::*;
use crate::core::material::*;
use crate::core::medium::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::primitive::*;
use crate::core::scene::*;
use crate::core::shape::*;
use crate::core::spectrum::*;
use crate::core::stats::*;
use crate::core::texture::*;
use crate::core::transform::transform_set::*;
use crate::core::transform::*;
use crate::filters::*;
use crate::integrators::*;
use crate::lights::*;
use crate::materials::*;
use crate::media::*;
use crate::samplers::*;
use crate::shapes::*;
use crate::textures::*;

use log::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::ops::Deref;
use std::path::Path;
use std::sync::Arc;
use std::sync::RwLock;

thread_local!(static N_SHAPES: StatCounter = StatCounter::new("Scene/Shapes")); //"Scene/Shapes created"
thread_local!(static N_LIGHTS: StatCounter = StatCounter::new("Scene/Lights"));
thread_local!(static N_AREA_LIGHTS: StatCounter = StatCounter::new("Scene/AreaLights"));
thread_local!(static N_MATERIALS: StatCounter = StatCounter::new("Scene/Materials")); //"Scene/Materials created"
thread_local!(static N_OBJECT_INSTANCES_CREATED: StatCounter = StatCounter::new("Scene/Object instances created")); //"Scene/Object instances created"
thread_local!(static N_OBJECT_INSTANCES_USED: StatCounter = StatCounter::new("Scene/Object instances used")); //"Scene/Object instances created"

#[derive(Debug, PartialEq)]
enum APIState {
    OptionsBlock,
    WorldBlock,
}

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

type FloatTextureMap = HashMap<String, Arc<dyn Texture<Float>>>;
//type SpectrumTextureMap = HashMap<String, Arc<dyn Texture<Spectrum>>>;

pub struct SceneContext {
    current_api_state: APIState,
    transforms: Vec<RefCell<TransformSet>>,
    transform_bits: Vec<u32>,
    graphics_states: Vec<RefCell<GraphicsState>>,
    render_options: RefCell<RenderOptions>,
    named_coordinate_systems: HashMap<String, RefCell<TransformSet>>,
    float_textures: RefCell<HashMap<String, Arc<dyn Texture<Float>>>>,
    spectrum_textures: RefCell<HashMap<String, Arc<dyn Texture<Spectrum>>>>,
    work_dirs: Vec<String>,
}

impl SceneContext {
    pub fn new() -> Self {
        let mut ctx = SceneContext {
            current_api_state: APIState::OptionsBlock,
            transforms: Vec::new(),
            transform_bits: Vec::new(),
            graphics_states: Vec::new(),
            render_options: RefCell::new(RenderOptions::new()),
            named_coordinate_systems: HashMap::new(),
            float_textures: RefCell::new(HashMap::new()),
            spectrum_textures: RefCell::new(HashMap::new()),
            work_dirs: Vec::new(),
        };
        ctx.initialize();
        return ctx;
    }

    pub fn initialize(&mut self) {
        self.transforms.clear();
        self.transform_bits.clear();
        self.graphics_states.clear();
        self.transforms.push(RefCell::new(TransformSet::new()));
        self.transform_bits.push(ALL_TRANSFORM_BITS);
        self.graphics_states
            .push(RefCell::new(GraphicsState::new()));
    }

    pub fn push_transform(&mut self) {
        self.transforms
            .push(self.transforms[self.transforms.len() - 1].clone());
        self.transform_bits
            .push(self.transform_bits[self.transform_bits.len() - 1]);
    }

    pub fn pop_transform(&mut self) {
        self.transform_bits.pop();
        self.transforms.pop();
    }

    pub fn push_graphics_state(&mut self) {
        self.graphics_states
            .push(self.graphics_states[self.transforms.len() - 1].clone());
    }

    pub fn pop_graphics_state(&mut self) {
        self.graphics_states.pop();
    }

    pub fn print_error(&self, e: &PbrtError) {
        match e.kind {
            PbrtErrorKind::Error => {
                error!("{}", e.msg);
            }
            PbrtErrorKind::Warning => {
                warn!("{}", e.msg);
            }
        }
    }

    pub fn verify_initialized(&mut self, _fun: &str) {
        //
    }

    pub fn verify_options(&mut self, fun: &str) {
        if self.current_api_state != APIState::OptionsBlock {
            error!(
                "Options cannot be set inside world block;\n\"{}\" not allowed. Ignoring.",
                fun
            );
        }
    }

    pub fn verify_world(&mut self, fun: &str) {
        if self.current_api_state == APIState::WorldBlock {
            error!(
                "Scene description must be inside world block;\n\"{}\" not allowed. Ignoring.",
                fun
            );
        }
    }

    pub fn replace_params_int(
        &mut self,
        name: &str,
        key: &str,
        value: i32,
    ) -> Result<(), PbrtError> {
        match name {
            "Sampler" => {
                let mut opts = self.render_options.borrow_mut();
                let params = &mut opts.sampler_params;
                params.replace_one_int(key, value);
                return Ok(());
            }
            _ => {
                let msg = format!("\"{}\" is unknown params name.", name);
                return Err(PbrtError::error(&msg));
            }
        }
    }

    pub fn replace_params_string(
        &mut self,
        name: &str,
        key: &str,
        value: &str,
    ) -> Result<(), PbrtError> {
        match name {
            "Film" => {
                let mut opts = self.render_options.borrow_mut();
                let params = &mut opts.film_params;
                params.replace_one_string(key, value);
                return Ok(());
            }
            _ => {
                let msg = format!("\"{}\" is unknown params name.", name);
                return Err(PbrtError::error(&msg));
            }
        }
    }

    pub fn create_medium_interface(&self) -> MediumInterface {
        let attr = self.graphics_states[self.graphics_states.len() - 1].borrow();
        let opts = self.render_options.borrow();
        let inside_name = attr.current_inside_medium.clone();
        let outside_name = attr.current_outside_medium.clone();
        let mut m = MediumInterface::new();
        if let Some(medium) = opts.named_media.get(&inside_name) {
            //println!("inside_name: {}", inside_name);
            m.set_inside(medium);
        }
        if let Some(medium) = opts.named_media.get(&outside_name) {
            //println!("outside_name: {}", outside_name);
            m.set_outside(medium);
        }
        return m;
    }

    pub fn shape_may_set_material_parameters(params: &ParamSet) -> bool {
        for name in params.strings.keys() {
            if name != "alpha" && name != "shadowalpha" {
                return true;
            }
        }
        for (name, p) in &params.floats {
            let param = p.borrow();
            if param.len() == 1 && name != "radius" {
                return true;
            }
        }
        for (name, p) in &params.strings {
            let param = p.borrow();
            if param.len() == 1 && name != "filename" && name != "type" && name != "scheme" {
                return true;
            }
        }
        for p in params.bools.values() {
            let param = p.borrow();
            if param.len() == 1 {
                return true;
            }
        }
        for p in params.ints.values() {
            let param = p.borrow();
            if param.len() == 1 {
                return true;
            }
        }
        for p in params.spectrums.values() {
            let param = p.borrow();
            if param.len() == 1 {
                return true;
            }
        }
        for p in params.points.values() {
            let param = p.borrow();
            if param.len() == 3 {
                return true;
            }
        }
        return false;
    }

    pub fn get_material_for_shape(&self, params: &ParamSet) -> Option<Arc<dyn Material>> {
        if Self::shape_may_set_material_parameters(params) {
            let attr = self.graphics_states[self.graphics_states.len() - 1].borrow();
            if let Some(mat_inst) = attr.current_material.as_ref() {
                let mat_inst = mat_inst.borrow();
                let f_tex = attr.float_textures.as_ref().borrow();
                let s_tex = attr.spectrum_textures.as_ref().borrow();
                let mp = TextureParams::new(params, &mat_inst.params, &f_tex, &s_tex);
                match self.make_material(&mat_inst.name, &mp) {
                    Ok(material) => {
                        return Some(material);
                    }
                    Err(e) => {
                        self.print_error(&e);
                        return None;
                    }
                }
            }
            return None;
        } else {
            let attr = self.graphics_states[self.graphics_states.len() - 1].borrow();
            if let Some(mat_inst) = attr.current_material.as_ref() {
                let mat_inst = mat_inst.borrow();
                return Some(Arc::clone(&mat_inst.material));
            }
        }
        return None;
    }

    pub fn get_named_material(&self, name: &str) -> Option<Arc<RefCell<MaterialInstance>>> {
        let attr = self.graphics_states[self.graphics_states.len() - 1].borrow();
        let materials = attr.named_materials.as_ref().borrow();
        let mat = materials.get(name);
        match mat {
            Some(m) => {
                return Some(Arc::clone(m));
            }
            None => None,
        }
    }

    pub fn get_current_instance_name(&self) -> Option<String> {
        let opts = self.render_options.borrow();
        match opts.current_instance_name.as_ref() {
            Some(name) => {
                return Some(name.clone());
            }
            None => None,
        }
    }

    pub fn get_filepath(&self, name: &str) -> Option<String> {
        let path = Path::new(name);
        if path.is_absolute() || self.work_dirs.is_empty() {
            return Some(String::from(path.to_str().unwrap()));
        } else {
            for d in self.work_dirs.iter().rev() {
                let dir = Path::new(d);
                let full_path = dir.join(name);
                if full_path.exists() {
                    return Some(String::from(full_path.to_str().unwrap()));
                }
            }
            return None;
        }
    }

    pub fn create_shapes(
        &self,
        name: &str,
        object2world: &Transform,
        reverse_orientation: bool,
        params: &ParamSet,
        float_textures: &FloatTextureMap,
    ) -> Result<Vec<Arc<dyn Shape>>, PbrtError> {
        match name {
            "trianglemesh" => {
                //Do not use make_absolute_path
                return create_shapes(
                    name,
                    object2world,
                    reverse_orientation,
                    params,
                    float_textures,
                );
            }
            "curve" => {
                //Do not use make_absolute_path
                return create_shapes(
                    name,
                    object2world,
                    reverse_orientation,
                    params,
                    float_textures,
                );
            }
            _ => {
                let params = self.make_absolute_path(params);
                return create_shapes(
                    name,
                    object2world,
                    reverse_orientation,
                    &params,
                    float_textures,
                );
            }
        }
    }

    pub fn make_shapes_core(
        &self,
        name: &str,
        object2world: &Transform,
        reverse_orientation: bool,
        params: &ParamSet,
    ) -> Vec<Arc<dyn Shape>> {
        let attr = self.graphics_states[self.graphics_states.len() - 1].borrow();
        let f_tex = attr.float_textures.borrow();
        match self.create_shapes(
            name,
            object2world,
            reverse_orientation,
            params,
            f_tex.deref(),
        ) {
            Ok(shapes) => {
                N_SHAPES.with(|c| c.add(shapes.len() as u64));
                return shapes;
            }
            Err(e) => {
                self.print_error(&e);
                return Vec::new();
            }
        }
    }

    pub fn make_shapes(
        &self,
        name: &str,
        object2world: &Transform,
        reverse_orientation: bool,
        params: &ParamSet,
    ) -> Vec<Arc<dyn Shape>> {
        let attr = self.graphics_states[self.graphics_states.len() - 1].borrow();
        if !attr.area_light_name.is_empty() {
            let two_sided = attr.area_light_params.find_one_bool("twosided", true);
            let two_sided = params.find_one_bool("twosided", two_sided);
            let mut params = params.clone();
            params.replace_one_bool("bool twosided", two_sided);
            return self.make_shapes_core(name, object2world, reverse_orientation, &params);
        }
        return self.make_shapes_core(name, object2world, reverse_orientation, params);
    }

    pub fn get_named_material_matte(&self, mp: &TextureParams, key: &str) -> Arc<dyn Material> {
        if let Some(mat) = self.get_named_material(key) {
            return Arc::clone(&mat.borrow().material);
        } else {
            //
            return self.make_material("matte", mp).unwrap();
        }
    }

    pub fn make_material(
        &self,
        name: &str,
        mp: &TextureParams,
    ) -> Result<Arc<dyn Material>, PbrtError> {
        if name == "subsurface" || name == "kdsubsurface" {
            let render_options = self.render_options.borrow();
            if render_options.integrator_name != "path"
                && render_options.integrator_name != "volpath"
            {
                warn!("Subsurface scattering material \"{}\" used, but \"{}\" integrator doesn't support subsurface scattering. Use \"path\" or \"volpath\".", name, render_options.integrator_name);
            }
        }

        let material = if name == "mix" {
            let m1 = mp.find_string("namedmaterial1", "");
            let m2 = mp.find_string("namedmaterial2", "");
            let mat1 = self.get_named_material_matte(mp, &m1);
            let mat2 = self.get_named_material_matte(mp, &m2);
            create_mix_material(mp, &mat1, &mat2)?
        } else {
            create_material(name, mp)?
        };
        N_MATERIALS.with(|c| c.inc());
        return Ok(material);
    }

    fn make_material_params(
        &self,
        name: &str,
        geom_params: &ParamSet,
        mat_params: &ParamSet,
    ) -> Option<Arc<dyn Material>> {
        let attr = self.graphics_states.last().unwrap().borrow();
        let f_tex = attr.float_textures.clone();
        let s_tex = attr.spectrum_textures.clone();
        let f = f_tex.as_ref().borrow();
        let s = s_tex.as_ref().borrow();
        let mp = TextureParams::new(geom_params, mat_params, f.deref(), s.deref());
        match self.make_material(name, &mp) {
            Ok(material) => {
                return Some(material);
            }
            Err(e) => {
                self.print_error(&e);
                return None;
            }
        }
    }

    pub fn make_float_texture(
        &self,
        name: &str,
        tex2world: &Transform,
        tp: &TextureParams,
    ) -> Option<Arc<dyn Texture<Float>>> {
        if let Some(texinfo) = create_texinfo(tp) {
            let key = format!("{:?}", texinfo);
            let textures = self.float_textures.borrow();
            if let Some(texture) = textures.get(&key) {
                return Some(Arc::clone(texture));
            }
        }

        match create_float_texture(name, tex2world, tp) {
            Ok(texture) => {
                if let Some(texinfo) = create_texinfo(tp) {
                    let key = format!("{:?}", texinfo);
                    let mut textures = self.float_textures.borrow_mut();
                    textures.insert(key, Arc::clone(&texture));
                }

                return Some(texture);
            }
            Err(e) => {
                self.print_error(&e);
                return None;
            }
        }
    }

    pub fn make_spectrum_texture(
        &self,
        name: &str,
        tex2world: &Transform,
        tp: &TextureParams,
    ) -> Option<Arc<dyn Texture<Spectrum>>> {
        if let Some(texinfo) = create_texinfo(tp) {
            let key = format!("{:?}", texinfo);
            let textures = self.spectrum_textures.borrow();
            if let Some(texture) = textures.get(&key) {
                return Some(Arc::clone(texture));
            }
        }

        match create_spectrum_texture(name, tex2world, tp) {
            Ok(texture) => {
                if let Some(texinfo) = create_texinfo(tp) {
                    let key = format!("{:?}", texinfo);
                    let mut textures = self.spectrum_textures.borrow_mut();
                    textures.insert(key, Arc::clone(&texture));
                }

                return Some(texture);
            }
            Err(e) => {
                self.print_error(&e);
                return None;
            }
        }
    }

    pub fn make_light(
        &self,
        name: &str,
        light2world: &Transform,
        medium_interface: &MediumInterface,
        params: &ParamSet,
    ) -> Option<Arc<dyn Light>> {
        match create_light(name, light2world, medium_interface, params) {
            Ok(light) => {
                N_LIGHTS.with(|c| c.inc());
                return Some(light);
            }
            Err(e) => {
                self.print_error(&e);
                return None;
            }
        }
    }

    pub fn make_area_light(
        &self,
        name: &str,
        light2world: &Transform,
        medium_interface: &MediumInterface,
        params: &ParamSet,
        shape: &Arc<dyn Shape>,
    ) -> Option<Arc<dyn Light>> {
        let attr = self.graphics_states[self.graphics_states.len() - 1].borrow();
        if attr.area_light_name.is_empty() {
            return None;
        }
        match create_area_light(name, light2world, medium_interface, params, shape) {
            Ok(light) => {
                N_LIGHTS.with(|c| c.inc());
                N_AREA_LIGHTS.with(|c| c.inc());
                return Some(light);
            }
            Err(e) => {
                self.print_error(&e);
                return None;
            }
        }
    }

    pub fn make_accelerator(
        &self,
        name: &str,
        prims: &[Arc<dyn Primitive>],
        params: &ParamSet,
    ) -> Option<Arc<dyn Primitive>> {
        match create_accelerator(name, prims, params) {
            Ok(accelerator) => {
                return Some(accelerator);
            }
            Err(_e) => {
                return None;
            }
        }
    }

    pub fn make_medium(
        &self,
        name: &str,
        params: &ParamSet,
        medium2world: &Transform,
    ) -> Option<Arc<dyn Medium>> {
        return create_medium(name, params, medium2world);
    }

    pub fn make_filter(&self) -> Option<Arc<dyn Filter>> {
        let opts = self.render_options.borrow();
        match create_filter(&opts.filter_name, &opts.filter_params) {
            Ok(filter) => {
                return Some(filter);
            }
            Err(e) => {
                self.print_error(&e);
                return None;
            }
        }
    }

    pub fn make_film(&self, filter: &Arc<dyn Filter>) -> Option<Arc<RwLock<Film>>> {
        let opts = self.render_options.borrow();
        match create_film(&opts.film_name, &opts.film_params, filter) {
            Ok(film) => {
                return Some(film);
            }
            Err(e) => {
                self.print_error(&e);
                return None;
            }
        }
    }

    pub fn make_camera(&self) -> Option<Arc<dyn Camera>> {
        if let Some(camera_transform) = self.named_coordinate_systems.get("camera") {
            let camera_to_world = &camera_transform.borrow();
            if let Some(filter) = self.make_filter() {
                if let Some(film) = self.make_film(&filter) {
                    let medium_interface = self.create_medium_interface();
                    let opts = self.render_options.borrow();
                    let start_time = opts.transform_start_time;
                    let end_time = opts.transform_end_time;
                    let animated_cam_2_world = AnimatedTransform::new(
                        &camera_to_world[0],
                        start_time,
                        &camera_to_world[1],
                        end_time,
                    );
                    match create_camera(
                        &opts.camera_name,
                        &opts.camera_params,
                        &animated_cam_2_world,
                        &film,
                        &medium_interface,
                    ) {
                        Ok(camara) => {
                            return Some(camara);
                        }
                        Err(e) => {
                            self.print_error(&e);
                            return None;
                        }
                    }
                }
            }
        }
        return None;
    }

    pub fn make_integrator(&self) -> Option<Arc<RwLock<dyn Integrator>>> {
        let camera = self.make_camera();
        if camera.is_none() {
            error!("Unable to create camera");
            return None;
        }
        let camera = camera.unwrap();
        let film = camera.as_ref().get_film();
        let opts = self.render_options.borrow();
        let sampler = create_sampler(&opts.sampler_name, &opts.sampler_params, &film);
        if sampler.is_err() {
            self.print_error(&sampler.err().unwrap());
            return None;
        }
        let sampler = sampler.unwrap();

        if opts.have_scattering_media
            && opts.integrator_name != "volpath"
            && opts.integrator_name != "bdpt"
            && opts.integrator_name != "mlt"
        {
            warn!("Scene has scattering media but \"{}\" integrator doesn't support volume scattering. Use \"volpath\", \"bdpt\" or \"mlt\".", opts.integrator_name);
        }
        // Warn if no light sources are defined
        if opts.lights.is_empty() {
            warn!("No light sources defined in scene; rendering a black image.");
        }
        match create_integrator(
            &opts.integrator_name,
            &opts.integrator_params,
            &sampler,
            &camera,
        ) {
            Ok(integrator) => {
                return Some(integrator);
            }
            Err(e) => {
                self.print_error(&e);
                return None;
            }
        }
    }

    pub fn make_scene(&self) -> Option<Arc<Scene>> {
        let opts = self.render_options.borrow();
        if let Some(accelerator) = self.make_accelerator(
            &opts.accelerator_name,
            &opts.primitives,
            &opts.accelerator_params,
        ) {
            let scene = Arc::new(Scene::new(&accelerator, &opts.lights));
            return Some(scene);
        }
        return None;
    }

    //AbsolutePath
    pub fn make_absolute_path(&self, params: &ParamSet) -> ParamSet {
        let mut n_params = params.clone();
        let keys = params.get_keys();
        //Spectrum
        {
            let target_keys: Vec<_> = keys
                .iter()
                .filter(|key| -> bool {
                    let (t, _) = get_param_type(key);
                    return t == "spectrum";
                })
                .collect();
            for key in target_keys {
                if let Some(names) = params.get_strings_ref(key) {
                    for name in names.iter() {
                        if let Some(path) = self.get_filepath(name) {
                            let ret = SampledSpectrum::load_sampled_spectrum_file(&path);
                            match ret {
                                Ok(spc) => {
                                    let spc = Spectrum::from(&spc);
                                    n_params.add_spectrum_no_key(name, &spc);
                                }
                                Err(e) => {
                                    self.print_error(&e);
                                }
                            }
                        }
                    }
                }
            }
        }
        // blackbody
        {
            let target_keys: Vec<_> = keys
                .iter()
                .filter(|key| -> bool {
                    let (t, _) = get_param_type(key);
                    return t == "blackbody";
                })
                .collect();
            for key in target_keys {
                if let Some(values) = params.get_floats_ref(key) {
                    let spc = SampledSpectrum::from_blackbody(&values);
                    let spc = Spectrum::from(&spc);
                    n_params.add_spectrum_no_key(key, &spc);
                }
            }
        }

        //filaname && bsdffile
        {
            let file_path_keys = ["filename", "mapname", "bsdffile", "lensfile"];

            let target_keys: Vec<_> = keys
                .iter()
                .filter(|key| -> bool {
                    let (_, key) = get_param_type(key);
                    for file_path_key in file_path_keys {
                        if file_path_key == key {
                            return true;
                        }
                    }
                    return false;
                })
                .collect();
            let mut replaces = Vec::new();
            for key in target_keys {
                if let Some(names) = params.get_strings_ref(key) {
                    for name in names.iter() {
                        if let Some(path) = self.get_filepath(name) {
                            replaces.push((key.clone(), path));
                        }
                    }
                }
            }
            for (key, value) in replaces.iter() {
                n_params.replace_one_string(key, value);
            }
        }

        return n_params;
    }
}

impl ParseContext for SceneContext {
    fn pbrt_cleanup(&mut self) {
        self.initialize();
    }

    fn pbrt_identity(&mut self) {
        self.verify_initialized("Identity");

        let t = Transform::identity();
        let mut trns = self.transforms[self.transforms.len() - 1].borrow_mut();
        trns.set_transform(&t, self.transform_bits[self.transform_bits.len() - 1]);
    }

    fn pbrt_translate(&mut self, dx: Float, dy: Float, dz: Float) {
        self.verify_initialized("Translate");

        let t = Transform::translate(dx, dy, dz);
        let mut trns = self.transforms[self.transforms.len() - 1].borrow_mut();
        trns.mul_transform(&t, self.transform_bits[self.transform_bits.len() - 1]);
    }

    fn pbrt_rotate(&mut self, angle: Float, ax: Float, ay: Float, az: Float) {
        self.verify_initialized("Rotate");

        let t = Transform::rotate(angle, ax, ay, az);
        let mut trns = self.transforms[self.transforms.len() - 1].borrow_mut();
        trns.mul_transform(&t, self.transform_bits[self.transform_bits.len() - 1]);
    }

    fn pbrt_scale(&mut self, sx: Float, sy: Float, sz: Float) {
        self.verify_initialized("Scale");

        let t = Transform::scale(sx, sy, sz);
        let mut trns = self.transforms[self.transforms.len() - 1].borrow_mut();
        trns.mul_transform(&t, self.transform_bits[self.transform_bits.len() - 1]);
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
        self.verify_initialized("LookAt");

        let t = Transform::look_at(ex, ey, ez, lx, ly, lz, ux, uy, uz); //w2c
        let mut trns = self.transforms[self.transforms.len() - 1].borrow_mut();
        trns.mul_transform(&t, self.transform_bits[self.transform_bits.len() - 1]);
    }

    fn pbrt_concat_transform(&mut self, t: &[Float]) {
        self.verify_initialized("ConcatTransform");

        #[rustfmt::skip]
        let m = Matrix4x4::from([
            t[0], t[4], t[8], t[12],
            t[1], t[5], t[9], t[13],
            t[2], t[6], t[10], t[14],
            t[3], t[7], t[11], t[15],
        ]);
        if let Some(im) = m.inverse() {
            let t = Transform::from((m, im));
            let mut trns = self.transforms[self.transforms.len() - 1].borrow_mut();
            trns.mul_transform(&t, self.transform_bits[self.transform_bits.len() - 1]);
        } else {
            error!("Singular matrix in MatrixInvert");
        }
    }

    fn pbrt_transform(&mut self, t: &[Float]) {
        self.verify_initialized("Transform");

        #[rustfmt::skip]
        let m = Matrix4x4::from([
            t[0], t[4], t[8], t[12],
            t[1], t[5], t[9], t[13],
            t[2], t[6], t[10], t[14],
            t[3], t[7], t[11], t[15],
        ]);
        if let Some(im) = m.inverse() {
            let t = Transform::from((m, im));
            let mut trns = self.transforms[self.transforms.len() - 1].borrow_mut();
            trns.set_transform(&t, self.transform_bits[self.transform_bits.len() - 1]);
        } else {
            error!("Singular matrix in MatrixInvert");
        }
    }

    fn pbrt_coordinate_system(&mut self, name: &str) {
        self.verify_initialized("CoordinateSystem");

        let trns = self.transforms[self.transforms.len() - 1].borrow();
        self.named_coordinate_systems
            .insert(String::from(name), RefCell::new(trns.clone()));
    }

    fn pbrt_coord_sys_transform(&mut self, name: &str) {
        self.verify_initialized("CoordSysTransform");

        if let Some(t) = self.named_coordinate_systems.get(&String::from(name)) {
            let t = t.borrow();
            let mut trns = self.transforms[self.transforms.len() - 1].borrow_mut();
            trns.set(&t, self.transform_bits[self.transform_bits.len() - 1]); //tmx
        } else {
            warn!("Couldn't find named coordinate system \"{}\"", name);
        }
    }

    fn pbrt_active_transform_all(&mut self) {
        let last = self.transform_bits.len() - 1;
        self.transform_bits[last] = ALL_TRANSFORM_BITS;
    }

    fn pbrt_active_transform_end_time(&mut self) {
        let last = self.transform_bits.len() - 1;
        self.transform_bits[last] = END_TRANSFORM_BITS;
    }

    fn pbrt_active_transform_start_time(&mut self) {
        let last = self.transform_bits.len() - 1;
        self.transform_bits[last] = START_TRANSFORM_BITS;
    }

    fn pbrt_transform_times(&mut self, start: Float, end: Float) {
        self.verify_options("TransformTimes");

        let mut opts = self.render_options.borrow_mut();
        opts.transform_start_time = start;
        opts.transform_end_time = end;
    }

    fn pbrt_pixel_filter(&mut self, name: &str, params: &ParamSet) {
        self.verify_options("PixelFilter");

        let mut opts = self.render_options.borrow_mut();
        opts.filter_name = String::from(name);
        opts.filter_params.set(params);
    }

    fn pbrt_film(&mut self, name: &str, params: &ParamSet) {
        self.verify_options("Film");
        let mut opts = self.render_options.borrow_mut();
        opts.film_name = String::from(name);
        opts.film_params.set(params);
    }

    fn pbrt_sampler(&mut self, name: &str, params: &ParamSet) {
        self.verify_options("Sampler");

        let mut opts = self.render_options.borrow_mut();
        opts.sampler_name = String::from(name);
        opts.sampler_params.set(params);
    }

    fn pbrt_accelerator(&mut self, name: &str, params: &ParamSet) {
        self.verify_options("Accelerator");

        let mut opts = self.render_options.borrow_mut();
        opts.accelerator_name = String::from(name);
        opts.accelerator_params.set(params);
    }

    fn pbrt_integrator(&mut self, name: &str, params: &ParamSet) {
        self.verify_options("Integrator");

        let mut opts = self.render_options.borrow_mut();
        opts.integrator_name = String::from(name);
        opts.integrator_params.set(params);
    }

    fn pbrt_camera(&mut self, name: &str, params: &ParamSet) {
        self.verify_options("Camera");

        let params = self.make_absolute_path(params);

        let trns = self.transforms[self.transforms.len() - 1].borrow();
        let mut opts = self.render_options.borrow_mut();
        opts.camera_name = String::from(name);
        opts.camera_params.set(&params);
        //Camera to World
        self.named_coordinate_systems.insert(
            String::from("camera"),
            RefCell::new(TransformSet::inverse(&trns)),
        );
    }

    fn pbrt_make_named_medium(&mut self, name: &str, params: &ParamSet) {
        self.verify_initialized("MakeNamedMedium");

        let t = params.find_one_string("type", "");
        if t == "" {
            error!("No parameter string \"type\" found in MakeNamedMedium");
        } else {
            let trns = self.transforms[self.transforms.len() - 1].borrow();
            if let Some(m) = self.make_medium(&t, params, &trns[0]) {
                let mut opts = self.render_options.borrow_mut();
                opts.named_media.insert(String::from(name), m);
            } else {
                error!("Unable to create medium \"{}\" \"{}\"", name, t);
            }
        }
    }

    fn pbrt_medium_interface(&mut self, inside_name: &str, outside_name: &str) {
        self.verify_initialized("MediumInterface");

        let mut atrs = self.graphics_states[self.graphics_states.len() - 1].borrow_mut();
        let mut opts = self.render_options.borrow_mut();
        atrs.current_inside_medium = String::from(inside_name);
        atrs.current_outside_medium = String::from(outside_name);
        opts.have_scattering_media = true;
    }

    fn pbrt_world_begin(&mut self) {
        self.verify_options("WorldBegin");

        self.current_api_state = APIState::OptionsBlock;
        let mut trns = self.transforms[self.transforms.len() - 1].borrow_mut();
        trns.set_transform(&Transform::identity(), ALL_TRANSFORM_BITS);
        let last = self.transform_bits.len() - 1;
        self.transform_bits[last] = ALL_TRANSFORM_BITS;
        self.named_coordinate_systems
            .insert(String::from("world"), RefCell::new(trns.clone()));
    }

    fn pbrt_attribute_begin(&mut self) {
        self.verify_world("AttributeBegin");

        self.push_graphics_state();
        self.push_transform();
        let mut attr = self.graphics_states[self.graphics_states.len() - 1].borrow_mut();
        attr.float_textures_shared = true;
        attr.spectrum_textures_shared = true;
        attr.named_materials_shared = true;
    }

    fn pbrt_attribute_end(&mut self) {
        self.verify_world("AttributeEnd");

        self.pop_transform();
        self.pop_graphics_state();
    }

    fn pbrt_transform_begin(&mut self) {
        self.verify_world("TransformBegin");

        self.push_transform();
    }

    fn pbrt_transform_end(&mut self) {
        self.verify_world("Transformend");

        self.pop_transform();
    }

    fn pbrt_texture(&mut self, name: &str, t: &str, tex_name: &str, params: &ParamSet) {
        self.verify_world("Texture");

        let params = self.make_absolute_path(params);

        // Check for filename errors
        if let Some(maps) = params.get_strings_ref("filename") {
            if maps.len() == 1 {
                let path = Path::new(&maps[0]);
                if !path.exists() {
                    warn!("Texture file \"{}\" does not exist.", path.display());
                    return;
                }
            }
        }

        let attr = self.graphics_states[self.graphics_states.len() - 1].borrow_mut();

        let mut f_tex = attr.float_textures.as_ref().borrow_mut();
        let mut s_tex = attr.spectrum_textures.as_ref().borrow_mut();
        let tp = TextureParams::new(&params, &params, &f_tex, &s_tex);

        if t == "float" {
            let cur_trans = self.transforms[self.transforms.len() - 1].borrow();
            if let Some(ft) = self.make_float_texture(tex_name, &cur_trans[0], &tp) {
                //warn
                f_tex.insert(String::from(name), ft);
            }
        } else if t == "color" || t == "rgb" || t == "spectrum" {
            let cur_trans = self.transforms[self.transforms.len() - 1].borrow();
            if let Some(st) = self.make_spectrum_texture(tex_name, &cur_trans[0], &tp) {
                //warn
                s_tex.insert(String::from(name), st);
            }
        } else {
            //error
        }
    }

    fn pbrt_material(&mut self, name: &str, params: &ParamSet) {
        self.verify_world("Material");

        if name == "" {
            let mut attr = self.graphics_states.last().unwrap().borrow_mut();
            attr.current_material = None;
        } else {
            let params = self.make_absolute_path(params);
            let empty_params = ParamSet::new();
            if let Some(mtl) = self.make_material_params(name, &params, &empty_params) {
                let mut attr = self.graphics_states.last().unwrap().borrow_mut();
                attr.current_material = Some(Arc::new(RefCell::new(MaterialInstance::new(
                    name, mtl, &params,
                ))));
            }
        }
    }

    fn pbrt_make_named_material(&mut self, name: &str, params: &ParamSet) {
        //self.verify_world("MakeNamedMaterial");

        let params = self.make_absolute_path(params);
        let empty_params = ParamSet::new();
        let mat_name = params.find_one_string("type", "");
        if mat_name.is_empty() {
            error!("No parameter string \"type\" found in MakeNamedMaterial");
        }
        if let Some(mat) = self.make_material_params(&mat_name, &params, &empty_params) {
            let mut attr = self.graphics_states.last().unwrap().borrow_mut();
            {
                let named_materials = attr.named_materials.as_ref().borrow();
                if named_materials.contains_key(name) {
                    warn!("Named material \"{}\" redefined.", name);
                }
            }

            if attr.named_materials_shared {
                let materials = attr.clone_named_materials();
                attr.named_materials = materials;
                attr.named_materials_shared = false;
            }
            {
                let mut named_materials = attr.named_materials.as_ref().borrow_mut();
                named_materials.insert(
                    String::from(name),
                    Arc::new(RefCell::new(MaterialInstance::new(&mat_name, mat, &params))),
                );
            }
        }
    }

    fn pbrt_named_material(&mut self, name: &str) {
        self.verify_world("NamedMaterial");

        if let Some(mtl) = self.get_named_material(name) {
            let mut attr = self.graphics_states[self.graphics_states.len() - 1].borrow_mut();
            attr.current_material = Some(mtl);
        } else {
            error!("NamedMaterial \"{}\" unknown.", name);
        }
    }

    fn pbrt_light_source(&mut self, name: &str, params: &ParamSet) {
        self.verify_world("LightSource");

        let params = self.make_absolute_path(params);

        let cur_trans = self.transforms[self.transforms.len() - 1].borrow()[0];
        let mi = self.create_medium_interface();
        if let Some(lt) = self.make_light(name, &cur_trans, &mi, &params) {
            let mut opts = self.render_options.borrow_mut();
            opts.lights.push(lt);
        }
    }

    fn pbrt_area_light_source(&mut self, name: &str, params: &ParamSet) {
        self.verify_world("AreaLightSource");

        let params = self.make_absolute_path(params);

        let mut attr = self.graphics_states[self.graphics_states.len() - 1].borrow_mut();
        attr.area_light_name = String::from(name);
        attr.area_light_params = params;
    }

    fn pbrt_shape(&mut self, name: &str, params: &ParamSet) {
        self.verify_world("Shape");

        let mut prims: Vec<Arc<dyn Primitive>> = Vec::new();
        let mut area_lights: Vec<Arc<dyn Light>> = Vec::new();

        let attr = self.graphics_states[self.graphics_states.len() - 1].borrow();
        //
        let cur_trans = self.transforms[self.transforms.len() - 1].borrow();
        if !cur_trans.is_animated() {
            let object2world = cur_trans[0].clone();
            let shapes = self.make_shapes(name, &object2world, attr.reverse_orientation, params);
            if !shapes.is_empty() {
                prims.reserve(shapes.len());
                let mtl = self.get_material_for_shape(params);
                let mi = self.create_medium_interface();

                for s in shapes {
                    let area = self.make_area_light(
                        &attr.area_light_name,
                        &object2world,
                        &mi,
                        &attr.area_light_params,
                        &s,
                    );
                    if let Some(a) = area.as_ref() {
                        area_lights.push(a.clone());
                    }
                    let prim = Arc::new(GeometricPrimitive::new(&s, &mtl, &area, &mi));
                    prims.push(prim);
                }
            }
        } else {
            // Initialize _prims_ and _areaLights_ for animated shape

            // Create initial shape or shapes for animated shape
            if attr.area_light_name != "" {
                warn!("Ignoring currently set area light when creating animated shape");
            }

            let identity = Transform::identity();
            let shapes = self.make_shapes(name, &identity, attr.reverse_orientation, params);
            if shapes.is_empty() {
                return;
            }

            // Create _GeometricPrimitive_(s) for animated shape
            let mtl = self.get_material_for_shape(params);
            let mi = self.create_medium_interface();
            prims.reserve(shapes.len());
            for s in shapes {
                let prim = Arc::new(GeometricPrimitive::new(&s, &mtl, &None, &mi));
                prims.push(prim);
            }
            // Create single _TransformedPrimitive_ for _prims_

            // Get _animatedObjectToWorld_ transform for shape
            let cur_trans = self.transforms[self.transforms.len() - 1].borrow();
            if cur_trans.len() < 2 {
                error!("Transformation matrix has missing components");
                return;
            }
            let obj_to_world = [cur_trans[0].clone(), cur_trans[1].clone()];
            let animated_transform = {
                let opts = self.render_options.borrow();
                //println!("start_time:{}, end_time:{}", opts.transform_start_time, opts.transform_end_time);
                AnimatedTransform::new(
                    &obj_to_world[0],
                    opts.transform_start_time,
                    &obj_to_world[1],
                    opts.transform_end_time,
                )
            };

            if prims.len() > 1 {
                let opts = self.render_options.borrow();
                let accelerator_name = opts.accelerator_name.clone();
                let accelerator_params = opts.accelerator_params.clone();
                if let Some(accel) =
                    self.make_accelerator(&accelerator_name, &prims, &accelerator_params)
                {
                    prims.clear();
                    prims.push(accel);
                }
            }

            if !prims.is_empty() {
                assert!(prims.len() == 1);
                let prim = Arc::new(TransformedPrimitive::new(&prims[0], &animated_transform));
                prims.clear();
                prims.push(prim);
            }
        }

        {
            if let Some(instance_name) = self.get_current_instance_name() {
                if !area_lights.is_empty() {
                    warn!("Area lights not supported with object instancing"); //ignored
                }
                let mut opts = self.render_options.borrow_mut();
                if let Some(out_primitives) = opts.instances.get_mut(&instance_name) {
                    for prim in prims {
                        out_primitives.push(prim);
                    }
                }
            } else {
                let mut opts = self.render_options.borrow_mut();
                for prim in prims {
                    opts.primitives.push(prim);
                }
                if !area_lights.is_empty() {
                    for l in area_lights {
                        opts.lights.push(l);
                    }
                }
            }
        }
    }

    fn pbrt_reverse_orientation(&mut self) {
        self.verify_world("ReverseOrientation");

        let mut atrs = self.graphics_states[self.graphics_states.len() - 1].borrow_mut();
        atrs.reverse_orientation = !atrs.reverse_orientation;
    }

    fn pbrt_object_begin(&mut self, name: &str) {
        self.verify_world("ObjectBegin");

        self.pbrt_attribute_begin();
        {
            let mut opts = self.render_options.borrow_mut();
            opts.instances.insert(String::from(name), Vec::new());
            opts.current_instance_name = Some(String::from(name));
        }
    }

    fn pbrt_object_end(&mut self) {
        self.verify_world("ObjectEnd");

        {
            let mut opts = self.render_options.borrow_mut();
            opts.current_instance_name = None;
        }
        self.pbrt_attribute_end();
        N_OBJECT_INSTANCES_CREATED.with(|c| c.inc());
    }

    fn pbrt_object_instance(&mut self, name: &str) {
        self.verify_world("ObjectInstance");

        {
            let opts = self.render_options.borrow();
            if opts.current_instance_name.as_ref().is_some() {
                return;
            }
        }

        {
            let mut opts = self.render_options.borrow_mut();
            let accelerator_name = opts.accelerator_name.clone();
            let accelerator_params = opts.accelerator_params.clone();
            let start_time = opts.transform_start_time;
            let end_time = opts.transform_end_time;
            let instances = &mut opts.instances;
            if let Some(prims) = instances.get_mut(name) {
                N_OBJECT_INSTANCES_USED.with(|c| c.inc());

                if prims.len() > 1 {
                    if let Some(accel) =
                        self.make_accelerator(&accelerator_name, prims, &accelerator_params)
                    {
                        prims.clear();
                        prims.push(accel);
                    }
                }
                if !prims.is_empty() {
                    assert_eq!(prims.len(), 1);
                    let cur_trans = self.transforms[self.transforms.len() - 1].borrow();
                    let animated_transform = AnimatedTransform::new(
                        &cur_trans.t[0],
                        start_time,
                        &cur_trans.t[1],
                        end_time,
                    );
                    let prim = Arc::new(TransformedPrimitive::new(&prims[0], &animated_transform));
                    opts.primitives.push(prim);
                }
            }
        }
    }

    fn pbrt_world_end(&mut self) {
        self.verify_world("WorldEnd");
        //Do not anything!
    }

    fn pbrt_parse_file(&mut self, _file_name: &str) {
        //Do not anything!
    }
    fn pbrt_parse_string(&mut self, _s: &str) {
        //Do not anything!
    }

    fn pbrt_work_dir_begin(&mut self, path: &str) {
        self.work_dirs.push(String::from(path));
    }

    fn pbrt_work_dir_end(&mut self) {
        self.work_dirs.pop();
    }
}
