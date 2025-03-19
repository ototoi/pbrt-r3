use crate::core::material::*;
use crate::core::param_set::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use crate::core::texture::*;
use crate::materials::*;

use std::cell::RefCell;
use std::collections::HashMap;
use std::ops::Deref;
use std::sync::Arc;

/*
std::string name;
    std::shared_ptr<Material> material;
    ParamSet params;
*/
pub struct MaterialInstance {
    pub name: String,
    pub material: Arc<dyn Material>,
    pub params: ParamSet,
}

impl MaterialInstance {
    pub fn new(name: &str, material: Arc<dyn Material>, params: &ParamSet) -> Self {
        MaterialInstance {
            name: String::from(name),
            material,
            params: params.clone(),
        }
    }
}

type FloatTextureMap = HashMap<String, Arc<dyn Texture<Float>>>;
type SpectrumTextureMap = HashMap<String, Arc<dyn Texture<Spectrum>>>;
type NamedMaterialMap = HashMap<String, Arc<RefCell<MaterialInstance>>>;

pub struct GraphicsState {
    pub current_inside_medium: String,
    pub current_outside_medium: String,
    pub float_textures: Arc<RefCell<FloatTextureMap>>,
    pub float_textures_shared: bool,
    pub spectrum_textures: Arc<RefCell<SpectrumTextureMap>>,
    pub spectrum_textures_shared: bool,
    pub named_materials: Arc<RefCell<NamedMaterialMap>>,
    pub named_materials_shared: bool,
    pub current_material: Option<Arc<RefCell<MaterialInstance>>>,

    pub area_light_params: ParamSet,
    pub area_light_name: String,
    pub reverse_orientation: bool,
}

fn create_default_material(
    f_tex: &Arc<RefCell<FloatTextureMap>>,
    s_tex: &Arc<RefCell<SpectrumTextureMap>>,
) -> Arc<dyn Material> {
    let params = ParamSet::new();
    let f = f_tex.borrow();
    let s = s_tex.borrow();
    let mp = TextureParams::new(&params, &params, f.deref(), s.deref());
    let mat = create_matte_material(&mp).unwrap();
    return mat;
}

impl GraphicsState {
    pub fn new() -> Self {
        let f_tex = Arc::new(RefCell::new(FloatTextureMap::new()));
        let s_tex = Arc::new(RefCell::new(SpectrumTextureMap::new()));
        let mat = create_default_material(&f_tex, &s_tex);

        GraphicsState {
            current_inside_medium: String::from(""),
            current_outside_medium: String::from(""),
            float_textures: f_tex,
            float_textures_shared: false,
            spectrum_textures: s_tex,
            spectrum_textures_shared: false,
            named_materials: Arc::new(RefCell::new(NamedMaterialMap::new())),
            named_materials_shared: false,
            current_material: Some(Arc::new(RefCell::new(MaterialInstance::new(
                "matte",
                mat,
                &ParamSet::new(),
            )))),

            area_light_params: ParamSet::new(),
            area_light_name: String::from(""),
            reverse_orientation: false,
        }
    }

    pub fn clone_named_materials(&self) -> Arc<RefCell<NamedMaterialMap>> {
        let materials = self.named_materials.borrow();
        return Arc::new(RefCell::new(materials.clone()));
    }
}

impl Clone for GraphicsState {
    fn clone(&self) -> Self {
        GraphicsState {
            current_inside_medium: self.current_inside_medium.clone(),
            current_outside_medium: self.current_outside_medium.clone(),
            float_textures: Arc::clone(&self.float_textures),
            float_textures_shared: self.float_textures_shared,
            spectrum_textures: Arc::clone(&self.spectrum_textures),
            spectrum_textures_shared: self.spectrum_textures_shared,
            named_materials: Arc::clone(&self.named_materials),
            named_materials_shared: self.named_materials_shared,
            current_material: self.current_material.clone(),
            area_light_params: self.area_light_params.clone(),
            area_light_name: self.area_light_name.clone(),
            reverse_orientation: self.reverse_orientation,
        }
    }
}
