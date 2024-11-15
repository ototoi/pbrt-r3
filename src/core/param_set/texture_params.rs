use super::param_set::ParamSet;
use crate::core::pbrt::*;
use log::*;
use std::collections::HashMap;
use std::sync::Arc;

type FloatTextureMap = HashMap<String, Arc<dyn Texture<Float>>>;
type SpectrumTextureMap = HashMap<String, Arc<dyn Texture<Spectrum>>>;

pub struct TextureParams<'a> {
    pub geom_params: &'a ParamSet,
    pub mat_params: &'a ParamSet,
    pub f_tex: &'a FloatTextureMap,
    pub s_tex: &'a SpectrumTextureMap,
}

impl<'a> TextureParams<'a> {
    pub fn new(
        geom_params: &'a ParamSet,
        mat_params: &'a ParamSet,
        f_tex: &'a FloatTextureMap,
        s_tex: &'a SpectrumTextureMap,
    ) -> Self {
        TextureParams::<'a> {
            geom_params,
            mat_params,
            f_tex,
            s_tex,
        }
    }

    fn get_float(&self, key: &str) -> Option<f32> {
        if let Some(c) = self.mat_params.get_floats_ref(key) {
            if c.len() >= 1 {
                if c.len() > 1 {
                    warn!("More than one texture present in parameter file for key \"{}\". Using first.", key);
                }
                return Some(c[0]);
            }
        }
        if let Some(c) = self.geom_params.get_floats_ref(key) {
            if c.len() >= 1 {
                if c.len() > 1 {
                    warn!("More than one texture present in parameter file for key \"{}\". Using first.", key);
                }
                return Some(c[0]);
            }
        }
        return None;
    }

    fn get_spectrum_from(params: &ParamSet, key: &str) -> Option<Spectrum> {
        if let Some(c) = params.get_points_ref(key) {
            if c.len() >= 3 {
                return Some(Spectrum::from_rgb(
                    &[c[0], c[1], c[2]],
                    SpectrumType::Reflectance,
                ));
            }
        } else if let Some(c) = params.get_spectrums_ref(key) {
            if c.len() >= 1 {
                if c.len() > 1 {
                    warn!("More than one texture present in parameter file for key \"{}\". Using first.", key);
                }
                return Some(c[0]);
            }
        }
        return None;
    }

    fn get_spectrum(&self, key: &str) -> Option<Spectrum> {
        if let Some(c) = Self::get_spectrum_from(self.mat_params, key) {
            return Some(c);
        }
        if let Some(c) = Self::get_spectrum_from(self.geom_params, key) {
            return Some(c);
        }
        return None;
    }

    fn get_texture_from(params: &ParamSet, key: &str) -> Option<String> {
        if let Some(s) = params.get_textures_ref(key) {
            if s.len() >= 1 {
                if s.len() > 1 {
                    warn!("More than one texture present in parameter file for key \"{}\". Using first.", key);
                }
                return Some(s[0].clone());
            }
        }
        return None;
    }

    fn get_texture(&self, key: &str) -> Option<String> {
        if let Some(c) = Self::get_texture_from(self.geom_params, key) {
            return Some(c);
        }
        if let Some(c) = Self::get_texture_from(self.mat_params, key) {
            return Some(c);
        }
        return None;
    }

    pub fn get_float_texture_or_null(&self, key: &str) -> Option<Arc<dyn Texture<Float>>> {
        if let Some(name) = self.get_texture(key) {
            if let Some(tex) = self.f_tex.get(&name) {
                return Some(Arc::clone(tex));
            }
        } else if let Some(s) = self.get_float(key) {
            return Some(Arc::new(ConstantTexture::new(&s)));
        }
        return None;
    }

    pub fn get_float_texture(&self, key: &str, value: Float) -> Arc<dyn Texture<Float>> {
        if let Some(tex) = self.get_float_texture_or_null(key) {
            return tex;
        } else {
            return Arc::new(ConstantTexture::new(&value));
        }
    }

    pub fn get_spectrum_texture_or_null(&self, key: &str) -> Option<Arc<dyn Texture<Spectrum>>> {
        if let Some(name) = self.get_texture(key) {
            //println!("get_spectrum_texture:key:{}", name);
            if let Some(tex) = self.s_tex.get(&name) {
                return Some(Arc::clone(tex));
            }
            if let Some(s) = self.get_spectrum(&name) {
                return Some(Arc::new(ConstantTexture::new(&s)));
            }
        } else if let Some(s) = self.get_spectrum(key) {
            return Some(Arc::new(ConstantTexture::new(&s)));
        }
        return None;
    }

    pub fn get_spectrum_texture(&self, key: &str, value: &Spectrum) -> Arc<dyn Texture<Spectrum>> {
        if let Some(tex) = self.get_spectrum_texture_or_null(key) {
            return tex;
        } else {
            return Arc::new(ConstantTexture::new(value));
        }
    }

    pub fn find_float(&self, key: &str, value: Float) -> Float {
        return self
            .geom_params
            .find_one_float(key, self.mat_params.find_one_float(key, value));
    }

    pub fn find_int(&self, key: &str, value: i32) -> i32 {
        return self
            .geom_params
            .find_one_int(key, self.mat_params.find_one_int(key, value));
    }

    pub fn find_bool(&self, key: &str, value: bool) -> bool {
        return self
            .geom_params
            .find_one_bool(key, self.mat_params.find_one_bool(key, value));
    }

    pub fn find_vector3f(&self, key: &str, value: &Vector3f) -> Vector3f {
        return self
            .geom_params
            .find_one_vector3f(key, &self.mat_params.find_one_vector3f(key, value));
    }

    pub fn find_spectrum(&self, key: &str, value: &Spectrum) -> Spectrum {
        return self
            .geom_params
            .find_one_spectrum(key, &self.mat_params.find_one_spectrum(key, value));
    }

    pub fn find_string(&self, key: &str, value: &str) -> String {
        return self
            .geom_params
            .find_one_string(key, &self.mat_params.find_one_string(key, value));
    }

    pub fn find_filename(&self, key: &str, value: &str) -> String {
        return self.find_string(key, value);
    }
}
