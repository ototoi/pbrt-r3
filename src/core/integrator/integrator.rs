use crate::core::camera::*;
use crate::core::scene::*;

use std::sync::Arc;

pub trait Integrator {
    fn render(&mut self, scene: &Scene);
    fn get_camera(&self) -> Arc<dyn Camera>;
}
