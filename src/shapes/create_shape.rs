use super::cone::create_cone_shape;
use super::curve::create_curve_shape;
use super::cylinder::create_cylinder_shape;
use super::disk::create_disk_shape;
use super::heightfield::create_heightfield;
use super::hyperboloid::create_hyperboloid_shape;
use super::loopsubdiv::create_loop_subdiv;
use super::nurbs::create_nurbs;
use super::paraboloid::create_paraboloid_shape;
use super::plymesh::create_ply_mesh;
use super::sphere::create_sphere_shape;
use super::triangle::create_triangle_mesh_shape;
use crate::core::prelude::*;

use std::collections::HashMap;
use std::sync::Arc;

type FloatTextureMap = HashMap<String, Arc<dyn Texture<Float>>>;

pub fn create_shapes(
    name: &str,
    object2world: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
    float_textures: &FloatTextureMap,
) -> Result<Vec<Arc<dyn Shape>>, PbrtError> {
    //println!("{} shape create", name);
    let mut shapes: Vec<Arc<dyn Shape>> = Vec::new();
    match name {
        "sphere" => {
            let s = create_sphere_shape(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            )?;
            shapes.push(Arc::new(s));
        }
        "cylinder" => {
            let s = create_cylinder_shape(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            )?;
            shapes.push(Arc::new(s));
        }
        "disk" => {
            let s = create_disk_shape(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            )?;
            shapes.push(Arc::new(s));
        }
        "cone" => {
            let s = create_cone_shape(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            )?;
            shapes.push(Arc::new(s));
        }
        "paraboloid" => {
            let s = create_paraboloid_shape(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            )?;
            shapes.push(Arc::new(s));
        }
        "hyperboloid" => {
            let s = create_hyperboloid_shape(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            )?;
            shapes.push(Arc::new(s));
        }
        "curve" => {
            return create_curve_shape(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            );
        }
        "trianglemesh" => {
            return create_triangle_mesh_shape(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
                float_textures,
            );
        }
        "plymesh" => {
            return create_ply_mesh(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
                float_textures,
            );
        }
        "heightfield" => {
            return create_heightfield(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            );
        }
        "loopsubdiv" => {
            return create_loop_subdiv(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            );
        }
        "nurbs" => {
            return create_nurbs(
                object2world,
                &object2world.inverse(),
                reverse_orientation,
                params,
            );
        }
        s => {
            return Err(PbrtError::warning(&format!("{} shape cannot create", s)));
        }
    }
    return Ok(shapes);
}
