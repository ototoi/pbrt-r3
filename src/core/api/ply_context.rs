use super::parse_context::*;
use crate::core::pbrt::*;
use std::cell::RefCell;
use std::path::Path;
use std::sync::Arc;

use log::*;

use ply_rs::ply::*;
use ply_rs::writer::Writer;

pub struct PlyContext {
    dir: String,
    context: Arc<RefCell<dyn ParseContext>>,
    index: usize,
}

impl PlyContext {
    pub fn new(dir: &str, context: Arc<RefCell<dyn ParseContext>>) -> Self {
        PlyContext {
            dir: dir.to_string(),
            context,
            index: 0,
        }
    }
}

#[derive(Debug, Default)]
struct Mesh {
    pub p: Vec<Point3f>,
    //pub s: Vec<Vector3f>,
    pub n: Vec<Vector3f>,
    pub uv: Vec<Point2f>,
    pub vertex_indices: Vec<u32>,
    pub face_indices: Vec<u32>,
}

fn convert_to_mesh(params: &ParamSet) -> Option<Mesh> {
    let mut vertex_indices = Vec::new();
    let mut p: Vec<Vector3f> = Vec::new();
    let mut s: Vec<Vector3f> = Vec::new();
    let mut n: Vec<Vector3f> = Vec::new();
    let mut uv: Vec<Vector2f> = Vec::new();
    let mut face_indices = Vec::new();

    if let Some(vi) = params.get_ints_ref("indices") {
        vertex_indices.resize(vi.len(), 0);
        for i in 0..vi.len() {
            vertex_indices[i] = vi[i] as u32;
        }
    }

    if let Some(ps) = params.get_points_ref("P") {
        let sz = ps.len() / 3;
        p.resize(sz, Point3f::zero());
        for i in 0..sz {
            p[i] = Point3f::new(ps[3 * i + 0], ps[3 * i + 1], ps[3 * i + 2]);
        }
    }

    if let Some(ps) = params.get_floats_ref("uv") {
        let sz = ps.len() / 2;
        uv.resize(sz, Vector2::zero());
        for i in 0..sz {
            uv[i] = Vector2::new(ps[2 * i + 0], ps[2 * i + 1]);
        }
    } else if let Some(ps) = params.get_floats_ref("st") {
        let sz = ps.len() / 2;
        uv.resize(sz, Vector2::zero());
        for i in 0..sz {
            uv[i] = Vector2::new(ps[2 * i + 0], ps[2 * i + 1]);
        }
    }

    if let Some(ps) = params.get_points_ref("S") {
        let sz = ps.len() / 3;
        s.resize(sz, Vector3::zero());
        for i in 0..sz {
            s[i] = Vector3f::new(ps[3 * i + 0], ps[3 * i + 1], ps[3 * i + 2]);
        }
    }

    if let Some(ps) = params.get_points_ref("N") {
        let sz = ps.len() / 3;
        n.resize(sz, Normal3f::zero());
        for i in 0..sz {
            n[i] = Normal3f::new(ps[3 * i + 0], ps[3 * i + 1], ps[3 * i + 2]);
        }
    }

    if let Some(vi) = params.get_ints_ref("faceIndices") {
        face_indices.resize(vi.len(), 0);
        for i in 0..vi.len() {
            face_indices[i] = vi[i] as u32;
        }
    }

    if !p.is_empty() && !vertex_indices.is_empty() {
        let mesh = Mesh {
            p,
            s,
            n,
            uv,
            vertex_indices,
            face_indices,
        };
        return Some(mesh);
    } else {
        return None;
    }
}

/*
example:
ply
format binary_little_endian 1.0
element vertex 33
property float x
property float y
property float z
property float nx
property float ny
property float nz
element face 56
property list uint8 int vertex_indices
*/
fn write_mesh_to_ply(mesh: &Mesh, file_name: &Path) -> Result<(), PbrtError> {
    let p = &mesh.p;
    let n = &mesh.n;
    let uv = &mesh.uv;

    let vertex_indices = &mesh.vertex_indices;
    let face_indices = &mesh.face_indices;

    let vertex_count = p.len();
    let face_count = vertex_indices.len() / 3;
    let mut header = Header::new();
    header.encoding = Encoding::BinaryLittleEndian;
    {
        let mut element = ElementDef::new("vertex".to_string());
        element.count = vertex_count;
        element.properties.insert(
            "x".to_string(),
            PropertyDef::new("x".to_string(), PropertyType::Scalar(ScalarType::Float)),
        );
        element.properties.insert(
            "y".to_string(),
            PropertyDef::new("y".to_string(), PropertyType::Scalar(ScalarType::Float)),
        );
        element.properties.insert(
            "z".to_string(),
            PropertyDef::new("z".to_string(), PropertyType::Scalar(ScalarType::Float)),
        );
        if n.len() > 0 {
            element.properties.insert(
                "nx".to_string(),
                PropertyDef::new("nx".to_string(), PropertyType::Scalar(ScalarType::Float)),
            );
            element.properties.insert(
                "ny".to_string(),
                PropertyDef::new("ny".to_string(), PropertyType::Scalar(ScalarType::Float)),
            );
            element.properties.insert(
                "nz".to_string(),
                PropertyDef::new("nz".to_string(), PropertyType::Scalar(ScalarType::Float)),
            );
        }
        if uv.len() > 0 {
            element.properties.insert(
                "u".to_string(),
                PropertyDef::new("u".to_string(), PropertyType::Scalar(ScalarType::Float)),
            );
            element.properties.insert(
                "v".to_string(),
                PropertyDef::new("v".to_string(), PropertyType::Scalar(ScalarType::Float)),
            );
        }
        header.elements.insert("vertex".to_string(), element);
    }
    {
        let mut element = ElementDef::new("face".to_string());
        element.count = face_count;
        element.properties.insert(
            "vertex_indices".to_string(),
            PropertyDef::new(
                "vertex_indices".to_string(),
                PropertyType::List(ScalarType::UChar, ScalarType::Int),
            ),
        );
        if face_indices.len() > 0 {
            element.properties.insert(
                "face_indices".to_string(),
                PropertyDef::new(
                    "face_indices".to_string(),
                    PropertyType::Scalar(ScalarType::Int),
                ),
            );
        }
        header.elements.insert("face".to_string(), element);
    }

    let mut ply = Ply::<DefaultElement>::new();
    ply.header = header;
    {
        let mut vertices = Vec::new();
        for i in 0..vertex_count {
            let mut vertex = DefaultElement::new();
            vertex.set_property("x".to_string(), Property::Float(p[i].x));
            vertex.set_property("y".to_string(), Property::Float(p[i].y));
            vertex.set_property("z".to_string(), Property::Float(p[i].z));
            if n.len() > 0 {
                vertex.set_property("nx".to_string(), Property::Float(n[i].x));
                vertex.set_property("ny".to_string(), Property::Float(n[i].y));
                vertex.set_property("nz".to_string(), Property::Float(n[i].z));
            }
            if uv.len() > 0 {
                vertex.set_property("u".to_string(), Property::Float(uv[i].x));
                vertex.set_property("v".to_string(), Property::Float(uv[i].y));
            }
            vertices.push(vertex);
        }
        ply.payload.insert("vertex".to_string(), vertices);
    }
    {
        let mut faces = Vec::new();
        for i in 0..face_count {
            let mut face = DefaultElement::new();
            face.set_property(
                "vertex_indices".to_string(),
                Property::ListInt(vec![
                    vertex_indices[3 * i + 0] as i32,
                    vertex_indices[3 * i + 1] as i32,
                    vertex_indices[3 * i + 2] as i32,
                ]),
            );
            if face_indices.len() > 0 {
                face.set_property(
                    "face_indices".to_string(),
                    Property::Int(face_indices[i] as i32),
                );
            }
            faces.push(face);
        }
        ply.payload.insert("face".to_string(), faces);
    }

    let mut buf = std::fs::File::create(file_name)?;
    let writer = Writer::new();
    writer.write_ply(&mut buf, &mut ply)?;
    Ok(())
}

fn create_plymesh_params(params: &ParamSet) -> ParamSet {
    let mut p = ParamSet::new();
    let keys = params.get_keys();
    for key in keys {
        let keyname = params.get_key_name(&key);
        match keyname.as_str() {
            "indices" | "P" | "uv" | "st" | "S" | "N" | "faceIndices" => {
                continue;
            }
            _ => {
                let key = key.as_str();
                if let Some(v) = params.get_bools_ref(key) {
                    p.add_bools(key, &v);
                } else if let Some(v) = params.get_ints_ref(key) {
                    p.add_ints(key, &v);
                } else if let Some(v) = params.get_floats_ref(key) {
                    p.add_floats(key, &v);
                } else if let Some(v) = params.get_points_ref(key) {
                    p.add_floats(key, &v);
                } else if let Some(v) = params.get_strings_ref(key) {
                    let v = v.iter().map(|x| x.as_str()).collect::<Vec<&str>>();
                    p.add_strings(key, &v);
                } else if let Some(v) = params.get_textures_ref(key) {
                    let v = v.iter().map(|x| x.as_str()).collect::<Vec<&str>>();
                    p.add_strings(key, &v);
                } else {
                    warn!("Unsupported type for key: {}", key);
                }
            }
        }
    }
    return p;
}

impl ParseContext for PlyContext {
    fn pbrt_cleanup(&mut self) {
        self.context.borrow_mut().pbrt_cleanup();
    }

    fn pbrt_identity(&mut self) {
        self.context.borrow_mut().pbrt_identity();
    }

    fn pbrt_translate(&mut self, dx: Float, dy: Float, dz: Float) {
        self.context.borrow_mut().pbrt_translate(dx, dy, dz);
    }

    fn pbrt_rotate(&mut self, angle: Float, ax: Float, ay: Float, az: Float) {
        self.context.borrow_mut().pbrt_rotate(angle, ax, ay, az);
    }

    fn pbrt_scale(&mut self, sx: Float, sy: Float, sz: Float) {
        self.context.borrow_mut().pbrt_scale(sx, sy, sz);
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
        self.context
            .borrow_mut()
            .pbrt_look_at(ex, ey, ez, lx, ly, lz, ux, uy, uz);
    }

    fn pbrt_concat_transform(&mut self, transform: &[Float]) {
        self.context.borrow_mut().pbrt_concat_transform(transform);
    }

    fn pbrt_transform(&mut self, transform: &[Float]) {
        self.context.borrow_mut().pbrt_transform(transform);
    }

    fn pbrt_coordinate_system(&mut self, name: &str) {
        self.context.borrow_mut().pbrt_coordinate_system(name);
    }

    fn pbrt_coord_sys_transform(&mut self, name: &str) {
        self.context.borrow_mut().pbrt_coord_sys_transform(name);
    }

    fn pbrt_active_transform_all(&mut self) {
        self.context.borrow_mut().pbrt_active_transform_all();
    }

    fn pbrt_active_transform_end_time(&mut self) {
        self.context.borrow_mut().pbrt_active_transform_end_time();
    }

    fn pbrt_active_transform_start_time(&mut self) {
        self.context.borrow_mut().pbrt_active_transform_start_time();
    }

    fn pbrt_transform_times(&mut self, start: Float, end: Float) {
        self.context.borrow_mut().pbrt_transform_times(start, end);
    }

    fn pbrt_pixel_filter(&mut self, name: &str, params: &ParamSet) {
        self.context.borrow_mut().pbrt_pixel_filter(name, params);
    }

    fn pbrt_film(&mut self, name: &str, params: &ParamSet) {
        self.context.borrow_mut().pbrt_film(name, params);
    }

    fn pbrt_sampler(&mut self, name: &str, params: &ParamSet) {
        self.context.borrow_mut().pbrt_sampler(name, params);
    }

    fn pbrt_accelerator(&mut self, name: &str, params: &ParamSet) {
        self.context.borrow_mut().pbrt_accelerator(name, params);
    }

    fn pbrt_integrator(&mut self, name: &str, params: &ParamSet) {
        self.context.borrow_mut().pbrt_integrator(name, params);
    }
    fn pbrt_camera(&mut self, name: &str, params: &ParamSet) {
        self.context.borrow_mut().pbrt_camera(name, params);
    }

    fn pbrt_make_named_medium(&mut self, name: &str, params: &ParamSet) {
        self.context
            .borrow_mut()
            .pbrt_make_named_medium(name, params);
    }

    fn pbrt_medium_interface(&mut self, inside_name: &str, outside_name: &str) {
        self.context
            .borrow_mut()
            .pbrt_medium_interface(inside_name, outside_name);
    }

    fn pbrt_world_begin(&mut self) {
        self.context.borrow_mut().pbrt_world_begin();
    }

    fn pbrt_attribute_begin(&mut self) {
        self.context.borrow_mut().pbrt_attribute_begin();
    }

    fn pbrt_attribute_end(&mut self) {
        self.context.borrow_mut().pbrt_attribute_end();
    }

    fn pbrt_transform_begin(&mut self) {
        self.context.borrow_mut().pbrt_transform_begin();
    }

    fn pbrt_transform_end(&mut self) {
        self.context.borrow_mut().pbrt_transform_end();
    }

    fn pbrt_texture(&mut self, name: &str, t: &str, tex_name: &str, params: &ParamSet) {
        self.context
            .borrow_mut()
            .pbrt_texture(name, t, tex_name, params);
    }

    fn pbrt_material(&mut self, name: &str, params: &ParamSet) {
        self.context.borrow_mut().pbrt_material(name, params);
    }

    fn pbrt_make_named_material(&mut self, name: &str, params: &ParamSet) {
        self.context
            .borrow_mut()
            .pbrt_make_named_material(name, params);
    }

    fn pbrt_named_material(&mut self, name: &str) {
        self.context.borrow_mut().pbrt_named_material(name);
    }

    fn pbrt_light_source(&mut self, name: &str, params: &ParamSet) {
        self.context.borrow_mut().pbrt_light_source(name, params);
    }

    fn pbrt_area_light_source(&mut self, name: &str, params: &ParamSet) {
        self.context
            .borrow_mut()
            .pbrt_area_light_source(name, params);
    }

    fn pbrt_shape(&mut self, name: &str, params: &ParamSet) {
        if name == "trianglemesh" {
            if let Some(mesh) = convert_to_mesh(params) {
                let dir = Path::new(&self.dir);
                let dir = dir.join("geometry");
                if let Err(e) = std::fs::create_dir_all(&dir) {
                    error!("Error: {}", e);
                    return;
                }
                let file_name = dir.join(format!("mesh_{}.ply", self.index));
                match write_mesh_to_ply(&mesh, &file_name) {
                    Ok(_) => {
                        let mut params = create_plymesh_params(params);
                        let filename = file_name.file_name().unwrap().to_str().unwrap();
                        let filename = Path::new("geometry").join(filename);
                        let filename = filename.to_str().unwrap();
                        params.add_string("string filename", filename);
                        self.context.borrow_mut().pbrt_shape("plymesh", &params);
                    }
                    Err(e) => {
                        error!("Error: {}", e);
                    }
                }
            }
            self.index += 1;
        } else {
            self.context.borrow_mut().pbrt_shape(name, params);
        }
    }

    fn pbrt_reverse_orientation(&mut self) {
        self.context.borrow_mut().pbrt_reverse_orientation();
    }

    fn pbrt_object_begin(&mut self, name: &str) {
        self.context.borrow_mut().pbrt_object_begin(name);
    }

    fn pbrt_object_end(&mut self) {
        self.context.borrow_mut().pbrt_object_end();
    }
    fn pbrt_object_instance(&mut self, name: &str) {
        self.context.borrow_mut().pbrt_object_instance(name);
    }

    fn pbrt_world_end(&mut self) {
        self.context.borrow_mut().pbrt_world_end();
    }

    fn pbrt_parse_file(&mut self, file_name: &str) {
        self.context.borrow_mut().pbrt_parse_file(file_name);
    }

    fn pbrt_parse_string(&mut self, s: &str) {
        self.context.borrow_mut().pbrt_parse_string(s);
    }

    fn pbrt_work_dir_begin(&mut self, path: &str) {
        self.context.borrow_mut().pbrt_work_dir_begin(path);
    }

    fn pbrt_work_dir_end(&mut self) {
        self.context.borrow_mut().pbrt_work_dir_end();
    }

    fn pbrt_include(&mut self, file_name: &str, params: &ParamSet) {
        self.context.borrow_mut().pbrt_include(file_name, params);
    }
}
