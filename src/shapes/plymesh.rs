use super::alphamask::AlphaMaskShape;
use super::triangle::*;
use crate::core::prelude::*;

use std::collections::HashMap;
use std::sync::Arc;

use ply_rs::parser;
use ply_rs::ply;

const VERTEX_P: u32 = 1;
const VERTEX_N: u32 = 2;
const VERTEX_UV: u32 = 8;

#[derive(Debug)]
struct Vertex {
    x: f32,
    y: f32,
    z: f32,
    nx: f32,
    ny: f32,
    nz: f32,
    u: f32,
    v: f32,
    flags: u32,
}

#[derive(Debug)]
struct Face {
    vertex_index: Vec<i32>,
    n: [f32; 3],
}

impl ply::PropertyAccess for Vertex {
    fn new() -> Self {
        Vertex {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            nx: 0.0,
            ny: 0.0,
            nz: 0.0,
            u: 0.0,
            v: 0.0,
            flags: 0,
        }
    }
    fn set_property(&mut self, key: String, property: ply::Property) {
        match property {
            ply::Property::Float(v) => match key.as_str() {
                "x" => {
                    self.x = v;
                    self.flags |= VERTEX_P;
                }
                "y" => {
                    self.y = v;
                    self.flags |= VERTEX_P;
                }
                "z" => {
                    self.z = v;
                    self.flags |= VERTEX_P;
                }
                "nx" => {
                    self.nx = v;
                    self.flags |= VERTEX_N;
                }
                "ny" => {
                    self.ny = v;
                    self.flags |= VERTEX_N;
                }
                "nz" => {
                    self.nz = v;
                    self.flags |= VERTEX_N;
                }
                "u" => {
                    self.u = v;
                    self.flags |= VERTEX_UV;
                }
                "v" => {
                    self.v = v;
                    self.flags |= VERTEX_UV;
                }
                "s" => {
                    self.u = v;
                    self.flags |= VERTEX_UV;
                }
                "t" => {
                    self.v = v;
                    self.flags |= VERTEX_UV;
                }
                "texture_u" => {
                    self.u = v;
                    self.flags |= VERTEX_UV;
                }
                "texture_v" => {
                    self.v = v;
                    self.flags |= VERTEX_UV;
                }
                "texture_s" => {
                    self.u = v;
                    self.flags |= VERTEX_UV;
                }
                "texture_t" => {
                    self.v = v;
                    self.flags |= VERTEX_UV;
                }
                k => panic!("Vertex: Unexpected key/value combination: key: {}", k),
            },
            _ => {
                panic!(
                    "Vertex: Unexpected key/value combination: key: {}, type: {:?}",
                    &key, property
                );
            }
        }
    }
}

// same thing for Face
impl ply::PropertyAccess for Face {
    fn new() -> Self {
        Face {
            vertex_index: Vec::new(),
            n: [0.0; 3],
        }
    }
    fn set_property(&mut self, key: String, property: ply::Property) {
        if key == "vertex_indices" || key == "vertex_index" {
            match property {
                ply::Property::ListInt(vec) => {
                    for i in 0..vec.len() {
                        self.vertex_index.push(vec[i]);
                    }
                }
                ply::Property::ListUInt(vec) => {
                    for i in 0..vec.len() {
                        self.vertex_index.push(vec[i] as i32);
                    }
                }
                _ => {
                    panic!(
                        "Face: Unexpected key/value combination: key: {}, type: {:?}",
                        &key, property
                    );
                }
            }
        } else if key == "nx" || key == "ny" || key == "nz" {
            match property {
                ply::Property::Float(v) => match key.as_str() {
                    "nx" => self.n[0] = v,
                    "ny" => self.n[1] = v,
                    "nz" => self.n[2] = v,
                    _ => {}
                },
                _ => {
                    panic!(
                        "Face: Unexpected key/value combination: key: {}, type: {:?}",
                        &key, property
                    );
                }
            }
        }
    }
}

/*
fn create_bound_mesh(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    _: Vec<u32>,
    p: Vec<Point3f>,
    _s: Vec<Vector3f>,
    _n: Vec<Vector3f>,
    _uv: Vec<Point2f>,
    params: &ParamSet,
) -> Vec<Arc<dyn Shape>> {
    //return create_triangle_mesh(o2w, w2o, reverse_orientation, vertex_indices, p, s, n, uv);
    let min = p[0]; //o2w.transform_point(&p[0]);
    let min = p.iter().fold(min, |acc, e| -> Vector3f {
        return Vector3f::new(
            Float::min(acc.x, e.x),
            Float::min(acc.y, e.y),
            Float::min(acc.z, e.z),
        );
    });
    let max = p[0]; //2w.transform_point(&p[0]);
    let max = p.iter().fold(max, |acc, e| -> Vector3f {
        return Vector3f::new(
            Float::max(acc.x, e.x),
            Float::max(acc.y, e.y),
            Float::max(acc.z, e.z),
        );
    });

    let bp = vec![
        Vector3f::new(min.x, min.y, min.z), //0
        Vector3f::new(min.x, min.y, max.z), //1
        Vector3f::new(min.x, max.y, min.z), //2
        Vector3f::new(min.x, max.y, max.z), //3
        Vector3f::new(max.x, min.y, min.z), //4
        Vector3f::new(max.x, min.y, max.z), //5
        Vector3f::new(max.x, max.y, min.z), //6
        Vector3f::new(max.x, max.y, max.z), //7
    ];

    #[rustfmt::skip]
    let bi = vec![
        0, 1, 2,
        1, 3, 2,
        4, 5, 0,
        5, 1, 0,
        6, 7, 4,
        7, 5, 4,
        2, 3, 6,
        3, 7, 6,
        1, 5, 3,
        5, 7, 3,
        4, 0, 6,
        0, 2, 6
    ];

    let bs = Vec::new();
    let bn = Vec::new();
    let buv = Vec::new();

    return create_triangle_mesh(o2w, w2o, reverse_orientation, bi, bp, bs, bn, buv, params);
}
*/

type FloatTextureMap = HashMap<String, Arc<dyn Texture<Float>>>;

pub fn create_ply_mesh(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
    float_textures: &FloatTextureMap,
) -> Result<Vec<Arc<dyn Shape>>, PbrtError> {
    let filename = params.find_one_string("filename", "");
    match std::fs::File::open(filename) {
        Ok(f) => {
            let mut reader = std::io::BufReader::new(f);
            let vertex_parser = parser::Parser::<Vertex>::new();
            let face_parser = parser::Parser::<Face>::new();
            let header = vertex_parser.read_header(&mut reader).unwrap();

            let mut p = Vec::new();
            let mut vertex_indices: Vec<u32> = Vec::new();
            //let mut face_list = Vec::new();
            let mut n = Vec::new();
            let s = Vec::new();
            let mut uv = Vec::new();
            for (_name, element) in header.elements.iter() {
                //println!("{:?}", name);
                // we could also just parse them in sequence, but the file format might change
                match element.name.as_ref() {
                    "vertex" => {
                        let r =
                            vertex_parser.read_payload_for_element(&mut reader, element, &header);
                        match r {
                            Ok(vertex_list) => {
                                if !vertex_list.is_empty() {
                                    let flags = vertex_list[0].flags;
                                    if (flags & VERTEX_P) != 0 {
                                        p.reserve(vertex_list.len());
                                        for v in vertex_list.iter() {
                                            p.push(Vector3f::new(
                                                v.x as Float,
                                                v.y as Float,
                                                v.z as Float,
                                            ));
                                        }
                                    }
                                    if (flags & VERTEX_N) != 0 {
                                        n.reserve(vertex_list.len());
                                        for v in vertex_list.iter() {
                                            n.push(Normal3f::new(
                                                v.nx as Float,
                                                v.ny as Float,
                                                v.nz as Float,
                                            ));
                                        }
                                    }
                                    if (flags & VERTEX_UV) != 0 {
                                        uv.reserve(vertex_list.len());
                                        for v in vertex_list.iter() {
                                            uv.push(Vector2f::new(v.u as Float, v.v as Float));
                                        }
                                    }
                                }
                            }
                            Err(e) => {
                                return Err(PbrtError::from(e));
                            }
                        }
                    }
                    "face" => {
                        let r = face_parser.read_payload_for_element(&mut reader, element, &header);
                        match r {
                            Ok(face_list) => {
                                vertex_indices.reserve(face_list.len() * 3);
                                for face in face_list {
                                    let n_vert = face.vertex_index.len();
                                    match n_vert {
                                        3 => {
                                            for idx in face.vertex_index {
                                                vertex_indices.push(idx as u32);
                                            }
                                        }
                                        4 => {
                                            let i0 = face.vertex_index[0] as u32;
                                            let i1 = face.vertex_index[1] as u32;
                                            let i2 = face.vertex_index[2] as u32;
                                            let i3 = face.vertex_index[3] as u32;
                                            vertex_indices.push(i0);
                                            vertex_indices.push(i1);
                                            vertex_indices.push(i2);
                                            vertex_indices.push(i3);
                                            vertex_indices.push(i0);
                                            vertex_indices.push(i2);
                                        }
                                        _ => {
                                            let msg = format!("plymesh: Ignoring face with {} vertices (only triangles and quads are supported!)", n_vert);
                                            return Err(PbrtError::error(&msg));
                                        }
                                    }
                                }
                            }
                            Err(e) => {
                                return Err(PbrtError::from(e));
                            }
                        }
                    }
                    _ => {}
                }
            }

            let mut mesh = create_triangle_mesh(
                o2w,
                w2o,
                reverse_orientation,
                vertex_indices,
                p,
                s,
                n,
                uv,
                params,
            );
            //let mesh =
            //    create_bound_mesh(o2w, w2o, reverse_orientation, vertex_indices, p, s, n, uv);
            let alpha_mask_info = get_alpha_texture(params, float_textures);
            let shadow_alpha_mask_info = get_shadow_alpha_texture(params, float_textures);
            if alpha_mask_info.is_some() || shadow_alpha_mask_info.is_some() {
                for i in 0..mesh.len() {
                    mesh[i] = Arc::new(AlphaMaskShape::new(
                        &mesh[i],
                        &alpha_mask_info,
                        &shadow_alpha_mask_info,
                    ));
                }
            }

            return Ok(mesh);
        }
        Err(e) => {
            let msg = e.to_string();
            return Err(PbrtError::error(&msg));
        }
    }
}
