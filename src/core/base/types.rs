use crate::core::geometry::*;

#[cfg(not(feature = "float-as-double"))]
pub type Float = f32;
#[cfg(feature = "float-as-double")]
pub type Float = f64;

pub type Vector2i = Vector2<i32>;
pub type Point2i = Vector2<i32>;
pub type Vector2f = Vector2<Float>;
pub type Point2f = Vector2<Float>;
pub type Normal2f = Vector2<Float>;

pub type Vector3i = Vector3<i32>;
pub type Point3i = Vector3<i32>;
pub type Vector3f = Vector3<Float>;
pub type Point3f = Vector3<Float>;
pub type Normal3f = Vector3<Float>;

//pub type Bounds2i = Bounds2<i32>;
//pub type Boundsd2f = Bounds2<Float>;

//pub type Bounds3i = Bounds3<i32>;
//pub type Bounds3f = Bounds3<Float>;
