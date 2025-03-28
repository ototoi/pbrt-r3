use crate::core::base::*;

pub trait Filter {
    fn evaluate(&self, p: &Point2f) -> Float;
    fn get_radius(&self) -> Vector2f;
}
