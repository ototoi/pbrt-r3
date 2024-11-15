use crate::core::pbrt::*;

pub struct ConstantTexture<T> {
    pub value: T,
}

impl<T: Copy> ConstantTexture<T> {
    pub fn new(value: &T) -> Self {
        return ConstantTexture::<T> { value: *value };
    }
}

impl<T: Copy> Texture<T> for ConstantTexture<T> {
    fn evaluate(&self, _si: &SurfaceInteraction) -> T {
        return self.value;
    }
}
