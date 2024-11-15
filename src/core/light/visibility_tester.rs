use crate::core::pbrt::*;

#[derive(Clone)]
pub struct VisibilityTester {
    pub p0: Interaction,
    pub p1: Interaction,
}

impl VisibilityTester {
    pub fn new() -> Self {
        VisibilityTester {
            p0: Interaction::zero(),
            p1: Interaction::zero(),
        }
    }

    pub fn unoccluded(&self, scene: &Scene) -> bool {
        return !scene.intersect_p(&self.p0.spawn_ray_to(&self.p1));
    }

    pub fn tr(&self, scene: &Scene, sampler: &mut dyn Sampler) -> Spectrum {
        let mut ray = self.p0.spawn_ray_to(&self.p1);
        let mut tr = Spectrum::one();
        loop {
            if let Some(isect) = scene.intersect(&ray) {
                if let Some(primitive) = isect.get_primitive() {
                    if primitive.get_material().is_some() {
                        return Spectrum::zero();
                    }
                }

                if let Some(medium) = ray.medium.as_ref() {
                    tr *= medium.as_ref().tr(&ray, sampler);
                }

                let tisect = Interaction::from(isect);
                ray = tisect.spawn_ray_to(&self.p1);
            } else {
                if let Some(medium) = ray.medium.as_ref() {
                    tr *= medium.as_ref().tr(&ray, sampler);
                }
                break;
            }
        }
        return tr;
    }
}

impl From<(&Interaction, &Interaction)> for VisibilityTester {
    fn from(value: (&Interaction, &Interaction)) -> Self {
        VisibilityTester {
            p0: value.0.clone(),
            p1: value.1.clone(),
        }
    }
}

impl From<(Interaction, Interaction)> for VisibilityTester {
    fn from(value: (Interaction, Interaction)) -> Self {
        VisibilityTester {
            p0: value.0,
            p1: value.1,
        }
    }
}
