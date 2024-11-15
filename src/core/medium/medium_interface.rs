use super::medium::Medium;

use std::sync::Arc;

#[derive(Clone, Default, Debug)]
pub struct MediumInterface {
    pub inside: Option<Arc<dyn Medium>>,
    pub outside: Option<Arc<dyn Medium>>,
}

impl MediumInterface {
    pub fn new() -> Self {
        MediumInterface {
            inside: None,
            outside: None,
        }
    }

    pub fn is_medium_transition(&self) -> bool {
        if let Some(inside) = self.inside.as_ref() {
            if let Some(outside) = self.outside.as_ref() {
                //print!("inside: {:p}\n", inside);
                let inside = inside.as_ref();
                let inside_ptr = inside as *const dyn Medium;
                let outside = outside.as_ref();
                let outside_ptr = outside as *const dyn Medium;
                return !std::ptr::eq(inside_ptr, outside_ptr);
                //return inside.as_ptr() != outside.as_ptr();
            } else {
                return true;
            }
        }
        return false;
    }
}

impl From<&Option<Arc<dyn Medium>>> for MediumInterface {
    fn from(medium: &Option<Arc<dyn Medium>>) -> Self {
        MediumInterface {
            inside: medium.clone(),
            outside: medium.clone(),
        }
    }
}

impl From<&Arc<dyn Medium>> for MediumInterface {
    fn from(medium: &Arc<dyn Medium>) -> Self {
        MediumInterface {
            inside: Some(Arc::clone(medium)),
            outside: Some(Arc::clone(medium)),
        }
    }
}
