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
        match (self.inside.as_ref(), self.outside.as_ref()) {
            (Some(inside), Some(outside)) => {
                return !std::ptr::eq(inside, outside);
            }
            (Some(_), None) => true,
            (None, Some(_)) => true,
            (None, None) => false,
        }
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
