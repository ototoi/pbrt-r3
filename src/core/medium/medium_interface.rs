use super::medium::Medium;

use std::sync::Arc;
use std::sync::Weak;

#[derive(Clone, Default, Debug)]
pub struct MediumInterface {
    pub inside: Option<Weak<dyn Medium>>,
    pub outside: Option<Weak<dyn Medium>>,
}

impl MediumInterface {
    pub fn new() -> Self {
        MediumInterface {
            inside: None,
            outside: None,
        }
    }

    pub fn set_inside(&mut self, medium: &Arc<dyn Medium>) {
        self.inside = Some(Arc::downgrade(medium));
    }

    pub fn set_outside(&mut self, medium: &Arc<dyn Medium>) {
        self.outside = Some(Arc::downgrade(medium));
    }

    pub fn get_inside(&self) -> Option<Arc<dyn Medium>> {
        match &self.inside {
            Some(inside) => inside.upgrade(),
            None => None,
        }
    }

    pub fn get_outside(&self) -> Option<Arc<dyn Medium>> {
        match &self.outside {
            Some(outside) => outside.upgrade(),
            None => None,
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
        match medium {
            Some(medium) => MediumInterface::from(medium),
            None => MediumInterface::new(),
        }
    }
}

impl From<&Arc<dyn Medium>> for MediumInterface {
    fn from(medium: &Arc<dyn Medium>) -> Self {
        let medium = Arc::downgrade(medium);
        MediumInterface {
            inside: Some(medium.clone()),
            outside: Some(medium.clone()),
        }
    }
}
