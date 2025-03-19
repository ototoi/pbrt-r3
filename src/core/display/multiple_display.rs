use super::display::{Display, DisplayTile};
use crate::core::error::*;
use std::sync::Arc;
use std::sync::RwLock;

pub struct MutipleDisplay {
    pub displays: Vec<Arc<RwLock<dyn Display>>>,
}

impl MutipleDisplay {
    pub fn new() -> Self {
        MutipleDisplay {
            displays: Vec::new(),
        }
    }

    pub fn add_display(&mut self, display: &Arc<RwLock<dyn Display>>) {
        self.displays.push(display.clone());
    }

    pub fn len(&self) -> usize {
        return self.displays.len();
    }

    pub fn is_empty(&self) -> bool {
        return self.displays.is_empty();
    }
}

unsafe impl Sync for MutipleDisplay {}

impl Display for MutipleDisplay {
    fn start(
        &mut self,
        title: &str,
        resolution: &[usize; 2],
        channel_names: &[&str],
    ) -> Result<(), PbrtError> {
        for d in self.displays.iter() {
            let mut display = d.write().unwrap();
            display.start(title, resolution, channel_names)?;
        }
        return Ok(());
    }

    fn update(&mut self, tile: &DisplayTile) -> Result<(), PbrtError> {
        //TODO:multi thread
        for d in self.displays.iter() {
            let mut display = d.write().unwrap();
            display.update(tile)?;
        }
        return Ok(());
    }

    fn end(&mut self) -> Result<(), PbrtError> {
        for d in self.displays.iter() {
            let mut display = d.write().unwrap();
            display.end()?;
        }
        return Ok(());
    }
}
