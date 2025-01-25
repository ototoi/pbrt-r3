use std::sync::Arc;
use std::sync::RwLock;

pub struct ScopedAssignment<T: Clone> {
    pub target: Arc<RwLock<T>>,
    pub backup: T,
}

impl<T: Clone> ScopedAssignment<T> {
    pub fn new(target: &Arc<RwLock<T>>, value: &T) -> Self {
        let target = target.clone();
        let backup;
        {
            let mut target = target.as_ref().write().unwrap();
            backup = target.clone();
            *target = value.clone();
        }
        Self {
            target: target,
            backup: backup,
        }
    }
}

impl<T: Clone> Drop for ScopedAssignment<T> {
    fn drop(&mut self) {
        let mut target = self.target.as_ref().write().unwrap();
        *target = self.backup.clone();
    }
}
