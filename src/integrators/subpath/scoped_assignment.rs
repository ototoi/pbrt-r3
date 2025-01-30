use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, Clone)]
pub struct ScopedValue<T: Clone> {
    pub value: Arc<RwLock<T>>,
}

impl<T: Clone> ScopedValue<T> {
    pub fn new(value: T) -> Self {
        Self {
            value: Arc::new(RwLock::new(value)),
        }
    }

    pub fn set(&self, value: T) {
        let mut target = self.value.write().unwrap();
        *target = value;
    }

    pub fn get(&self) -> T {
        let value = self.value.read().unwrap();
        return value.clone();
    }
}

impl<T: Clone> From<T> for ScopedValue<T> {
    fn from(value: T) -> Self {
        Self::new(value)
    }
}

pub struct ScopedAssignment<T: Clone> {
    pub target: Arc<RwLock<ScopedValue<T>>>,
    pub backup: Arc<RwLock<T>>,
}

impl<T: Clone> ScopedAssignment<T> {
    pub fn new(target: &Arc<RwLock<ScopedValue<T>>>, src: &ScopedValue<T>) -> Self {
        let target = target.clone();
        let backup = target.read().unwrap().value.clone();
        {
            let mut target = target.write().unwrap();
            target.value = src.value.clone();
        }
        Self { target, backup }
    }
    pub fn new_value(target: &Arc<RwLock<ScopedValue<T>>>, value: T) -> Self {
        let target = target.clone();
        let backup = target.read().unwrap().value.clone();
        {
            let mut target = target.write().unwrap();
            target.value = Arc::new(RwLock::new(value));
        }
        Self { target, backup }
    }
}

impl<T: Clone> Drop for ScopedAssignment<T> {
    fn drop(&mut self) {
        let mut target = self.target.write().unwrap();
        target.value = self.backup.clone();
    }
}
