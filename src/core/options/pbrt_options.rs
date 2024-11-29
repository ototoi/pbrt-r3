use std::sync::LazyLock;
use std::sync::Mutex;

#[derive(Debug, Copy, Clone)]
pub struct PbrtOptions {
    pub quick_render: bool,
    pub quick_render_full_resolution: bool,
    //pub crop_window: [f32; 4],
}

impl Default for PbrtOptions {
    fn default() -> Self {
        PbrtOptions {
            quick_render: false,
            quick_render_full_resolution: false,
            //crop_window: [0.0, 1.0, 0.0, 1.0],
        }
    }
}

static PBRT_OPTIONS: LazyLock<Mutex<PbrtOptions>> =
    LazyLock::new(|| Mutex::new(PbrtOptions::new()));

impl PbrtOptions {
    pub fn new() -> Self {
        PbrtOptions::default()
    }

    pub fn set(opt: PbrtOptions) {
        let mut options = PBRT_OPTIONS.lock().unwrap();
        *options = opt;
    }

    pub fn get() -> PbrtOptions {
        let options = PBRT_OPTIONS.lock().unwrap();
        return *options;
    }
}
