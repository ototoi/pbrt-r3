use std::cell::LazyCell;
use std::sync::LazyLock;
use std::sync::Mutex;

#[derive(Debug, Copy, Clone)]
pub struct PbrtOptions {
    pub quick_render: bool,
    pub quick_render_full_resolution: bool,
    //pub crop_window: [f32; 4],
    pub no_stats: bool,
    pub no_profile: bool,
}

impl Default for PbrtOptions {
    fn default() -> Self {
        PbrtOptions {
            quick_render: false,
            quick_render_full_resolution: false,
            no_stats: false,
            no_profile: false,
            //crop_window: [0.0, 1.0, 0.0, 1.0],
        }
    }
}

static PBRT_OPTIONS: LazyLock<Mutex<PbrtOptions>> =
    LazyLock::new(|| Mutex::new(PbrtOptions::new()));

thread_local!(
    static LOCAL_OPTIONS: LazyCell<PbrtOptions> = LazyCell::new(|| PbrtOptions::get_lock());
);

impl PbrtOptions {
    pub fn new() -> Self {
        PbrtOptions::default()
    }

    pub fn set(opt: PbrtOptions) {
        let mut options = PBRT_OPTIONS.lock().unwrap();
        *options = opt;
    }

    pub fn get() -> PbrtOptions {
        LOCAL_OPTIONS.with(|opt| **opt)
    }

    pub fn get_lock() -> PbrtOptions {
        let options = PBRT_OPTIONS.lock().unwrap();
        return *options;
    }
}
