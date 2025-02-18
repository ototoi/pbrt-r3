use super::profile_category::ProfileCategory;

#[derive(Debug)]
pub struct ProfilePhase {
    #[allow(dead_code)]
    reset: bool,
    pub category_bit: u64,
}

#[cfg(feature = "profile")]
mod _impl {
    use super::super::profiler::set_profiler_state;
    use super::*;
    use crate::core::options::PbrtOptions;
    use std::cell::Cell;

    thread_local!(static PROFILER_STATE: Cell<u64> = Cell::new(0));

    impl ProfilePhase {
        pub fn new(category: ProfileCategory) -> Self {
            let options = PbrtOptions::get();
            if !options.profile {
                ProfilePhase {
                    reset: false,
                    category_bit: 0,
                }
            } else {
                let category_bit = category.to_bits();
                let reset = (PROFILER_STATE.get() & category_bit) == 0;
                PROFILER_STATE.with(|state| {
                    let bit = state.get() | category_bit;
                    state.set(bit);
                    if reset {
                        set_profiler_state(bit);
                    }
                });
                ProfilePhase {
                    reset,
                    category_bit,
                }
            }
        }
    }

    impl Drop for ProfilePhase {
        fn drop(&mut self) {
            if self.reset {
                PROFILER_STATE.with(|state| {
                    let bit = state.get() & !self.category_bit;
                    state.set(bit);
                    set_profiler_state(bit);
                });
            }
        }
    }
}

#[cfg(not(feature = "profile"))]
mod _impl {
    use super::*;

    impl ProfilePhase {
        pub fn new(_category: ProfileCategory) -> Self {
            ProfilePhase {
                reset: false,
                category_bit: 0,
            }
        }
    }
}
