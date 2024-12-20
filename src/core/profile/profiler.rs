// For a given profiler state (i.e., a set of "on" bits corresponding to
// profiling categories that are active), ProfileSample stores a count of
// the number of times that state has been active when the timer interrupt
// to record a profiling sample has fired.
struct ProfileSample {
    pub profile_state: u64,
    pub count: u64,
}

pub fn init_profiler() {
    //
}

pub fn suspend_profiler() {
    //
}
