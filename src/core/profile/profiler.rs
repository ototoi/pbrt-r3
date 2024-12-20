use std::collections::HashMap;
use std::sync::LazyLock;
use std::sync::RwLock;
use std::thread::ThreadId;

type ThreadStateMap = HashMap<ThreadId, u64>;
static THREAD_STATE_MAP: LazyLock<RwLock<ThreadStateMap>> =
    LazyLock::new(|| RwLock::new(ThreadStateMap::new()));

pub fn set_profiler_state(state_bit: u64) {
    let tid = std::thread::current().id();
    let mut map = THREAD_STATE_MAP.write().unwrap();
    map.insert(tid, state_bit);
}

// For a given profiler state (i.e., a set of "on" bits corresponding to
// profiling categories that are active), ProfileSample stores a count of
// the number of times that state has been active when the timer interrupt
// to record a profiling sample has fired.
//struct ProfileSample {
//    pub profile_state: u64,
//    pub count: u64,
//}
type StateCountMap = HashMap<u64, u64>;
static STATE_COUNT_MAP: LazyLock<RwLock<StateCountMap>> =
    LazyLock::new(|| RwLock::new(StateCountMap::new()));

// SampleStateCount() increments the count of the number of times each
fn sample_state_count() {
    let mut count_map = STATE_COUNT_MAP.write().unwrap();
    let state_map = THREAD_STATE_MAP.read().unwrap();
    for (_, state) in state_map.iter() {
        let state = *state;
        if state == 0 {
            continue;
        }
        let count = count_map.entry(state).or_insert(0);
        *count += 1;
    }
}

#[derive(Debug, Clone, Copy)]
enum SamplerMessage {
    Stop,
}

type Sender = std::sync::mpsc::Sender<SamplerMessage>;
//type Receiver = std::sync::mpsc::Receiver<SamplerMessage>;

struct ProfileSampler {
    handles: Vec<(Sender, std::thread::JoinHandle<()>)>,
}

impl ProfileSampler {
    pub fn new() -> Self {
        ProfileSampler {
            handles: Vec::new(),
        }
    }

    pub fn start_sample(&mut self, interval: std::time::Duration) {
        if self.handles.len() > 0 {
            return;
        }
        let (tx, rx) = std::sync::mpsc::channel::<SamplerMessage>();
        let handle = std::thread::spawn(move || loop {
            match rx.try_recv() {
                Ok(SamplerMessage::Stop) => {
                    break;
                }
                _ => {}
            }
            sample_state_count();
            std::thread::sleep(interval);
        });
        self.handles.push((tx, handle));
    }

    pub fn stop_sample(&mut self) {
        for (tx, _) in self.handles.iter() {
            tx.send(SamplerMessage::Stop).unwrap();
        }
        for (_, handle) in self.handles.drain(..) {
            handle.join().unwrap();
        }
        self.handles.clear();
    }
}

impl Drop for ProfileSampler {
    fn drop(&mut self) {
        self.stop_sample();
    }
}

static PROFILE_SAMPLER: LazyLock<RwLock<ProfileSampler>> =
    LazyLock::new(|| RwLock::new(ProfileSampler::new()));

//-----------------------------------------------------------------------
pub fn start_profiler() {
    clear_profiler();
    let interval = std::time::Duration::from_millis(1000 / 100); // 100 Hz sampling
    let mut sampler = PROFILE_SAMPLER.write().unwrap();
    sampler.start_sample(interval);
}

pub fn stop_profiler() {
    let mut sampler = PROFILE_SAMPLER.write().unwrap();
    sampler.stop_sample();
}

pub fn clear_profiler() {
    THREAD_STATE_MAP.write().unwrap().clear();
    STATE_COUNT_MAP.write().unwrap().clear();
}

pub fn init_profiler() {
    //
}
pub fn cleanup_profiler() {
    //
}

pub fn suspend_profiler() {
    //
}

pub fn resume_profiler() {
    //
}

pub fn print_profiler() {
    report_profiler_results();
}

pub fn report_profiler_results() {
    let count_map = STATE_COUNT_MAP.read().unwrap();
    for (state, count) in count_map.iter() {
        println!("Profile state {:b}: {}", state, count);
    }
}
