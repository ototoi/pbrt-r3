use super::stats_accumlator::*;

use std::sync::Arc;
use std::sync::Mutex;
use std::sync::RwLock;

pub trait StatReporter: Send + Sync {
    fn report(&self, accum: &mut StatsAccumulator);
    fn clear(&mut self);
}

pub type StatReporterRef = Arc<RwLock<dyn StatReporter>>;

static REGISTER_REPORTERS: Mutex<Vec<StatReporterRef>> = Mutex::new(Vec::new());

pub fn register_stat_reporter(reporter: StatReporterRef) {
    let mut reporters = REGISTER_REPORTERS.lock().unwrap();
    reporters.push(reporter);
}

pub fn report_stats(accum: &mut StatsAccumulator) {
    let reporters = REGISTER_REPORTERS.lock().unwrap();
    for reporter in reporters.iter() {
        let reporter = reporter.read().unwrap();
        reporter.report(accum);
    }
}

pub fn print_stats() {
    let mut accum = StatsAccumulator::new();
    report_stats(&mut accum);
    println!("{}", accum);
}

pub fn clear_stats() {
    let reporters = REGISTER_REPORTERS.lock().unwrap();
    for reporter in reporters.iter() {
        let mut reporter = reporter.write().unwrap();
        reporter.clear();
    }
}
