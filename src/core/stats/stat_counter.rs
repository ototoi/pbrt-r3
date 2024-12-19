use super::stat_reporter::*;
use super::stats_accumlator::StatsAccumulator;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, Clone)]
struct BaseReporter {
    pub name: String,
    pub value: u64,
}

impl BaseReporter {
    pub fn new(name: &str) -> Self {
        BaseReporter {
            name: name.to_string(),
            value: 0,
        }
    }
    pub fn add(&mut self, val: u64) {
        self.value += val;
    }
}

#[derive(Debug, Clone)]
struct CountReporter(BaseReporter);
impl CountReporter {
    pub fn new(name: &str) -> Self {
        CountReporter(BaseReporter::new(name))
    }
}
impl StatReporter for CountReporter {
    fn report(&self, accum: &mut StatsAccumulator) {
        accum.report_counter(&self.0.name, self.0.value);
    }
    fn clear(&mut self) {
        self.0.value = 0;
    }
    fn add_int(&mut self, val: u64) {
        self.0.add(val);
    }
}

#[derive(Debug, Clone)]
struct MemoryReporter(BaseReporter);
impl MemoryReporter {
    pub fn new(name: &str) -> Self {
        MemoryReporter(BaseReporter::new(name))
    }
}
impl StatReporter for MemoryReporter {
    fn report(&self, accum: &mut StatsAccumulator) {
        accum.report_memory_counter(&self.0.name, self.0.value);
    }
    fn clear(&mut self) {
        self.0.value = 0;
    }
    fn add_int(&mut self, val: u64) {
        self.0.add(val);
    }
}

#[derive(Debug, Clone)]
pub struct IntDistributionReporter {
    pub name: String,
    pub sum: u64,
    pub count: u64,
    pub min: u64,
    pub max: u64,
}

impl IntDistributionReporter {
    pub fn new(name: &str) -> Self {
        IntDistributionReporter {
            name: name.to_string(),
            sum: 0,
            count: 0,
            min: std::u64::MAX,
            max: std::u64::MIN,
        }
    }

    pub fn add(&mut self, val: u64) {
        self.sum += val;
        self.count += 1;
        self.min = self.min.min(val);
        self.max = self.max.max(val);
    }
}

impl StatReporter for IntDistributionReporter {
    fn report(&self, accum: &mut StatsAccumulator) {
        accum.report_int_distribution(&self.name, self.sum, self.count, self.min, self.max);
    }
    fn clear(&mut self) {
        self.sum = 0;
        self.count = 0;
        self.min = std::u64::MAX;
        self.max = std::u64::MIN;
    }
    fn add_int(&mut self, val: u64) {
        self.add(val);
    }
}

#[derive(Debug, Clone)]
pub struct FloatDistributionReporter {
    pub name: String,
    pub sum: f64,
    pub count: u64,
    pub min: f64,
    pub max: f64,
}

impl FloatDistributionReporter {
    pub fn new(name: &str) -> Self {
        FloatDistributionReporter {
            name: name.to_string(),
            sum: 0.0,
            count: 0,
            min: std::f64::MAX,
            max: std::f64::MIN,
        }
    }

    pub fn add(&mut self, val: f64) {
        self.sum += val;
        self.count += 1;
        self.min = self.min.min(val);
        self.max = self.max.max(val);
    }
}

impl StatReporter for FloatDistributionReporter {
    fn report(&self, accum: &mut StatsAccumulator) {
        accum.report_float_distribution(&self.name, self.sum);
    }
    fn clear(&mut self) {
        self.sum = 0.0;
        self.count = 0;
        self.min = std::f64::MAX;
        self.max = std::f64::MIN;
    }
    fn add_float(&mut self, val: f64) {
        self.add(val);
    }
}

//-----------------------------------------------------------------------

pub struct StatCounter {
    reporter: Arc<RwLock<dyn StatReporter>>,
}

impl StatCounter {
    pub fn new(name: &str) -> Self {
        let reporter = Arc::new(RwLock::new(CountReporter::new(name)));
        register_stat_reporter(reporter.clone());
        StatCounter { reporter }
    }
    pub fn inc(&self) {
        self.add(1);
    }
    pub fn add(&self, val: u64) {
        let mut reporter = self.reporter.write().unwrap();
        reporter.add_int(val);
    }
}

pub struct StatMemoryCounter {
    reporter: Arc<RwLock<dyn StatReporter>>,
}

impl StatMemoryCounter {
    pub fn new(name: &str) -> Self {
        let reporter = Arc::new(RwLock::new(MemoryReporter::new(name)));
        register_stat_reporter(reporter.clone());
        StatMemoryCounter { reporter }
    }
    pub fn inc(&self) {
        self.add(1);
    }
    pub fn add(&self, val: usize) {
        let mut reporter = self.reporter.write().unwrap();
        reporter.add_int(val as u64);
    }
}

pub struct StatIntDistribution {
    reporter: Arc<RwLock<dyn StatReporter>>,
}

impl StatIntDistribution {
    pub fn new(name: &str) -> Self {
        let reporter = Arc::new(RwLock::new(IntDistributionReporter::new(name)));
        register_stat_reporter(reporter.clone());
        StatIntDistribution { reporter }
    }
    pub fn add(&self, val: u64) {
        let mut reporter = self.reporter.write().unwrap();
        reporter.add_int(val);
    }
}

pub struct StatFloatDistribution {
    reporter: Arc<RwLock<dyn StatReporter>>,
}

impl StatFloatDistribution {
    pub fn new(name: &str) -> Self {
        let reporter = Arc::new(RwLock::new(FloatDistributionReporter::new(name)));
        register_stat_reporter(reporter.clone());
        StatFloatDistribution { reporter }
    }
    pub fn add(&self, val: f64) {
        let mut reporter = self.reporter.write().unwrap();
        reporter.add_float(val);
    }
}
