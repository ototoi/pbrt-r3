use super::stat_reporter::*;
use super::stats_accumlator::StatsAccumulator;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Debug, Clone)]
struct BaseCountReporter {
    pub name: String,
    pub value: u64,
}

impl BaseCountReporter {
    pub fn new(name: &str) -> Self {
        BaseCountReporter {
            name: name.to_string(),
            value: 0,
        }
    }
    pub fn add(&mut self, val: u64) {
        self.value += val;
    }
}

#[derive(Debug, Clone)]
struct CountReporter(BaseCountReporter);
impl CountReporter {
    pub fn new(name: &str) -> Self {
        CountReporter(BaseCountReporter::new(name))
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
struct MemoryReporter(BaseCountReporter);
impl MemoryReporter {
    pub fn new(name: &str) -> Self {
        MemoryReporter(BaseCountReporter::new(name))
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

#[derive(Debug, Clone)]
pub struct BaseFractionReporter {
    pub name: String,
    pub num: u64,
    pub denom: u64,
}

impl BaseFractionReporter {
    pub fn new(name: &str) -> Self {
        BaseFractionReporter {
            name: name.to_string(),
            num: 0,
            denom: 0,
        }
    }
    fn clear(&mut self) {
        self.num = 0;
        self.denom = 0;
    }
    pub fn add_num(&mut self, val: u64) {
        self.num += val;
    }
    pub fn add_denom(&mut self, val: u64) {
        self.denom += val;
    }
}

#[derive(Debug, Clone)]
struct PercentageReporter(BaseFractionReporter);

impl PercentageReporter {
    pub fn new(name: &str) -> Self {
        PercentageReporter(BaseFractionReporter::new(name))
    }
}

impl StatReporter for PercentageReporter {
    fn report(&self, accum: &mut StatsAccumulator) {
        accum.report_percentage(&self.0.name, self.0.num, self.0.denom);
    }
    fn clear(&mut self) {
        self.0.clear();
    }
    fn add_num(&mut self, val: u64) {
        self.0.add_num(val);
    }
    fn add_denom(&mut self, val: u64) {
        self.0.add_denom(val);
    }
}

#[derive(Debug, Clone)]
struct RatioReporter(BaseFractionReporter);

impl RatioReporter {
    pub fn new(name: &str) -> Self {
        RatioReporter(BaseFractionReporter::new(name))
    }
}

impl StatReporter for RatioReporter {
    fn report(&self, accum: &mut StatsAccumulator) {
        accum.report_percentage(&self.0.name, self.0.num, self.0.denom);
    }
    fn clear(&mut self) {
        self.0.clear();
    }
    fn add_num(&mut self, val: u64) {
        self.0.add_num(val);
    }
    fn add_denom(&mut self, val: u64) {
        self.0.add_denom(val);
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

pub struct StatPercent {
    reporter: Arc<RwLock<dyn StatReporter>>,
}

impl StatPercent {
    pub fn new(name: &str) -> Self {
        let reporter = Arc::new(RwLock::new(PercentageReporter::new(name)));
        register_stat_reporter(reporter.clone());
        StatPercent { reporter }
    }
    pub fn add_num(&self, val: u64) {
        let mut reporter = self.reporter.write().unwrap();
        reporter.add_num(val);
    }
    pub fn add_denom(&self, val: u64) {
        let mut reporter = self.reporter.write().unwrap();
        reporter.add_denom(val);
    }
}

pub struct StatRatio {
    reporter: Arc<RwLock<dyn StatReporter>>,
}

impl StatRatio {
    pub fn new(name: &str) -> Self {
        let reporter = Arc::new(RwLock::new(PercentageReporter::new(name)));
        register_stat_reporter(reporter.clone());
        StatRatio { reporter }
    }
    pub fn add_num(&self, val: u64) {
        let mut reporter = self.reporter.write().unwrap();
        reporter.add_num(val);
    }
    pub fn add_denom(&self, val: u64) {
        let mut reporter = self.reporter.write().unwrap();
        reporter.add_denom(val);
    }
}
