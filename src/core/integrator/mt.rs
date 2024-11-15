use std::sync::Arc;
use std::sync::Mutex;
use std::sync::RwLock;

pub trait ThreadContainer: Send + Sync {
    fn func(&self);
}

pub struct AContainer {}
impl ThreadContainer for AContainer {
    fn func(&self) {
        let _idx = rayon::current_thread_index().unwrap();
        //println!("{}:AContainer::func",idx);
    }
}

pub struct ThreadExecutor {
    pub container: Arc<Mutex<dyn ThreadContainer>>, //Arc<Mutex<dyn ThreadContainer>>,
    //pub c1: Box<dyn ThreadContainer>,
    //pub c2: Arc<RefCell<dyn ThreadContainer>>,
    pub c3: Arc<RwLock<dyn ThreadContainer>>,
}

impl ThreadExecutor {
    pub fn new() -> Self {
        ThreadExecutor {
            container: Arc::new(Mutex::new(AContainer {})), //Arc::new(Mutex::new(AContainer{})),
            //c1: Box::new(AContainer{}),
            //c2: Arc::new(RwLock::new(AContainer{}))
            c3: Arc::new(RwLock::new(AContainer {})),
        }
    }
    pub fn test_task(&self, _x: i32, _y: i32) {
        {
            let c3 = self.c3.read().unwrap();
            c3.func();
        }
        {
            let _weak = Arc::downgrade(&self.c3);
            //weak.
        }

        {
            let c = self.container.lock().unwrap();
            c.func();
        }

        //println!("{}, {}", x, y);
        //self.container.borrow().func();
    }
    pub fn execute(&self) {
        let inputs = Arc::new(Mutex::new(Vec::new()));
        {
            let mut queue = inputs.lock().unwrap();
            for y in 0..100 {
                for x in 0..100 {
                    queue.push((x, y));
                }
            }
        }
        {
            let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
            //println!("{}", pool.current_num_threads());
            loop {
                let mut queue = inputs.lock().unwrap();
                if queue.is_empty() {
                    break;
                } else {
                    if let Some((x, y)) = queue.pop() {
                        pool.install(|| {
                            self.test_task(x, y);
                        });
                    }
                }
            }
        }
    }
}
