use std::env;
use std::fs;
use std::fs::{File, Metadata};
use std::path::Path;

/*
macro_rules! p {
    ($($tokens: tt)*) => {
        println!("cargo:warning={}", format!($($tokens)*))
    }
}
*/

pub fn get_meta_data(path: &str) -> Result<Metadata, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let meta = file.metadata()?;
    return Ok(meta);
}

pub fn should_build(path: &str, depends: &[String]) -> bool {
    match get_meta_data(path) {
        Ok(meta_a) => {
            let updated_a = meta_a.modified().unwrap();
            for dep_path in depends {
                match get_meta_data(dep_path) {
                    Ok(meta_b) => {
                        let updated_b = meta_b.modified().unwrap();
                        if updated_a < updated_b {
                            return true;
                        }
                    }
                    Err(_) => {
                        return true;
                    }
                }
            }
            return false;
        }
        Err(_) => {
            return true;
        }
    }
}

pub fn make_depends_path(depends: &[&str]) -> Vec<String> {
    let mut paths = Vec::new();
    let root_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    for path in depends {
        let path = Path::new(&root_dir).join(*path);
        let path = String::from(path.to_str().unwrap());
        paths.push(path);
    }
    return paths;
}

pub fn copy_if_modified(src: &str, dst: &str) -> Result<u64, String> {
    if !Path::new(dst).exists() {
        match fs::copy(src, dst) {
            Ok(size) => return Ok(size),
            Err(e) => return Err(e.to_string()),
        }
    } else {
        if let Ok(meta_src) = get_meta_data(src) {
            if let Ok(meta_dst) = get_meta_data(dst) {
                let modified_src = meta_src.modified().unwrap();
                let modified_dst = meta_dst.modified().unwrap();
                if modified_src > modified_dst {
                    match fs::copy(src, dst) {
                        Ok(size) => return Ok(size),
                        Err(e) => return Err(e.to_string()),
                    }
                }
            }
        }
        return Err(String::from(""));
    }
}
