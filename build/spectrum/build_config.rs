use super::super::utils::*;
use std::env;
use std::path::Path;

pub fn build() {
    println!("cargo:rerun-if-changed=build/spectrum/spectrum_config.rs");
    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("spectrum_config.rs");
    let src_path = Path::new("build/spectrum/spectrum_config.rs");
    let _ret = copy_if_modified(src_path.to_str().unwrap(), dest_path.to_str().unwrap());
}
