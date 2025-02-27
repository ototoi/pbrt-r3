use super::super::utils::*;
use super::cie_data;
use super::config::*;
use super::to_string::spectrum_to_string;
use super::utils::sample_spectrum;
use std::env;
use std::fs;
use std::path::Path;

pub fn build_core(path: &str) {
    let mut m: Vec<(&str, [Float; SPECTRAL_SAMPLES])> = Vec::new();
    m.push((
        "CIE_X",
        sample_spectrum(&cie_data::CIE_LAMBDA, &cie_data::CIE_X),
    ));

    m.push((
        "CIE_Y",
        sample_spectrum(&cie_data::CIE_LAMBDA, &cie_data::CIE_Y),
    ));

    m.push((
        "CIE_Z",
        sample_spectrum(&cie_data::CIE_LAMBDA, &cie_data::CIE_Z),
    ));

    let mut contents = String::from("");
    contents += "use crate::core::pbrt::Float;\n";
    contents += "\n";

    contents += &format!("const SPECTRAL_SAMPLES: usize = {};\n", SPECTRAL_SAMPLES);
    contents += "\n";
    for (key, v) in m.iter() {
        let name = format!("ARRAY_{}", *key);
        contents += &spectrum_to_string(&name, v);
        contents += "\n";
    }

    fs::write(path, &contents).unwrap();
}

pub fn build() {
    println!("cargo:rerun-if-changed=build/spectrum/cie_data.rs;build/spectrum/build_xyz.rs");
    let depends = ["build/spectrum/cie_data.rs", "build/spectrum/build_xyz.rs"];
    let target = "spectrum_data_xyz.rs";
    let out_dir = env::var("OUT_DIR").unwrap();
    let path = Path::new(&out_dir).join(target);
    let path = String::from(path.to_str().unwrap());
    let depends = make_depends_path(&depends);
    if should_build(&path, &depends) {
        build_core(&path);
    }
}
