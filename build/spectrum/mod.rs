mod build_cie;
mod build_config;
mod build_rgb_refl;
mod build_utils;
mod build_xyz;
mod cie_data;
mod config;
mod spectrum_config;
mod rgb_data;
mod to_string;
mod utils;

pub fn build() {
    build_config::build();
    build_utils::build();
    build_cie::build();
    build_xyz::build();
    build_rgb_refl::build();
}
