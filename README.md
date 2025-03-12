# pbrt-r3
[![Rust](https://github.com/ototoi/pbrt-r3/actions/workflows/rust.yml/badge.svg)](https://github.com/ototoi/pbrt-r3/actions/workflows/rust.yml)
[![License](https://img.shields.io/github/license/ototoi/pbrt-r3)](LICENSE)
[![GitHub Release](https://img.shields.io/github/v/release/ototoi/pbrt-r3)](https://github.com/ototoi/pbrt-r3/releases/latest)



## What is pbrt-r3
pbrt-r3 is a rust implementation version of the [pbrt-v3](https://github.com/mmp/pbrt-v3).

## License
pbrt-r3 is distributed under [the BSD license](LICENSE) based on [the original pbrt-v3](https://github.com/mmp/pbrt-v3/blob/master/LICENSE.txt).

## Build
You can build this using cargo.
```
cargo build --release
```

## How to use
You can render using pbrt-r3 with the following command.
```
./target/release/pbrt-r3 -i <example.pbrt>
```
## Tev display
pbrt-r3 supports the [tev](https://github.com/Tom94/tev) display implemented in pbrt-v4.
After starting tev, you can display your rendering progress by adding the following option.
```
./target/release/pbrt-r3 -i <example.pbrt> --display-server localhost:14158
```
## Stats and Profile
If you want to use the stats and profile features, please do as follows.
### Build
```
cargo build --release --features stats --features profile
```
### Use
```
./target/release/pbrt-r3 -i <example.pbrt> --stats --profile
```

## Build Options
pbrt-r3 has several options as Rust features.
| Feature | Description |
|---------|-------------|
| `profile` | Enables the `--profile` option. |
| `stats` | Enables the `--stats` option. |
| `float-as-double` | Uses double precision (64-bit) for floating-point calculations. This increases precision but also increases execution time and memory usage. |
| `sampled-spectrum` | Uses SampledSpectrum to represent the Spectrum type for colors. SampledSpectrum represents visible light with more samples (60 samples) instead of just three colors. |

## Example scenes
pbrt-r3 can take pbrt-v3 scene files as input.
See the official [pbrt-v3 scenes page](http://pbrt.org/scenes-v3.html) on the pbrt website for information about how to download them.

![images](https://github.com/user-attachments/assets/ce1bebc6-8377-4da7-8b49-38e5073a397e)
Rendered images are stored at [pbrt-r3-devkit](https://github.com/ototoi/pbrt-r3-devkit).


## Differences between pbrt-v3 and pbrt-r3

- **Language**: pbrt-v3 is implemented in C++, while pbrt-r3 is a re-implementation in Rust.
  - **Memory Safety**: Rust's ownership model in pbrt-r3 helps prevent common memory safety issues such as null pointer dereferencing and buffer overflows.
  - **Build System**: pbrt-r3 uses Cargo, Rust's package manager and build system, simplifying dependency management and build processes.
  - **Syntax**: Differences arise due to the syntax differences between C++ and Rust.
    - **Function Return Values**: When it is necessary to return multiple values in a function, C++ achieves this by passing pointers to arguments and return values. In Rust, this is directly achieved using tuples.
    - **Function Overloading**: Function overloading exists in C++ but not in Rust. In Rust, this is achieved by using different names for each function.
    - **Class Inheritance**: Class inheritance exists in C++, but not in Rust. In Rust, composition is used. In other words, a class member holds an instance of the parent class as `base`, and this is called as needed.
- **External Libraries**: Some features use external libraries.
  - **Parallel Processing**: In pbrt-v3, `OpenMP` was used for parallel processing. On the other hand, pbrt-r3 uses `rayon` crate for parallel processing. 
  - **Parser**: In pbrt-v3, a custom parser was implemented. On the other hand, in Rust, the `nom` crate was used to implement the parser.
  - **Progress Bar**: In pbrt-v3, a custom progress bar was implemented, while in pbrt-r3, the `indicatif` crate is used.
  - **Command Line Options**: In pbrt-v3, command line options were implemented independently, while in pbrt-r3, the `clap` crate is used.
  - **Logging**: In pbrt-v3, Google glog was used, while in pbrt-r3, the `log` and `env_logger` crates are used.
  - **Others**: In addition to the above, pbrt-r3 uses the following crates:
    - Image loading and saving: `image`
    - Hash functions: `rust-crypto`
    - PLY file loading and saving: `ply-rs`
    - JSON loading and saving: `serde`, `serde_json`
- **Additional Features**: pbrt-r3 implements some additional features from pbrt-v3.
  - **Tev Display**: You can display the image during rendering using the external image viewer Tev.
  - **QBVH Accelerator**: QBVH is implemented as an accelerator, enabling fast traversal using SIMD instructions.
  - **AOV Integrator**: When AOV is specified as an Integrator, it can render some attribute information of the scene (e.g., position, normal, etc.). This is useful for compositing and debugging rendering results.
- **Unimplemented Features**: Some features are not implemented in pbrt-r3.
  - **PTex texture**: PTex texturing is not implemented in pbrt-r3.
  - **Other Commands**: Several commands are implemented under `src/tools` in pbrt-v3. However, they are not necessarily needed in pbrt-r3, so their priority is low.
- **Bug Fixes**: pbrt-r3 has fixed several bugs from pbrt-v3.
  - Do not make `pdf` negative: There are places where pdf becomes negative due to calculation errors during the pdf calculation process. If left as is, it can cause infinity or negative colors. For example, this was fixed by using `pdf = pdf.max(0.0)`.
  - Direction of `shading.dpdu/dpdv`: There was a bug where the directions of `shading.dpdu` and `shading.dpdv` were reversed.

## Future Plans
There are some remaining tasks:
- Benchmark: Rendering speed tends to be slower than pbrt-v3 implemented in C++. Create and run benchmarks to identify bottlenecks.
- Register on Crates.io: Register pbrt-r3 on Crates.io.
- Copy Comments: To respect the original implementation of pbrt-v3, we want to transplant as many original comments as possible to pbrt-r3.
- Organize License: Correctly describe the license terms.
- Create Documentation: Create documentation.

