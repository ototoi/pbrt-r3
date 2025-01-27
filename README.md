# pbrt-r3
[![Rust](https://github.com/ototoi/pbrt-r3/actions/workflows/rust.yml/badge.svg)](https://github.com/ototoi/pbrt-r3/actions/workflows/rust.yml)
![License](https://img.shields.io/github/license/ototoi/pbrt-r3)
![GitHub Release](https://img.shields.io/github/v/release/ototoi/pbrt-r3)



## What is pbrt-r3
pbrt-r3 is a rust implementation version of the [pbrt-v3](https://github.com/mmp/pbrt-v3).

## License
pbrt-r3 is distributed under [the BSD license based on the original pbrt-v3](https://github.com/mmp/pbrt-v3/blob/master/LICENSE.txt).

## Build
You can build this using cargo.
```
cargo build --release
```

## How to use
You can render using pbrt-r3 with the following command.
```
./target/release/pbrt-r3 -i <exmaple.pbrt>
```
## Tev display
pbrt-r3 supports the [tev](https://github.com/Tom94/tev) display implemented in pbrt-v4.
After starting tev, you can display your rendering progress by adding the following option.
```
./target/release/pbrt-r3 -i <exmaple.pbrt> --display-server localhost:14158
```
## Stats and Profile
If you want to use the stats and profile features, please do as follows.
### Build
```
cargo build --release --features stats --features profile
```
### Use
```
./target/release/pbrt-r3 -i <exmaple.pbrt> --stats --profile
```

## Example scenes
pbrt-r3 can take pbrt-v3 scene files as input.
See the official [pbrt-v3 scenes page](http://pbrt.org/scenes-v3.html) on the pbrt website for information about how to download them.

![images](https://github.com/user-attachments/assets/ce1bebc6-8377-4da7-8b49-38e5073a397e)
Rendered images are stored at [pbrt-r3-devkit](https://github.com/ototoi/pbrt-r3-devkit).




## Differences

## Feature works
