# pbrt-r3

## What is pbrt-r3
pbrt-r3 is the rust implementation version of the pbrt-v3.

## License
pbrt-r3 is distributed under the BSD license based on the original pbrt-v3.

## Build
You can build it using cargo.
```
cargo build --release
```

## How to use
You can render using pbrt-r3 with the following command.
```
./target/release/pbrt-r3 -i <exmaple.pbrt>
```
### Tev display
pbrt-r3 supports the [tev](https://github.com/Tom94/tev) display implemented in pbrt-v4.
After starting tev, you can display the rendering progress by adding the following option.
```
./target/release/pbrt-r3 -i <exmaple.pbrt> --display-server localhost:14158
```


## Example scenes
pbrt-r3 can take pbrt-v3 scene files as input.
See the official [pbrt-v3 scenes page](http://pbrt.org/scenes-v3.html) on the pbrt website for information about how to download them.
### pavilion-night
![pavilion-night](https://github.com/user-attachments/assets/aa492fa3-493c-45f1-94e3-5c52a60086d6)
### bmw-m6
![bmw-m6](https://github.com/user-attachments/assets/33195ce5-f638-4e33-8414-e1d6e2c28485)





## Results
Rendered images are stored at [pbrt-r3-devkit](https://github.com/ototoi/pbrt-r3-devkit).

## Differences

## Feature works
