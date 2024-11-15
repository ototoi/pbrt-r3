use super::super::display::MutipleDisplay;
use super::super::imageio::*;
use super::film_tile::*;
use super::splat_tile::*;
use crate::core::pbrt::*;
//use atomic_float::AtomicF32;
use log::*;
use std::env;
use std::path::Path;
use std::sync::{Arc, Mutex, RwLock};

//pub const FILTER_TABLE_WIDTH: usize = 16;
//pub const FT_W: usize = FILTER_TABLE_WIDTH;
//pub const FT_SZ: usize = FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH;

#[derive(Debug, Default, Copy, Clone)]
struct Pixel {
    pub xyz: [Float; 3],
    pub filter_weight_sum: Float,
}

impl Pixel {
    pub fn zero() -> Self {
        Pixel {
            xyz: [0.0; 3],
            filter_weight_sum: 0.0,
        }
    }
}

pub struct Film {
    pub full_resolution: Point2i,
    pub diagonal: Float,
    pub filter: Arc<dyn Filter>,
    pub filename: String,
    pub cropped_pixel_bounds: Bounds2i,
    pub display: MutipleDisplay,

    pixels: Mutex<Vec<Pixel>>,

    splat_tiles: Vec<Arc<RwLock<SplatTile>>>, //pixels_atomic: Vec<PixelAtomic>,
    splat_size: Vector2i,

    filter_table: Arc<RwLock<[Float; FT_SZ]>>,

    scale: Float,
    max_sample_luminance: Float,
}

impl Film {
    pub fn new(
        resolution: &Point2i,
        crop_window: &Bounds2f,
        filter: &Arc<dyn Filter>,
        diagonal: Float,
        filename: &str,
        scale: Float,
        max_sample_luminance: Float,
    ) -> Self {
        let x0 = i32::max(
            0,
            Float::floor(resolution.x as f32 * crop_window.min.x) as i32,
        ); //ceil -> floor
        let y0 = i32::max(
            0,
            Float::floor(resolution.y as f32 * crop_window.min.y) as i32,
        ); //ceil -> floor
        let x1 = i32::min(
            Float::ceil(resolution.x as f32 * crop_window.max.x) as i32,
            resolution[0],
        );
        let y1 = i32::min(
            Float::ceil(resolution.y as f32 * crop_window.max.y) as i32,
            resolution[1],
        );

        let cropped_pixel_bounds = Bounds2i::from(((x0, y0), (x1, y1)));

        info!(
            "Created film with full resolution {:?}. Crop window of {:?} -> CroppedPixelBounds {:?}",
            resolution, crop_window, cropped_pixel_bounds
        );

        let pixels: Vec<Pixel> = vec![Pixel::zero(); cropped_pixel_bounds.area() as usize];
        let mut filter_table: [Float; FT_SZ] = [1.0; FT_SZ];
        {
            {
                let f = filter.as_ref();
                let radius = f.get_radius();
                for y in 0..FT_W {
                    for x in 0..FT_W {
                        //original
                        //let xx = (x as f32 + 0.5) * (radius.x / FT_W as f32);
                        //let yy = (y as f32 + 0.5) * (radius.y / FT_W as f32);

                        //fixed
                        let xx = (x as f32) * (radius.x / (FT_W - 1) as f32);
                        let yy = (y as f32) * (radius.y / (FT_W - 1) as f32);
                        filter_table[y * FT_W + x] = f.evaluate(&Vector2f::new(xx, yy));
                    }
                }
            }
            /*
            {
                let mut total = 0.0;
                for y in 0..FT_W {
                    for x in 0..FT_W {
                        total += filter_table[y * FT_W + x];
                    }
                }
                let inv_total = 1.0 / total;
                for y in 0..FT_W {
                    for x in 0..FT_W {
                        filter_table[y * FT_W + x] *= inv_total;
                    }
                }
                //println!("{:?}", filter_table);
            }
            */
        }

        let mut splat_tiles = Vec::new();
        let splat_size;
        {
            let w = (cropped_pixel_bounds.max.x - cropped_pixel_bounds.min.x) as usize;
            let h = (cropped_pixel_bounds.max.y - cropped_pixel_bounds.min.y) as usize;
            let nx = (w + ST_W - 1) / ST_W;
            let ny = (h + ST_W - 1) / ST_W;
            splat_size = Vector2i::new(nx as i32, ny as i32);

            for ty in 0..ny {
                for tx in 0..nx {
                    let x0 = cropped_pixel_bounds.min.x as usize + tx * ST_W;
                    let x1 = usize::min(x0 + ST_W, cropped_pixel_bounds.max.x as usize);
                    let y0 = cropped_pixel_bounds.min.y as usize + ty * ST_W;
                    let y1 = usize::min(y0 + ST_W, cropped_pixel_bounds.max.y as usize);
                    let pixel_bounds = Bounds2i::new(
                        &Vector2i::new(x0 as i32, y0 as i32),
                        &Vector2i::new(x1 as i32, y1 as i32),
                    );
                    splat_tiles.push(Arc::new(RwLock::new(SplatTile::new(&pixel_bounds))));
                }
            }
        }
        let diagonal = 0.001 * diagonal;
        Film {
            full_resolution: *resolution,
            diagonal,
            filter: Arc::clone(filter),
            filename: String::from(filename),
            cropped_pixel_bounds,
            display: MutipleDisplay::new(),

            pixels: Mutex::new(pixels),
            splat_tiles: splat_tiles,
            splat_size: splat_size,
            filter_table: Arc::new(RwLock::new(filter_table)),
            scale,
            max_sample_luminance,
        }
    }

    pub fn get_sample_bounds(&self) -> Bounds2i {
        let radius = self.filter.as_ref().get_radius();
        /*
        let x0 = f32::floor(self.cropped_pixel_bounds.min.x as f32 + 0.5 - radius.x) as i32;
        let y0 = f32::floor(self.cropped_pixel_bounds.min.y as f32 + 0.5 - radius.y) as i32;
        let x1 = f32::ceil(self.cropped_pixel_bounds.max.x as f32 + 0.5 + radius.x) as i32;
        let y1 = f32::ceil(self.cropped_pixel_bounds.max.y as f32 + 0.5 + radius.y) as i32;
        */
        let x0 = f32::floor(self.cropped_pixel_bounds.min.x as f32 - radius.x) as i32;
        let y0 = f32::floor(self.cropped_pixel_bounds.min.y as f32 - radius.y) as i32;
        let x1 = f32::ceil(self.cropped_pixel_bounds.max.x as f32 + radius.x) as i32;
        let y1 = f32::ceil(self.cropped_pixel_bounds.max.y as f32 + radius.y) as i32;
        Bounds2i::new(&Vector2i::new(x0, y0), &Vector2i::new(x1, y1))
    }

    pub fn get_physical_extent(&self) -> Bounds2f {
        let aspect = self.full_resolution.y as f32 / self.full_resolution.x as f32;
        let x = f32::sqrt(self.diagonal * self.diagonal / (1.0 + aspect * aspect));
        let y = aspect * x;
        Bounds2f::new(
            &Point2f::new(-x / 2.0, -y / 2.0),
            &Point2f::new(x / 2.0, y / 2.0),
        )
    }

    fn i2f(v: &Vector2i) -> Vector2f {
        Vector2f::new(v.x as Float, v.y as Float)
    }

    fn ceil(v: &Vector2f) -> Vector2i {
        Vector2i::new(v.x.ceil() as i32, v.y.ceil() as i32)
    }

    fn floor(v: &Vector2f) -> Vector2i {
        Vector2i::new(v.x.floor() as i32, v.y.floor() as i32)
    }

    pub fn get_film_tile(&self, sample_bounds: &Bounds2i) -> FilmTile {
        let radius = self.filter.as_ref().get_radius();
        //let half_pixel = Vector2f::new(0.5, 0.5);
        //let p0 = Self::ceil(&(Self::i2f(&sample_bounds.min) - half_pixel - radius));
        //let p1 = Self::floor(&(Self::i2f(&sample_bounds.max) - half_pixel + radius)) + Vector2i::new(1,1);
        let p0 = Self::floor(&(Self::i2f(&sample_bounds.min) - radius));
        let p1 = Self::ceil(&(Self::i2f(&sample_bounds.max) + radius));
        let tile_pixel_bounds = Bounds2i::new(&p0, &p1).intersect(&self.cropped_pixel_bounds);
        FilmTile::new(
            &tile_pixel_bounds,
            &radius,
            &self.filter_table,
            self.max_sample_luminance,
        )
    }

    pub fn merge_film_tile(&mut self, tile: &FilmTile) {
        let bounds = tile.get_pixel_bounds();
        {
            let mut pixels = self.pixels.lock().unwrap();
            for y in bounds.min.y..bounds.max.y {
                assert!(y >= 0);
                for x in bounds.min.x..bounds.max.x {
                    assert!(x >= 0);
                    let p = Vector2i::from((x, y));
                    let src_index = tile.get_pixel_index(&p);
                    let dst_index = self.get_pixel_index(&p);
                    let xyz = tile.pixels[src_index].contrib_sum.to_xyz();
                    let filter_weight_sum = tile.pixels[src_index].filter_weight_sum;
                    for i in 0..3 {
                        pixels[dst_index].xyz[i] += xyz[i];
                    }
                    pixels[dst_index].filter_weight_sum += filter_weight_sum;
                }
            }
        }
        {
            self.update_display(&bounds);
        }
    }

    pub fn merge_splats(&mut self, splat_scale: Float) {
        for tile in self.splat_tiles.iter() {
            let mut tile = tile.write().unwrap();
            if !tile.dirty {
                continue;
            }

            let bounds = tile.pixel_bounds;
            let x0 = bounds.min.x as usize;
            let x1 = bounds.max.x as usize;
            let y0 = bounds.min.y as usize;
            let y1 = bounds.max.y as usize;
            {
                let mut pixels = self.pixels.lock().unwrap();
                for y in y0..y1 {
                    for x in x0..x1 {
                        //
                        let p = Vector2i::from((x as i32, y as i32));
                        let sx = x - x0;
                        let sy = y - y0;
                        let src_index = sy * ST_W + sx;
                        let dst_index = self.get_pixel_index(&p);
                        for i in 0..3 {
                            pixels[dst_index].xyz[i] += splat_scale * tile.pixels[src_index][i];
                        }
                    }
                }
            }

            for i in 0..tile.pixels.len() {
                tile.pixels[i] = [0.0; 3];
            }
            tile.dirty = false;
        }
    }

    pub fn update_display(&mut self, bounds: &Bounds2i) {
        self.update_display_scale(bounds, 1.0);
    }

    pub fn update_display_scale(&mut self, bounds: &Bounds2i, scale: Float) {
        if self.display.is_empty() {
            return;
        }

        let scale = self.scale * scale;

        let tx0 = i32::max(self.cropped_pixel_bounds.min.x, bounds.min.x) as usize;
        let ty0 = i32::max(self.cropped_pixel_bounds.min.y, bounds.min.y) as usize;
        let tx1 = i32::min(self.cropped_pixel_bounds.max.x, bounds.max.x) as usize;
        let ty1 = i32::min(self.cropped_pixel_bounds.max.y, bounds.max.y) as usize;

        if tx0 >= tx1 || ty0 >= ty1 {
            return;
        }

        let twidth = tx1 - tx0;
        let theight = ty1 - ty0;

        let mut buffer: Vec<f32> = vec![0.0; (3 * twidth * theight) as usize];
        {
            let pixels = self.pixels.lock().unwrap();
            for y in ty0..ty1 {
                for x in tx0..tx1 {
                    let p = Vector2i::from((x as i32, y as i32));
                    let src_index = self.get_pixel_index(&p);
                    let pixel = &pixels[src_index];
                    let mut c = xyz_to_rgb(&pixel.xyz);
                    let filter_weight_sum = pixel.filter_weight_sum;
                    if filter_weight_sum > 0.0 {
                        let inv_wt = 1.0 / filter_weight_sum;
                        c[0] = f32::max(0.0, c[0] * inv_wt);
                        c[1] = f32::max(0.0, c[1] * inv_wt);
                        c[2] = f32::max(0.0, c[2] * inv_wt);
                    }

                    let by = y - ty0;
                    let bx = x - tx0;
                    let idx = by * twidth + bx;

                    buffer[3 * idx + 0] = c[0] * scale;
                    buffer[3 * idx + 1] = c[1] * scale;
                    buffer[3 * idx + 2] = c[2] * scale;
                }
            }
        }

        {
            let display_tile = DisplayTile {
                x: tx0,
                y: ty0,
                width: twidth,
                height: theight,
                buffer,
            };
            self.display
                .update(&display_tile)
                .or_else(|e| -> Result<(), PbrtError> {
                    warn!("{:?}", e);
                    return Ok(());
                })
                .unwrap();
        }
    }

    pub fn set_image(&mut self, img: &[Spectrum]) {
        let mut pixels = self.pixels.lock().unwrap();
        let n_pixels = self.cropped_pixel_bounds.area() as usize;
        if n_pixels <= img.len() {
            for i in 0..n_pixels {
                let xyz = img[i].to_xyz();
                pixels[i].xyz = xyz;
                pixels[i].filter_weight_sum = 1.0;
            }
        }
    }

    pub fn add_splat(&mut self, p: &Vector2f, v: &Spectrum) {
        let pi = Point2i::new(p.x.floor() as i32, p.y.floor() as i32);
        if !self.cropped_pixel_bounds.inside_exclusive(&pi) {
            return;
        }

        let pi = Point2i::from((
            pi.x - self.cropped_pixel_bounds.min.x,
            pi.y - self.cropped_pixel_bounds.min.y,
        ));
        let tx = pi.x as usize / ST_W;
        let ty = pi.y as usize / ST_W;
        let xyz = v.to_xyz();
        let x = pi.x as usize - tx * ST_W;
        let y = pi.y as usize - ty * ST_W;

        {
            let tindex = ty * self.splat_size.x as usize + tx;
            if tindex >= self.splat_tiles.len() {
                println!("pi: {:?}, tx: {}, ty: {}", pi, tx, ty);
                println!("cropped_pixel_bounds: {:?}", self.cropped_pixel_bounds);
                println!(
                    "tindex: {}, splat_tiles.len(): {}",
                    tindex,
                    self.splat_tiles.len()
                );

                panic!("tindex >= self.splat_tiles.len()");
            }
            let mut splat_tile = self.splat_tiles[ty * self.splat_size.x as usize + tx]
                .write()
                .unwrap();
            //if xyz[0] != 0.0 || xyz[1] != 0.0 || xyz[2] != 0.0 {
            //println!("Splatting at ({}, {})", x, y);
            //}
            splat_tile.pixels[y * ST_W + x][0] += xyz[0];
            splat_tile.pixels[y * ST_W + x][1] += xyz[1];
            splat_tile.pixels[y * ST_W + x][2] += xyz[2];
            splat_tile.dirty = true;
        }
    }

    pub fn clear(&mut self) {
        let mut pixels = self.pixels.lock().unwrap();
        for index in 0..pixels.len() {
            for i in 0..3 {
                pixels[index].xyz[i] = 0.0;
            }
            pixels[index].filter_weight_sum = 0.0;
        }
    }

    pub fn render_start(&mut self) {
        let resolution = [
            self.full_resolution[0] as usize,
            self.full_resolution[1] as usize,
        ];
        let channel_names = ["R", "G", "B"];
        let _ = self
            .display
            .start(&self.filename, &resolution, &channel_names);
    }

    pub fn render_end(&mut self) {
        let _ = self.display.end();
    }

    pub fn write_image(&self) {
        info!("Converting image to RGB and computing final weighted pixel values");
        // _splat_scale: Float
        let mut rgb = vec![0.0; 3 * self.cropped_pixel_bounds.area() as usize];
        let pixels = self.pixels.lock().unwrap();
        {
            for offset in 0..pixels.len() {
                let pixel = &pixels[offset];
                let mut c = xyz_to_rgb(&pixel.xyz);
                let filter_weight_sum = pixel.filter_weight_sum;
                if filter_weight_sum > 0.0 {
                    let inv_wt = 1.0 / filter_weight_sum;
                    c[0] = f32::max(0.0, c[0] * inv_wt);
                    c[1] = f32::max(0.0, c[1] * inv_wt);
                    c[2] = f32::max(0.0, c[2] * inv_wt);
                }
                rgb[3 * offset + 0] = c[0];
                rgb[3 * offset + 1] = c[1];
                rgb[3 * offset + 2] = c[2];

                let scale = self.scale;
                rgb[3 * offset + 0] *= scale;
                rgb[3 * offset + 1] *= scale;
                rgb[3 * offset + 2] *= scale;
            }
        }
        info!(
            "Writing image {} with bounds {:?}",
            self.filename, self.cropped_pixel_bounds
        );
        let _ = write_image(
            &self.filename,
            &rgb,
            &self.cropped_pixel_bounds,
            &self.full_resolution,
        );
    }

    fn get_pixel_index(&self, p: &Vector2i) -> usize {
        let cmaxx = self.cropped_pixel_bounds.max.x as usize;
        let cminx = self.cropped_pixel_bounds.min.x as usize;
        let cminy = self.cropped_pixel_bounds.min.y as usize;
        let x = p.x as usize;
        let y = p.y as usize;
        let width: usize = cmaxx - cminx;
        let offset: usize = (x - cminx) + (y - cminy) * width;
        return offset;
    }

    pub fn add_display(&mut self, display: &Arc<RwLock<dyn Display>>) {
        self.display.add_display(display);
    }
}

unsafe impl Send for Film {}
unsafe impl Sync for Film {}

fn get_crop_window(params: &ParamSet) -> Result<Bounds2f, PbrtError> {
    if let Some(cropwindow) = params.get_floats_ref("cropwindow") {
        if cropwindow.len() == 4 {
            let crop_window = Bounds2f::from((
                (cropwindow[0], cropwindow[2]),
                (cropwindow[1], cropwindow[3]),
            ));
            Ok(crop_window)
        } else {
            let msg = format!(
                "{} values supplied for \"cropwindow\". Expected 4.",
                cropwindow.len()
            );
            return Err(PbrtError::error(&msg));
        }
    } else {
        let crop_window = Bounds2f::from(((0.0, 0.0), (1.0, 1.0)));
        Ok(crop_window)
    }
}

fn get_full_path(filename: &str) -> String {
    let path = Path::new(filename);
    if path.is_absolute() {
        return String::from(path.to_str().unwrap());
    } else {
        let dir = env::current_dir().unwrap();
        let full_path = dir.join(path);
        return String::from(full_path.to_str().unwrap());
    }
}

pub fn create_film(
    _name: &str,
    params: &ParamSet,
    filter: &Arc<dyn Filter>,
) -> Result<Arc<RwLock<Film>>, PbrtError> {
    let filename = params.find_one_string("filename", "pbrt.exr");
    let filepath = get_full_path(&filename);
    //let filepath = params.find_one_string("filepath", &filename);
    let xres = params.find_one_int("xresolution", 1280);
    let yres = params.find_one_int("yresolution", 720);
    let crop = get_crop_window(params)?;
    let scale = params.find_one_float("scale", 1.0);
    let diagonal = params.find_one_float("diagonal", 35.0);
    let maxsampleluminance = params.find_one_float("maxsampleluminance", f32::INFINITY);
    let resolution = Point2i::from((xres, yres));
    let film = Film::new(
        &resolution,
        &crop,
        filter,
        diagonal,
        &filepath,
        scale,
        maxsampleluminance,
    );
    return Ok(Arc::new(RwLock::new(film)));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filters::*;

    fn create_film_test(t: (i32, i32)) -> Film {
        let resolution = Vector2i::from(t);
        let crop_window = Bounds2f::new(&Vector2f::new(0.0, 0.0), &Vector2f::new(1.0, 1.0));
        let filter: Arc<dyn Filter> = Arc::new(BoxFilter::new(&Vector2f::new(1.0, 1.0)));
        let film = Film::new(&resolution, &crop_window, &filter, 1.0, "test", 1.0, 1000.0);
        return film;
    }

    #[test]
    fn test_001() {
        let mut film = create_film_test((3, 3));
        assert_eq!(film.filename, "test");
        assert_eq!(film.diagonal, 1.0);
        assert_eq!(film.filter.as_ref().get_radius(), Vector2f::new(1.0, 1.0));
        let mut tile =
            film.get_film_tile(&Bounds2i::new(&Vector2i::new(1, 1), &Vector2i::new(2, 2)));
        assert_eq!(tile.pixel_bounds.area(), 9);
        assert_eq!(tile.pixels.len(), 9);

        tile.pixels[4].contrib_sum = Spectrum::from([1.0, 0.0, 0.0]);
        tile.pixels[4].filter_weight_sum = 1.0;
        tile.pixels[5].filter_weight_sum = 2.0;

        film.merge_film_tile(&tile);
        {
            let pixels = film.pixels.lock().unwrap();
            assert_eq!(pixels[4].filter_weight_sum, 1.0);
            assert_ne!(pixels[5].filter_weight_sum, 1.0);
            assert_eq!(pixels[4].xyz, tile.pixels[4].contrib_sum.to_xyz());
        }
    }

    #[test]
    fn test_002() {
        let col1 = Spectrum::from([1.0, 0.0, 0.0]);
        let col2 = Spectrum::from([1.0, 1.0, 0.0]);
        let mut film = create_film_test((5, 5));
        let mut img: Vec<Spectrum> = vec![col1; 5 * 5];
        img[1] = col2;
        film.set_image(&img);
        {
            let pixels = film.pixels.lock().unwrap();
            assert_eq!(pixels[0].xyz, col1.to_xyz());
            assert_eq!(pixels[1].xyz, col2.to_xyz());
        }
        film.clear();
        {
            let pixels = film.pixels.lock().unwrap();
            assert_eq!(pixels[0].xyz[0], 0.0);
            assert_eq!(pixels[0].filter_weight_sum, 0.0);
        }
    }

    #[test]
    fn test_003() {
        /*
        pub full_resolution: Point2i,
        pub diagonal: Float,
        pub filter: Arc<dyn Filter>,
        pub filename: String,
        pub cropped_pixel_bounds: Bounds2i,
        */
        let mut params = ParamSet::new();

        let filename = "pbrt.exr";
        params.add_string("filename", filename);
        params.add_int("xresolution", 1280);
        params.add_int("yresolution", 720);
        params.add_floats("cropwindow", &[0.0, 1.0, 0.0, 1.0]);

        params.add_float("diagonal", 2.0);
        params.add_float("scale", 3.0);

        params.add_float("maxsampleluminance", 12345.0);

        let filter: Arc<dyn Filter> = Arc::new(BoxFilter::new(&Vector2::new(1.0, 1.0)));
        let opt_film = create_film(filename, &params, &filter);
        assert!(opt_film.is_ok());
        let opt_film = opt_film.unwrap();
        let film = opt_film.read().unwrap();
        //assert_eq!(film.filename, filename);
        assert_eq!(film.full_resolution, Vector2i::new(1280, 720));
        assert_eq!(
            film.cropped_pixel_bounds,
            Bounds2i::new(&Vector2i::new(0, 0), &Vector2i::new(1280, 720))
        );
        assert_eq!(film.diagonal, 2.0);
        assert_eq!(film.scale, 3.0);
        assert_eq!(film.max_sample_luminance, 12345.0);
    }
}
