use crate::core::pbrt::*;
use std::sync::Arc;
use std::sync::RwLock;

#[derive(Clone)]
pub struct ProxyFilm {
    film: Arc<RwLock<Film>>,
}

impl ProxyFilm {
    pub fn new(film: &Arc<RwLock<Film>>) -> Self {
        ProxyFilm {
            film: Arc::clone(film),
        }
    }

    pub fn get_film_tile(&self, sample_bounds: &Bounds2i) -> FilmTile {
        return self
            .film
            .as_ref()
            .read()
            .unwrap()
            .get_film_tile(sample_bounds);
    }

    pub fn merge_film_tile(&mut self, tile: &FilmTile) {
        let mut film = self.film.as_ref().write().unwrap();
        film.merge_film_tile(tile);
    }

    pub fn merge_splats(&mut self, splat_scale: Float) {
        let mut film = self.film.as_ref().write().unwrap();
        film.merge_splats(splat_scale);
    }

    pub fn add_splat(&mut self, p: &Vector2f, v: &Spectrum) {
        let mut film = self.film.as_ref().write().unwrap();
        film.add_splat(p, v);
    }

    pub fn get_filename(&self) -> String {
        let film = self.film.as_ref().read().unwrap();
        return film.filename.clone();
    }
}

unsafe impl Send for ProxyFilm {}
unsafe impl Sync for ProxyFilm {}
