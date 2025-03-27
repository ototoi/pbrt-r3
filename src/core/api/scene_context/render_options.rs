use crate::core::base::*;
use crate::core::light::*;
use crate::core::medium::*;
use crate::core::param_set::*;
use crate::core::primitive::*;

use std::collections::HashMap;
use std::sync::Arc;

/*
Float transformStartTime = 0, transformEndTime = 1;
    std::string FilterName = "box";
    ParamSet FilterParams;
    std::string FilmName = "image";
    ParamSet FilmParams;
    std::string SamplerName = "halton";
    ParamSet SamplerParams;
    std::string AcceleratorName = "bvh";
    ParamSet AcceleratorParams;
    std::string IntegratorName = "path";
    ParamSet IntegratorParams;
    std::string CameraName = "perspective";
    ParamSet CameraParams;
    TransformSet CameraToWorld;
    std::map<std::string, std::shared_ptr<Medium>> namedMedia;
    std::vector<std::shared_ptr<Light>> lights;
    std::vector<std::shared_ptr<Primitive>> primitives;
    std::map<std::string, std::vector<std::shared_ptr<Primitive>>> instances;
    std::vector<std::shared_ptr<Primitive>> *currentInstance = nullptr;
    bool haveScatteringMedia = false;
*/

type PrimitivesMap = HashMap<String, Vec<Arc<dyn Primitive>>>;
type NamedMediumMap = HashMap<String, Arc<dyn Medium>>;

pub struct RenderOptions {
    pub transform_start_time: Float,
    pub transform_end_time: Float,
    pub filter_name: String,
    pub filter_params: ParamSet,
    pub film_name: String,
    pub film_params: ParamSet,
    pub sampler_name: String,
    pub sampler_params: ParamSet,
    pub accelerator_name: String,
    pub accelerator_params: ParamSet,
    pub integrator_name: String,
    pub integrator_params: ParamSet,
    pub camera_name: String,
    pub camera_params: ParamSet,
    //CameraToWorld
    pub named_media: NamedMediumMap,
    pub lights: Vec<Arc<dyn Light>>,
    pub primitives: Vec<Arc<dyn Primitive>>,
    pub instances: PrimitivesMap,
    pub current_instance_name: Option<String>,
    pub have_scattering_media: bool,
}

impl RenderOptions {
    pub fn new() -> Self {
        let accelerator_name =
            std::env::var("PBRT_ACCELERATOR").unwrap_or_else(|_| "bvh".to_string());
        RenderOptions {
            transform_start_time: 0.0,
            transform_end_time: 1.0,
            filter_name: String::from("box"),
            filter_params: ParamSet::new(),
            film_name: String::from("image"),
            film_params: ParamSet::new(),
            sampler_name: String::from("halton"),
            sampler_params: ParamSet::new(),
            accelerator_name: accelerator_name,
            accelerator_params: ParamSet::new(),
            integrator_name: String::from("path"),
            integrator_params: ParamSet::new(),
            camera_name: String::from("perspective"),
            camera_params: ParamSet::new(),
            //
            named_media: NamedMediumMap::new(),
            lights: Vec::new(),
            primitives: Vec::new(),
            instances: PrimitivesMap::new(),
            current_instance_name: None,
            have_scattering_media: false,
        }
    }
}
