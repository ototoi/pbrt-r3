use super::grid::GridDensityMedium;
use super::homogeneous::HomogeneousMedium;
use crate::core::prelude::*;

use std::sync::Arc;

use log::*;

const SUBSURFACE_PARAMETER_TABLE: [(&str, [Float; 3], [Float; 3]); 47] = [
    // From "A Practical Model for Subsurface Light Transport"
    // Jensen, Marschner, Levoy, Hanrahan
    // Proc SIGGRAPH 2001
    ("Apple", [2.29, 2.39, 1.97], [0.0030, 0.0034, 0.046]),
    ("Chicken1", [0.15, 0.21, 0.38], [0.015, 0.077, 0.19]),
    ("Chicken2", [0.19, 0.25, 0.32], [0.018, 0.088, 0.20]),
    ("Cream", [7.38, 5.47, 3.15], [0.0002, 0.0028, 0.0163]),
    ("Ketchup", [0.18, 0.07, 0.03], [0.061, 0.97, 1.45]),
    ("Marble", [2.19, 2.62, 3.0], [0.0021, 0.0041, 0.0071]),
    ("Potato", [0.68, 0.70, 0.55], [0.0024, 0.0090, 0.12]),
    ("Skimmilk", [0.70, 1.22, 1.90], [0.0014, 0.0025, 0.0142]),
    ("Skin1", [0.74, 0.88, 1.01], [0.032, 0.17, 0.48]),
    ("Skin2", [1.09, 1.59, 1.79], [0.013, 0.070, 0.145]),
    ("Spectralon", [11.6, 20.4, 14.9], [0.0, 0.0, 0.0]),
    ("Wholemilk", [2.55, 3.21, 3.77], [0.0011, 0.0024, 0.014]),
    // From "Acquiring Scattering Properties of Participating Media by
    // Dilution",
    // Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen
    // Proc SIGGRAPH 2006
    (
        "Lowfat Milk",
        [0.89187, 1.5136, 2.532],
        [0.002875, 0.00575, 0.0115],
    ),
    (
        "Reduced Milk",
        [2.4858, 3.1669, 4.5214],
        [0.0025556, 0.0051111, 0.012778],
    ),
    (
        "Regular Milk",
        [4.5513, 5.8294, 7.136],
        [0.0015333, 0.0046, 0.019933],
    ),
    (
        "Espresso",
        [0.72378, 0.84557, 1.0247],
        [4.7984, 6.5751, 8.8493],
    ),
    (
        "Mint Mocha Coffee",
        [0.31602, 0.38538, 0.48131],
        [3.772, 5.8228, 7.82],
    ),
    (
        "Lowfat Soy Milk",
        [0.30576, 0.34233, 0.61664],
        [0.0014375, 0.0071875, 0.035937],
    ),
    (
        "Regular Soy Milk",
        [0.59223, 0.73866, 1.4693],
        [0.0019167, 0.0095833, 0.065167],
    ),
    (
        "Lowfat Chocolate Milk",
        [0.64925, 0.83916, 1.1057],
        [0.0115, 0.0368, 0.1564],
    ),
    (
        "Regular Chocolate Milk",
        [1.4585, 2.1289, 2.9527],
        [0.010063, 0.043125, 0.14375],
    ),
    (
        "Coke",
        [8.9053e-05, 8.372e-05, 0.0],
        [0.10014, 0.16503, 0.2468],
    ),
    (
        "Pepsi",
        [6.1697e-05, 4.2564e-05, 0.0],
        [0.091641, 0.14158, 0.20729],
    ),
    (
        "Sprite",
        [6.0306e-06, 6.4139e-06, 6.5504e-06],
        [0.001886, 0.0018308, 0.0020025],
    ),
    (
        "Gatorade",
        [0.0024574, 0.003007, 0.0037325],
        [0.024794, 0.019289, 0.008878],
    ),
    (
        "Chardonnay",
        [1.7982e-05, 1.3758e-05, 1.2023e-05],
        [0.010782, 0.011855, 0.023997],
    ),
    (
        "White Zinfandel",
        [1.7501e-05, 1.9069e-05, 1.288e-05],
        [0.012072, 0.016184, 0.019843],
    ),
    (
        "Merlot",
        [2.1129e-05, 0.0, 0.0],
        [0.11632, 0.25191, 0.29434],
    ),
    (
        "Budweiser Beer",
        [2.4356e-05, 2.4079e-05, 1.0564e-05],
        [0.011492, 0.024911, 0.057786],
    ),
    (
        "Coors Light Beer",
        [5.0922e-05, 4.301e-05, 0.0],
        [0.006164, 0.013984, 0.034983],
    ),
    (
        "Clorox",
        [0.0024035, 0.0031373, 0.003991],
        [0.0033542, 0.014892, 0.026297],
    ),
    (
        "Apple Juice",
        [0.00013612, 0.00015836, 0.000227],
        [0.012957, 0.023741, 0.052184],
    ),
    (
        "Cranberry Juice",
        [0.00010402, 0.00011646, 7.8139e-05],
        [0.039437, 0.094223, 0.12426],
    ),
    (
        "Grape Juice",
        [5.382e-05, 0.0, 0.0],
        [0.10404, 0.23958, 0.29325],
    ),
    (
        "Ruby Grapefruit Juice",
        [0.011002, 0.010927, 0.011036],
        [0.085867, 0.18314, 0.25262],
    ),
    (
        "White Grapefruit Juice",
        [0.22826, 0.23998, 0.32748],
        [0.0138, 0.018831, 0.056781],
    ),
    (
        "Shampoo",
        [0.0007176, 0.0008303, 0.0009016],
        [0.014107, 0.045693, 0.061717],
    ),
    (
        "Strawberry Shampoo",
        [0.00015671, 0.00015947, 1.518e-05],
        [0.01449, 0.05796, 0.075823],
    ),
    (
        "Head & Shoulders Shampoo",
        [0.023805, 0.028804, 0.034306],
        [0.084621, 0.15688, 0.20365],
    ),
    (
        "Lemon Tea Powder",
        [0.040224, 0.045264, 0.051081],
        [2.4288, 4.5757, 7.2127],
    ),
    (
        "Orange Powder",
        [0.00015617, 0.00017482, 0.0001762],
        [0.001449, 0.003441, 0.007863],
    ),
    (
        "Pink Lemonade Powder",
        [0.00012103, 0.00013073, 0.00012528],
        [0.001165, 0.002366, 0.003195],
    ),
    (
        "Cappuccino Powder",
        [1.8436, 2.5851, 2.1662],
        [35.844, 49.547, 61.084],
    ),
    (
        "Salt Powder",
        [0.027333, 0.032451, 0.031979],
        [0.28415, 0.3257, 0.34148],
    ),
    (
        "Sugar Powder",
        [0.00022272, 0.00025513, 0.000271],
        [0.012638, 0.031051, 0.050124],
    ),
    (
        "Suisse Mocha Powder",
        [2.7979, 3.5452, 4.3365],
        [17.502, 27.004, 35.433],
    ),
    (
        "Pacific Ocean Surface Water",
        [0.0001764, 0.00032095, 0.00019617],
        [0.031845, 0.031324, 0.030147],
    ),
];

pub fn get_medium_scattering_properties(name: &str) -> Option<(Spectrum, Spectrum)> {
    if !name.is_empty() {
        for mss in SUBSURFACE_PARAMETER_TABLE {
            if name == mss.0 {
                return Some((Spectrum::from(mss.2), Spectrum::from(mss.1)));
            }
        }
    }
    return None;
}

pub fn create_medium(
    name: &str,
    params: &ParamSet,
    medium_2_world: &Transform,
) -> Option<Arc<dyn Medium>> {
    let mut sig_a = Spectrum::from([0.0011, 0.0024, 0.014]);
    let mut sig_s = Spectrum::from([2.55, 3.21, 3.77]);
    let preset = params.find_one_string("preset", "");
    if let Some((a, b)) = get_medium_scattering_properties(&preset) {
        sig_a = params.find_one_spectrum("sigma_a", &a);
        sig_s = params.find_one_spectrum("sigma_s", &b);
    } else {
        if !preset.is_empty() {
            warn!("Material preset \"{}\" not found.  Using defaults.", preset);
        }
    }
    let scale = params.find_one_float("scale", 1.0);
    let g = params.find_one_float("g", 0.0);
    sig_a = params.find_one_spectrum("sigma_a", &sig_a) * scale;
    sig_s = params.find_one_spectrum("sigma_s", &sig_s) * scale;

    match name {
        "homogeneous" => Some(Arc::new(HomogeneousMedium::new(&sig_a, &sig_s, g))),
        "heterogeneous" => {
            if let Some(data) = params.get_floats_ref("density") {
                let nx = params.find_one_int("nx", 1) as u32;
                let ny = params.find_one_int("ny", 1) as u32;
                let nz = params.find_one_int("nz", 1) as u32;
                let p0 = params.find_one_point3f("p0", &Point3f::new(0.0, 0.0, 0.0));
                let p1 = params.find_one_point3f("p1", &Point3f::new(1.0, 1.0, 1.0));
                if data.len() == (nx * ny * nz) as usize {
                    let data_2_medium = Transform::translate(p0.x, p0.y, p0.z)
                        * Transform::scale(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
                    let t = (*medium_2_world) * data_2_medium;
                    return Some(Arc::new(GridDensityMedium::new(
                        &sig_a, &sig_s, g, nx, ny, nz, &t, &data,
                    )));
                }
            }
            return None;
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_001() {
        let r = get_medium_scattering_properties("Apple");
        assert!(r.is_some());
        let r2 = r.unwrap();
        assert_eq!(r2.0, Spectrum::from([0.0030, 0.0034, 0.046]));
        assert_eq!(r2.1, Spectrum::from([2.29, 2.39, 1.97]));
    }

    #[test]
    fn test_002() {
        let params = ParamSet::new();
        let t = Transform::identity();
        let m = create_medium("hoge", &params, &t);
        assert!(m.is_none());
    }

    #[test]
    fn test_003() {
        /*
        let preset = params.find_one_string("preset", "");
        let g = params.find_one_float("g", 0.0);
        let scale = params.find_one_float("scale", 1.0);
        */
        let mut params = ParamSet::new();
        params.add_floats("g", &[0.1]);
        params.add_floats("scale", &[1.0]);
        params.add_string("preset", "Cream");

        let t = Transform::identity();
        let m = create_medium("homogeneous", &params, &t);
        assert!(m.is_some());
    }
}
