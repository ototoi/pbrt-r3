#[inline]
pub fn float_to_bits(f: f32) -> u32 {
    f.to_bits()
}

#[inline]
pub fn bits_to_float(ui: u32) -> f32 {
    f32::from_bits(ui)
}

pub fn next_float_down(x: f32) -> f32 {
    let mut v = x;
    if f32::is_infinite(v) && v < 0.0 {
        return v;
    }
    if v == 0.0 {
        v = -0.0;
    }

    let mut ui = float_to_bits(v);
    if v > 0.0 {
        ui = ui.saturating_sub(1);
    } else {
        ui = ui.saturating_add(1);
    }
    return bits_to_float(ui);
}

#[inline]
pub fn next_float_up(x: f32) -> f32 {
    let mut v = x;
    if f32::is_infinite(v) && v > 0.0 {
        return v;
    }
    if v == -0.0 {
        v = 0.0;
    }
    let mut ui = float_to_bits(v);
    if v >= 0.0 {
        ui = ui.saturating_add(1);
    } else {
        ui = ui.saturating_sub(1);
    }
    return bits_to_float(ui);
}
