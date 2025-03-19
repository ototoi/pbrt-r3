use crate::core::misc::*;
use crate::core::pbrt::*;

#[inline]
pub fn offset_ray_origin(p: &Point3f, p_error: &Vector3f, n: &Normal3f, w: &Vector3f) -> Point3f {
    let d = Vector3f::dot(&n.abs(), p_error);

    let mut offset = d * n.clone();
    if Vector3f::dot(w, n) < 0.0 {
        offset = -offset;
    }
    let po = *p + offset;
    let mut a = [po.x, po.y, po.z];

    for i in 0..3 {
        if offset[i] > 0.0 {
            a[i] = next_float_up(a[i]);
        } else if offset[i] < 0.0 {
            a[i] = next_float_down(a[i]);
        }
    }
    return Point3f::new(a[0], a[1], a[2]);
}

#[inline]
pub fn face_forward(n: &Vector3f, v: &Vector3f) -> Vector3f {
    if Vector3f::dot(n, v) < 0.0 {
        return *n * -1.0;
    } else {
        return *n;
    }
}

#[inline]
pub fn max_dimension(v: &Vector3f) -> i32 {
    if v.x > v.y {
        if v.x > v.z {
            return 0;
        } else {
            return 2;
        }
    } else {
        if v.y > v.z {
            return 1;
        } else {
            return 2;
        }
    }
}

#[inline]
pub fn permute(v: &Vector3f, x: usize, y: usize, z: usize) -> Vector3f {
    return Vector3f::new(v[x], v[y], v[z]);
}

#[inline]
pub fn max_component(v: &Vector3f) -> Float {
    return Float::max(v.x, Float::max(v.y, v.z));
}

#[inline]
pub fn coordinate_system(v1: &Vector3f) -> (Vector3f, Vector3f) {
    let v2 = if Float::abs(v1.x) > Float::abs(v1.y) {
        Vector3f::new(-v1.z, 0.0, v1.x) / Float::sqrt(v1.x * v1.x + v1.z * v1.z)
    } else {
        Vector3f::new(0.0, v1.z, -v1.y) / Float::sqrt(v1.y * v1.y + v1.z * v1.z)
    };
    let v3 = Vector3f::cross(v1, &v2).normalize();
    return (v2, v3);
}

#[inline]
pub fn spherical_direction(sin_theta: Float, cos_theta: Float, phi: Float) -> Vector3f {
    return Vector3f::new(
        sin_theta * Float::cos(phi),
        sin_theta * Float::sin(phi),
        cos_theta,
    );
}

#[inline]
pub fn spherical_direction_axes(
    sin_theta: Float,
    cos_theta: Float,
    phi: Float,
    x: &Vector3f,
    y: &Vector3f,
    z: &Vector3f,
) -> Vector3f {
    return (sin_theta * Float::cos(phi) * *x)
        + (sin_theta * Float::sin(phi) * *y)
        + (cos_theta * *z);
}

#[inline]
pub fn spherical_theta(v: &Vector3f) -> Float {
    return Float::acos(Float::clamp(v.z, -1.0, 1.0));
}

#[inline]
pub fn spherical_phi(v: &Vector3f) -> Float {
    let p = Float::atan2(v.y, v.x);
    return if p < 0.0 { p + 2.0 * PI } else { p };
}
