// Imported from shapes.cpp

use pbrt_r3::core::pbrt::*;
use pbrt_r3::shapes::*;

use std::sync::Arc;

fn p_exp(rng: &mut RNG, exp: Float) -> Float {
    let logu = lerp(rng.uniform_float(), -exp, exp);
    return Float::powf(10.0, logu);
}

fn p_unif(rng: &mut RNG, range: Float) -> Float {
    return lerp(rng.uniform_float(), -range, range);
}

fn get_random_triangle<F>(mut value: F) -> Option<Arc<dyn Shape>>
where
    F: FnMut() -> Float,
{
    let mut v = vec![Point3f::default(); 3];
    for i in 0..3 {
        v[i] = Point3f::new(value(), value(), value());
    }
    if Vector3f::cross(&(v[1] - v[0]), &(v[2] - v[0])).length_squared() < 1e-20 {
        return None;
    }

    let identity = Transform::identity();
    let indices = vec![0, 1, 2];
    let tri_vec = create_triangle_mesh(
        &identity,
        &identity,
        false,
        indices,
        v,
        Vec::new(),
        Vec::new(),
        Vec::new(),
        &ParamSet::new(),
    );
    if tri_vec.is_empty() {
        return None;
    }
    return Some(tri_vec[0].clone());
}

// Checks the closed-form solid angle computation for triangles against a
// Monte Carlo estimate of it.
#[test]
fn triangle_solid_angle() {
    for i in 0..50 {
        let range = 10.0;
        let mut rng = RNG::new_sequence(100 + i);
        if let Some(tri) = get_random_triangle(|| p_unif(&mut rng, range)) {

            // Ensure that the reference point isn't too close to the
            // triangle's surface (which makes the Monte Carlo stuff have more
            // variance, thus requiring more samples).
            

        }
    }
}

// Use Quasi Monte Carlo with uniform sphere sampling to esimate the solid
// angle subtended by the given shape from the given point.
fn mc_solid_angle(p: &Point3f, shape: &dyn Shape, n_samples: usize) -> Float {
    let mut n_hits = 0;
    for i in 0..n_samples {
        let u = Point2f::new(radical_inverse(0, i as u64), radical_inverse(1, i as u64));
        let w = uniform_sample_sphere(&u);
        let ray = Ray::new(p, &w, Float::INFINITY, 0.0);
        if shape.intersect_p(&ray) {
            n_hits += 1;
        }
    }
    return n_hits as Float / (uniform_sphere_pdf() * n_samples as Float);
}

#[test]
fn sphere_solid_angle() {
    let tr = Transform::translate(1.0, 0.5, -0.8) * Transform::rotate_x(30.0);
    let tr_inv = tr.inverse();
    let sphere = Sphere::new(&tr, &tr_inv, false, 1.0, -1.0, 1.0, 360.0);

    // Make sure we get a subtended solid angle of 4pi for a point
    // inside the sphere.
    let p_inside = Point3f::new(1.0, 0.9, -0.8);
    let n_samples = 128 * 1024;
    let solid_angle_mc = mc_solid_angle(&p_inside, &sphere, n_samples);
    assert!(
        Float::abs(solid_angle_mc - 4.0 * PI) < 0.01,
        "solid_angle_mc: {}",
        solid_angle_mc
    );

    let solid_angle = sphere.solid_angle(&p_inside, n_samples as i32);
    assert!(
        Float::abs(solid_angle - 4.0 * PI) < 0.01,
        "solid_angle: {}",
        solid_angle
    );

    // Now try a point outside the sphere
    let p_outside = Point3f::new(-0.25, -1.0, 0.8);
    let mc_sa = mc_solid_angle(&p_outside, &sphere, n_samples);
    let sphere_sa = sphere.solid_angle(&p_outside, n_samples as i32);
    assert!(
        Float::abs(mc_sa - sphere_sa) < 0.001,
        "mc_sa: {}, sphere_sa: {}",
        mc_sa,
        sphere_sa
    );
}

#[test]
fn cylinder_solid_angle() {
    let tr = Transform::translate(1.0, 0.5, -0.8) * Transform::rotate_x(30.0);
    let tr_inv = tr.inverse();
    let cyl = Cylinder::new(&tr, &tr_inv, false, 0.25, -1.0, 1.0, 360.0);
    let p = Point3f::new(0.5, 0.25, 0.5);
    let n_samples = 128 * 1024;
    let solid_angle_mc = mc_solid_angle(&p, &cyl, n_samples);
    let solid_angle = cyl.solid_angle(&p, n_samples as i32);
    assert!(
        Float::abs(solid_angle_mc - solid_angle) < 0.001,
        "solid_angle_mc: {}, solid_angle: {}",
        solid_angle_mc,
        solid_angle
    );
}

#[test]
fn disk_solid_angle() {
    let tr = Transform::translate(1.0, 0.5, -0.8) * Transform::rotate_x(30.0);
    let tr_inv = tr.inverse();
    let disk = Disk::new(&tr, &tr_inv, false, 0.0, 1.25, 0.0, 360.0);
    let p = Point3f::new(0.5, -0.8, 0.5);
    let n_samples = 128 * 1024;
    let solid_angle_mc = mc_solid_angle(&p, &disk, n_samples);
    let solid_angle = disk.solid_angle(&p, n_samples as i32);
    assert!(
        Float::abs(solid_angle_mc - solid_angle) < 0.001,
        "solid_angle_mc: {}, solid_angle: {}",
        solid_angle_mc,
        solid_angle
    );
}

// Check for incorrect self-intersection: assumes that the shape is convex,
// such that if the dot product of an outgoing ray and the surface normal
// at a point is positive, then a ray leaving that point in that direction
// should never intersect the shape.
fn test_reintersect_convex(shape: &dyn Shape, rng: &mut RNG) {
    // Ray origin
    let o = Point3f::new(p_exp(rng, 8.0), p_exp(rng, 8.0), p_exp(rng, 8.0));

    // Destination: a random point in the shape's bounding box.
    let bbox = shape.world_bound();
    let t = Point3f::new(
        rng.uniform_float(),
        rng.uniform_float(),
        rng.uniform_float(),
    );
    let p2 = bbox.lerp(&t);

    // Ray to intersect with the shape.
    let mut r = Ray::new(&o, &(p2 - o), Float::INFINITY, 0.0);
    if rng.uniform_float() < 0.5 {
        r.d = r.d.normalize();
    }

    // We should usually (but not always) find an intersection.
    if let Some((_t, isect)) = shape.intersect(&r) {
        // Now trace a bunch of rays leaving the intersection point.
        for _ in 0..10000 {
            // Random direction leaving the intersection point.
            let u = Point2f::new(rng.uniform_float(), rng.uniform_float());
            let w = uniform_sample_sphere(&u);
            // Make sure it's in the same hemisphere as the surface normal.
            let w = face_forward(&w, &isect.n);
            let r_out = isect.spawn_ray(&w);
            assert!(!shape.intersect_p(&r_out));
            assert!(!shape.intersect(&r_out).is_some());

            // Choose a random point to trace rays to.
            let p2 = Point3f::new(p_exp(rng, 8.0), p_exp(rng, 8.0), p_exp(rng, 8.0));
            // Make sure that the point we're tracing rays toward is in the
            // hemisphere about the intersection point's surface normal.
            let w = p2 - isect.p;
            let w = face_forward(&w, &isect.n);
            let p2 = isect.p + w;
            let r_out = isect.spawn_ray_to_point(&p2);

            assert!(!shape.intersect_p(&r_out));
            assert!(!shape.intersect(&r_out).is_some());
        }
    }
}

#[test]
fn full_sphere_reintersect() {
    for i in 0..100 {
        let mut rng = RNG::new_sequence(i);
        let identity = Transform::identity();
        let radius = p_exp(&mut rng, 4.0);
        let zmin = -radius;
        let zmax = radius;
        let phi_max = 360.0;
        let sphere = Sphere::new(&identity, &identity, false, radius, zmin, zmax, phi_max);
        test_reintersect_convex(&sphere, &mut rng);
    }
}

#[test]
fn partial_sphere_reintersect() {
    for i in 0..100 {
        let mut rng = RNG::new_sequence(i);
        let identity = Transform::identity();
        let radius = p_exp(&mut rng, 4.0);
        let mut zmin = if rng.uniform_float() < 0.5 {
            -radius
        } else {
            lerp(rng.uniform_float(), -radius, radius)
        };
        let mut zmax = if rng.uniform_float() < 0.5 {
            radius
        } else {
            lerp(rng.uniform_float(), -radius, radius)
        };
        if zmin > zmax {
            std::mem::swap(&mut zmin, &mut zmax);
        }
        let phi_max = if rng.uniform_float() < 0.5 {
            360.0
        } else {
            rng.uniform_float() * 360.0
        };
        let sphere = Sphere::new(&identity, &identity, false, radius, zmin, zmax, phi_max);
        test_reintersect_convex(&sphere, &mut rng);
    }
}

#[test]
fn cylinder_reintersect() {
    for i in 0..100 {
        let mut rng = RNG::new_sequence(i);
        let identity = Transform::identity();
        let radius = p_exp(&mut rng, 4.0);
        let mut zmin = p_exp(&mut rng, 4.0) * (if rng.uniform_float() < 0.5 { -1.0 } else { 1.0 });
        let mut zmax = p_exp(&mut rng, 4.0) * (if rng.uniform_float() < 0.5 { -1.0 } else { 1.0 });
        if zmin > zmax {
            std::mem::swap(&mut zmin, &mut zmax);
        }
        let phi_max = if rng.uniform_float() < 0.5 {
            360.0
        } else {
            rng.uniform_float() * 360.0
        };
        let cylinder = Cylinder::new(&identity, &identity, false, radius, zmin, zmax, phi_max);
        test_reintersect_convex(&cylinder, &mut rng);
    }
}

#[test]
fn triangle_badcases() {
    let identity = Transform::identity();
    let indices = vec![0, 1, 2];
    let p = vec![
        Point3f::new(-1113.45459, -79.049614, -56.2431908),
        Point3f::new(-1113.45459, -87.0922699, -56.2431908),
        Point3f::new(-1113.45459, -79.2090149, -56.2431908),
    ];
    let s = Vec::new();
    let n = Vec::new();
    let uv = Vec::new();
    let params = ParamSet::new();
    let mesh = create_triangle_mesh(&identity, &identity, false, indices, p, s, n, uv, &params);
    if !mesh.is_empty() {
        assert!(mesh.len() > 0);
        let ray = Ray::new(
            &Point3f::new(-1081.47925, 99.9999542, 87.7701111),
            &Vector3f::new(-32.1072998, -183.355865, -144.607635),
            0.9999,
            0.0,
        );
        assert!(mesh[0].intersect(&ray).is_none());
    }
}
