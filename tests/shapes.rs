// Imported from shapes.cpp

use pbrt_r3::core::prelude::*;
use pbrt_r3::shapes::*;

use std::sync::Arc;

fn p_exp(rng: &mut RNG, exp: Float) -> Float {
    let logu = lerp(rng.uniform_float(), -exp, exp);
    return Float::powf(10.0, logu);
}

fn p_unif(rng: &mut RNG, range: Float) -> Float {
    return lerp(rng.uniform_float(), -range, range);
}

#[test]
fn triangle_watertight() {
    let mut rng = RNG::new_sequence(12111);
    let n_theta = 16 as u32;
    let n_phi = 8 as u32;
    assert!(n_theta >= 3);
    assert!(n_phi >= 4);

    // Make a triangle mesh representing a triangulated sphere (with
    // vertices randomly offset along their normal), centered at the
    // origin.
    let n_vertices = (n_theta * n_phi) as usize;
    let mut vertices = Vec::new();
    for t in 0..n_theta {
        let theta = PI * (t as Float) / ((n_theta - 1) as Float);
        let cos_theta = theta.cos();
        let sin_theta = theta.sin();
        for p in 0..n_phi {
            let phi = 2.0 * PI * (p as Float) / ((n_phi - 1) as Float);
            let mut radius = 1.0;
            // Make sure all of the top and bottom vertices are coincident.
            if t == 0 {
                vertices.push(Point3f::new(0.0, 0.0, radius));
            } else if t == n_theta - 1 {
                vertices.push(Point3f::new(0.0, 0.0, -radius));
            } else if p == n_phi - 1 {
                // Close it up exactly at the end
                let v = vertices[vertices.len() - (n_phi - 1) as usize];
                vertices.push(v);
            } else {
                radius += 5.0 * rng.uniform_float();
                let v = Point3f::zero() + radius * spherical_direction(sin_theta, cos_theta, phi);
                vertices.push(v);
            }
        }
    }
    assert_eq!(vertices.len(), n_vertices);

    let mut indices: Vec<u32> = Vec::new();
    // fan at the top
    let offset = |t: u32, p: u32| -> u32 { t * n_phi + p };
    for p in 0..n_phi - 1 {
        indices.push(offset(0, 0));
        indices.push(offset(1, p));
        indices.push(offset(1, p + 1));
    }

    // quads in the middle rows
    for t in 1..n_theta - 2 {
        for p in 0..n_phi - 1 {
            indices.push(offset(t, p));
            indices.push(offset(t + 1, p));
            indices.push(offset(t + 1, p + 1));

            indices.push(offset(t, p));
            indices.push(offset(t + 1, p + 1));
            indices.push(offset(t, p + 1));
        }
    }

    // fan at bottom
    for p in 0..n_phi - 1 {
        indices.push(offset(n_theta - 1, 0));
        indices.push(offset(n_theta - 2, p));
        indices.push(offset(n_theta - 2, p + 1));
    }

    let identity = Transform::identity();
    let tris = create_triangle_mesh(
        &identity,
        &identity,
        false,
        indices,
        vertices.clone(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        &ParamSet::new(),
    );

    for _ in 0..100000 {
        let u = Point2f::new(rng.uniform_float(), rng.uniform_float());
        let p = Point3f::zero() + uniform_sample_sphere(&u) * 0.5;

        let u = Point2f::new(rng.uniform_float(), rng.uniform_float());
        let ray = Ray::new(&p, &uniform_sample_sphere(&u), Float::INFINITY, 0.0);
        let mut n_hits = 0;
        for tri in tris.iter() {
            if tri.intersect(&ray).is_some() {
                n_hits += 1;
            }
        }
        assert!(n_hits >= 1);

        // Now tougher: shoot directly at a vertex.
        let p_vertex = vertices[rng.uniform_uint32_threshold(vertices.len() as u32) as usize];
        let ray = Ray::new(&p, &(p_vertex - p), Float::INFINITY, 0.0);
        let mut n_hits = 0;
        for tri in tris.iter() {
            if tri.intersect(&ray).is_some() {
                n_hits += 1;
            }
        }
        assert!(n_hits >= 1, "p_vertex: {:?}", p_vertex);
    }
}

fn get_random_triangle<F>(mut value: F) -> Arc<dyn Shape>
where
    F: FnMut() -> Float,
{
    loop {
        let mut v = vec![Point3f::default(); 3];
        for i in 0..3 {
            v[i] = Point3f::new(value(), value(), value());
        }
        if Vector3f::cross(&(v[1] - v[0]), &(v[2] - v[0])).length_squared() < 1e-20 {
            continue;
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
        if !tri_vec.is_empty() {
            return tri_vec[0].clone();
        }
    }
}

#[test]
fn triangle_reintersect() {
    for i in 0..1000 {
        let mut rng = RNG::new_sequence(i);
        let tri = get_random_triangle(|| p_unif(&mut rng, 10.0));
        test_reintersect_convex(tri.as_ref(), &mut rng);

        // Sample a point on the triangle surface to shoot the ray toward.
        let u = Point2f::new(rng.uniform_float(), rng.uniform_float());
        if let Some((p_tri, _pdf)) = tri.sample(&u) {
            // Choose a ray origin.
            let o = Point3f::new(
                p_exp(&mut rng, 8.0),
                p_exp(&mut rng, 8.0),
                p_exp(&mut rng, 8.0),
            );
            let r = Ray::new(&o, &(p_tri.get_p() - o), Float::INFINITY, 0.0);
            if let Some((_t_hit, isect)) = tri.intersect(&r) {
                // Now trace a bunch of rays leaving the intersection point.
                for _ in 0..10000 {
                    // Random direction leaving the intersection point.
                    let u = Point2f::new(rng.uniform_float(), rng.uniform_float());
                    let w = uniform_sample_sphere(&u);
                    let r_out = isect.spawn_ray(&w);
                    assert!(!tri.intersect_p(&r_out));

                    // Choose a random point to trace rays to.
                    let p2 = Point3f::new(
                        p_exp(&mut rng, 8.0),
                        p_exp(&mut rng, 8.0),
                        p_exp(&mut rng, 8.0),
                    );
                    let r_out = isect.spawn_ray_to_point(&p2);

                    assert!(!tri.intersect_p(&r_out));
                    assert!(!tri.intersect(&r_out).is_some());
                }
            }
        }
    }
}

// Checks the closed-form solid angle computation for triangles against a
// Monte Carlo estimate of it.
#[test]
fn triangle_solid_angle() {
    for i in 0..50 {
        let range = 10.0;
        let mut rng = RNG::new_sequence(100 + i);
        let tri = get_random_triangle(|| p_unif(&mut rng, range));

        // Ensure that the reference point isn't too close to the
        // triangle's surface (which makes the Monte Carlo stuff have more
        // variance, thus requiring more samples).
        let mut pc = Point3f::new(
            p_unif(&mut rng, range),
            p_unif(&mut rng, range),
            p_unif(&mut rng, range),
        );
        pc[(rng.uniform_uint32() % 3) as usize] = if rng.uniform_float() > 0.5 {
            -range - 3.0
        } else {
            range + 3.0
        };

        // Compute reference value using Monte Carlo with uniform spherical
        // sampling.
        let count = 512 * 1024;
        let mut hits = 0;
        for j in 0..count {
            let u = Point2f::new(radical_inverse(0, j as u64), radical_inverse(1, j as u64));
            let w = uniform_sample_sphere(&u);
            let ray = Ray::new(&pc, &w, Float::INFINITY, 0.0);
            if tri.intersect_p(&ray) {
                hits += 1;
            }
        }
        let unif_estimate = hits as Float / (count as Float * uniform_sphere_pdf());

        // Now use Triangle::Sample()...
        let inter = Interaction::from((
            pc,
            Normal3f::default(),
            Vector3f::default(),
            Vector3f::new(0.0, 0.0, 1.0),
            0.0,
            MediumInterface::default(),
        ));
        let mut tri_sample_estimate = 0.0;
        for j in 0..count {
            let u = Point2f::new(radical_inverse(0, j as u64), radical_inverse(1, j as u64));
            if let Some((_p, pdf)) = tri.sample_from(&inter, &u) {
                assert!(pdf > 0.0);
                tri_sample_estimate += 1.0 / (count as Float * pdf);
            } else {
                assert!(false);
            }
        }

        // Now make sure that the two computed solid angle values are
        // fairly close.
        // Absolute error for small solid angles, relative for large.
        let error = |a: Float, b: Float| {
            if a.abs() < 1e-4 || b.abs() < 1e-4 {
                return Float::abs(a - b);
            } else {
                return Float::abs((a - b) / b);
            }
        };

        // Don't compare really small triangles, since uniform sampling
        // doesn't get a good estimate for them.
        if tri_sample_estimate > 1e-3 {
            assert!(
                error(tri_sample_estimate, unif_estimate) < 0.1,
                "unif_estimate: {}, tri_sample_estimate: {}, tri index: {}",
                unif_estimate,
                tri_sample_estimate,
                i
            );
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
        assert!(!mesh[0].intersect(&ray).is_some());
    }
}
