// Imported from animatedtransform.cpp

use pbrt_r3::core::pbrt::*;

use pbrt_r3::core::transform::animated_transform::decompose::decompose;

fn nearly_equal(a: Float, b: Float, e: Float) -> bool {
    Float::abs(a - b) < e
}

fn nearly_equal_matrix(a: &Matrix4x4, b: &Matrix4x4, e: Float) -> bool {
    for i in 0..16 {
        if !nearly_equal(a.m[i], b.m[i], e) {
            return false;
        }
    }
    return true;
}

fn invertible(m: &Matrix4x4) -> bool {
    let m = *m;
    if let Some(im) = m.inverse() {
        let ident = im * m;
        if !nearly_equal_matrix(&ident, &Matrix4x4::identity(), 1e-3) {
            println!("ident: {:?}", ident);
            println!("Matrix4x4::identity(): {:?}", Matrix4x4::identity());
            return false;
        } else {
            return true;
        }
    } else {
        false
    }
}

fn random_transform(rng: &mut RNG) -> Transform {
    let mut t = Transform::new();
    let r = |rng: &mut RNG| -10.0 + 20.0 * rng.uniform_float();
    for _ in 0..10 {
        match rng.uniform_uint32_threshold(3) {
            0 => {
                t = t * Transform::scale(r(rng).abs(), r(rng).abs(), r(rng).abs());
                //t = t * Transform::scale(1.0 + rng.uniform_float(), 1.0 + rng.uniform_float(), 1.0 + rng.uniform_float());
            }
            1 => {
                t = t * Transform::translate(r(rng), r(rng), r(rng));
            }
            2 => {
                let theta = r(rng) * 20.0;
                let uv = Point2f::new(rng.uniform_float(), rng.uniform_float());
                let axis = uniform_sample_sphere(&uv);
                //println!("axis: {:?}", axis);
                assert!((axis.length() - 1.0).abs() < 1e-3);
                t = t * Transform::rotate(theta, axis.x, axis.y, axis.z);
            }
            _ => unreachable!(),
        }
    }

    {
        let r = |rng: &mut RNG| -10.0 + 20.0 * rng.uniform_float();

        let ss = Transform::scale(r(rng).abs(), r(rng).abs(), r(rng).abs());
        let tt = Transform::translate(r(rng), r(rng), r(rng));

        let theta = r(rng) * 20.0;
        let uv = Point2f::new(rng.uniform_float(), rng.uniform_float());
        let axis = uniform_sample_sphere(&uv).normalize();
        //println!("axis: {:?}", axis);
        assert!((axis.length() - 1.0).abs() < 1e-3);
        //let qq = Quaternion::from_axis_angle(&axis, theta);
        let q0 = Quaternion::from_angle_axis(theta, &axis);
        let rr0 = Transform::from(q0.to_matrix());
        let rr = Transform::rotate(theta, axis.x, axis.y, axis.z);
        let q1 = Quaternion::from(rr.m);
        let rr2 = Transform::from(q1.to_matrix());
        let q2 = Quaternion::from(rr2.m);
        println!("theta: {}", theta);
        println!("axis: {:?}", axis);
        println!("rr0: {:?}", rr0.m);
        println!("rr1: {:?}", rr.m);
        println!("rr2: {:?}", rr2.m);
        println!("q0: {:?}, {:?}", q0, q0.to_angle_axis());
        println!("q1: {:?}, {:?}", q1, q1.to_angle_axis());
        println!("q2: {:?}, {:?}", q2, q2.to_angle_axis());

        assert!(invertible(&ss.m));
        assert!(invertible(&tt.m));
        assert!(invertible(&rr.m));

        let m = tt * rr * ss;
        let (nt, nr, ns) = decompose(&m.m, 1e-6, 100).unwrap();
        let mt = Matrix4x4::translate(nt.x, nt.y, nt.z);
        let mr = nr.to_matrix();
        let ms = ns;
        println!("mt: {:?}", mt);
        println!("tt: {:?}", tt.m);
        println!("mr: {:?}", mr);
        println!("rr: {:?}", rr.m);
        println!("ms: {:?}", ms);
        println!("ss: {:?}", ss.m);
        //assert!(nearly_equal_matrix(&mt, &tt.m, 1e-3));
        //assert!(nearly_equal_matrix(&mr, &rr.m, 1e-3));
        //assert!(nearly_equal_matrix(&ms, &ss.m, 1e-3));
    }
    return t;
}

#[test]
fn animated_transform_randoms() {
    let mut rng = RNG::new();
    let r = |rng: &mut RNG| -10.0 + 20.0 * rng.uniform_float();

    for _ in 0..200 {
        // Generate a pair of random transformation matrices.
        let t0 = random_transform(&mut rng);
        let t1 = random_transform(&mut rng);
        let at = AnimatedTransform::new(&t0, 0.0, &t1, 1.0);

        for _ in 0..5 {
            // Generate a pair of random transformation matrices.
            let bounds = Bounds3f::from((
                (r(&mut rng), r(&mut rng), r(&mut rng)),
                (r(&mut rng), r(&mut rng), r(&mut rng)),
            ));
            //let bounds = Bounds3f::from(((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)));
            assert!(bounds.diagonal().x >= 0.0);
            assert!(bounds.diagonal().y >= 0.0);
            assert!(bounds.diagonal().z >= 0.0);

            let motion_bounds = at.motion_bounds(&bounds);

            let mut t = 0.0;
            while t <= 1.0 {
                // Now, interpolate the transformations at a bunch of times
                // along the time range and then transform the bounding box
                // with the result.
                assert!(t >= 0.0 && t <= 1.0);
                let tr = at.interpolate(t);
                let mut tb = tr.transform_bounds(&bounds);

                // Add a little slop to allow for floating-point round-off
                // error in computing the motion extrema times.
                let diagonal = tb.diagonal();
                tb.min += 1e-2 * diagonal;
                tb.max -= 1e-2 * diagonal;

                println!("has_rotation: {}, t: {}", at.has_rotation, t);

                println!(
                    "{}/{} .. {}/{}",
                    motion_bounds.min.x, tb.min.x, tb.max.x, motion_bounds.max.x
                );
                println!(
                    "{}/{} .. {}/{}",
                    motion_bounds.min.y, tb.min.y, tb.max.y, motion_bounds.max.y
                );
                println!(
                    "{}/{} .. {}/{}",
                    motion_bounds.min.z, tb.min.z, tb.max.z, motion_bounds.max.z
                );

                println!("t0: {:?}", t0.m);
                println!("t1: {:?}", t1.m);
                println!("tr: {:?}", tr.m);
                println!();

                // Now, the transformed bounds should be inside the motion
                // bounds.
                assert!(tb.min.x >= motion_bounds.min.x);
                assert!(tb.min.y >= motion_bounds.min.y);
                assert!(tb.min.z >= motion_bounds.min.z);

                assert!(tb.max.x <= motion_bounds.max.x);
                assert!(tb.max.y <= motion_bounds.max.y);
                assert!(tb.max.z <= motion_bounds.max.z);

                t += 1e-3 * rng.uniform_float();
                //t += 0.999999;
            }
        }
    }
}
