use super::decompose::*;
use super::derivatives::*;
use super::interval::*;
use crate::core::pbrt::*;

#[derive(Debug, PartialEq, Clone)]
pub struct AnimatedTransform {
    pub transforms: [Transform; 2],
    pub times: [Float; 2],
    pub actually_animated: bool,
    pub has_rotation: bool,
    pub t: [Vector3f; 2],
    pub r: [Quaternion; 2],
    pub s: [Vector3f; 2],
    pub derivatives: Option<[[DerivativeTerm; 3]; 5]>,
}

impl AnimatedTransform {
    pub fn new(
        start_transform: &Transform,
        start_time: Float,
        end_transform: &Transform,
        end_time: Float,
    ) -> Self {
        const EPS: f32 = f32::EPSILON * 1e+2;
        let transforms = [*start_transform, *end_transform];
        let times = [start_time, end_time];
        let actually_animated = start_transform != end_transform;
        let (t0, r0, s0) = decompose(&start_transform.m, EPS, 100).unwrap();
        let (t1, r1, s1) = decompose(&end_transform.m, EPS, 100).unwrap();
        let t = [t0, t1];
        let mut r = [r0, r1];
        let s = [s0, s1];
        if Quaternion::dot(&r[0], &r[1]) < 0.0 {
            r[1] = -r[1];
        }
        let has_rotation = Quaternion::dot(&r[0], &r[1]) < 0.9995;
        let derivatives = if has_rotation {
            let s = [
                Matrix4x4::scale(s0.x, s0.y, s0.z),
                Matrix4x4::scale(s1.x, s1.y, s1.z),
            ];
            Some(get_derivatives(&t, &r, &s))
        } else {
            None
        };
        AnimatedTransform {
            transforms,
            times,
            actually_animated,
            has_rotation,
            t,
            r,
            s,
            derivatives,
        }
    }

    pub fn interpolate(&self, time: Float) -> Transform {
        if !self.actually_animated || time <= self.times[0] {
            return self.transforms[0];
        }

        if self.times[1] <= time {
            return self.transforms[1];
        }

        let dt = (time - self.times[0]) / (self.times[1] - self.times[0]);
        //println!("{:?}",dt);
        let trans = (1.0 - dt) * self.t[0] + dt * self.t[1];
        let rotate = Quaternion::slerp(dt, &self.r[0], &self.r[1]);
        let scale = (1.0 - dt) * self.s[0] + dt * self.s[1];
        let m = Matrix4x4::translate(trans.x, trans.y, trans.z)
            * rotate.to_matrix()
            * Matrix4x4::scale(scale.x, scale.y, scale.z);
        return Transform::from(m);
    }

    pub fn transform_point(&self, time: Float, p: &Point3f) -> Point3f {
        if !self.actually_animated || time <= self.times[0] {
            return self.transforms[0].transform_point(p);
        }

        if self.times[1] <= time {
            return self.transforms[1].transform_point(p);
        }

        let m = self.interpolate(time);
        return m.transform_point(p);
    }

    pub fn transform_vector(&self, time: Float, v: &Vector3f) -> Point3f {
        if !self.actually_animated || time <= self.times[0] {
            return self.transforms[0].transform_vector(v);
        }

        if self.times[1] <= time {
            return self.transforms[1].transform_vector(v);
        }

        let m = self.interpolate(time);
        return m.transform_vector(v);
    }

    pub fn transform_normal(&self, time: Float, n: &Normal3f) -> Point3f {
        if !self.actually_animated || time <= self.times[0] {
            return self.transforms[0].transform_normal(n);
        }

        if self.times[1] <= time {
            return self.transforms[1].transform_normal(n);
        }

        let m = self.interpolate(time);
        return m.transform_normal(n);
    }

    pub fn transform_ray(&self, r: &Ray) -> (Ray, Vector3f, Vector3f) {
        if !self.actually_animated || r.time <= self.times[0] {
            return self.transforms[0].transform_ray(r);
        } else if self.times[1] <= r.time {
            return self.transforms[1].transform_ray(r);
        } else {
            let t = self.interpolate(r.time);
            return t.transform_ray(r);
        }
    }

    pub fn transform_ray_differential(
        &self,
        r: &RayDifferential,
    ) -> (RayDifferential, Vector3f, Vector3f) {
        if !self.actually_animated || r.ray.time <= self.times[0] {
            return self.transforms[0].transform_ray_differential(r);
        } else if r.ray.time >= self.times[1] {
            return self.transforms[1].transform_ray_differential(r);
        } else {
            let t = self.interpolate(r.ray.time);
            return t.transform_ray_differential(r);
        }
    }

    pub fn motion_bounds(&self, b: &Bounds3f) -> Bounds3f {
        if !self.actually_animated {
            return self.transforms[0].transform_bounds(b);
        }
        let b0 = self.transforms[0].transform_bounds(b);
        let b1 = self.transforms[1].transform_bounds(b);
        let mut bounds = Bounds3f::union(&b0, &b1);
        if !self.has_rotation {
            return bounds;
        } else {
            if true {
                let count = 64;
                (1..count).into_iter().map(|i| {
                    let t = lerp(i as Float / count as Float, self.times[0], self.times[1]);
                    let tr = self.interpolate(t);
                    tr.transform_bounds(b)
                }).fold(bounds, |acc, x| Bounds3f::union(&acc, &x));
            } else {
                // Return motion bounds accounting for animated rotation
                for corner in 0..8 {
                    bounds = bounds.union(&self.bound_point_motion(&b.corner(corner)));
                }
            }
            return bounds;
        }
    }

    pub fn bound_point_motion(&self, p: &Point3f) -> Bounds3f {
        if !self.actually_animated {
            let p = self.transforms[0].transform_point(p);
            return Bounds3f::from((p.x, p.y, p.z));
        }
        if !self.has_rotation {
            let p0 = self.transforms[0].transform_point(p);
            let p1 = self.transforms[1].transform_point(p);
            let b0 = Bounds3f::from((p0.x, p0.y, p0.z));
            let b1 = Bounds3f::from((p1.x, p1.y, p1.z));
            return Bounds3f::union(&b0, &b1);
        } else {
            let derivatives = self.derivatives.as_ref().unwrap();
            let p0 = self.transforms[0].transform_point(p);
            let p1 = self.transforms[1].transform_point(p);
            let mut bounds = Bounds3f::new(&p0, &p1);
            let cos_theta = Quaternion::dot(&self.r[0], &self.r[1]);
            let theta = Float::acos(Float::clamp(cos_theta, -1.0, 1.0));
            let max_depth = 8;
            for c in 0..3 {
                let c1 = derivatives[0][c].eval(p);
                let c2 = derivatives[1][c].eval(p);
                let c3 = derivatives[2][c].eval(p);
                let c4 = derivatives[3][c].eval(p);
                let c5 = derivatives[4][c].eval(p);
                let zeros = interval_find_zeros(
                    c1,
                    c2,
                    c3,
                    c4,
                    c5,
                    theta,
                    Interval::new(0.0, 1.0),
                    max_depth,
                );
                for zero in zeros {
                    let t = lerp(zero, self.times[0], self.times[1]);
                    let pz = self.transform_point(t, p);
                    bounds = bounds.union_p(&pz);
                }
            }
            return bounds;
        }
    }

    pub fn is_animated(&self) -> bool {
        return self.actually_animated;
    }

    pub fn has_scale(&self) -> bool {
        return self.transforms[0].has_scale() || self.transforms[1].has_scale();
    }
}
