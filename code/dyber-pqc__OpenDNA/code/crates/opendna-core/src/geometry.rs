use serde::{Deserialize, Serialize};
use std::ops::{Add, Mul, Neg, Sub};

/// A 3D vector / point in Ångström.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct Vec3 {
    /// x coordinate
    pub x: f64,
    /// y coordinate
    pub y: f64,
    /// z coordinate
    pub z: f64,
}

impl Vec3 {
    /// Construct a vector.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Dot product.
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product.
    #[must_use]
    pub fn cross(self, o: Self) -> Self {
        Self::new(self.y * o.z - self.z * o.y, self.z * o.x - self.x * o.z, self.x * o.y - self.y * o.x)
    }

    /// Euclidean length.
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// Unit vector in the same direction. Returns the zero vector unchanged.
    #[must_use]
    pub fn normalized(self) -> Self {
        let n = self.norm();
        if n == 0.0 {
            self
        } else {
            self * (1.0 / n)
        }
    }

    /// Distance to another point.
    pub fn distance(self, o: Self) -> f64 {
        (self - o).norm()
    }

    /// Squared distance (avoids a square root in neighbour searches).
    pub fn distance_squared(self, o: Self) -> f64 {
        let d = self - o;
        d.dot(d)
    }

    /// Angle at `self` formed by `a`-`self`-`b`, in degrees.
    pub fn angle(a: Self, vertex: Self, b: Self) -> f64 {
        let u = (a - vertex).normalized();
        let v = (b - vertex).normalized();
        u.dot(v).clamp(-1.0, 1.0).acos().to_degrees()
    }

    /// IUPAC dihedral angle p0-p1-p2-p3 in degrees, in the range (-180, 180].
    pub fn dihedral(p0: Self, p1: Self, p2: Self, p3: Self) -> f64 {
        let b0 = p0 - p1;
        let b1 = (p2 - p1).normalized();
        let b2 = p3 - p2;
        // Project b0 and b2 onto the plane perpendicular to b1.
        let v = b0 - b1 * b0.dot(b1);
        let w = b2 - b1 * b2.dot(b1);
        let x = v.dot(w);
        let y = b1.cross(v).dot(w);
        y.atan2(x).to_degrees()
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn right_angle() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let o = Vec3::default();
        let b = Vec3::new(0.0, 1.0, 0.0);
        assert!(close(Vec3::angle(a, o, b), 90.0));
    }

    #[test]
    fn dihedral_signs_follow_iupac_convention() {
        let p0 = Vec3::new(1.0, 0.0, 0.0);
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(0.0, 1.0, 0.0);
        assert!(close(Vec3::dihedral(p0, p1, p2, Vec3::new(1.0, 1.0, 0.0)), 0.0));
        assert!(close(Vec3::dihedral(p0, p1, p2, Vec3::new(-1.0, 1.0, 0.0)), 180.0));
        assert!(close(Vec3::dihedral(p0, p1, p2, Vec3::new(0.0, 1.0, -1.0)), 90.0));
        assert!(close(Vec3::dihedral(p0, p1, p2, Vec3::new(0.0, 1.0, 1.0)), -90.0));
    }
}
