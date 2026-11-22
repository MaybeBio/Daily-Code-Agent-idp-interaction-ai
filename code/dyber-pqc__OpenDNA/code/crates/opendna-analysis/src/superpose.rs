//! Optimal rigid-body superposition (Kabsch 1976).

use crate::{shape::centroid, AnalysisError};
use nalgebra::{Matrix3, Vector3};
use opendna_core::Vec3;
use serde::Serialize;

/// A rotation followed by a translation: `x' = R·x + t`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Transform {
    /// Row-major 3×3 rotation matrix.
    pub rotation: [[f64; 3]; 3],
    /// Translation in Å.
    pub translation: Vec3,
}

impl Transform {
    /// Apply to a point.
    pub fn apply(&self, p: Vec3) -> Vec3 {
        let r = &self.rotation;
        Vec3::new(
            r[0][0] * p.x + r[0][1] * p.y + r[0][2] * p.z,
            r[1][0] * p.x + r[1][1] * p.y + r[1][2] * p.z,
            r[2][0] * p.x + r[2][1] * p.y + r[2][2] * p.z,
        ) + self.translation
    }
}

/// Result of superposing `mobile` onto `target`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Superposition {
    /// Transform that maps `mobile` onto `target`.
    pub transform: Transform,
    /// RMSD after superposition, in Å.
    pub rmsd: f64,
}

/// Find the proper rotation and translation minimising RMSD between paired
/// points. Reflections are excluded.
///
/// # Errors
/// [`AnalysisError::LengthMismatch`] or [`AnalysisError::Empty`].
pub fn superpose(mobile: &[Vec3], target: &[Vec3]) -> Result<Superposition, AnalysisError> {
    if mobile.len() != target.len() {
        return Err(AnalysisError::LengthMismatch(mobile.len(), target.len()));
    }
    let cm = centroid(mobile)?;
    let ct = centroid(target)?;
    let v = |p: Vec3| Vector3::new(p.x, p.y, p.z);

    let mut h = Matrix3::<f64>::zeros();
    for (&m, &t) in mobile.iter().zip(target) {
        h += v(m - cm) * v(t - ct).transpose();
    }
    let svd = h.svd(true, true);
    let (Some(u), Some(vt)) = (svd.u, svd.v_t) else {
        return Err(AnalysisError::InvalidParameter("SVD did not converge"));
    };
    let d = (vt.transpose() * u.transpose()).determinant().signum();
    let correction = Matrix3::from_diagonal(&Vector3::new(1.0, 1.0, d));
    let r = vt.transpose() * correction * u.transpose();

    let rotation = [
        [r[(0, 0)], r[(0, 1)], r[(0, 2)]],
        [r[(1, 0)], r[(1, 1)], r[(1, 2)]],
        [r[(2, 0)], r[(2, 1)], r[(2, 2)]],
    ];
    let rotated_cm = Transform { rotation, translation: Vec3::default() }.apply(cm);
    let transform = Transform { rotation, translation: ct - rotated_cm };
    let rmsd = rmsd(&mobile.iter().map(|&p| transform.apply(p)).collect::<Vec<_>>(), target)?;
    Ok(Superposition { transform, rmsd })
}

/// RMSD between paired points without superposition.
///
/// # Errors
/// [`AnalysisError::LengthMismatch`] or [`AnalysisError::Empty`].
pub fn rmsd(a: &[Vec3], b: &[Vec3]) -> Result<f64, AnalysisError> {
    if a.len() != b.len() {
        return Err(AnalysisError::LengthMismatch(a.len(), b.len()));
    }
    if a.is_empty() {
        return Err(AnalysisError::Empty);
    }
    let sum: f64 = a.iter().zip(b).map(|(p, q)| p.distance_squared(*q)).sum();
    Ok((sum / a.len() as f64).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn rotate_z(p: Vec3, theta: f64) -> Vec3 {
        let (s, c) = theta.sin_cos();
        Vec3::new(c * p.x - s * p.y, s * p.x + c * p.y, p.z)
    }

    proptest! {
        #[test]
        fn recovers_rigid_motion(
            pts in prop::collection::vec((-20.0f64..20.0, -20.0f64..20.0, -20.0f64..20.0), 4..40),
            theta in -3.1f64..3.1,
            shift in (-10.0f64..10.0, -10.0f64..10.0, -10.0f64..10.0),
        ) {
            let mobile: Vec<Vec3> = pts.iter().map(|&(x, y, z)| Vec3::new(x, y, z)).collect();
            let t = Vec3::new(shift.0, shift.1, shift.2);
            let target: Vec<Vec3> = mobile.iter().map(|&p| rotate_z(p, theta) + t).collect();
            let s = superpose(&mobile, &target).unwrap();
            prop_assert!(s.rmsd < 1e-6, "rmsd {}", s.rmsd);
        }
    }

    #[test]
    fn does_not_reflect() {
        // A mirror image cannot be superposed exactly by a proper rotation.
        let a = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let mirrored: Vec<Vec3> = a.iter().map(|p| Vec3::new(p.x, p.y, -p.z)).collect();
        assert!(superpose(&a, &mirrored).unwrap().rmsd > 0.1);
    }
}
