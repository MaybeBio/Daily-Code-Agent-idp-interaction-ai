//! Global shape descriptors.

use crate::AnalysisError;
use opendna_core::Vec3;

/// Arithmetic mean of a set of points.
///
/// # Errors
/// [`AnalysisError::Empty`] for no points.
pub fn centroid(points: &[Vec3]) -> Result<Vec3, AnalysisError> {
    if points.is_empty() {
        return Err(AnalysisError::Empty);
    }
    let sum = points.iter().fold(Vec3::default(), |acc, &p| acc + p);
    Ok(sum * (1.0 / points.len() as f64))
}

/// Unweighted radius of gyration, in Å. Pass Cα positions for the usual
/// protein Rg, or all heavy atoms for a mass-agnostic variant.
///
/// # Errors
/// [`AnalysisError::Empty`] for no points.
pub fn radius_of_gyration(points: &[Vec3]) -> Result<f64, AnalysisError> {
    let c = centroid(points)?;
    let msd = points.iter().map(|p| p.distance_squared(c)).sum::<f64>() / points.len() as f64;
    Ok(msd.sqrt())
}
