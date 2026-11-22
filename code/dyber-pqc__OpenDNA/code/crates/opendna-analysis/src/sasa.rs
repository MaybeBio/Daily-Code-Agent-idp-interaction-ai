//! Solvent-accessible surface area by the Shrake–Rupley method (1973).
//!
//! Each atom is a sphere of its van der Waals radius plus the probe radius,
//! sampled with a Fibonacci-lattice point set. A point is accessible if no
//! neighbouring sphere covers it. Hydrogens and waters are excluded.

use crate::AnalysisError;
use opendna_core::{AminoAcid, ResidueId, Structure, Vec3};
use serde::Serialize;
use std::collections::HashMap;

/// Parameters for [`compute`].
#[derive(Debug, Clone, Copy)]
pub struct SasaParams {
    /// Probe radius in Å (water ≈ 1.4).
    pub probe_radius: f64,
    /// Sample points per atom. 960 gives errors well under 1 %.
    pub points: usize,
    /// Include HETATM residues other than water.
    pub include_hetero: bool,
}

impl Default for SasaParams {
    fn default() -> Self {
        Self { probe_radius: 1.4, points: 960, include_hetero: false }
    }
}

/// Accessible area of one residue.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ResidueSasa {
    /// Chain id.
    pub chain: String,
    /// Depositor numbering.
    pub id: ResidueId,
    /// Residue name.
    pub name: String,
    /// Absolute area in Å².
    pub area: f64,
    /// Area divided by the theoretical maximum of Tien et al. (2013), when
    /// the residue is a standard amino acid.
    pub relative: Option<f64>,
}

/// Total and per-residue SASA.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Sasa {
    /// Sum over all atoms, Å².
    pub total: f64,
    /// Per-residue breakdown in structure order.
    pub residues: Vec<ResidueSasa>,
}

/// Bondi (1964) van der Waals radii in Å; unknown elements default to 1.8.
fn vdw_radius(element: &str) -> f64 {
    match element {
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "SE" => 1.90,
        "F" => 1.47,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        // S, P and anything not listed.
        _ => 1.80,
    }
}

/// Maximum accessible area per residue (Tien et al. 2013, theoretical).
pub fn max_asa(aa: AminoAcid) -> Option<f64> {
    Some(match aa {
        AminoAcid::Ala => 129.0,
        AminoAcid::Arg => 274.0,
        AminoAcid::Asn => 195.0,
        AminoAcid::Asp => 193.0,
        AminoAcid::Cys => 167.0,
        AminoAcid::Gln => 225.0,
        AminoAcid::Glu => 223.0,
        AminoAcid::Gly => 104.0,
        AminoAcid::His | AminoAcid::Met => 224.0,
        AminoAcid::Ile => 197.0,
        AminoAcid::Leu => 201.0,
        AminoAcid::Lys => 236.0,
        AminoAcid::Phe => 240.0,
        AminoAcid::Pro => 159.0,
        AminoAcid::Ser => 155.0,
        AminoAcid::Thr => 172.0,
        AminoAcid::Trp => 285.0,
        AminoAcid::Tyr => 263.0,
        AminoAcid::Val => 174.0,
        AminoAcid::Unknown => return None,
    })
}

fn sphere_points(n: usize) -> Vec<Vec3> {
    let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    (0..n)
        .map(|i| {
            let y = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
            let r = (1.0 - y * y).sqrt();
            let theta = golden * i as f64;
            Vec3::new(r * theta.cos(), y, r * theta.sin())
        })
        .collect()
}

/// Compute SASA.
///
/// # Errors
/// [`AnalysisError::Empty`] if no atoms qualify, or
/// [`AnalysisError::InvalidParameter`] for a zero point count or negative probe.
pub fn compute(structure: &Structure, params: SasaParams) -> Result<Sasa, AnalysisError> {
    if params.points == 0 || params.probe_radius < 0.0 {
        return Err(AnalysisError::InvalidParameter("points must be > 0 and probe radius >= 0"));
    }

    // Flatten qualifying atoms, remembering which output residue each belongs to.
    let mut centers = Vec::new();
    let mut radii = Vec::new();
    let mut owner = Vec::new();
    let mut residues = Vec::new();
    for (chain, residue) in structure.residues() {
        if residue.is_water() || (residue.hetero && !params.include_hetero && residue.amino_acid().is_none())
        {
            continue;
        }
        let index = residues.len();
        let mut any = false;
        for atom in residue.atoms.iter().filter(|a| a.element != "H" && a.element != "D") {
            centers.push(atom.position);
            radii.push(vdw_radius(&atom.element) + params.probe_radius);
            owner.push(index);
            any = true;
        }
        if any {
            residues.push(ResidueSasa {
                chain: chain.id.clone(),
                id: residue.id,
                name: residue.name.clone(),
                area: 0.0,
                relative: None,
            });
        }
    }
    if centers.is_empty() {
        return Err(AnalysisError::Empty);
    }

    // Uniform grid for neighbour lookup; cell edge = largest possible contact distance.
    let max_r = radii.iter().copied().fold(0.0, f64::max);
    let cell = 2.0 * max_r;
    let key =
        |p: Vec3| ((p.x / cell).floor() as i64, (p.y / cell).floor() as i64, (p.z / cell).floor() as i64);
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, &c) in centers.iter().enumerate() {
        grid.entry(key(c)).or_default().push(i);
    }

    let unit = sphere_points(params.points);
    let mut neighbours = Vec::new();
    let mut total = 0.0;
    for i in 0..centers.len() {
        let (ci, ri) = (centers[i], radii[i]);
        neighbours.clear();
        let (kx, ky, kz) = key(ci);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(cands) = grid.get(&(kx + dx, ky + dy, kz + dz)) {
                        for &j in cands {
                            let reach = ri + radii[j];
                            if j != i && ci.distance_squared(centers[j]) < reach * reach {
                                neighbours.push(j);
                            }
                        }
                    }
                }
            }
        }
        let accessible = unit
            .iter()
            .filter(|&&u| {
                let p = ci + u * ri;
                neighbours.iter().all(|&j| p.distance_squared(centers[j]) >= radii[j] * radii[j])
            })
            .count();
        let area = 4.0 * std::f64::consts::PI * ri * ri * accessible as f64 / params.points as f64;
        residues[owner[i]].area += area;
        total += area;
    }

    for r in &mut residues {
        r.relative = AminoAcid::from_three_letter(&r.name).and_then(max_asa).map(|m| r.area / m);
    }
    Ok(Sasa { total, residues })
}

#[cfg(test)]
mod tests {
    use super::*;
    use opendna_core::{Atom, Chain, Residue};

    fn single_atom(element: &str) -> Structure {
        Structure {
            chains: vec![Chain {
                id: "A".into(),
                residues: vec![Residue {
                    name: "LIG".into(),
                    id: ResidueId { number: 1, insertion: None },
                    hetero: true,
                    atoms: vec![Atom {
                        serial: 1,
                        name: element.into(),
                        element: element.into(),
                        position: Vec3::default(),
                        occupancy: 1.0,
                        b_factor: 0.0,
                    }],
                }],
            }],
        }
    }

    #[test]
    fn isolated_atom_has_full_sphere_area() {
        let params = SasaParams { include_hetero: true, ..SasaParams::default() };
        let s = compute(&single_atom("C"), params).unwrap();
        let r: f64 = 1.70 + 1.4;
        assert!((s.total - 4.0 * std::f64::consts::PI * r * r).abs() < 1e-6);
    }

    #[test]
    fn ligands_excluded_by_default() {
        assert_eq!(compute(&single_atom("C"), SasaParams::default()), Err(AnalysisError::Empty));
    }
}
