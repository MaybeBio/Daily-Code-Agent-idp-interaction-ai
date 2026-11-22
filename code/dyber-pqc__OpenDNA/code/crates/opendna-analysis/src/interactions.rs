//! Pairwise residue interactions derived from geometry.

use opendna_core::{AminoAcid, ResidueId, Structure, Vec3};
use serde::Serialize;

/// Identifies a residue in a structure.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct ResidueRef {
    /// Chain id.
    pub chain: String,
    /// Depositor numbering.
    pub id: ResidueId,
    /// Amino acid.
    pub amino_acid: AminoAcid,
}

/// A pair of residues with the shortest qualifying atom–atom distance.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ResiduePair {
    /// First residue.
    pub a: ResidueRef,
    /// Second residue.
    pub b: ResidueRef,
    /// Distance in Å.
    pub distance: f64,
}

/// SG–SG cutoff for a disulfide bond, in Å (observed range 2.0–2.1; 2.5 tolerates low resolution).
pub const DISULFIDE_MAX: f64 = 2.5;
/// Charged N–O cutoff for a salt bridge, in Å (Barlow & Thornton 1983).
pub const SALT_BRIDGE_MAX: f64 = 4.0;

fn atoms_of(structure: &Structure, pick: impl Fn(AminoAcid, &str) -> bool) -> Vec<(ResidueRef, Vec3)> {
    let mut out = Vec::new();
    for (chain, residue) in structure.residues() {
        let Some(aa) = residue.amino_acid() else { continue };
        for atom in residue.atoms.iter().filter(|a| pick(aa, &a.name)) {
            out.push((ResidueRef { chain: chain.id.clone(), id: residue.id, amino_acid: aa }, atom.position));
        }
    }
    out
}

fn closest_pairs(
    left: &[(ResidueRef, Vec3)],
    right: &[(ResidueRef, Vec3)],
    max: f64,
    same_set: bool,
) -> Vec<ResiduePair> {
    let mut out: Vec<ResiduePair> = Vec::new();
    for (i, (ra, pa)) in left.iter().enumerate() {
        let start = if same_set { i + 1 } else { 0 };
        for (rb, pb) in &right[start..] {
            if ra == rb {
                continue;
            }
            let d = pa.distance(*pb);
            if d > max {
                continue;
            }
            if let Some(existing) = out.iter_mut().find(|p| &p.a == ra && &p.b == rb) {
                existing.distance = existing.distance.min(d);
            } else {
                out.push(ResiduePair { a: ra.clone(), b: rb.clone(), distance: d });
            }
        }
    }
    out
}

/// Cys–Cys pairs whose SG atoms are within [`DISULFIDE_MAX`].
pub fn disulfides(structure: &Structure) -> Vec<ResiduePair> {
    let sg = atoms_of(structure, |aa, name| aa == AminoAcid::Cys && name == "SG");
    closest_pairs(&sg, &sg, DISULFIDE_MAX, true)
}

/// Acid–base residue pairs with a side-chain O–N distance within [`SALT_BRIDGE_MAX`].
pub fn salt_bridges(structure: &Structure) -> Vec<ResiduePair> {
    let acid = atoms_of(structure, |aa, name| {
        matches!((aa, name), (AminoAcid::Asp, "OD1" | "OD2") | (AminoAcid::Glu, "OE1" | "OE2"))
    });
    let base = atoms_of(structure, |aa, name| {
        matches!(
            (aa, name),
            (AminoAcid::Lys, "NZ") | (AminoAcid::Arg, "NE" | "NH1" | "NH2") | (AminoAcid::His, "ND1" | "NE2")
        )
    });
    closest_pairs(&acid, &base, SALT_BRIDGE_MAX, false)
}

/// Symmetric Cα contact map: `map[i][j]` is `true` when residues `i` and `j`
/// (in [`crate::backbone::extract`] order) have Cα atoms within `cutoff` Å.
pub fn contact_map(backbone: &[crate::backbone::BackboneResidue], cutoff: f64) -> Vec<Vec<bool>> {
    let c2 = cutoff * cutoff;
    backbone.iter().map(|a| backbone.iter().map(|b| a.ca.distance_squared(b.ca) <= c2).collect()).collect()
}
