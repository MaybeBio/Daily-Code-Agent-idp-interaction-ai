//! Backbone extraction and torsion angles.

use opendna_core::{AminoAcid, ResidueId, Structure, Vec3};
use serde::Serialize;

/// Maximum C(i)–N(i+1) distance, in Å, for two residues to count as
/// covalently linked. Larger gaps are treated as chain breaks (DSSP uses 2.5).
pub const PEPTIDE_BOND_MAX: f64 = 2.5;

/// Backbone atoms of one amino-acid residue.
#[derive(Debug, Clone, Serialize)]
pub struct BackboneResidue {
    /// Chain identifier.
    pub chain: String,
    /// Depositor numbering.
    pub id: ResidueId,
    /// Residue identity.
    pub amino_acid: AminoAcid,
    /// Amide nitrogen.
    pub n: Vec3,
    /// Alpha carbon.
    pub ca: Vec3,
    /// Carbonyl carbon.
    pub c: Vec3,
    /// Carbonyl oxygen.
    pub o: Vec3,
    /// `true` if this residue is covalently bonded to the previous one in the list.
    pub linked_to_previous: bool,
}

/// All amino-acid residues that have N, CA, C and O, in chain order.
/// Residues missing any backbone atom are skipped and cause a chain break.
pub fn extract(structure: &Structure) -> Vec<BackboneResidue> {
    let mut out: Vec<BackboneResidue> = Vec::new();
    for chain in &structure.chains {
        let mut previous_in_chain = false;
        for residue in &chain.residues {
            let (Some(aa), Some(n), Some(ca), Some(c), Some(o)) = (
                residue.amino_acid(),
                residue.atom("N"),
                residue.atom("CA"),
                residue.atom("C"),
                residue.atom("O").or_else(|| residue.atom("OXT")),
            ) else {
                previous_in_chain = false;
                continue;
            };
            let linked =
                previous_in_chain && out.last().is_some_and(|p| p.c.distance(n.position) <= PEPTIDE_BOND_MAX);
            out.push(BackboneResidue {
                chain: chain.id.clone(),
                id: residue.id,
                amino_acid: aa,
                n: n.position,
                ca: ca.position,
                c: c.position,
                o: o.position,
                linked_to_previous: linked,
            });
            previous_in_chain = true;
        }
    }
    out
}

/// Backbone torsions of one residue, in degrees. `None` at chain termini and breaks.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Torsions {
    /// φ: C(i-1)–N–CA–C
    pub phi: Option<f64>,
    /// ψ: N–CA–C–N(i+1)
    pub psi: Option<f64>,
    /// ω: CA–C–N(i+1)–CA(i+1)
    pub omega: Option<f64>,
}

/// φ/ψ/ω for each residue returned by [`extract`].
pub fn torsions(backbone: &[BackboneResidue]) -> Vec<Torsions> {
    (0..backbone.len())
        .map(|i| {
            let r = &backbone[i];
            let prev = (i > 0 && r.linked_to_previous).then(|| &backbone[i - 1]);
            let next = backbone.get(i + 1).filter(|n| n.linked_to_previous);
            Torsions {
                phi: prev.map(|p| Vec3::dihedral(p.c, r.n, r.ca, r.c)),
                psi: next.map(|n| Vec3::dihedral(r.n, r.ca, r.c, n.n)),
                omega: next.map(|n| Vec3::dihedral(r.ca, r.c, n.n, n.ca)),
            }
        })
        .collect()
}
