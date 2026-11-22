//! Checks against published properties of well-characterised PDB entries.
//!
//! Fixtures are unmodified wwPDB files (CC0): 1UBQ (ubiquitin, 1.8 Å) and
//! 1CRN (crambin, 1.5 Å).

#![allow(clippy::unwrap_used, missing_docs)]

use opendna_analysis::{backbone, dssp, interactions, sasa, shape, superpose};
use opendna_core::Structure;

fn load(name: &str) -> Structure {
    let path = format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    opendna_io::read_structure(&text).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// Residue numbers (depositor numbering) whose DSSP state is `code`.
fn residues_with(assignment: &dssp::Assignment, code: char) -> Vec<i32> {
    assignment
        .residues
        .iter()
        .zip(&assignment.states)
        .filter(|(_, s)| s.code() == code)
        .map(|(r, _)| r.id.number)
        .collect()
}

fn assert_all_in(found: &[i32], range: std::ops::RangeInclusive<i32>, what: &str) {
    for r in range.clone() {
        assert!(found.contains(&r), "{what}: residue {r} not assigned; got {found:?}");
    }
}

const UBIQUITIN: &str = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

#[test]
fn ubiquitin_pdb_and_mmcif_agree() {
    let pdb = load("1ubq.pdb");
    let cif = load("1ubq.cif");
    let a = pdb.chain("A").unwrap();
    let b = cif.chain("A").unwrap();
    assert_eq!(a.sequence().to_string(), UBIQUITIN);
    assert_eq!(b.sequence().to_string(), UBIQUITIN);

    let ca = |s: &Structure| backbone::extract(s).iter().map(|r| r.ca).collect::<Vec<_>>();
    assert!(superpose::rmsd(&ca(&pdb), &ca(&cif)).unwrap() < 1e-3);
    assert_eq!(pdb.atom_count(), cif.atom_count());
}

#[test]
fn ubiquitin_secondary_structure() {
    // Published DSSP for 1UBQ: β1 1-7, β2 10-17, α1 23-34, β3 40-45, β5 66-72.
    let a = dssp::assign(&load("1ubq.pdb"));
    assert_eq!(a.residues.len(), 76);
    assert_all_in(&residues_with(&a, 'H'), 24..=33, "ubiquitin α1");
    let strands = residues_with(&a, 'E');
    assert_all_in(&strands, 2..=6, "ubiquitin β1");
    assert_all_in(&strands, 12..=16, "ubiquitin β2");
    assert_all_in(&strands, 41..=44, "ubiquitin β3");
    assert_all_in(&strands, 67..=71, "ubiquitin β5");
}

#[test]
fn crambin_secondary_structure_and_disulfides() {
    // Crambin: author HELIX records 7-19 and 23-30; disulfides 3-40, 4-32, 16-26.
    // Under DSSP rules α1 ends at 17: extending it would need an i→i+4 bond
    // to the NH of Pro19, which has none.
    let s = load("1crn.pdb");
    let a = dssp::assign(&s);
    let helix = residues_with(&a, 'H');
    assert_all_in(&helix, 7..=17, "crambin α1");
    assert!(!helix.contains(&18) && !helix.contains(&19));
    assert_all_in(&helix, 24..=29, "crambin α2");

    let mut pairs: Vec<(i32, i32)> =
        interactions::disulfides(&s).iter().map(|p| (p.a.id.number, p.b.id.number)).collect();
    pairs.sort_unstable();
    assert_eq!(pairs, [(3, 40), (4, 32), (16, 26)]);
}

#[test]
fn helix_torsions_are_alpha_region() {
    let s = load("1ubq.pdb");
    let bb = backbone::extract(&s);
    let t = backbone::torsions(&bb);
    let helix: Vec<_> =
        bb.iter().zip(&t).filter(|(r, _)| (24..=33).contains(&r.id.number)).map(|(_, t)| *t).collect();
    let mean = |f: fn(&backbone::Torsions) -> Option<f64>| {
        helix.iter().filter_map(f).sum::<f64>() / helix.len() as f64
    };
    let (phi, psi) = (mean(|t| t.phi), mean(|t| t.psi));
    assert!((-75.0..-55.0).contains(&phi), "mean φ {phi}");
    assert!((-50.0..-30.0).contains(&psi), "mean ψ {psi}");
    // Termini have undefined torsions.
    assert!(t[0].phi.is_none() && t[75].psi.is_none());
    // Trans peptide bonds throughout.
    assert!(t.iter().filter_map(|t| t.omega).all(|w| w.abs() > 150.0));
}

#[test]
fn ubiquitin_salt_bridge_k11_e34() {
    let pairs = interactions::salt_bridges(&load("1ubq.pdb"));
    assert!(
        pairs.iter().any(|p| p.a.id.number == 34 && p.b.id.number == 11),
        "expected the E34–K11 salt bridge; got {:?}",
        pairs.iter().map(|p| (p.a.id.number, p.b.id.number)).collect::<Vec<_>>()
    );
}

#[test]
fn ubiquitin_size_and_surface() {
    let s = load("1ubq.pdb");
    let ca: Vec<_> = backbone::extract(&s).iter().map(|r| r.ca).collect();
    let rg = shape::radius_of_gyration(&ca).unwrap();
    assert!((11.0..12.5).contains(&rg), "Cα Rg {rg}");

    // Lee–Richards/Shrake–Rupley totals for 1UBQ (heavy atoms, 1.4 Å probe)
    // fall around 4,800–5,000 Å² depending on the radius set.
    let area = sasa::compute(&s, sasa::SasaParams::default()).unwrap();
    assert!((4600.0..5300.0).contains(&area.total), "SASA {}", area.total);
    assert_eq!(area.residues.len(), 76);
    let relative = |n: i32| area.residues.iter().find(|r| r.id.number == n).unwrap().relative.unwrap();
    // Hydrophobic core residues packed between the helix and the sheet are buried...
    for core in [26, 30, 67] {
        assert!(relative(core) < 0.05, "residue {core} relative SASA {}", relative(core));
    }
    // ...and the flexible C-terminal tail is fully exposed.
    assert!(relative(75) > 0.5 && relative(76) > 0.5);
}
