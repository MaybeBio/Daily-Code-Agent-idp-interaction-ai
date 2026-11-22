//! Sequence-derived physicochemical properties.
//!
//! Values and methods match `ExPASy` `ProtParam` / Biopython `ProteinAnalysis`:
//! average isotopic residue masses, Kyte–Doolittle hydropathy (1982),
//! extinction coefficients of Pace et al. (1995), aliphatic index of Ikai
//! (1980), and the pKa set used by Biopython's `IsoelectricPoint`.

use crate::AnalysisError;
use opendna_core::{AminoAcid, Sequence};
use serde::Serialize;

/// Physicochemical summary of a sequence.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Properties {
    /// Number of residues.
    pub length: usize,
    /// Average molecular weight in Da.
    pub molecular_weight: f64,
    /// Isoelectric point.
    pub isoelectric_point: f64,
    /// Net charge at pH 7.0.
    pub charge_at_ph7: f64,
    /// Grand average of hydropathy (Kyte–Doolittle).
    pub gravy: f64,
    /// Fraction of F + W + Y.
    pub aromaticity: f64,
    /// Ikai aliphatic index.
    pub aliphatic_index: f64,
    /// Molar extinction coefficient at 280 nm (M⁻¹cm⁻¹) assuming all Cys reduced.
    pub extinction_reduced: f64,
    /// Same, assuming all Cys pairs form cystines.
    pub extinction_cystines: f64,
}

/// Compute [`Properties`].
///
/// # Errors
/// [`AnalysisError::Empty`] or [`AnalysisError::NonStandardResidue`].
pub fn compute(seq: &Sequence) -> Result<Properties, AnalysisError> {
    require_standard(seq)?;
    let n = seq.len() as f64;
    let count = |aa| seq.count(aa) as f64;
    let (w, y, c) = (count(AminoAcid::Trp), count(AminoAcid::Tyr), count(AminoAcid::Cys));
    let extinction_reduced = 5500.0 * w + 1490.0 * y;

    Ok(Properties {
        length: seq.len(),
        molecular_weight: molecular_weight(seq)?,
        isoelectric_point: isoelectric_point(seq)?,
        charge_at_ph7: net_charge(seq, 7.0)?,
        gravy: seq.residues().iter().map(|&aa| kyte_doolittle(aa)).sum::<f64>() / n,
        aromaticity: (count(AminoAcid::Phe) + w + y) / n,
        aliphatic_index: 100.0
            * (count(AminoAcid::Ala)
                + 2.9 * count(AminoAcid::Val)
                + 3.9 * (count(AminoAcid::Ile) + count(AminoAcid::Leu)))
            / n,
        extinction_reduced,
        extinction_cystines: extinction_reduced + 125.0 * (c / 2.0).floor(),
    })
}

fn require_standard(seq: &Sequence) -> Result<(), AnalysisError> {
    if seq.is_empty() {
        return Err(AnalysisError::Empty);
    }
    match seq.residues().iter().position(|aa| !aa.is_standard()) {
        Some(position) => Err(AnalysisError::NonStandardResidue { position }),
        None => Ok(()),
    }
}

/// Average mass of the free amino acid in Da (IUPAC average isotopic masses).
fn amino_acid_mass(aa: AminoAcid) -> f64 {
    match aa {
        AminoAcid::Ala => 89.0932,
        AminoAcid::Arg => 174.2010,
        AminoAcid::Asn => 132.1179,
        AminoAcid::Asp => 133.1027,
        AminoAcid::Cys => 121.1582,
        AminoAcid::Gln => 146.1445,
        AminoAcid::Glu => 147.1293,
        AminoAcid::Gly => 75.0666,
        AminoAcid::His => 155.1546,
        AminoAcid::Ile | AminoAcid::Leu => 131.1729,
        AminoAcid::Lys => 146.1876,
        AminoAcid::Met => 149.2113,
        AminoAcid::Phe => 165.1891,
        AminoAcid::Pro => 115.1305,
        AminoAcid::Ser => 105.0926,
        AminoAcid::Thr => 119.1192,
        AminoAcid::Trp => 204.2252,
        AminoAcid::Tyr => 181.1885,
        AminoAcid::Val => 117.1463,
        AminoAcid::Unknown => f64::NAN,
    }
}

const WATER_MASS: f64 = 18.015_28;

/// Average molecular weight of the linear peptide in Da.
///
/// # Errors
/// See [`compute`].
pub fn molecular_weight(seq: &Sequence) -> Result<f64, AnalysisError> {
    require_standard(seq)?;
    let sum: f64 = seq.residues().iter().map(|&aa| amino_acid_mass(aa)).sum();
    Ok(sum - WATER_MASS * (seq.len() as f64 - 1.0))
}

/// Kyte–Doolittle hydropathy value.
pub fn kyte_doolittle(aa: AminoAcid) -> f64 {
    match aa {
        AminoAcid::Ala => 1.8,
        AminoAcid::Arg => -4.5,
        AminoAcid::Asn | AminoAcid::Asp | AminoAcid::Gln | AminoAcid::Glu => -3.5,
        AminoAcid::Cys => 2.5,
        AminoAcid::Gly => -0.4,
        AminoAcid::His => -3.2,
        AminoAcid::Ile => 4.5,
        AminoAcid::Leu => 3.8,
        AminoAcid::Lys => -3.9,
        AminoAcid::Met => 1.9,
        AminoAcid::Phe => 2.8,
        AminoAcid::Pro => -1.6,
        AminoAcid::Ser => -0.8,
        AminoAcid::Thr => -0.7,
        AminoAcid::Trp => -0.9,
        AminoAcid::Tyr => -1.3,
        AminoAcid::Val => 4.2,
        AminoAcid::Unknown => f64::NAN,
    }
}

/// Sliding-window Kyte–Doolittle profile. Element `k` is the mean over
/// residues `k..k + window`, so the output has `len - window + 1` values.
///
/// # Errors
/// [`AnalysisError::InvalidParameter`] if `window` is zero or longer than the
/// sequence, plus the errors of [`compute`].
pub fn hydropathy_profile(seq: &Sequence, window: usize) -> Result<Vec<f64>, AnalysisError> {
    require_standard(seq)?;
    if window == 0 || window > seq.len() {
        return Err(AnalysisError::InvalidParameter("window must be between 1 and the sequence length"));
    }
    let values: Vec<f64> = seq.residues().iter().map(|&aa| kyte_doolittle(aa)).collect();
    Ok(values.windows(window).map(|w| w.iter().sum::<f64>() / window as f64).collect())
}

// pKa values as used by Biopython's Bio.SeqUtils.IsoelectricPoint.
const PK_N_TERM_DEFAULT: f64 = 9.0;
const PK_C_TERM_DEFAULT: f64 = 2.0;

fn pk_n_terminal(aa: AminoAcid) -> f64 {
    match aa {
        AminoAcid::Ala => 7.59,
        AminoAcid::Met => 7.0,
        AminoAcid::Ser => 6.93,
        AminoAcid::Pro => 8.36,
        AminoAcid::Thr => 6.82,
        AminoAcid::Val => 7.44,
        AminoAcid::Glu => 7.7,
        _ => PK_N_TERM_DEFAULT,
    }
}

fn pk_c_terminal(aa: AminoAcid) -> f64 {
    match aa {
        AminoAcid::Asp => 4.55,
        AminoAcid::Glu => 4.75,
        _ => PK_C_TERM_DEFAULT,
    }
}

/// Net charge at the given pH (Henderson–Hasselbalch, independent sites).
///
/// # Errors
/// See [`compute`].
pub fn net_charge(seq: &Sequence, ph: f64) -> Result<f64, AnalysisError> {
    require_standard(seq)?;
    let r = seq.residues();
    let positive = |pk: f64| {
        let cr = 10f64.powf(pk - ph);
        cr / (cr + 1.0)
    };
    let negative = |pk: f64| {
        let cr = 10f64.powf(ph - pk);
        cr / (cr + 1.0)
    };
    let count = |aa| seq.count(aa) as f64;

    let pos = positive(pk_n_terminal(r[0]))
        + count(AminoAcid::Lys) * positive(10.0)
        + count(AminoAcid::Arg) * positive(12.0)
        + count(AminoAcid::His) * positive(5.98);
    let neg = negative(pk_c_terminal(r[r.len() - 1]))
        + count(AminoAcid::Asp) * negative(4.05)
        + count(AminoAcid::Glu) * negative(4.45)
        + count(AminoAcid::Cys) * negative(9.0)
        + count(AminoAcid::Tyr) * negative(10.0);
    Ok(pos - neg)
}

/// Isoelectric point by bisection on [`net_charge`] to 1e-4 pH units.
///
/// # Errors
/// See [`compute`].
pub fn isoelectric_point(seq: &Sequence) -> Result<f64, AnalysisError> {
    let (mut lo, mut hi) = (0.0_f64, 14.0_f64);
    while hi - lo > 1e-4 {
        // `midpoint` rather than `0.5 * (lo + hi)`: the same value here, where
        // the bracket is [0, 14], and correct in general because it cannot
        // overflow on the way to halving.
        let mid = f64::midpoint(lo, hi);
        if net_charge(seq, mid)? > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(f64::midpoint(lo, hi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    /// Human ubiquitin (`UniProt` P0CG48, residues 1-76).
    const UBIQUITIN: &str = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

    fn seq(s: &str) -> Sequence {
        Sequence::parse(s).unwrap()
    }

    #[test]
    fn ubiquitin_matches_protparam() {
        // Reference: ExPASy ProtParam for P0CG48 1-76.
        let p = compute(&seq(UBIQUITIN)).unwrap();
        assert_eq!(p.length, 76);
        assert_abs_diff_eq!(p.molecular_weight, 8564.8, epsilon = 0.1);
        assert_abs_diff_eq!(p.gravy, -0.489, epsilon = 0.001);
        assert_abs_diff_eq!(p.extinction_reduced, 1490.0);
        assert!((6.0..7.2).contains(&p.isoelectric_point), "pI {}", p.isoelectric_point);
    }

    #[test]
    fn dipeptide_mass_subtracts_one_water() {
        assert_abs_diff_eq!(
            molecular_weight(&seq("GG")).unwrap(),
            2.0 * 75.0666 - WATER_MASS,
            epsilon = 1e-9
        );
    }

    #[test]
    fn basic_and_acidic_peptides_have_expected_pi() {
        assert!(isoelectric_point(&seq("KKKKKKKK")).unwrap() > 10.0);
        assert!(isoelectric_point(&seq("DDDDDDDD")).unwrap() < 4.0);
    }

    #[test]
    fn unknown_residue_is_rejected_with_position() {
        let s = Sequence::parse("MKXV").unwrap();
        assert_eq!(compute(&s), Err(AnalysisError::NonStandardResidue { position: 2 }));
    }

    #[test]
    fn profile_has_expected_length() {
        assert_eq!(hydropathy_profile(&seq(UBIQUITIN), 9).unwrap().len(), 68);
        assert!(hydropathy_profile(&seq("MK"), 3).is_err());
    }

    proptest! {
        #[test]
        fn charge_decreases_with_ph(s in "[ACDEFGHIKLMNPQRSTVWY]{1,100}", a in 0.0f64..13.0, d in 0.1f64..1.0) {
            let s = seq(&s);
            prop_assert!(net_charge(&s, a).unwrap() > net_charge(&s, a + d).unwrap());
        }

        #[test]
        fn charge_at_pi_is_near_zero(s in "[ACDEFGHIKLMNPQRSTVWY]{1,100}") {
            let s = seq(&s);
            let pi = isoelectric_point(&s).unwrap();
            if pi > 0.01 && pi < 13.99 {
                prop_assert!(net_charge(&s, pi).unwrap().abs() < 1e-2);
            }
        }
    }
}
