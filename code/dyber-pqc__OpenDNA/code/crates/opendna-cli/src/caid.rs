//! The CAID prediction format, for submitting the disorder head to the
//! Critical Assessment of Intrinsic Disorder.
//!
//! CAID participants hand the organisers software that is run on blind
//! targets, which is what makes a CAID ranking worth more than any number we
//! compute ourselves: we never see the answers. This module writes the file
//! their scorer reads.
//!
//! # The format, as their parser defines it
//!
//! Verified against `vectorized_metrics/parsers.py` in `BioComputingUP/CAID`
//! and against `AIUPred.caid`, a real CAID3 entry, rather than inferred from
//! prose:
//!
//! ```text
//! >DP02732
//! 1\tM\t0.868\t1
//! 2\tE\t0.870\t1
//! ```
//!
//! One `>accession` block per target, then one tab-separated line per
//! residue: 1-based position, residue, score, and the binary call. Their
//! parser accepts three columns and derives the call by thresholding at 0.5
//! — the same threshold the head was evaluated at — but we emit four, so the
//! call is ours on the record rather than inferred.
//!
//! # Why there are no comments in this format
//!
//! Every other output this CLI writes carries a `#` provenance block. This
//! one cannot. Their reader groups lines on whether they begin with `>`, so a
//! leading comment becomes a block of its own, is read as an accession, and
//! shifts every subsequent target's body onto the wrong protein. The result
//! would parse without complaint and score nonsense. Provenance therefore
//! goes to a sidecar file instead, and the reason is written here because the
//! omission otherwise looks like an oversight.

use anyhow::{bail, Result};
use std::fmt::Write as _;

/// Decimal places used for the score column, matching the published entries
/// and the default their parser rounds to.
pub const DECIMALS: usize = 3;

/// One residue's prediction.
pub struct Residue {
    /// 1-based position.
    pub position: usize,
    /// Amino acid, one letter.
    pub code: char,
    /// Disorder probability.
    pub score: f64,
    /// Whether the score is at or above the calling threshold, at full
    /// precision. [`called`] is what a written file should use.
    pub disordered: bool,
}

/// Whether a score counts as disordered once written to `decimals` places.
///
/// Thresholding the full-precision score and then printing it rounded leaves
/// a file that contradicts itself: a probability of 0.49951 prints as `0.500`
/// and is called ordered, so anyone who thresholds the published number gets
/// the opposite answer. CAID's own parser prefers the state column and would
/// score it faithfully, but a file whose two columns disagree reads as
/// malformed to anyone checking it, and that is not a thing to hand an
/// assessment.
///
/// So the call follows the number actually published. It changes the answer
/// only for residues within half a printed unit of the threshold, where
/// either call is equally defensible, and in exchange every written file is
/// internally consistent.
pub fn called(score: f64, threshold: f64, decimals: usize) -> bool {
    let scale = 10f64.powi(decimals as i32);
    (score * scale).round() / scale >= threshold
}

/// Render per-target predictions in CAID format.
///
/// `targets` is `(accession, residues)` in the order they should appear.
///
/// # Errors
/// If a target has no residues, or if its positions are not the contiguous
/// 1..=n a per-residue prediction must cover. Their scorer joins on position,
/// so a gap would silently align our scores to the wrong residues.
pub fn render(targets: &[(String, Vec<Residue>)], threshold: f64) -> Result<String> {
    if targets.is_empty() {
        bail!("no targets to write");
    }
    let mut out = String::new();
    for (accession, residues) in targets {
        if residues.is_empty() {
            bail!("{accession}: no residues predicted");
        }
        for (expected, residue) in (1..=residues.len()).zip(residues) {
            if residue.position != expected {
                bail!(
                    "{accession}: expected position {expected}, found {}. CAID joins on position, \
                     so a gap or a repeat would align these scores to the wrong residues.",
                    residue.position
                );
            }
        }
        let _ = writeln!(out, ">{accession}");
        for r in residues {
            // Three decimals is what their parser rounds to and what the
            // published entries use.
            // The call follows the score as published, not the full-precision
            // one, so the two columns of this file cannot disagree.
            let _ = writeln!(
                out,
                "{}\t{}\t{:.*}\t{}",
                r.position,
                r.code,
                DECIMALS,
                r.score,
                u8::from(called(r.score, threshold, DECIMALS))
            );
        }
    }
    Ok(out)
}

/// The accession CAID will match against its reference.
///
/// Their reference parser takes `acc, *desc = header[1:].split()`, so the
/// accession is the first whitespace-delimited token of the FASTA header and
/// anything after it is a description. Emitting the whole header would
/// produce an accession that matches nothing, and their parser skips
/// unmatched targets silently — so the submission would score zero targets
/// and look like it simply failed.
pub fn accession(header: &str) -> &str {
    header.split_whitespace().next().unwrap_or(header)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn residues(scores: &[f64]) -> Vec<Residue> {
        scores
            .iter()
            .enumerate()
            .map(|(i, s)| Residue {
                position: i + 1,
                code: "MKVLA".chars().nth(i).unwrap_or('X'),
                score: *s,
                disordered: *s >= 0.5,
            })
            .collect()
    }

    #[test]
    fn the_layout_matches_a_published_entry() {
        let out = render(&[("DP02732".to_string(), residues(&[0.868, 0.870]))], 0.5).unwrap();
        assert_eq!(out, ">DP02732\n1\tM\t0.868\t1\n2\tK\t0.870\t1\n");
    }

    #[test]
    fn the_binary_call_follows_the_threshold() {
        let out = render(&[("A".to_string(), residues(&[0.49, 0.5, 0.51]))], 0.5).unwrap();
        assert!(out.contains("1\tM\t0.490\t0"), "{out}");
        assert!(out.contains("2\tK\t0.500\t1"), "{out}");
        assert!(out.contains("3\tV\t0.510\t1"), "{out}");
    }

    #[test]
    fn several_targets_are_separate_blocks() {
        let out =
            render(&[("A".to_string(), residues(&[0.9])), ("B".to_string(), residues(&[0.1]))], 0.5).unwrap();
        assert_eq!(out.matches('>').count(), 2);
        assert!(out.starts_with(">A\n1\tM\t0.900\t1\n>B\n"), "{out}");
    }

    /// Their scorer joins on position. A gap would align our scores to the
    /// wrong residues and still parse.
    #[test]
    fn a_gap_in_the_positions_is_refused() {
        let mut r = residues(&[0.9, 0.8, 0.7]);
        r[1].position = 5;
        let err = render(&[("A".to_string(), r)], 0.5).unwrap_err().to_string();
        assert!(err.contains("expected position 2, found 5"), "{err}");
        assert!(err.contains("wrong residues"), "{err}");
    }

    #[test]
    fn an_empty_target_is_refused() {
        assert!(render(&[("A".to_string(), vec![])], 0.5).is_err());
        assert!(render(&[], 0.5).is_err());
    }

    /// Found by validating a real 319-target file: a probability of 0.49951
    /// printed as `0.500` while the call said ordered, so the two columns of
    /// the file contradicted each other. CAID's parser prefers the state
    /// column and would have scored it faithfully, but a self-contradicting
    /// file reads as malformed to anyone who checks it.
    #[test]
    fn the_call_agrees_with_the_score_as_published() {
        let r = vec![Residue { position: 1, code: 'M', score: 0.499_51, disordered: false }];
        let out = render(&[("A".to_string(), r)], 0.5).unwrap();
        assert_eq!(out, ">A\n1\tM\t0.500\t1\n", "the printed score and the call must agree");
    }

    #[test]
    fn calling_follows_the_published_precision() {
        // Rounds up to the threshold, so it is called.
        assert!(called(0.499_51, 0.5, 3));
        // Rounds below it, so it is not.
        assert!(!called(0.499_49, 0.5, 3));
        assert!(called(0.5, 0.5, 3));
        assert!(!called(0.0, 0.5, 3));
        assert!(called(1.0, 0.5, 3));
    }

    /// A header with a description must yield only the accession, or the
    /// submission matches no targets at all.
    #[test]
    fn the_accession_is_the_first_token_of_the_header() {
        assert_eq!(accession("DP02732 some description here"), "DP02732");
        assert_eq!(accession("DP02732"), "DP02732");
        assert_eq!(accession("sp|P04637|P53_HUMAN Cellular tumor antigen"), "sp|P04637|P53_HUMAN");
    }

    /// No line may start with '#': their reader groups on '>' and would read
    /// a comment as an accession, shifting every later target's body onto the
    /// wrong protein.
    #[test]
    fn the_output_carries_no_comments() {
        let out = render(&[("A".to_string(), residues(&[0.9, 0.1]))], 0.5).unwrap();
        assert!(!out.lines().any(|l| l.starts_with('#')), "{out}");
    }
}
