//! `opendna humanness`: how often each peptide of a sequence occurs in human
//! antibody repertoires.
//!
//! No model runs here. Given a prevalence table, the answer is an exact count:
//! cut the sequence into overlapping 9-mers and look each one up. That is why
//! there is no confidence interval to report and no accuracy to quote -- the
//! number is a property of the table, and what the table is built from decides
//! what the number means.
//!
//! Which is the reason this prints the per-window detail by default rather
//! than a single score. "81% human" invites a comparison with a published
//! figure computed from different repertoires, and the comparison would be
//! wrong. The list of rare windows invites the question the measure can
//! actually answer: which positions are unlike anything seen in a human, and
//! therefore where to look.

use crate::exit::{Classify as _, Failure, Fault};
use anyhow::{bail, Context, Result};
use opendna_analysis::humanness::{PrevalenceTable, Score, Stringency, Window};
use std::path::PathBuf;

/// What `humanness` was asked to do.
pub struct Request {
    /// FASTA file of sequences to score.
    pub input: PathBuf,
    /// The prevalence table to score against.
    pub table: PathBuf,
    /// Which threshold the summary identity is computed at.
    pub stringency: Stringency,
    /// Show every window, not just the rare ones.
    pub all_windows: bool,
    /// How many of the rarest windows to list per sequence.
    pub worst: usize,
}

/// Score every sequence in the input.
///
/// # Errors
/// A missing or malformed table is the installation's problem; unreadable
/// FASTA or a sequence too short to hold a single window is the input's.
pub fn run(args: &Request) -> Result<(), Failure> {
    let file = std::fs::File::open(&args.table)
        .with_context(|| {
            format!(
                "opening the prevalence table {}. Build one with `opendna-data build humanness` \
                 from the dataset pipeline.",
                args.table.display()
            )
        })
        .blame(Fault::Install)?;
    let table = PrevalenceTable::read(std::io::BufReader::with_capacity(1 << 20, file))
        .with_context(|| format!("reading {}", args.table.display()))
        .blame(Fault::Install)?;

    let text = std::fs::read_to_string(&args.input)
        .with_context(|| format!("reading {}", args.input.display()))
        .blame(Fault::Input)?;
    // `read_raw` already refuses a file with no records, so there is no empty
    // case to handle here.
    let sequences = opendna_io::fasta::read_raw(&text)
        .with_context(|| format!("parsing {} as FASTA", args.input.display()))
        .blame(Fault::Input)?;

    // Said before the numbers, not after them. A table's donor count decides
    // which thresholds mean anything, and a threshold finer than the table
    // can express is the one mistake this measure invites.
    eprintln!("table: {} donors, {} peptides, from {}", table.donors(), table.len(), args.table.display());
    if !args.stringency.expressible_with(table.donors()) {
        return Err(Failure {
            fault: Fault::Checks,
            error: anyhow::anyhow!(
                "this table has {} donors, so its smallest non-zero prevalence is {:.2}%. \
                 The {} threshold cannot be evaluated against it, and reporting an identity \
                 at that threshold would be reporting a number nothing measured.",
                table.donors(),
                100.0 / f64::from(table.donors().max(1)),
                args.stringency.label()
            ),
        });
    }

    let mut scores = Vec::new();
    for (name, sequence) in &sequences {
        let score =
            score_one(&table, sequence).with_context(|| format!("scoring {name}")).blame(Fault::Input)?;
        scores.push((name.clone(), score));
    }

    report(&scores, args);
    Ok(())
}

/// Score one sequence, turning "nothing scorable" into an error rather than a
/// zero.
fn score_one(table: &PrevalenceTable, sequence: &str) -> Result<Score> {
    let score = opendna_analysis::humanness::score(table, sequence)?;
    if score.scored == 0 {
        bail!(
            "no scorable window: the sequence has {} residues and {} windows were skipped for \
             holding a non-standard residue. A score of zero would say this looks unlike a human \
             antibody, which is not what was measured.",
            sequence.len(),
            score.skipped
        );
    }
    Ok(score)
}

fn report(scores: &[(String, Score)], args: &Request) {
    println!("{:<24} {:>8} {:>8} {:>9}  rarest window", "name", "identity", "windows", "skipped");
    for (name, score) in scores {
        let identity = score
            .identity(args.stringency)
            .map_or_else(|| "n/a".to_owned(), |value| format!("{:.1}%", value * 100.0));
        let rarest = score.rarest().map_or_else(String::new, |window| {
            format!(
                "{} at {} ({})",
                String::from_utf8_lossy(&window.peptide),
                window.start,
                window
                    .prevalence
                    .map_or_else(|| "not scored".to_owned(), |p| format!("{:.2}% of donors", p * 100.0))
            )
        });
        println!("{name:<24} {identity:>8} {:>8} {:>9}  {rarest}", score.scored, score.skipped);
    }

    for (name, score) in scores {
        if args.all_windows {
            println!("\n{name}: every scored window, in order");
            for window in score.windows.iter().filter(|w| w.prevalence.is_some()) {
                println!(
                    "  {:>5}  {}  {:>7.2}% of donors",
                    window.start,
                    window.as_str(),
                    window.prevalence.unwrap_or(0.0) * 100.0
                );
            }
            continue;
        }

        // Runs, not the n rarest windows. Rare windows cluster -- nine
        // overlapping windows over one unusual CDR3 are one finding -- so a
        // list of the rarest came back as the same stretch reported nine
        // times, every entry at the same floor. A run says which stretch of
        // the sequence is unlike human repertoires, and how long it is.
        let mut runs = score.rare_runs(args.stringency);
        if runs.is_empty() {
            println!("\n{name}: no stretch below the {} threshold.", args.stringency.label());
            continue;
        }
        let total = runs.len();
        runs.sort_by(|a, b| b.windows.cmp(&a.windows).then(a.first.cmp(&b.first)));
        runs.truncate(args.worst);
        runs.sort_by_key(|run| run.first);
        println!(
            "\n{name}: {total} stretch{} below the {}{}",
            if total == 1 { "" } else { "es" },
            args.stringency.label(),
            if runs.len() < total { format!(", {} longest shown", runs.len()) } else { String::new() }
        );
        println!("  {:>13}  {:>7}  {:>8}  rarest peptide in the stretch", "residues", "windows", "lowest");
        for run in runs {
            let rarest = score
                .windows
                .iter()
                .filter(|w| w.start >= run.first && w.start <= run.last && w.prevalence.is_some())
                .min_by(|a, b| a.prevalence.unwrap_or(0.0).total_cmp(&b.prevalence.unwrap_or(0.0)));
            let span = run.residues();
            println!(
                "  {:>13}  {:>7}  {:>7.2}%  {}",
                format!("{}-{}", span.start, span.end),
                run.windows,
                run.lowest * 100.0,
                rarest.map_or_else(String::new, Window::as_str)
            );
        }
    }

    if !args.all_windows {
        println!(
            "\nA stretch is consecutive windows that all fall below the threshold. Residue spans \
             can\noverlap, because a window is nine residues wide and one window that clears the \
             threshold\nis enough to separate two stretches."
        );
    }
    println!(
        "\nPrevalence is the fraction of the table's donors whose repertoire contains that exact \
         peptide.\nIt is a count, not a prediction: there is no model here and no uncertainty to \
         report. A low\nprevalence means the peptide was not seen in those donors' sampled \
         sequences -- not that no\nhuman has it -- and a high one means it occurs in human \
         repertoires, not that an antibody\ncarrying it is safe or non-immunogenic. The table's \
         own card states what it was built from."
    );
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use opendna_analysis::humanness::PEPTIDE_LEN;

    fn table() -> PrevalenceTable {
        // Three donors, so a 50% threshold is expressible and a 1% one is not.
        PrevalenceTable::build(
            vec![vec!["ACDEFGHIKLMNPQRSTVWY"], vec!["ACDEFGHIKLMNPQRSTVWY"], vec!["ACDEFGHIKLMNPQRSTVWYY"]],
            1,
        )
    }

    /// A sequence of nothing but `X` has no scorable window. Returning zero
    /// would assert it is unlike a human antibody, which is the one thing the
    /// measure must not say when it has measured nothing.
    #[test]
    fn a_sequence_with_no_scorable_window_is_an_error_not_a_zero() {
        let t = table();
        let error = score_one(&t, &"X".repeat(PEPTIDE_LEN * 2)).unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("no scorable window"), "{message}");
        assert!(message.contains("would say this looks unlike"), "{message}");
    }

    #[test]
    fn a_sequence_shorter_than_one_window_is_an_error() {
        let t = table();
        assert!(score_one(&t, "ACDEF").is_err());
    }

    #[test]
    fn a_peptide_every_donor_has_scores_as_universal() {
        let t = table();
        let score = score_one(&t, "ACDEFGHIKLMNPQRSTVWY").unwrap();
        assert_eq!(score.identity(Stringency::Strict), Some(1.0));
        // Every window of this sequence is in every donor, so even the rarest
        // is at 100%.
        assert!((score.rarest().unwrap().prevalence.unwrap() - 1.0).abs() < 1e-9);
    }
}
