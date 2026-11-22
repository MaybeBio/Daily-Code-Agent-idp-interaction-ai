//! Peptide humanness: how much of a sequence occurs in human antibody
//! repertoires (Prihoda et al. 2022, the `OASis` score).
//!
//! The question a humanisation campaign actually asks is not "is this human",
//! which the V gene answers by eye, but "which parts of this are not". So
//! nothing here produces a single verdict and stops: a sequence is cut into
//! overlapping 9-residue windows, each window is looked up in a table of how
//! many human donors have that exact peptide somewhere in their repertoire,
//! and what comes back is every window with its prevalence. The summary score
//! is the fraction of windows that clear a prevalence threshold, and it is the
//! least interesting part of the output.
//!
//! # Why this is in this crate
//!
//! It computes a definition rather than estimating a quantity. "The fraction
//! of this sequence's 9-mers that occur in at least a tenth of human donors"
//! is a count, exact given the table, with no model and no training. A head
//! that predicted humanness would belong in `opendna-ml`; this does not.
//!
//! # Nine residues, and why
//!
//! Nine is the length of a peptide presented by MHC class II, which is the
//! route by which a therapeutic antibody is seen by a T cell. It is also what
//! `OASis` uses, so scores are comparable with published ones. It is not tuned
//! here, and tuning it would make them incomparable.
//!
//! # What a score is not
//!
//! Not an immunogenicity prediction. A sequence made entirely of peptides
//! common in human repertoires can still provoke a response, and a rare
//! peptide is not an epitope. The published correlation with clinical
//! immunogenicity is a correlation across a set of antibodies, not a
//! guarantee about one.
//!
//! Not a measure of anything if the table is thin. Prevalence from `n` donors
//! takes `n + 1` distinct values, so a table built from twenty donors cannot
//! express "one percent of the population" at all. [`PrevalenceTable::donors`]
//! is reported with every score for that reason, and
//! [`Stringency::expressible_with`] says whether a threshold means anything
//! for a given table.
//!
//! > Prihoda, D. et al. (2022). BioPhi: A platform for antibody design,
//! > humanization, and humanness evaluation based on natural antibody
//! > repertoires and deep learning. *mAbs* 14(1), 2020203.
//! > <https://doi.org/10.1080/19420862.2021.2020203>

use crate::AnalysisError;
use std::io::{Read, Write};

/// Residues per window. See the module documentation.
pub const PEPTIDE_LEN: usize = 9;

/// One window's worth of sequence.
pub type Peptide = [u8; PEPTIDE_LEN];

/// File magic, so a wrong file is refused rather than misread.
const MAGIC: &[u8; 16] = b"OPENDNA-HUMAN\0\0\0";

/// Format version of the table file.
pub const FORMAT: u32 = 1;

/// How common a peptide must be to count as human.
///
/// Three named thresholds, matching the stringencies `OASis` reports, so a
/// score here can be compared with a published one. A threshold is a fraction
/// of donors, not of sequences: one deeply sequenced donor must not be able
/// to make a peptide look widespread.
/// The order of the variants is the order of the thresholds, and `Ord` is
/// derived from it, so reordering them would silently reorder every report
/// that sorts by stringency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Stringency {
    /// Present in at least 1% of donors.
    Loose,
    /// At least 10%.
    Medium,
    /// At least 50%.
    Strict,
}

impl Stringency {
    /// Every threshold, for a report that gives all three.
    pub const ALL: [Self; 3] = [Self::Loose, Self::Medium, Self::Strict];

    /// The threshold as a fraction of donors.
    #[must_use]
    pub fn fraction(self) -> f64 {
        match self {
            Self::Loose => 0.01,
            Self::Medium => 0.10,
            Self::Strict => 0.50,
        }
    }

    /// How the threshold is named in output a person reads.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Loose => "loose (>=1% of donors)",
            Self::Medium => "medium (>=10% of donors)",
            Self::Strict => "strict (>=50% of donors)",
        }
    }

    /// Whether a table of `donors` donors can express this threshold.
    ///
    /// With `n` donors the smallest non-zero prevalence is `1/n`. A threshold
    /// below that is not strict at all -- it means "at least one donor" -- and
    /// reporting it as "1% of donors" would overstate what was measured.
    #[must_use]
    pub fn expressible_with(self, donors: u32) -> bool {
        donors > 0 && self.fraction() >= 1.0 / f64::from(donors)
    }
}

/// How many human donors have each peptide.
///
/// Sorted by peptide and searched by bisection: the table runs to millions of
/// entries, and a hash map of them costs several times the memory for a
/// lookup that is not the bottleneck.
#[derive(Debug, Clone, Default)]
pub struct PrevalenceTable {
    donors: u32,
    /// `(peptide, donors containing it)`, sorted by peptide.
    entries: Vec<(Peptide, u32)>,
}

impl PrevalenceTable {
    /// How many donors the table was built from.
    #[must_use]
    pub fn donors(&self) -> u32 {
        self.donors
    }

    /// How many distinct peptides it holds.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the table holds no peptides at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Every peptide and the number of donors holding it, in peptide order.
    ///
    /// For reporting what a table contains -- how much of it is near-universal
    /// and how much is one donor's idiosyncrasy -- which is the honest way to
    /// describe a table whose coverage depends entirely on how it was built.
    pub fn entries(&self) -> impl ExactSizeIterator<Item = (&Peptide, u32)> {
        self.entries.iter().map(|(peptide, donors)| (peptide, *donors))
    }

    /// Drop peptides held by fewer than `min_donors` donors.
    ///
    /// Separate from [`PrevalenceTable::build`] so a builder can count the
    /// whole distribution, report it, and only then discard the tail. The
    /// donor total is unchanged: it is the denominator of the measure, not a
    /// count of what survived.
    pub fn restrict(&mut self, min_donors: u32) {
        self.entries.retain(|(_, donors)| *donors >= min_donors);
    }

    /// The fraction of donors holding `peptide`; zero if it is not in the
    /// table.
    ///
    /// Absent means "not seen in any donor's repertoire", which is what a
    /// prevalence of zero says. It does not mean the peptide cannot occur in
    /// a human: see the module documentation on thin tables.
    #[must_use]
    pub fn prevalence(&self, peptide: &Peptide) -> f64 {
        if self.donors == 0 {
            return 0.0;
        }
        match self.entries.binary_search_by(|(p, _)| p.cmp(peptide)) {
            Ok(i) => f64::from(self.entries[i].1) / f64::from(self.donors),
            Err(_) => 0.0,
        }
    }

    /// Build a table from one repertoire per donor.
    ///
    /// A peptide counts once per donor however many of that donor's sequences
    /// contain it, which is the whole point of the measure.
    ///
    /// `min_donors` drops peptides seen in fewer than that many donors. Two is
    /// a reasonable floor: a peptide in exactly one donor's sample is not
    /// evidence of prevalence, and the single-donor tail is most of the table.
    /// Dropping it is recorded in the file rather than left implicit.
    pub fn build<I, S>(repertoires: impl IntoIterator<Item = I>, min_donors: u32) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut counts: std::collections::HashMap<Peptide, u32> = std::collections::HashMap::new();
        let mut donors = 0_u32;
        for repertoire in repertoires {
            donors += 1;
            // Per donor, so a peptide in a thousand of their sequences counts
            // once.
            let mut seen: std::collections::HashSet<Peptide> = std::collections::HashSet::new();
            for sequence in repertoire {
                for peptide in peptides(sequence.as_ref()) {
                    seen.insert(peptide);
                }
            }
            for peptide in seen {
                *counts.entry(peptide).or_default() += 1;
            }
        }
        let mut entries: Vec<(Peptide, u32)> = counts.into_iter().filter(|(_, n)| *n >= min_donors).collect();
        entries.sort_unstable_by_key(|(peptide, _)| *peptide);
        Self { donors, entries }
    }

    /// Write the table.
    ///
    /// # Errors
    /// If the writer fails.
    pub fn write(&self, mut out: impl Write) -> Result<(), AnalysisError> {
        let io = |e: std::io::Error| AnalysisError::Malformed(format!("writing humanness table: {e}"));
        out.write_all(MAGIC).map_err(io)?;
        out.write_all(&FORMAT.to_le_bytes()).map_err(io)?;
        out.write_all(&(PEPTIDE_LEN as u32).to_le_bytes()).map_err(io)?;
        out.write_all(&self.donors.to_le_bytes()).map_err(io)?;
        out.write_all(&(self.entries.len() as u64).to_le_bytes()).map_err(io)?;
        for (peptide, count) in &self.entries {
            out.write_all(peptide).map_err(io)?;
            out.write_all(&count.to_le_bytes()).map_err(io)?;
        }
        Ok(())
    }

    /// Read a table.
    ///
    /// # Errors
    /// [`AnalysisError::Malformed`] if the magic, version or peptide length do
    /// not match, if a count exceeds the donor total, or if the entries are
    /// not sorted -- the last two because a table that lies about itself
    /// would silently produce prevalences above one or miss lookups.
    pub fn read(mut input: impl Read) -> Result<Self, AnalysisError> {
        let io = |e: std::io::Error| AnalysisError::Malformed(format!("reading humanness table: {e}"));
        let mut magic = [0_u8; 16];
        input.read_exact(&mut magic).map_err(io)?;
        if &magic != MAGIC {
            return Err(AnalysisError::Malformed(
                "not an OpenDNA humanness table: the file does not start with the expected magic".into(),
            ));
        }
        let mut four = [0_u8; 4];
        input.read_exact(&mut four).map_err(io)?;
        let format = u32::from_le_bytes(four);
        if format != FORMAT {
            return Err(AnalysisError::Malformed(format!(
                "humanness table is format {format}, this build reads {FORMAT}"
            )));
        }
        input.read_exact(&mut four).map_err(io)?;
        let length = u32::from_le_bytes(four) as usize;
        if length != PEPTIDE_LEN {
            return Err(AnalysisError::Malformed(format!(
                "humanness table uses {length}-residue peptides, this build uses {PEPTIDE_LEN}"
            )));
        }
        input.read_exact(&mut four).map_err(io)?;
        let donors = u32::from_le_bytes(four);
        let mut eight = [0_u8; 8];
        input.read_exact(&mut eight).map_err(io)?;
        let count = usize::try_from(u64::from_le_bytes(eight))
            .map_err(|_| AnalysisError::Malformed("humanness table is too large to read".into()))?;

        let mut entries = Vec::with_capacity(count.min(1 << 20));
        let mut buffer = [0_u8; PEPTIDE_LEN + 4];
        for _ in 0..count {
            input.read_exact(&mut buffer).map_err(io)?;
            let mut peptide = [0_u8; PEPTIDE_LEN];
            peptide.copy_from_slice(&buffer[..PEPTIDE_LEN]);
            let donors_with = u32::from_le_bytes([
                buffer[PEPTIDE_LEN],
                buffer[PEPTIDE_LEN + 1],
                buffer[PEPTIDE_LEN + 2],
                buffer[PEPTIDE_LEN + 3],
            ]);
            if donors_with == 0 || donors_with > donors {
                return Err(AnalysisError::Malformed(format!(
                    "humanness table says a peptide occurs in {donors_with} of {donors} donors"
                )));
            }
            if entries.last().is_some_and(|(p, _): &(Peptide, u32)| p >= &peptide) {
                return Err(AnalysisError::Malformed(
                    "humanness table is not sorted by peptide, so lookups would miss entries".into(),
                ));
            }
            entries.push((peptide, donors_with));
        }
        Ok(Self { donors, entries })
    }
}

/// Whether a byte is one of the twenty standard amino acids.
///
/// Through `opendna_core` rather than a local table, so there is one answer
/// to this question in the project.
fn is_standard_residue(byte: u8) -> bool {
    opendna_core::AminoAcid::from_one_letter(byte as char).is_some_and(opendna_core::AminoAcid::is_standard)
}

/// Every 9-residue window of `sequence`, in order.
///
/// Residues are upper-cased. Windows containing anything that is not one of
/// the twenty standard amino acids are still produced: whether to score them
/// is [`score`]'s decision, not this function's.
pub fn peptides(sequence: &str) -> impl Iterator<Item = Peptide> + '_ {
    let bytes = sequence.as_bytes();
    (0..bytes.len().saturating_sub(PEPTIDE_LEN - 1)).map(move |i| {
        let mut peptide = [0_u8; PEPTIDE_LEN];
        for (slot, byte) in peptide.iter_mut().zip(&bytes[i..i + PEPTIDE_LEN]) {
            *slot = byte.to_ascii_uppercase();
        }
        peptide
    })
}

/// A stretch of consecutive windows that all fail a threshold.
///
/// What a humanisation campaign acts on. See [`Score::rare_runs`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RareRun {
    /// Start position of the first window in the run.
    pub first: usize,
    /// Start position of the last window in the run.
    pub last: usize,
    /// How many windows it covers.
    pub windows: usize,
    /// The lowest prevalence in the run.
    pub lowest: f64,
}

impl RareRun {
    /// The residue range the run covers, as a half-open interval.
    ///
    /// A window starting at `last` extends [`PEPTIDE_LEN`] residues, so the
    /// stretch of sequence implicated is wider than the window positions
    /// suggest.
    #[must_use]
    pub fn residues(&self) -> std::ops::Range<usize> {
        self.first..self.last + PEPTIDE_LEN
    }
}

/// One window of a scored sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct Window {
    /// Zero-based offset of the window in the sequence.
    pub start: usize,
    /// The peptide, upper-cased.
    pub peptide: Peptide,
    /// Fraction of donors holding it, or `None` when it could not be looked
    /// up because it contains a residue the table cannot hold.
    pub prevalence: Option<f64>,
}

impl Window {
    /// The peptide as text.
    #[must_use]
    pub fn as_str(&self) -> String {
        String::from_utf8_lossy(&self.peptide).into_owned()
    }

    /// Whether this window clears `stringency`.
    #[must_use]
    pub fn is_human(&self, stringency: Stringency) -> bool {
        self.prevalence.is_some_and(|p| p >= stringency.fraction())
    }
}

/// A scored sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct Score {
    /// Every window, so the reader can see which parts are rare. This is the
    /// output that is actually useful; `identity` is a summary of it.
    pub windows: Vec<Window>,
    /// Windows that could be looked up.
    pub scored: usize,
    /// Windows skipped because they contain a non-standard residue. Excluded
    /// from `identity` rather than counted as non-human: an `X` is not
    /// evidence either way.
    pub skipped: usize,
    /// Donors behind the table, carried so a score is never read without it.
    pub donors: u32,
}

impl Score {
    /// Fraction of scorable windows that clear `stringency`.
    ///
    /// `None` when nothing was scorable -- a sequence shorter than nine
    /// residues, or one made entirely of ambiguity codes. Zero would be a
    /// claim that the sequence is not human, which is not what no measurement
    /// means.
    #[must_use]
    pub fn identity(&self, stringency: Stringency) -> Option<f64> {
        if self.scored == 0 {
            return None;
        }
        let human = self.windows.iter().filter(|w| w.is_human(stringency)).count();
        Some(human as f64 / self.scored as f64)
    }

    /// Stretches of consecutive windows that all fail `stringency`.
    ///
    /// Listing the *n* rarest windows is close to useless in practice: the
    /// rare ones cluster, and the whole list comes back as one CDR3 reported
    /// nine times over, every entry at the same floor. A run says the thing
    /// worth knowing -- which stretch of the sequence is unlike human
    /// repertoires, and how long it is.
    ///
    /// Two windows are consecutive when their start positions differ by one.
    /// An unscorable window between them therefore ends a run, because
    /// bridging it would assert a continuity across a position nothing
    /// measured.
    #[must_use]
    pub fn rare_runs(&self, stringency: Stringency) -> Vec<RareRun> {
        let mut runs: Vec<RareRun> = Vec::new();
        for window in &self.windows {
            let Some(prevalence) = window.prevalence else { continue };
            if window.is_human(stringency) {
                continue;
            }
            match runs.last_mut() {
                Some(run) if run.last + 1 == window.start => {
                    run.last = window.start;
                    run.windows += 1;
                    run.lowest = run.lowest.min(prevalence);
                }
                _ => runs.push(RareRun {
                    first: window.start,
                    last: window.start,
                    windows: 1,
                    lowest: prevalence,
                }),
            }
        }
        runs
    }

    /// The rarest scorable window, which is where a humanisation campaign
    /// starts.
    #[must_use]
    pub fn rarest(&self) -> Option<&Window> {
        self.windows
            .iter()
            .filter(|w| w.prevalence.is_some())
            .min_by(|a, b| a.prevalence.partial_cmp(&b.prevalence).unwrap_or(std::cmp::Ordering::Equal))
    }
}

/// Score a sequence against a table.
///
/// # Errors
/// [`AnalysisError::Malformed`] if the table holds no donors, since every
/// prevalence would be zero and the output would look like a sequence with
/// nothing human in it.
pub fn score(table: &PrevalenceTable, sequence: &str) -> Result<Score, AnalysisError> {
    if table.donors() == 0 {
        return Err(AnalysisError::Malformed(
            "humanness table was built from no donors, so it cannot score anything".into(),
        ));
    }
    let mut windows = Vec::new();
    let mut scored = 0;
    let mut skipped = 0;
    for (start, peptide) in peptides(sequence).enumerate() {
        let standard = peptide.iter().all(|c| is_standard_residue(*c));
        let prevalence = if standard {
            scored += 1;
            Some(table.prevalence(&peptide))
        } else {
            skipped += 1;
            None
        };
        windows.push(Window { start, peptide, prevalence });
    }
    Ok(Score { windows, scored, skipped, donors: table.donors() })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two donors, with one peptide in common.
    fn table() -> PrevalenceTable {
        PrevalenceTable::build(
            vec![
                vec!["AAAAAAAAACCCCC"], // AAAAAAAAA, AAAAAAAAC, ...
                vec!["AAAAAAAAADDDDD"], // shares AAAAAAAAA
            ],
            1,
        )
    }

    /// A run is the useful unit precisely because rare windows cluster: nine
    /// overlapping windows covering one unusual stretch are one finding, not
    /// nine.
    #[test]
    fn consecutive_rare_windows_collapse_into_one_run() {
        // Two donors share the leading A-run; the tail is unique to one.
        let table =
            PrevalenceTable::build(vec![vec!["AAAAAAAAAAAAAAAAA"], vec!["AAAAAAAAAAAAAAAAAWWWWWWWW"]], 1);
        // Scored against the first donor's own sequence plus a foreign tail.
        let score = score(&table, "AAAAAAAAAAAAAAAAAQQQQQQQQ").unwrap();
        let runs = score.rare_runs(Stringency::Strict);
        assert_eq!(runs.len(), 1, "one unusual stretch is one finding: {runs:?}");
        let run = runs[0];
        assert!(run.windows > 1, "the Q tail spans several windows: {run:?}");
        assert_eq!(run.lowest, 0.0, "those peptides are in no donor");
        // The run's last window still covers nine residues, so the stretch
        // implicated reaches the end of the sequence.
        assert_eq!(run.residues().end, 25);
        assert!(run.residues().start > 0, "the shared A-run is not rare");
    }

    /// A sequence entirely of common peptides has no rare run, and that has to
    /// be an empty list rather than one run covering everything.
    #[test]
    fn a_sequence_of_common_peptides_has_no_rare_run() {
        let table =
            PrevalenceTable::build(vec![vec!["ACDEFGHIKLMNPQRSTVWY"], vec!["ACDEFGHIKLMNPQRSTVWY"]], 1);
        let score = score(&table, "ACDEFGHIKLMNPQRSTVWY").unwrap();
        assert_eq!(score.rare_runs(Stringency::Strict), vec![]);
    }

    #[test]
    fn restricting_drops_the_rare_tail_and_keeps_the_denominator() {
        let mut t = table();
        let before = t.len();
        let shared = *b"AAAAAAAAA";
        assert_eq!(t.prevalence(&shared), 1.0, "both donors hold it");

        t.restrict(2);
        assert!(t.len() < before, "the single-donor peptides should be gone");
        // The denominator is how many donors the table was built from, not how
        // many peptides survived, so prevalence does not move.
        assert_eq!(t.donors(), 2);
        assert_eq!(t.prevalence(&shared), 1.0);
        assert_eq!(t.prevalence(b"AAAAAAAAC"), 0.0, "one donor only, now absent");
        assert!(t.entries().all(|(_, donors)| donors >= 2));
    }

    #[test]
    fn entries_are_in_peptide_order_and_count_donors() {
        let t = table();
        let listed: Vec<(String, u32)> =
            t.entries().map(|(p, n)| (String::from_utf8_lossy(p).into_owned(), n)).collect();
        assert_eq!(listed.len(), t.len());
        assert!(listed.windows(2).all(|w| w[0].0 < w[1].0), "{listed:?}");
        assert_eq!(listed.iter().find(|(p, _)| p == "AAAAAAAAA").map(|(_, n)| *n), Some(2));
        assert!(listed.iter().all(|(_, n)| *n <= t.donors()));
    }

    #[test]
    fn windows_are_every_overlapping_nine_mer() {
        let found: Vec<String> =
            peptides("ACDEFGHIKLM").map(|p| String::from_utf8_lossy(&p).into_owned()).collect();
        assert_eq!(found, ["ACDEFGHIK", "CDEFGHIKL", "DEFGHIKLM"]);
        // Shorter than a window yields none rather than a padded one.
        assert_eq!(peptides("ACDEFGHI").count(), 0);
        assert_eq!(peptides("").count(), 0);
        assert_eq!(peptides("ACDEFGHIK").count(), 1);
    }

    #[test]
    fn lower_case_input_is_scored_the_same() {
        let one: Vec<Peptide> = peptides("acdefghik").collect();
        let two: Vec<Peptide> = peptides("ACDEFGHIK").collect();
        assert_eq!(one, two);
    }

    /// A peptide counts once per donor however often it occurs in that
    /// donor's sequences. Otherwise a deeply sequenced donor would make its
    /// own peptides look universal.
    #[test]
    fn prevalence_counts_donors_not_sequences() {
        let deep = vec!["AAAAAAAAA"; 1000];
        let other = vec!["CCCCCCCCC"];
        let built = PrevalenceTable::build(vec![deep, other], 1);
        assert_eq!(built.donors(), 2);
        assert!((built.prevalence(b"AAAAAAAAA") - 0.5).abs() < 1e-12);
        assert!((built.prevalence(b"CCCCCCCCC") - 0.5).abs() < 1e-12);
    }

    #[test]
    fn a_shared_peptide_is_in_every_donor() {
        let built = table();
        assert!((built.prevalence(b"AAAAAAAAA") - 1.0).abs() < 1e-12);
        // Present in one of two.
        assert!((built.prevalence(b"AAAAAAAAC") - 0.5).abs() < 1e-12);
        // Never seen.
        assert!((built.prevalence(b"WWWWWWWWW")).abs() < 1e-12);
    }

    #[test]
    fn a_rare_tail_can_be_dropped_at_build_time() {
        let both = table();
        let common_only = PrevalenceTable::build(vec![vec!["AAAAAAAAACCCCC"], vec!["AAAAAAAAADDDDD"]], 2);
        assert!(common_only.len() < both.len());
        assert!((common_only.prevalence(b"AAAAAAAAA") - 1.0).abs() < 1e-12);
        // The single-donor peptide is gone, which reads as prevalence zero.
        assert!(common_only.prevalence(b"AAAAAAAAC").abs() < 1e-12);
    }

    #[test]
    fn a_table_round_trips_through_its_file_format() {
        let built = table();
        let mut bytes = Vec::new();
        built.write(&mut bytes).unwrap();
        let read = PrevalenceTable::read(bytes.as_slice()).unwrap();
        assert_eq!(read.donors(), built.donors());
        assert_eq!(read.len(), built.len());
        assert!((read.prevalence(b"AAAAAAAAA") - 1.0).abs() < 1e-12);
    }

    #[test]
    fn a_file_that_is_not_a_table_is_refused() {
        let error = PrevalenceTable::read(b"not a table at all, really".as_slice()).unwrap_err();
        assert!(error.to_string().contains("does not start with the expected magic"), "{error}");
    }

    /// A table claiming a peptide is in more donors than it has would produce
    /// a prevalence above one, which would then pass any threshold.
    #[test]
    fn a_table_that_contradicts_itself_is_refused() {
        let mut bytes = Vec::new();
        PrevalenceTable { donors: 2, entries: vec![(*b"AAAAAAAAA", 5)] }.write(&mut bytes).unwrap();
        let error = PrevalenceTable::read(bytes.as_slice()).unwrap_err();
        assert!(error.to_string().contains("occurs in 5 of 2 donors"), "{error}");
    }

    /// Unsorted entries would make bisection miss them, so a score would
    /// quietly read lower than it should.
    #[test]
    fn an_unsorted_table_is_refused_rather_than_searched() {
        let mut bytes = Vec::new();
        PrevalenceTable { donors: 2, entries: vec![(*b"CCCCCCCCC", 1), (*b"AAAAAAAAA", 1)] }
            .write(&mut bytes)
            .unwrap();
        let error = PrevalenceTable::read(bytes.as_slice()).unwrap_err();
        assert!(error.to_string().contains("not sorted"), "{error}");
    }

    #[test]
    fn scoring_reports_every_window_and_the_rarest() {
        let built = table();
        let scored = score(&built, "AAAAAAAAACCCCC").unwrap();
        assert_eq!(scored.windows.len(), 6);
        assert_eq!(scored.scored, 6);
        assert_eq!(scored.skipped, 0);
        assert_eq!(scored.windows[0].as_str(), "AAAAAAAAA");
        assert_eq!(scored.windows[0].prevalence, Some(1.0));
        // The rarest window is one of the ones only this donor had.
        assert_eq!(scored.rarest().unwrap().prevalence, Some(0.5));
    }

    /// An ambiguity code is not evidence of anything, so its window is left
    /// out of the denominator rather than counted against the sequence.
    #[test]
    fn a_window_with_a_non_standard_residue_is_skipped_not_counted_against() {
        let built = table();
        let scored = score(&built, "AAAAAAAAAXAAAAAAAAA").unwrap();
        assert!(scored.skipped > 0, "{scored:?}");
        assert_eq!(scored.scored + scored.skipped, scored.windows.len());
        // Both all-A windows are present and human, so identity is 1 over the
        // windows that could be looked up.
        assert_eq!(scored.identity(Stringency::Strict), Some(1.0));
    }

    /// Nothing scorable must not read as zero, which would say "not human".
    #[test]
    fn a_sequence_too_short_to_score_has_no_identity() {
        let built = table();
        let scored = score(&built, "ACDEF").unwrap();
        assert_eq!(scored.scored, 0);
        assert_eq!(scored.identity(Stringency::Loose), None);
        assert_eq!(scored.rarest(), None);
    }

    #[test]
    fn an_empty_table_refuses_to_score() {
        let empty = PrevalenceTable::default();
        let error = score(&empty, "AAAAAAAAA").unwrap_err();
        assert!(error.to_string().contains("no donors"), "{error}");
    }

    /// The honesty guard on a thin table. With twenty donors the smallest
    /// non-zero prevalence is 5%, so "at least 1% of donors" is really "at
    /// least one donor" and must not be reported as the former.
    #[test]
    fn a_threshold_finer_than_the_table_is_not_expressible() {
        assert!(!Stringency::Loose.expressible_with(20));
        assert!(Stringency::Medium.expressible_with(20));
        assert!(Stringency::Strict.expressible_with(20));
        // `OASis` was built from 266 subjects, where 1% is expressible.
        assert!(Stringency::Loose.expressible_with(266));
        assert!(!Stringency::Loose.expressible_with(0));
    }

    #[test]
    fn the_thresholds_are_the_published_ones() {
        assert!((Stringency::Loose.fraction() - 0.01).abs() < 1e-12);
        assert!((Stringency::Medium.fraction() - 0.10).abs() < 1e-12);
        assert!((Stringency::Strict.fraction() - 0.50).abs() < 1e-12);
    }

    /// A stricter threshold can never call more windows human.
    #[test]
    fn identity_is_monotone_in_stringency() {
        let built = PrevalenceTable::build(
            vec![
                vec!["AAAAAAAAACCCCCCCCC"],
                vec!["AAAAAAAAAGGGGGGGGG"],
                vec!["AAAAAAAAATTTTTTTTT"],
                vec!["AAAAAAAAAYYYYYYYYY"],
            ],
            1,
        );
        let scored = score(&built, "AAAAAAAAACCCCCCCCC").unwrap();
        let loose = scored.identity(Stringency::Loose).unwrap();
        let medium = scored.identity(Stringency::Medium).unwrap();
        let strict = scored.identity(Stringency::Strict).unwrap();
        assert!(loose >= medium, "{loose} < {medium}");
        assert!(medium >= strict, "{medium} < {strict}");
    }
}
