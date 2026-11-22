//! Whether abstaining would help, measured before anything abstains.
//!
//! An abstention rule is easy to write and hard to justify. The failure mode
//! is not that it refuses too much — that is visible — but that it refuses the
//! wrong residues, leaving the errors in and manufacturing confidence in what
//! it kept. A rule that does that is worse than no rule, because it is
//! advertised as a safeguard.
//!
//! So this measures the thing a rule would depend on: does some quantity
//! available at inference time, with no labels, actually order residues by how
//! likely the prediction is to be wrong?
//!
//! The measurement is a risk–coverage curve. Order residues by confidence,
//! keep the most confident fraction, and measure the error rate on what was
//! kept. Three curves make the result readable:
//!
//! - **the signal under test**, whatever quantity is proposed;
//! - **random**, discarding residues without looking at anything, which is the
//!   null: if the signal does not beat this, it carries no information;
//! - **oracle**, discarding exactly the errors, which is the ceiling nothing
//!   can pass and which says how much there was to win.
//!
//! A signal that sits on the random curve is not a confidence measure, however
//! plausible its story. Reporting that is the point of measuring first.

use crate::calibrate::{Mapping, Protein};
use rand::prelude::*;
use rand::rngs::StdRng;
use std::collections::HashSet;
use std::fmt::Write as _;

/// One residue, reduced to what the measurement needs.
struct Scored {
    /// Confidence in this prediction. Higher means keep it for longer.
    confidence: f64,
    /// Whether the call at the operating threshold was wrong.
    wrong: bool,
}

/// One point on a risk–coverage curve.
///
/// A point sits at a boundary between tie groups, so it is a coverage a
/// threshold on the signal can actually produce.
#[derive(Debug, Clone, Copy)]
pub struct Point {
    /// Error rate among those kept.
    pub risk: f64,
    /// How many residues were kept. The coverage is this over the total,
    /// which the curve's owner knows; storing it again would be a second
    /// field that can disagree with this one.
    pub kept: usize,
    /// How many of those were wrong, which is what interpolation needs.
    pub wrong: usize,
}

/// Risk at an arbitrary coverage, by interpolating the counts.
///
/// Between two achievable points the operating point is reached by keeping a
/// random fraction of the next tie group, so both the kept count and the wrong
/// count grow linearly and the risk is their ratio. Interpolating the risks
/// directly would be wrong, and wrong in the flattering direction for a signal
/// with few distinct values.
fn risk_at(points: &[Point], total: usize, coverage: f64) -> Option<f64> {
    if points.is_empty() || total == 0 || coverage <= 0.0 {
        return None;
    }
    #[allow(clippy::cast_precision_loss)]
    let wanted = coverage * total as f64;
    let mut previous = (0.0_f64, 0.0_f64);
    for point in points {
        #[allow(clippy::cast_precision_loss)]
        let (kept, wrong) = (point.kept as f64, point.wrong as f64);
        if wanted <= kept {
            let span = kept - previous.0;
            if span <= 0.0 {
                return Some(wrong / kept.max(1.0));
            }
            let fraction = (wanted - previous.0) / span;
            let at_kept = previous.0 + fraction * span;
            let at_wrong = previous.1 + fraction * (wrong - previous.1);
            return Some(at_wrong / at_kept.max(1.0));
        }
        previous = (kept, wrong);
    }
    points.last().map(|p| p.risk)
}

/// What orders the residues.
pub enum Signal<'a> {
    /// Distance of the calibrated probability from the decision threshold.
    ///
    /// The natural confidence for a thresholded probability, and meaningful
    /// only because the probability is calibrated: on the uncalibrated head a
    /// margin of 0.2 meant nothing in particular.
    Margin,
    /// Whether the k-mer centred on the residue was seen in training.
    ///
    /// A proxy for "this is unlike what the head was trained on". A crude one:
    /// real novelty would be measured by alignment, and a k-mer is only a
    /// local view. If a crude signal works it is worth refining; if it carries
    /// nothing, a refined version probably does not either.
    Novelty(&'a [Vec<bool>]),
    /// No signal at all. The null hypothesis.
    Random(u64),
    /// Perfect knowledge of which predictions are wrong. The ceiling.
    Oracle,
}

impl Signal<'_> {
    /// Name for a report.
    pub fn name(&self) -> &'static str {
        match self {
            Self::Margin => "calibrated margin",
            Self::Novelty(_) => "k-mer novelty",
            Self::Random(_) => "random (null)",
            Self::Oracle => "oracle (ceiling)",
        }
    }
}

/// Reduce proteins to scored residues under a signal.
fn score(proteins: &[Protein], mapping: Mapping, threshold: f64, signal: &Signal<'_>) -> Vec<Scored> {
    let mut rng = match signal {
        Signal::Random(seed) => Some(StdRng::seed_from_u64(*seed)),
        _ => None,
    };
    let mut out = Vec::new();
    for (index, protein) in proteins.iter().enumerate() {
        for (position, (logit, label)) in protein.logits.iter().zip(&protein.labels).enumerate() {
            let p = mapping.apply(*logit);
            let wrong = (p >= threshold) != *label;
            let confidence = match signal {
                Signal::Margin => (p - threshold).abs(),
                Signal::Novelty(seen) => {
                    // Familiar is confident; novel is not.
                    let familiar = seen.get(index).and_then(|s| s.get(position)).copied().unwrap_or(true);
                    f64::from(u8::from(familiar))
                }
                Signal::Random(_) => rng.as_mut().map_or(0.0, |r| r.random_range(0.0..1.0)),
                Signal::Oracle => f64::from(u8::from(!wrong)),
            };
            out.push(Scored { confidence, wrong });
        }
    }
    out
}

/// Risk at each coverage level, most confident residues kept first.
///
/// Ties are broken by taking whole groups, so a signal with few distinct
/// values — novelty has two — produces the handful of honest points it can
/// support rather than a smooth curve interpolated out of nothing.
pub fn risk_coverage(
    proteins: &[Protein],
    mapping: Mapping,
    threshold: f64,
    signal: &Signal<'_>,
    steps: usize,
) -> Vec<Point> {
    let mut scored = score(proteins, mapping, threshold, signal);
    if scored.is_empty() {
        return Vec::new();
    }
    scored.sort_by(|a, b| b.confidence.total_cmp(&a.confidence));
    let total = scored.len();

    let mut points = Vec::new();
    let mut wrong_so_far = 0usize;
    let mut next_step = 1usize;
    for (index, residue) in scored.iter().enumerate() {
        wrong_so_far += usize::from(residue.wrong);
        let kept = index + 1;
        // Only emit a point where the next residue is less confident, so a
        // reported coverage is one a threshold on this signal can actually
        // produce.
        let boundary = scored.get(kept).is_none_or(|n| n.confidence < residue.confidence);
        #[allow(clippy::cast_precision_loss)]
        let reached = kept as f64 / total as f64 >= next_step as f64 / steps as f64;
        if boundary && reached {
            next_step += 1;
            #[allow(clippy::cast_precision_loss)]
            points.push(Point { risk: wrong_so_far as f64 / kept as f64, kept, wrong: wrong_so_far });
        }
    }
    // Full coverage always, as the reference every other point is read against.
    #[allow(clippy::cast_precision_loss)]
    if points.last().is_none_or(|p| p.kept < total) {
        let wrong = scored.iter().filter(|r| r.wrong).count();
        points.push(Point { risk: wrong as f64 / total as f64, kept: total, wrong });
    }
    points
}

/// Lowest coverage the area is integrated from.
///
/// Not zero: at a handful of residues the error rate is a coin flip, and a
/// signal would be ranked on noise from the far tail of its own ordering. The
/// floor is the same for every signal, which is what makes them comparable.
const COVERAGE_FLOOR: f64 = 0.05;

/// Steps used for integrating the area.
const INTEGRATION_STEPS: usize = 1000;

/// Area under the risk–coverage curve over `[COVERAGE_FLOOR, 1]`.
///
/// Lower is better. Every signal is integrated over the same coverage range on
/// the same grid, so the numbers are comparable even when one signal has two
/// distinct values and another has fifty thousand. It summarises the curve into
/// one number and therefore hides its shape: read it beside the curve.
pub fn aurc(points: &[Point], total: usize) -> f64 {
    if points.is_empty() || total == 0 {
        return 0.0;
    }
    let mut area = 0.0;
    let mut previous: Option<(f64, f64)> = None;
    for step in 0..=INTEGRATION_STEPS {
        #[allow(clippy::cast_precision_loss)]
        let coverage = COVERAGE_FLOOR + (1.0 - COVERAGE_FLOOR) * step as f64 / INTEGRATION_STEPS as f64;
        let Some(risk) = risk_at(points, total, coverage) else { continue };
        if let Some((last_coverage, last_risk)) = previous {
            area += (coverage - last_coverage) * (risk + last_risk) / 2.0;
        }
        previous = Some((coverage, risk));
    }
    area / (1.0 - COVERAGE_FLOOR)
}

/// Mark, per residue, whether its surrounding k-mer occurs in `training`.
///
/// Centred windows; residues too close to either terminus to have a full
/// window count as familiar, because refusing a prediction for being near a
/// terminus would be refusing it for the wrong reason.
pub fn novelty(training: &[String], sequences: &[String], k: usize) -> Vec<Vec<bool>> {
    let mut seen: HashSet<&[u8]> = HashSet::new();
    for sequence in training {
        let bytes = sequence.as_bytes();
        if bytes.len() >= k {
            for window in bytes.windows(k) {
                seen.insert(window);
            }
        }
    }
    sequences
        .iter()
        .map(|sequence| {
            let bytes = sequence.as_bytes();
            (0..bytes.len())
                .map(|position| {
                    let half = k / 2;
                    if bytes.len() < k || position < half || position + k - half > bytes.len() {
                        return true;
                    }
                    seen.contains(&bytes[position - half..position - half + k])
                })
                .collect()
        })
        .collect()
}

/// What the measurement concluded for one signal.
pub struct Measured {
    /// Which signal.
    pub name: &'static str,
    /// The curve.
    pub points: Vec<Point>,
    /// Area under it, over the shared coverage range.
    pub aurc: f64,
    /// Residues the curve was built from, needed to read it at a coverage.
    pub total: usize,
}

/// Compare every signal on one split.
pub struct Report {
    /// One entry per signal, in the order measured.
    pub signals: Vec<Measured>,
    /// Error rate at full coverage, which every curve ends at.
    pub base_risk: f64,
    /// Residues the comparison covers.
    pub residues: usize,
    /// Fraction of residues whose k-mer was novel, for context: a signal that
    /// fires on almost nothing cannot help much however well it orders.
    pub novel_fraction: f64,
}

/// Steps along the coverage axis.
const STEPS: usize = 20;

/// Measure every candidate signal against the null and the ceiling.
pub fn run(proteins: &[Protein], mapping: Mapping, threshold: f64, novel: &[Vec<bool>], seed: u64) -> Report {
    let signals = [Signal::Margin, Signal::Novelty(novel), Signal::Random(seed), Signal::Oracle];
    let measured: Vec<Measured> = signals
        .iter()
        .map(|signal| {
            let points = risk_coverage(proteins, mapping, threshold, signal, STEPS);
            let total = proteins.iter().map(|p| p.labels.len()).sum();
            Measured { name: signal.name(), aurc: aurc(&points, total), total, points }
        })
        .collect();

    let residues: usize = proteins.iter().map(|p| p.labels.len()).sum();
    let novel_count: usize = novel.iter().flatten().filter(|familiar| !**familiar).count();
    #[allow(clippy::cast_precision_loss)]
    let novel_fraction = if residues == 0 { 0.0 } else { novel_count as f64 / residues as f64 };

    Report {
        base_risk: measured.first().and_then(|m| m.points.last()).map_or(0.0, |p| p.risk),
        signals: measured,
        residues,
        novel_fraction,
    }
}

impl Report {
    /// The margin curve as confidence cut-offs, for recording beside the head.
    ///
    /// A coverage is not something a caller can apply -- they have one protein,
    /// not a distribution -- but a confidence cut-off is. So each point is
    /// reported as the margin a call must clear, with the coverage and error
    /// that were measured at it.
    pub fn trade_off(proteins: &[Protein], mapping: Mapping, threshold: f64) -> Vec<(f64, f64, f64)> {
        let mut scored: Vec<Scored> = score(proteins, mapping, threshold, &Signal::Margin);
        scored.sort_by(|a, b| b.confidence.total_cmp(&a.confidence));
        let total = scored.len();
        if total == 0 {
            return Vec::new();
        }
        // Round cut-offs so a recorded number is one somebody would type.
        let mut out = Vec::new();
        for cut in [0.0_f64, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30] {
            let kept: Vec<&Scored> = scored.iter().filter(|s| s.confidence >= cut).collect();
            if kept.is_empty() {
                continue;
            }
            #[allow(clippy::cast_precision_loss)]
            let coverage = kept.len() as f64 / total as f64;
            #[allow(clippy::cast_precision_loss)]
            let error = kept.iter().filter(|s| s.wrong).count() as f64 / kept.len() as f64;
            out.push((cut, coverage, error));
        }
        out
    }

    /// Share of the available reduction the margin captures, 0 to 1.
    pub fn captured(&self) -> f64 {
        let (Some(margin), Some(null), Some(oracle)) =
            (self.signal("calibrated margin"), self.signal("random (null)"), self.signal("oracle (ceiling)"))
        else {
            return 0.0;
        };
        let available = null.aurc - oracle.aurc;
        if available <= 0.0 {
            return 0.0;
        }
        ((null.aurc - margin.aurc) / available).clamp(0.0, 1.0)
    }

    /// Signals measured and found to carry nothing.
    pub fn rejected(&self) -> Vec<String> {
        self.signals
            .iter()
            .filter(|s| !matches!(s.name, "calibrated margin" | "random (null)" | "oracle (ceiling)"))
            .filter(|s| self.signal("random (null)").is_some_and(|null| s.aurc >= null.aurc))
            .map(|s| s.name.to_string())
            .collect()
    }

    /// The signal under test, by name.
    fn signal(&self, name: &str) -> Option<&Measured> {
        self.signals.iter().find(|s| s.name == name)
    }

    /// Whether a signal beats the null by enough to be worth shipping.
    ///
    /// "Beats the null" is the minimum bar, not the bar. The comparison is
    /// against how much the oracle shows was available: a signal that captures
    /// a tenth of the achievable reduction is a rule that mostly refuses
    /// correct predictions, and the right answer for one of those is to say so
    /// rather than to ship it behind a flag.
    pub fn verdict(&self, name: &str) -> String {
        let (Some(signal), Some(null), Some(oracle)) =
            (self.signal(name), self.signal("random (null)"), self.signal("oracle (ceiling)"))
        else {
            return format!("{name}: not measured");
        };
        let available = null.aurc - oracle.aurc;
        let captured = null.aurc - signal.aurc;
        if available <= 0.0 {
            return format!("{name}: nothing to win on this split");
        }
        let share = captured / available;
        if captured <= 0.0 {
            format!("{name}: no better than discarding residues at random -- carries no usable signal")
        } else if share < 0.25 {
            format!(
                "{name}: beats random but captures only {:.0}% of what the oracle shows is \
                 available, so a rule on it would refuse mostly correct predictions",
                100.0 * share
            )
        } else {
            format!("{name}: captures {:.0}% of the available reduction -- usable", 100.0 * share)
        }
    }

    /// A human-readable account.
    pub fn summary(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(
            out,
            "Selective prediction over {} residues, {:.1}% error at full coverage",
            self.residues,
            100.0 * self.base_risk
        );
        let _ = writeln!(
            out,
            "{:.1}% of residues sit on a k-mer unseen in training\n",
            100.0 * self.novel_fraction
        );

        let _ = writeln!(out, "  AURC (lower is better):");
        for signal in &self.signals {
            let _ = writeln!(out, "    {:<20} {:.4}", signal.name, signal.aurc);
        }
        let _ = writeln!(out);

        let _ = writeln!(out, "  Error rate among residues kept, by coverage:");
        let _ = write!(out, "    coverage ");
        for signal in &self.signals {
            let _ = write!(out, "{:>20}", signal.name);
        }
        let _ = writeln!(out);
        for step in [0.5_f64, 0.7, 0.8, 0.9, 0.95, 1.0] {
            let _ = write!(out, "    {:>7.0}% ", 100.0 * step);
            for signal in &self.signals {
                // Read at exactly the coverage asked for, interpolating within
                // a tie group. Showing a nearby point's value under a column
                // header it does not belong to would misreport a signal with
                // few distinct values, which is every signal worth doubting.
                match risk_at(&signal.points, signal.total, step) {
                    Some(risk) => {
                        let _ = write!(out, "{:>19.2}%", 100.0 * risk);
                    }
                    None => {
                        let _ = write!(out, "{:>20}", "-");
                    }
                }
            }
            let _ = writeln!(out);
        }
        let _ = writeln!(out);
        for name in ["calibrated margin", "k-mer novelty"] {
            let _ = writeln!(out, "  {}", self.verdict(name));
        }
        out
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::calibrate::logit;

    /// Proteins where confidence genuinely predicts error: confident residues
    /// are right, uncertain ones are coin flips.
    fn informative(count: usize, seed: u64) -> Vec<Protein> {
        let mut rng = StdRng::seed_from_u64(seed);
        (0..count)
            .map(|_| {
                let mut logits = Vec::new();
                let mut labels = Vec::new();
                for _ in 0..50 {
                    if rng.random_bool(0.5) {
                        // Confident and correct.
                        let positive = rng.random_bool(0.5);
                        logits.push(logit(if positive { 0.97 } else { 0.03 }));
                        labels.push(positive);
                    } else {
                        // Near the threshold and a coin flip.
                        logits.push(logit(rng.random_range(0.45..0.55)));
                        labels.push(rng.random_bool(0.5));
                    }
                }
                Protein { logits, labels }
            })
            .collect()
    }

    #[test]
    fn the_oracle_reaches_zero_risk_and_random_does_not_move() {
        let data = informative(40, 1);
        let novel = vec![vec![true; 50]; 40];
        let report = run(&data, Mapping::default(), 0.5, &novel, 7);

        let oracle = report.signal("oracle (ceiling)").unwrap();
        assert!(oracle.points[0].risk < 1e-9, "the oracle should discard every error first");

        // Random cannot change the expected error rate, so its curve is flat
        // at the base rate.
        let null = report.signal("random (null)").unwrap();
        for point in &null.points {
            assert!(
                (point.risk - report.base_risk).abs() < 0.05,
                "random risk {:.4} drifted from the base {:.4}",
                point.risk,
                report.base_risk
            );
        }
    }

    /// The case the whole module exists to detect: a confidence signal that
    /// works has to be recognised as working.
    #[test]
    fn a_real_confidence_signal_beats_the_null() {
        let data = informative(40, 2);
        let novel = vec![vec![true; 50]; 40];
        let report = run(&data, Mapping::default(), 0.5, &novel, 7);
        let margin = report.signal("calibrated margin").unwrap();
        let null = report.signal("random (null)").unwrap();
        assert!(margin.aurc < null.aurc, "margin {:.4} vs null {:.4}", margin.aurc, null.aurc);
        assert!(report.verdict("calibrated margin").contains("usable"), "{}", report.summary());
    }

    /// And the case that matters more: a signal carrying nothing must be
    /// reported as carrying nothing, not dressed up because it is plausible.
    #[test]
    fn a_useless_signal_is_called_useless() {
        let mut rng = StdRng::seed_from_u64(3);
        // Confidence unrelated to correctness.
        let data: Vec<Protein> = (0..40)
            .map(|_| Protein {
                logits: (0..50).map(|_| logit(rng.random_range(0.05..0.95))).collect(),
                labels: (0..50).map(|_| rng.random_bool(0.5)).collect(),
            })
            .collect();
        let novel = vec![vec![true; 50]; 40];
        let report = run(&data, Mapping::default(), 0.5, &novel, 7);
        let verdict = report.verdict("calibrated margin");
        assert!(
            verdict.contains("no better than") || verdict.contains("only"),
            "a signal unrelated to correctness must not read as usable: {verdict}"
        );
    }

    #[test]
    fn novelty_marks_unseen_kmers_and_spares_the_termini() {
        let training = vec!["AAAAAAAAAA".to_string()];
        let queries = vec!["AAAAAWWWWW".to_string()];
        let marked = novelty(&training, &queries, 5);
        assert_eq!(marked.len(), 1);
        assert_eq!(marked[0].len(), 10);
        // The first two positions have no full centred window, so they are
        // spared rather than refused for being near an end.
        assert!(marked[0][0] && marked[0][1]);
        // The all-A interior is familiar; the all-W tail is not.
        assert!(marked[0][2], "an AAAAA window was in training");
        assert!(!marked[0][7], "a WWWWW window was not");
    }

    #[test]
    fn a_sequence_shorter_than_the_window_is_all_familiar() {
        let marked = novelty(&["AAAAAAAA".to_string()], &["AW".to_string()], 5);
        assert_eq!(marked[0], vec![true, true]);
    }

    #[test]
    fn an_empty_split_measures_nothing_and_does_not_panic() {
        let report = run(&[], Mapping::default(), 0.5, &[], 1);
        assert_eq!(report.residues, 0);
        assert!(report.signals.iter().all(|s| s.points.is_empty()));
    }
}
