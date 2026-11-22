//! Whether a predicted probability means what it says.
//!
//! The disorder head emits a number between 0 and 1 and the CLI thresholds it
//! at 0.5. Nothing so far established that 0.7 corresponds to a 70% chance of
//! disorder, and there is a specific reason to expect it does not: the head is
//! trained with weighted cross-entropy, upweighting the minority class so it
//! does not predict "ordered" everywhere. That weight is deliberate and it
//! helps the ranking, but it systematically inflates the probabilities, which
//! is a cost nobody had measured.
//!
//! Two things live here. Measuring calibration — reliability bins, expected
//! calibration error, the Brier score — and correcting it by fitting a map
//! from the head's logit to a calibrated probability.
//!
//! # Which correction
//!
//! Temperature scaling, `sigmoid(z / T)`, is the usual first choice and is the
//! wrong shape here. It can only sharpen or flatten; it cannot move the
//! centre. Class weighting shifts the centre, so the fit needs an intercept
//! too: `sigmoid(a·z + b)`, which is Platt scaling. Both are implemented, the
//! choice is made on the validation split, and the test split is read once
//! afterwards. Choosing the method by test score would make the test number
//! meaningless.
//!
//! # Correlated residues
//!
//! Neighbouring residues almost always share a label, so a hundred thousand
//! residues are nothing like a hundred thousand independent observations. The
//! effective sample size is closer to the few hundred proteins in a split.
//! Every interval here is therefore bootstrapped over **proteins**, resampling
//! whole label arrays rather than individual residues. Bootstrapping over
//! residues would produce an interval several times too narrow, which is worse
//! than reporting none.
//!
//! # What calibration is against
//!
//! The labels. `DisProt` records disorder that has been demonstrated
//! experimentally, and residues that are simply unannotated are treated as
//! ordered during training. If that label set under-annotates disorder — and
//! it certainly does, by an unknown amount — then a probability calibrated
//! against it is calibrated against the annotation, not against nature. This
//! module cannot fix that and does not claim to.

use rand::prelude::*;
use rand::rngs::StdRng;
use std::fmt::Write as _;

/// How far a sigmoid output is kept from 0 and 1 before taking a logarithm.
const EPSILON: f64 = 1e-12;

/// One protein's predictions, kept together because the bootstrap resamples
/// proteins rather than residues.
#[derive(Debug, Clone)]
pub struct Protein {
    /// The head's raw logits, one per assessed residue.
    pub logits: Vec<f64>,
    /// The labels for the same residues.
    pub labels: Vec<bool>,
}

/// One bin of a reliability diagram.
#[derive(Debug, Clone, Copy)]
pub struct Bin {
    /// Lower edge of the predicted-probability range.
    pub lower: f64,
    /// Upper edge.
    pub upper: f64,
    /// Residues falling in it.
    pub count: usize,
    /// Mean predicted probability within the bin: what the model claimed.
    pub predicted: f64,
    /// Fraction actually labelled positive: what happened.
    pub observed: f64,
}

/// The map from a head's logit to a calibrated probability.
///
/// `probability = sigmoid(slope * logit + intercept)`. An identity map is
/// `slope = 1, intercept = 0`, which is what an uncalibrated head is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mapping {
    /// Multiplier on the logit.
    pub slope: f64,
    /// Added to the scaled logit.
    pub intercept: f64,
}

impl Default for Mapping {
    fn default() -> Self {
        Self { slope: 1.0, intercept: 0.0 }
    }
}

impl Mapping {
    /// Apply the map.
    pub fn apply(self, logit: f64) -> f64 {
        sigmoid(self.slope * logit + self.intercept)
    }

    /// The equivalent temperature, for readers who expect that vocabulary.
    ///
    /// `None` when the slope is zero, which would be a degenerate fit rather
    /// than an infinite temperature worth printing.
    pub fn temperature(self) -> Option<f64> {
        (self.slope.abs() > 1e-9).then(|| 1.0 / self.slope)
    }

    /// Whether this is the identity, within printing precision.
    #[cfg(test)]
    pub fn is_identity(self) -> bool {
        (self.slope - 1.0).abs() < 1e-9 && self.intercept.abs() < 1e-9
    }
}

/// Numerically stable logistic function.
pub fn sigmoid(z: f64) -> f64 {
    if z >= 0.0 {
        1.0 / (1.0 + (-z).exp())
    } else {
        let e = z.exp();
        e / (1.0 + e)
    }
}

/// Inverse of [`sigmoid`], for recovering a logit from a stored probability.
///
/// Clamped, because a probability that has been rounded to 1.0 would otherwise
/// give an infinite logit and poison every sum it enters.
///
/// The trainer reads logits from the head directly rather than using this,
/// because recovering one from a rounded probability loses precision exactly
/// where the model is most confident. It exists to build test fixtures.
#[cfg(test)]
pub fn logit(probability: f64) -> f64 {
    let p = probability.clamp(EPSILON, 1.0 - EPSILON);
    (p / (1.0 - p)).ln()
}

/// Flatten proteins into residue-level predictions under `mapping`.
fn flatten(proteins: &[Protein], mapping: Mapping) -> (Vec<f64>, Vec<bool>) {
    let n: usize = proteins.iter().map(|p| p.logits.len()).sum();
    let mut scores = Vec::with_capacity(n);
    let mut labels = Vec::with_capacity(n);
    for protein in proteins {
        scores.extend(protein.logits.iter().map(|z| mapping.apply(*z)));
        labels.extend(protein.labels.iter().copied());
    }
    (scores, labels)
}

/// Equal-width reliability bins.
///
/// Equal-width rather than equal-count because the question is "when the model
/// says 0.7, what happens", and that question is asked about a probability
/// range. Empty bins are returned as empty rather than dropped, since where a
/// model never predicts is itself worth seeing.
pub fn reliability(proteins: &[Protein], mapping: Mapping, bin_count: usize) -> Vec<Bin> {
    let (scores, labels) = flatten(proteins, mapping);
    let bins = bin_count.max(1);
    #[allow(clippy::cast_precision_loss)]
    let width = 1.0 / bins as f64;
    let mut sums = vec![(0usize, 0.0f64, 0usize); bins];
    for (score, label) in scores.iter().zip(&labels) {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let index = ((score / width) as usize).min(bins - 1);
        sums[index].0 += 1;
        sums[index].1 += score;
        sums[index].2 += usize::from(*label);
    }
    sums.iter()
        .enumerate()
        .map(|(index, &(count, total, positives))| {
            #[allow(clippy::cast_precision_loss)]
            let denominator = count.max(1) as f64;
            Bin {
                #[allow(clippy::cast_precision_loss)]
                lower: index as f64 * width,
                #[allow(clippy::cast_precision_loss)]
                upper: (index + 1) as f64 * width,
                count,
                predicted: if count == 0 { 0.0 } else { total / denominator },
                #[allow(clippy::cast_precision_loss)]
                observed: if count == 0 { 0.0 } else { positives as f64 / denominator },
            }
        })
        .collect()
}

/// Expected calibration error: the bin-count-weighted mean gap between what
/// was predicted and what happened.
///
/// Zero is perfect. It is a summary and hides shape — a model that is
/// overconfident at one end and underconfident at the other can score well —
/// so read it beside the reliability table rather than instead of it.
pub fn ece(proteins: &[Protein], mapping: Mapping, bin_count: usize) -> f64 {
    let bins = reliability(proteins, mapping, bin_count);
    let total: usize = bins.iter().map(|b| b.count).sum();
    if total == 0 {
        return 0.0;
    }
    #[allow(clippy::cast_precision_loss)]
    bins.iter()
        .filter(|b| b.count > 0)
        .map(|b| (b.count as f64 / total as f64) * (b.predicted - b.observed).abs())
        .sum()
}

/// Brier score: mean squared error of the probabilities.
///
/// Unlike the ECE it is a proper scoring rule, so it cannot be improved by
/// discarding information. A calibration map that lowered the ECE while
/// raising the Brier score would have thrown something away.
pub fn brier(proteins: &[Protein], mapping: Mapping) -> f64 {
    let (scores, labels) = flatten(proteins, mapping);
    if scores.is_empty() {
        return 0.0;
    }
    #[allow(clippy::cast_precision_loss)]
    let n = scores.len() as f64;
    scores.iter().zip(&labels).map(|(s, l)| (s - f64::from(*l)).powi(2)).sum::<f64>() / n
}

/// Mean negative log likelihood, which is what the fit minimises.
pub fn nll(proteins: &[Protein], mapping: Mapping) -> f64 {
    let (scores, labels) = flatten(proteins, mapping);
    if scores.is_empty() {
        return 0.0;
    }
    #[allow(clippy::cast_precision_loss)]
    let n = scores.len() as f64;
    -scores
        .iter()
        .zip(&labels)
        .map(|(s, l)| {
            let p = s.clamp(EPSILON, 1.0 - EPSILON);
            if *l {
                p.ln()
            } else {
                (1.0 - p).ln()
            }
        })
        .sum::<f64>()
        / n
}

/// Which parameters a fit is allowed to move.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Form {
    /// `sigmoid(z / T)`: slope only. Can sharpen or flatten, cannot shift.
    Temperature,
    /// `sigmoid(a·z + b)`: slope and intercept.
    Platt,
}

impl Form {
    /// Name for a report.
    pub fn name(self) -> &'static str {
        match self {
            Self::Temperature => "temperature",
            Self::Platt => "Platt",
        }
    }
}

/// Fit a calibration map by Newton's method on the log likelihood.
///
/// Two parameters at most, so the Hessian is 2×2 and inverted directly. The
/// objective is convex in both, so this converges from the identity in a few
/// steps and there is no local minimum to worry about.
///
/// # Panics
/// Does not. A singular Hessian or a non-finite step abandons the iteration
/// and returns the best map found so far, which in the worst case is the
/// identity — refusing to calibrate rather than returning a wild fit.
pub fn fit(proteins: &[Protein], form: Form) -> Mapping {
    let (logits, labels): (Vec<f64>, Vec<bool>) = {
        let mut z = Vec::new();
        let mut y = Vec::new();
        for protein in proteins {
            z.extend(protein.logits.iter().copied());
            y.extend(protein.labels.iter().copied());
        }
        (z, y)
    };
    if logits.is_empty() {
        return Mapping::default();
    }

    let mut best = Mapping::default();
    for _ in 0..100 {
        // Gradient and Hessian of the negative log likelihood.
        let (mut ga, mut gb) = (0.0f64, 0.0f64);
        let (mut haa, mut hab, mut hbb) = (0.0f64, 0.0f64, 0.0f64);
        for (z, y) in logits.iter().zip(&labels) {
            let s = best.apply(*z);
            let residual = s - f64::from(*y);
            let w = s * (1.0 - s);
            ga += residual * z;
            gb += residual;
            haa += w * z * z;
            hab += w * z;
            hbb += w;
        }
        // Ridge on the diagonal: the Hessian goes singular once every
        // prediction saturates, and a tiny damping term costs nothing while
        // keeping the step finite.
        haa += 1e-9;
        hbb += 1e-9;

        let (step_a, step_b) = match form {
            Form::Temperature => {
                if haa.abs() < 1e-12 {
                    break;
                }
                (-ga / haa, 0.0)
            }
            Form::Platt => {
                let determinant = haa * hbb - hab * hab;
                if determinant.abs() < 1e-12 {
                    break;
                }
                ((-ga * hbb + gb * hab) / determinant, (-hab * -ga + -haa * gb) / determinant)
            }
        };
        if !step_a.is_finite() || !step_b.is_finite() {
            break;
        }
        let next = Mapping { slope: best.slope + step_a, intercept: best.intercept + step_b };
        // Only accept a step that improves the objective. Newton can overshoot
        // on saturated data, and a worse map is not a fit.
        if nll(proteins, next) > nll(proteins, best) {
            break;
        }
        let moved = step_a.abs() + step_b.abs();
        best = next;
        if moved < 1e-10 {
            break;
        }
    }
    best
}

/// A measurement with a bootstrap interval over proteins.
#[derive(Debug, Clone, Copy)]
pub struct Interval {
    /// The value on the data as it is.
    pub point: f64,
    /// 2.5th percentile over resamples.
    pub low: f64,
    /// 97.5th percentile.
    pub high: f64,
}

impl std::fmt::Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.4} ({:.4}-{:.4})", self.point, self.low, self.high)
    }
}

/// Bootstrap a statistic by resampling whole proteins.
///
/// Resampling residues would treat neighbours as independent observations when
/// they almost always share a label, and would report an interval several
/// times too narrow. The unit of resampling has to be the unit of correlation.
pub fn bootstrap(
    proteins: &[Protein],
    resamples: usize,
    seed: u64,
    statistic: impl Fn(&[Protein]) -> f64,
) -> Interval {
    let point = statistic(proteins);
    if proteins.len() < 2 || resamples == 0 {
        return Interval { point, low: point, high: point };
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let mut values = Vec::with_capacity(resamples);
    for _ in 0..resamples {
        let sample: Vec<Protein> =
            (0..proteins.len()).map(|_| proteins[rng.random_range(0..proteins.len())].clone()).collect();
        values.push(statistic(&sample));
    }
    values.sort_by(f64::total_cmp);
    #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let at = |q: f64| values[((values.len() - 1) as f64 * q).round() as usize];
    Interval { point, low: at(0.025), high: at(0.975) }
}

/// What a fixed threshold does on a split.
#[derive(Debug, Clone, Copy)]
pub struct Operating {
    /// The threshold applied.
    pub threshold: f64,
    /// Residues called positive.
    pub called: usize,
    /// Matthews correlation at that threshold.
    pub mcc: f64,
}

/// Confusion-based summary at `threshold`.
fn operating(proteins: &[Protein], mapping: Mapping, threshold: f64) -> Operating {
    let (scores, labels) = flatten(proteins, mapping);
    let called = scores.iter().filter(|s| **s >= threshold).count();
    let confusion = crate::metrics::Confusion::at(&scores, &labels, threshold);
    Operating { threshold, called, mcc: confusion.mcc().unwrap_or(0.0) }
}

/// The threshold maximising Matthews correlation, searched on a grid.
///
/// Reported so a reader can see whether 0.5 is anywhere near the best
/// operating point, which on an uncalibrated head it usually is not.
fn best_threshold(proteins: &[Protein], mapping: Mapping) -> Operating {
    (1..100)
        .map(|i| operating(proteins, mapping, f64::from(i) / 100.0))
        .max_by(|a, b| a.mcc.total_cmp(&b.mcc))
        .unwrap_or(Operating { threshold: 0.5, called: 0, mcc: 0.0 })
}

/// Everything one calibration run concluded.
pub struct Report {
    /// The map chosen, and the form it came from.
    pub mapping: Mapping,
    /// Which form won on the validation split.
    pub form: Form,
    /// Validation ECE before and after, with intervals.
    pub validation: (Interval, Interval),
    /// Test ECE before and after.
    pub test: (Interval, Interval),
    /// Test Brier before and after, to confirm nothing was thrown away.
    pub test_brier: (f64, f64),
    /// Reliability on the test split after calibration.
    pub test_bins: Vec<Bin>,
    /// Reliability on the test split before.
    pub test_bins_raw: Vec<Bin>,
    /// How many proteins each interval was bootstrapped over.
    pub proteins: (usize, usize),
    /// AUROC on test before and after. Must be equal: the map is monotonic,
    /// so a difference here would mean something is wrong.
    pub test_auroc: (Option<f64>, Option<f64>),
    /// What the CLI's 0.5 threshold does on test, before and after.
    pub at_half: (Operating, Operating),
    /// The best threshold available on test, before and after.
    pub best: (Operating, Operating),
    /// Residues labelled disordered in the test split, for scale.
    pub positives: usize,
    /// Residues assessed in the test split.
    pub residues: usize,
    /// Threshold chosen on the validation split under the fitted map, and
    /// what it achieves on test. This is the one worth shipping: the `best`
    /// pair above is an upper bound nobody can claim, because finding it
    /// needed the test labels.
    pub recommended: (Operating, Operating),
}

/// Bins used throughout. Ten is the convention and keeps every bin populated
/// enough to mean something at these sample sizes.
pub const BINS: usize = 10;

/// Resamples for every interval.
const RESAMPLES: usize = 1000;

/// Measure calibration, fit a correction on validation, and report it on test.
///
/// The form is chosen by validation ECE and the test split is then read once.
/// Selecting on test would make the test number a training score.
pub fn run(validation: &[Protein], test: &[Protein], seed: u64) -> Report {
    let identity = Mapping::default();

    let candidates = [Form::Temperature, Form::Platt];
    let fitted: Vec<(Form, Mapping)> = candidates.iter().map(|f| (*f, fit(validation, *f))).collect();
    let (form, mapping) = fitted
        .iter()
        .min_by(|a, b| ece(validation, a.1, BINS).total_cmp(&ece(validation, b.1, BINS)))
        .copied()
        .unwrap_or((Form::Platt, identity));

    let interval = |set: &[Protein], map: Mapping, salt: u64| {
        bootstrap(set, RESAMPLES, seed.wrapping_add(salt), |sample| ece(sample, map, BINS))
    };

    let (test_scores_raw, test_labels) = flatten(test, identity);
    let (test_scores_cal, _) = flatten(test, mapping);

    Report {
        test_auroc: (
            crate::metrics::auroc(&test_scores_raw, &test_labels),
            crate::metrics::auroc(&test_scores_cal, &test_labels),
        ),
        at_half: (operating(test, identity, 0.5), operating(test, mapping, 0.5)),
        best: (best_threshold(test, identity), best_threshold(test, mapping)),
        recommended: {
            // Chosen without ever looking at the test labels.
            let chosen = best_threshold(validation, mapping);
            (chosen, operating(test, mapping, chosen.threshold))
        },
        positives: test_labels.iter().filter(|l| **l).count(),
        residues: test_labels.len(),
        mapping,
        form,
        validation: (interval(validation, identity, 1), interval(validation, mapping, 2)),
        test: (interval(test, identity, 3), interval(test, mapping, 4)),
        test_brier: (brier(test, identity), brier(test, mapping)),
        test_bins: reliability(test, mapping, BINS),
        test_bins_raw: reliability(test, identity, BINS),
        proteins: (validation.len(), test.len()),
    }
}

impl Report {
    /// Fraction of assessed test residues that are labelled disordered.
    fn fraction_positive(&self) -> f64 {
        if self.residues == 0 {
            return 0.0;
        }
        #[allow(clippy::cast_precision_loss)]
        let f = self.positives as f64 / self.residues as f64;
        f
    }

    /// Fraction of assessed test residues an operating point calls.
    fn called_fraction(&self, op: Operating) -> f64 {
        if self.residues == 0 {
            return 0.0;
        }
        #[allow(clippy::cast_precision_loss)]
        let f = op.called as f64 / self.residues as f64;
        f
    }

    /// Whether the fit is worth applying.
    ///
    /// Two conditions, both necessary. The test ECE has to improve by more
    /// than the width of its own interval — an improvement inside the noise is
    /// not an improvement — and the Brier score must not get worse, because a
    /// map that lowered the ECE while raising the Brier score bought
    /// calibration by discarding information.
    pub fn worth_applying(&self) -> bool {
        let (before, after) = self.test;
        let improvement = before.point - after.point;
        let noise = before.high - before.low;
        improvement > noise && self.test_brier.1 <= self.test_brier.0 + 1e-6
    }

    /// A human-readable account, for the terminal and the model card.
    pub fn summary(&self) -> String {
        let mut out = String::new();
        let _ =
            writeln!(out, "Calibration ({} proteins validation, {} test)", self.proteins.0, self.proteins.1);
        let _ = writeln!(out);
        let _ = writeln!(out, "  form fitted on validation: {}", self.form.name());
        let _ = writeln!(
            out,
            "  mapping: p = sigmoid({:.4} * logit + {:.4}){}",
            self.mapping.slope,
            self.mapping.intercept,
            self.mapping.temperature().map_or(String::new(), |t| format!("   (T = {t:.3})"))
        );
        let _ = writeln!(out);
        let _ = writeln!(out, "  ECE, 95% interval bootstrapped over proteins:");
        let _ = writeln!(out, "    validation  before {}  after {}", self.validation.0, self.validation.1);
        let _ = writeln!(out, "    test        before {}  after {}", self.test.0, self.test.1);
        let _ =
            writeln!(out, "  Brier (test): before {:.4}  after {:.4}", self.test_brier.0, self.test_brier.1);
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "  Ranking is untouched -- the map is monotonic. AUROC on test: before {}  after {}",
            self.test_auroc.0.map_or("n/a".into(), |v| format!("{v:.4}")),
            self.test_auroc.1.map_or("n/a".into(), |v| format!("{v:.4}"))
        );
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "  What it does to the calls. {} of {} test residues are labelled disordered ({:.1}%):",
            self.positives,
            self.residues,
            100.0 * self.fraction_positive()
        );
        for (label, op) in [("before", self.at_half.0), ("after ", self.at_half.1)] {
            let _ = writeln!(
                out,
                "    at 0.5, {label}: {:>7} residues called ({:>5.1}%)  MCC {:.4}",
                op.called,
                100.0 * self.called_fraction(op),
                op.mcc
            );
        }
        for (label, op) in [("before", self.best.0), ("after ", self.best.1)] {
            let _ = writeln!(
                out,
                "    best on test {label}: {:.2}  {:>7} residues called  MCC {:.4}  (upper bound, \
                 not claimable -- it used the test labels)",
                op.threshold, op.called, op.mcc
            );
        }
        let _ = writeln!(
            out,
            "    chosen on validation: {:.2}  ->  on test: {:>7} residues called ({:>5.1}%)  MCC {:.4}",
            self.recommended.0.threshold,
            self.recommended.1.called,
            100.0 * self.called_fraction(self.recommended.1),
            self.recommended.1.mcc
        );
        let _ = writeln!(out);
        let _ = writeln!(out, "  Reliability on the test split:");
        let _ = writeln!(out, "    predicted range      n   claimed  observed  (raw -> calibrated)");
        for (raw, cal) in self.test_bins_raw.iter().zip(&self.test_bins) {
            if raw.count == 0 && cal.count == 0 {
                continue;
            }
            let _ = writeln!(
                out,
                "    {:.1}-{:.1}  {:>8}   {:.3}     {:.3}      |  {:>8}   {:.3}     {:.3}",
                raw.lower,
                raw.upper,
                raw.count,
                raw.predicted,
                raw.observed,
                cal.count,
                cal.predicted,
                cal.observed
            );
        }
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "  verdict: {}",
            if self.worth_applying() {
                "apply it -- the test ECE improves by more than its own interval and the Brier score does not worsen"
            } else {
                "do not apply it -- the improvement is not larger than the noise, or the Brier score worsened"
            }
        );
        out
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// Build proteins whose true positive rate is `rate` and whose logits are
    /// offset from the truth by `shift`, which is what class weighting does.
    fn synthetic(proteins: usize, per_protein: usize, rate: f64, shift: f64, seed: u64) -> Vec<Protein> {
        let mut rng = StdRng::seed_from_u64(seed);
        (0..proteins)
            .map(|_| {
                let mut logits = Vec::with_capacity(per_protein);
                let mut labels = Vec::with_capacity(per_protein);
                for _ in 0..per_protein {
                    let truth = logit(rate);
                    let noise = rng.random_range(-2.0..2.0);
                    let true_p = sigmoid(truth + noise);
                    labels.push(rng.random_bool(true_p));
                    logits.push(truth + noise + shift);
                }
                Protein { logits, labels }
            })
            .collect()
    }

    #[test]
    fn a_perfectly_calibrated_model_has_almost_no_error() {
        let data = synthetic(200, 50, 0.3, 0.0, 1);
        let error = ece(&data, Mapping::default(), BINS);
        assert!(error < 0.02, "a well-calibrated model should not score {error:.4}");
    }

    /// The case this module exists for: a systematic shift, which is what
    /// upweighting the minority class produces.
    #[test]
    fn a_shifted_model_is_detected_and_a_shift_is_what_corrects_it() {
        let data = synthetic(200, 50, 0.3, 1.5, 2);
        let before = ece(&data, Mapping::default(), BINS);
        assert!(before > 0.1, "a 1.5-logit shift should be obvious, not {before:.4}");

        // Temperature alone cannot move the centre, so it should not fix this.
        let temperature = fit(&data, Form::Temperature);
        let after_temperature = ece(&data, temperature, BINS);

        let platt = fit(&data, Form::Platt);
        let after_platt = ece(&data, platt, BINS);

        assert!(
            after_platt < before / 3.0,
            "Platt should largely fix a shift: {before:.4} -> {after_platt:.4}"
        );
        assert!(
            after_platt < after_temperature,
            "Platt ({after_platt:.4}) should beat temperature ({after_temperature:.4}) on a pure shift"
        );
        // And it should recover roughly the shift that was applied.
        assert!(
            (platt.intercept + 1.5).abs() < 0.6,
            "expected an intercept near -1.5, got {:.3}",
            platt.intercept
        );
    }

    /// Temperature is the right tool for overconfidence, and should be allowed
    /// to win when that is the problem.
    #[test]
    fn temperature_corrects_pure_overconfidence() {
        // Logits scaled up by 2 with no shift: too sharp, correctly centred.
        let mut data = synthetic(200, 50, 0.5, 0.0, 3);
        for protein in &mut data {
            for z in &mut protein.logits {
                *z *= 2.0;
            }
        }
        let before = ece(&data, Mapping::default(), BINS);
        let fitted = fit(&data, Form::Temperature);
        let after = ece(&data, fitted, BINS);
        assert!(after < before, "{before:.4} -> {after:.4}");
        assert!(fitted.slope < 1.0, "it should flatten, not sharpen: slope {:.3}", fitted.slope);
    }

    /// The whole point of bootstrapping over proteins. Residue-level
    /// resampling of correlated data gives a far narrower interval, and the
    /// narrow one is the wrong one.
    #[test]
    fn the_interval_widens_when_residues_are_correlated() {
        // Each protein is entirely positive or entirely negative, which is the
        // extreme of the correlation real proteins have.
        let mut rng = StdRng::seed_from_u64(7);
        let blocks: Vec<Protein> = (0..40)
            .map(|_| {
                let positive = rng.random_bool(0.5);
                Protein { logits: vec![logit(0.5); 100], labels: vec![positive; 100] }
            })
            .collect();
        let by_protein = bootstrap(&blocks, 500, 1, |s| ece(s, Mapping::default(), BINS));

        // The same residues, each its own "protein", which is what resampling
        // residues amounts to.
        let by_residue: Vec<Protein> = blocks
            .iter()
            .flat_map(|p| {
                p.logits.iter().zip(&p.labels).map(|(z, l)| Protein { logits: vec![*z], labels: vec![*l] })
            })
            .collect();
        let residue_interval = bootstrap(&by_residue, 500, 1, |s| ece(s, Mapping::default(), BINS));

        let protein_width = by_protein.high - by_protein.low;
        let residue_width = residue_interval.high - residue_interval.low;
        assert!(
            protein_width > residue_width * 2.0,
            "protein-level interval {protein_width:.4} should be much wider than residue-level \
             {residue_width:.4}; if it is not, the bootstrap unit is not doing anything"
        );
    }

    #[test]
    fn the_identity_mapping_changes_nothing() {
        let identity = Mapping::default();
        assert!(identity.is_identity());
        for p in [0.01, 0.3, 0.5, 0.9, 0.99] {
            assert!((identity.apply(logit(p)) - p).abs() < 1e-9, "{p}");
        }
        assert_eq!(identity.temperature(), Some(1.0));
    }

    #[test]
    fn a_rounded_probability_does_not_produce_an_infinite_logit() {
        assert!(logit(1.0).is_finite());
        assert!(logit(0.0).is_finite());
        assert!(logit(1.0) > 0.0 && logit(0.0) < 0.0);
    }

    /// An improvement smaller than the noise is not an improvement, and the
    /// verdict has to say so rather than shipping it.
    #[test]
    fn a_fit_that_barely_helps_is_not_worth_applying() {
        let data = synthetic(150, 40, 0.3, 0.0, 11);
        let report = run(&data, &data, 99);
        assert!(
            !report.worth_applying(),
            "an already-calibrated model needs no correction:\n{}",
            report.summary()
        );
    }

    /// And one that clearly helps has to be recognised.
    #[test]
    fn a_fit_that_fixes_a_real_shift_is_worth_applying() {
        let validation = synthetic(200, 50, 0.3, 1.5, 21);
        let test = synthetic(200, 50, 0.3, 1.5, 22);
        let report = run(&validation, &test, 23);
        assert!(report.worth_applying(), "{}", report.summary());
        assert_eq!(report.form, Form::Platt, "a shift calls for Platt");
        assert!(report.test.1.point < report.test.0.point / 2.0, "{}", report.summary());
        assert!(report.summary().contains("apply it"));
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn empty_input_does_not_panic_and_fits_the_identity() {
        assert_eq!(fit(&[], Form::Platt), Mapping::default());
        assert_eq!(ece(&[], Mapping::default(), BINS), 0.0);
        assert_eq!(brier(&[], Mapping::default()), 0.0);
        let interval = bootstrap(&[], 100, 1, |s| ece(s, Mapping::default(), BINS));
        assert_eq!(interval.point, 0.0);
    }
}
