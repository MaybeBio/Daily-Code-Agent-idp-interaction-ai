//! How wrong a regressor's number is likely to be.
//!
//! The three regressors report a point: a ddG in kcal/mol, a solubility
//! percentage, a relative accessibility. A point with no spread invites being
//! read as a measurement, and the stability head in particular is a ranking
//! signal at about 0.44 Spearman out of domain — a number somebody should not
//! act on as though it were within 0.1 kcal/mol of the truth.
//!
//! Split conformal prediction turns the point into an interval. Take the
//! absolute residuals on a calibration split, take their `1 - alpha` quantile
//! with a finite-sample correction, and report `prediction ± that`. On
//! exchangeable data this covers the truth at least `1 - alpha` of the time,
//! with no assumption about the shape of the error — which is what makes it
//! usable on a model whose errors are nothing like Gaussian.
//!
//! # The two things it does not promise
//!
//! **Coverage is marginal, not conditional.** A constant half-width covers
//! `1 - alpha` of cases *on average*. It is too wide where the model is
//! reliable and too narrow where it is not, and nothing in the method says
//! which case you are looking at. So conditional coverage is reported across
//! bins of the predicted value: if it is badly uneven, the interval means
//! much less than its headline number suggests.
//!
//! **Exchangeability is assumed, and this project's splits deliberately break
//! it.** The dataset builder clusters by sequence identity so that near-
//! identical proteins cannot straddle a split, which makes the test split
//! systematically more distant from calibration than a random sample would
//! be. That is the right way to build the splits and it is exactly the
//! condition conformal needs. So the coverage actually achieved on the test
//! split is measured rather than assumed, and where it falls short of nominal
//! the shortfall is the result.

use std::fmt::Write as _;

/// One prediction against its truth.
#[derive(Debug, Clone, Copy)]
pub struct Prediction {
    /// What the head said.
    pub predicted: f64,
    /// What was measured.
    pub actual: f64,
}

impl Prediction {
    /// The absolute residual, which is the conformity score.
    fn residual(self) -> f64 {
        (self.predicted - self.actual).abs()
    }
}

/// Predictions that are not independent of each other.
///
/// Residues within one protein, or mutations within one domain, move together.
/// Grouping them keeps the unit of resampling honest when intervals are put
/// on the coverage figures.
#[derive(Debug, Clone)]
pub struct Group(pub Vec<Prediction>);

/// Flatten groups into residuals.
fn residuals(groups: &[Group]) -> Vec<f64> {
    groups.iter().flat_map(|g| g.0.iter().map(|p| p.residual())).collect()
}

/// The conformal half-width at level `alpha`.
///
/// The rank used is `ceil((n + 1) * (1 - alpha))`, which is the finite-sample
/// correction that makes the guarantee exact rather than asymptotic. When that
/// rank exceeds the number of calibration points there is no finite interval
/// with the requested coverage, and this returns `None` rather than quietly
/// handing back the largest residual it happens to have seen — a 99% interval
/// from 50 calibration points is not a 99% interval.
pub fn half_width(calibration: &[Group], alpha: f64) -> Option<f64> {
    let mut scores = residuals(calibration);
    if scores.is_empty() {
        return None;
    }
    scores.sort_by(f64::total_cmp);
    let n = scores.len();
    #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let rank = (((n + 1) as f64) * (1.0 - alpha)).ceil() as usize;
    if rank == 0 || rank > n {
        return None;
    }
    Some(scores[rank - 1])
}

/// Fraction of predictions whose interval contains the truth.
pub fn coverage(groups: &[Group], half_width: f64) -> f64 {
    let scores = residuals(groups);
    if scores.is_empty() {
        return 0.0;
    }
    #[allow(clippy::cast_precision_loss)]
    let covered = scores.iter().filter(|r| **r <= half_width).count() as f64;
    #[allow(clippy::cast_precision_loss)]
    let total = scores.len() as f64;
    covered / total
}

/// Coverage within one band of the predicted value.
#[derive(Debug, Clone, Copy)]
pub struct Band {
    /// Lower edge of the predicted value.
    pub low: f64,
    /// Upper edge.
    pub high: f64,
    /// Predictions in it.
    pub count: usize,
    /// Fraction of them the interval covered.
    pub coverage: f64,
}

/// Coverage across bands of the predicted value.
///
/// The check on the marginal figure. An interval covering 90% overall while
/// covering 99% of easy cases and 60% of hard ones is not a 90% interval for
/// anybody with a hard case, and the marginal number alone hides that.
pub fn conditional(groups: &[Group], half_width: f64, bands: usize) -> Vec<Band> {
    let mut all: Vec<Prediction> = groups.iter().flat_map(|g| g.0.iter().copied()).collect();
    if all.is_empty() || bands == 0 {
        return Vec::new();
    }
    all.sort_by(|a, b| a.predicted.total_cmp(&b.predicted));
    // Equal-count bands, so every reported figure rests on the same number of
    // observations and none is a single outlier's opinion.
    let per = all.len().div_ceil(bands);
    all.chunks(per)
        .map(|chunk| {
            #[allow(clippy::cast_precision_loss)]
            let covered = chunk.iter().filter(|p| p.residual() <= half_width).count() as f64;
            #[allow(clippy::cast_precision_loss)]
            let count = chunk.len() as f64;
            Band {
                low: chunk.first().map_or(0.0, |p| p.predicted),
                high: chunk.last().map_or(0.0, |p| p.predicted),
                count: chunk.len(),
                coverage: covered / count,
            }
        })
        .collect()
}

/// One confidence level, calibrated and then measured.
#[derive(Debug, Clone)]
pub struct Level {
    /// Nominal coverage, `1 - alpha`.
    pub nominal: f64,
    /// Half-width from the calibration split, in the head's own units.
    pub half_width: f64,
    /// Coverage actually achieved on the held-out split.
    pub achieved: f64,
    /// Coverage across bands of the predicted value on the held-out split.
    pub bands: Vec<Band>,
}

impl Level {
    /// Whether the held-out coverage reaches what was asked for.
    ///
    /// A small shortfall is sampling noise; a large one means exchangeability
    /// failed, which for splits clustered by sequence identity is the
    /// expected direction rather than a surprise.
    pub fn holds(&self) -> bool {
        self.achieved >= self.nominal - 0.02
    }

    /// The widest and narrowest conditional coverage.
    pub fn spread(&self) -> (f64, f64) {
        let lows = self.bands.iter().map(|b| b.coverage);
        let highs = self.bands.iter().map(|b| b.coverage);
        (lows.fold(f64::INFINITY, f64::min).min(1.0), highs.fold(f64::NEG_INFINITY, f64::max).max(0.0))
    }

    /// Bands that fall meaningfully short of what the level promised.
    ///
    /// The question is not whether the bands agree with each other but
    /// whether each delivers the headline. A 90% interval covering one band
    /// at 79% has failed the people holding those predictions, and that some
    /// other band reached 97% is no comfort to them.
    pub fn short_bands(&self) -> Vec<Band> {
        self.bands.iter().filter(|b| b.coverage < self.nominal - SHORTFALL).copied().collect()
    }
}

/// What a conformal run concluded for one head.
pub struct Report {
    /// The head's units, for printing.
    pub units: &'static str,
    /// One entry per requested level.
    pub levels: Vec<Level>,
    /// Calibration and held-out counts.
    pub counts: (usize, usize),
    /// Mean absolute error on the held-out split, for scale.
    pub mae: f64,
}

/// Levels measured by default.
pub const LEVELS: [f64; 3] = [0.5, 0.8, 0.9];

/// Bands used for the conditional check.
const BANDS: usize = 5;

/// How far below nominal a band may fall before it is called short.
///
/// Five points. Smaller than that is sampling noise at these band sizes;
/// larger means the interval does not do for that part of the range what its
/// label says.
const SHORTFALL: f64 = 0.05;

/// Calibrate on one split and measure on another.
pub fn run(calibration: &[Group], held_out: &[Group], units: &'static str) -> Report {
    let levels = LEVELS
        .iter()
        .filter_map(|nominal| {
            let half_width = half_width(calibration, 1.0 - nominal)?;
            Some(Level {
                nominal: *nominal,
                half_width,
                achieved: coverage(held_out, half_width),
                bands: conditional(held_out, half_width, BANDS),
            })
        })
        .collect();
    let scores = residuals(held_out);
    #[allow(clippy::cast_precision_loss)]
    let mae = if scores.is_empty() { 0.0 } else { scores.iter().sum::<f64>() / scores.len() as f64 };
    Report { units, levels, counts: (residuals(calibration).len(), scores.len()), mae }
}

impl Report {
    /// A human-readable account.
    pub fn summary(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(
            out,
            "Conformal intervals: {} calibration, {} held-out predictions, MAE {:.3} {}",
            self.counts.0, self.counts.1, self.mae, self.units
        );
        let _ = writeln!(out);
        let _ = writeln!(out, "  nominal   half-width   achieved on held-out");
        for level in &self.levels {
            let _ = writeln!(
                out,
                "  {:>5.0}%   +/- {:>8.3} {}   {:>5.1}%{}",
                100.0 * level.nominal,
                level.half_width,
                self.units,
                100.0 * level.achieved,
                if level.holds() { "" } else { "   <- short of nominal" }
            );
        }
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "  Coverage across bands of the predicted value, at the {:.0}% interval:",
            100.0 * self.levels.last().map_or(0.0, |l| l.nominal)
        );
        if let Some(level) = self.levels.last() {
            let _ = writeln!(out, "    predicted range            n   covered");
            for band in &level.bands {
                let _ = writeln!(
                    out,
                    "    {:>8.3} to {:>8.3}  {:>7}   {:>5.1}%",
                    band.low,
                    band.high,
                    band.count,
                    100.0 * band.coverage
                );
            }
            let (low, high) = level.spread();
            let short = level.short_bands();
            let _ = writeln!(out);
            let _ =
                writeln!(out, "  Conditional coverage ranges {:.0}% to {:.0}%.", 100.0 * low, 100.0 * high);
            if short.is_empty() {
                let _ = writeln!(
                    out,
                    "  Every band reaches within {:.0} points of the {:.0}% it promised.",
                    100.0 * SHORTFALL,
                    100.0 * level.nominal
                );
            } else {
                let _ = writeln!(
                    out,
                    "  {} of {} bands fall more than {:.0} points short of the {:.0}% promised:",
                    short.len(),
                    level.bands.len(),
                    100.0 * SHORTFALL,
                    100.0 * level.nominal
                );
                for band in &short {
                    let _ = writeln!(
                        out,
                        "    predicted {:.3} to {:.3}: {:.0}%",
                        band.low,
                        band.high,
                        100.0 * band.coverage
                    );
                }
                let _ = writeln!(
                    out,
                    "  A constant half-width cannot fix that: it is too wide where the head is \
                     reliable and too narrow where it is not, and the headline figure is an \
                     average of the two that describes neither."
                );
            }
        }
        out
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use rand::prelude::*;
    use rand::rngs::StdRng;

    fn gaussian_like(count: usize, spread: f64, seed: u64) -> Vec<Group> {
        let mut rng = StdRng::seed_from_u64(seed);
        (0..count)
            .map(|_| {
                let predicted = rng.random_range(-3.0..3.0);
                // Sum of uniforms: roughly bell-shaped, no normal needed.
                let noise: f64 = (0..6).map(|_| rng.random_range(-spread..spread)).sum::<f64>() / 3.0;
                Group(vec![Prediction { predicted, actual: predicted + noise }])
            })
            .collect()
    }

    /// The guarantee, which is the whole reason for using this method.
    ///
    /// Stated carefully, because a careless version of this test fails for
    /// the wrong reason. Conformal coverage holds *in expectation over the
    /// calibration draw*: for one calibration set the achieved coverage is
    /// itself a random quantity, and at these sizes it moves by a couple of
    /// percentage points either way. So the mean over several independent
    /// draws is what must reach the nominal level, and a single draw is
    /// allowed its sampling error.
    #[test]
    fn coverage_reaches_the_nominal_level_on_exchangeable_data() {
        const DRAWS: u64 = 12;
        let mut totals = [0.0_f64; LEVELS.len()];
        for draw in 0..DRAWS {
            let calibration = gaussian_like(2000, 1.0, 100 + draw * 2);
            let held_out = gaussian_like(2000, 1.0, 101 + draw * 2);
            let report = run(&calibration, &held_out, "units");
            for (index, level) in report.levels.iter().enumerate() {
                totals[index] += level.achieved;
                // Three standard errors on 2000 held-out points, doubled to
                // allow for the calibration draw's own variability.
                let tolerance = 3.0 * (0.25_f64 / 2000.0).sqrt() * 2.0;
                assert!(
                    level.achieved >= level.nominal - tolerance,
                    "draw {draw}: the {:.0}% interval covered {:.1}%, further below than \
                     sampling error explains",
                    100.0 * level.nominal,
                    100.0 * level.achieved
                );
            }
        }
        #[allow(clippy::cast_precision_loss)]
        for (index, nominal) in LEVELS.iter().enumerate() {
            let mean = totals[index] / DRAWS as f64;
            assert!(
                mean >= nominal - 0.005,
                "mean coverage over {DRAWS} draws was {mean:.3} against a nominal {nominal:.2}; \
                 the guarantee is about this average, so falling below it is a real failure"
            );
        }
    }

    /// And the failure it must not hide: a held-out set that is harder than
    /// calibration breaks the guarantee, which has to show up as a shortfall.
    #[test]
    fn a_harder_held_out_set_undercovers_and_says_so() {
        let calibration = gaussian_like(2000, 1.0, 3);
        let held_out = gaussian_like(2000, 2.5, 4);
        let report = run(&calibration, &held_out, "units");
        let ninety = report.levels.iter().find(|l| (l.nominal - 0.9).abs() < 1e-9).unwrap();
        assert!(ninety.achieved < 0.9, "coverage {:.3} should fall short", ninety.achieved);
        assert!(!ninety.holds());
        assert!(report.summary().contains("short of nominal"), "{}", report.summary());
    }

    /// Asking for more confidence than the calibration set can support must
    /// return nothing rather than the largest residual seen.
    #[test]
    fn too_few_calibration_points_cannot_support_a_tight_level() {
        let small = gaussian_like(10, 1.0, 5);
        // 99% from 10 points needs rank ceil(11 * 0.99) = 11 > 10.
        assert!(half_width(&small, 0.01).is_none());
        // 90% needs rank ceil(11 * 0.9) = 10, which it has.
        assert!(half_width(&small, 0.10).is_some());
    }

    #[test]
    fn the_half_width_is_the_corrected_quantile_not_the_plain_one() {
        // Residuals 1..=9 against a prediction of 0.
        let groups: Vec<Group> =
            (1..=9).map(|i| Group(vec![Prediction { predicted: 0.0, actual: f64::from(i) }])).collect();
        // n = 9, alpha = 0.1: rank = ceil(10 * 0.9) = 9, so the largest.
        assert!((half_width(&groups, 0.1).unwrap() - 9.0).abs() < 1e-9);
        // alpha = 0.5: rank = ceil(10 * 0.5) = 5.
        assert!((half_width(&groups, 0.5).unwrap() - 5.0).abs() < 1e-9);
    }

    /// Uneven conditional coverage has to be visible, because a constant
    /// half-width hides it behind the marginal figure.
    #[test]
    fn uneven_conditional_coverage_is_reported() {
        let mut rng = StdRng::seed_from_u64(7);
        // Error grows with the predicted value: easy at the bottom, hard at
        // the top, and the marginal coverage says nothing about either.
        let make = |rng: &mut StdRng, n: usize| -> Vec<Group> {
            (0..n)
                .map(|_| {
                    let predicted: f64 = rng.random_range(0.0..10.0);
                    let noise = rng.random_range(-1.0..1.0) * predicted;
                    Group(vec![Prediction { predicted, actual: predicted + noise }])
                })
                .collect()
        };
        let calibration = make(&mut rng, 2000);
        let held_out = make(&mut rng, 2000);
        let report = run(&calibration, &held_out, "units");
        let level = report.levels.last().unwrap();
        let (low, high) = level.spread();
        assert!(high - low > 0.2, "expected uneven coverage, got {low:.2} to {high:.2}");
        assert!(
            !level.short_bands().is_empty(),
            "a band should fall short of nominal when the error grows with the prediction"
        );
        assert!(report.summary().contains("short of the"), "{}", report.summary());
        assert!(report.summary().contains("constant half-width cannot fix"), "{}", report.summary());
    }

    /// A band short of its promise has to be named, even when the bands'
    /// range looks unremarkable.
    #[test]
    fn a_band_short_of_nominal_is_named() {
        let level = Level {
            nominal: 0.9,
            half_width: 1.0,
            achieved: 0.9,
            bands: vec![
                Band { low: 0.0, high: 1.0, count: 100, coverage: 0.97 },
                Band { low: 1.0, high: 2.0, count: 100, coverage: 0.83 },
                Band { low: 2.0, high: 3.0, count: 100, coverage: 0.90 },
            ],
        };
        let short = level.short_bands();
        assert_eq!(short.len(), 1, "the 83% band is more than five points short");
        assert!((short[0].coverage - 0.83).abs() < 1e-9);
        // The marginal figure is fine, which is exactly why the band matters.
        assert!(level.holds());
    }

    #[test]
    fn an_empty_split_yields_no_levels() {
        let report = run(&[], &[], "units");
        assert!(report.levels.is_empty());
        assert!(half_width(&[], 0.1).is_none());
    }
}
