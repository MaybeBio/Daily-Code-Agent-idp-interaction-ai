//! Evaluation metrics.
//!
//! These decide what the model cards report, so each is tested against a
//! worked example rather than trusted.

/// Pearson correlation coefficient. Returns `None` if either input has zero
/// variance, where the coefficient is undefined.
pub fn pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    assert_eq!(a.len(), b.len(), "pearson: length mismatch");
    let n = a.len() as f64;
    if a.len() < 2 {
        return None;
    }
    let (mean_a, mean_b) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let mut cov = 0.0;
    let mut var_a = 0.0;
    let mut var_b = 0.0;
    for (x, y) in a.iter().zip(b) {
        let (dx, dy) = (x - mean_a, y - mean_b);
        cov += dx * dy;
        var_a += dx * dx;
        var_b += dy * dy;
    }
    if var_a <= 0.0 || var_b <= 0.0 {
        return None;
    }
    Some(cov / (var_a * var_b).sqrt())
}

/// Ranks, averaging ties — the standard "fractional ranking" used by
/// Spearman's rho.
///
/// Ties are found by exact equality, which is what ranking means here: two
/// scores either are the same value or they are not, and a tolerance would
/// make the ranking depend on the order of comparisons.
#[allow(clippy::float_cmp, reason = "exact equality is the definition of a tie")]
fn ranks(values: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&i, &j| values[i].partial_cmp(&values[j]).unwrap_or(std::cmp::Ordering::Equal));
    let mut out = vec![0.0; values.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i + 1;
        while j < order.len() && values[order[j]] == values[order[i]] {
            j += 1;
        }
        // Ranks are 1-based; tied entries share their mean rank.
        let mean = ((i + 1 + j) as f64) / 2.0;
        for &k in &order[i..j] {
            out[k] = mean;
        }
        i = j;
    }
    out
}

/// Spearman rank correlation.
pub fn spearman(a: &[f64], b: &[f64]) -> Option<f64> {
    pearson(&ranks(a), &ranks(b))
}

/// Root mean squared error.
pub fn rmse(predicted: &[f64], actual: &[f64]) -> f64 {
    assert_eq!(predicted.len(), actual.len(), "rmse: length mismatch");
    let sum: f64 = predicted.iter().zip(actual).map(|(p, a)| (p - a).powi(2)).sum();
    (sum / predicted.len() as f64).sqrt()
}

/// Mean absolute error.
pub fn mae(predicted: &[f64], actual: &[f64]) -> f64 {
    assert_eq!(predicted.len(), actual.len(), "mae: length mismatch");
    predicted.iter().zip(actual).map(|(p, a)| (p - a).abs()).sum::<f64>() / predicted.len() as f64
}

/// Area under the ROC curve, computed from ranks so it handles ties.
/// Returns `None` unless both classes are present.
pub fn auroc(scores: &[f64], labels: &[bool]) -> Option<f64> {
    assert_eq!(scores.len(), labels.len(), "auroc: length mismatch");
    let positives = labels.iter().filter(|l| **l).count();
    let negatives = labels.len() - positives;
    if positives == 0 || negatives == 0 {
        return None;
    }
    // Mann-Whitney U via rank sums.
    let r = ranks(scores);
    let rank_sum: f64 = r.iter().zip(labels).filter(|(_, l)| **l).map(|(x, _)| *x).sum();
    let u = rank_sum - (positives * (positives + 1)) as f64 / 2.0;
    Some(u / (positives as f64 * negatives as f64))
}

/// Area under the precision-recall curve, by the trapezoidal rule over
/// thresholds at each distinct score. Returns `None` without any positives.
///
/// Scores are compared exactly: a threshold cannot separate two identical
/// scores, so they must be consumed together.
#[allow(clippy::float_cmp, reason = "identical scores must share a threshold")]
pub fn auprc(scores: &[f64], labels: &[bool]) -> Option<f64> {
    assert_eq!(scores.len(), labels.len(), "auprc: length mismatch");
    let total_positives = labels.iter().filter(|l| **l).count();
    if total_positives == 0 {
        return None;
    }
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&i, &j| scores[j].partial_cmp(&scores[i]).unwrap_or(std::cmp::Ordering::Equal));

    let (mut tp, mut fp) = (0.0, 0.0);
    let mut area = 0.0;
    let mut prev_recall = 0.0;
    let mut i = 0;
    while i < order.len() {
        // Advance past all entries sharing this score: a threshold cannot
        // separate them.
        let score = scores[order[i]];
        while i < order.len() && scores[order[i]] == score {
            if labels[order[i]] {
                tp += 1.0;
            } else {
                fp += 1.0;
            }
            i += 1;
        }
        let recall = tp / total_positives as f64;
        let precision = tp / (tp + fp);
        area += (recall - prev_recall) * precision;
        prev_recall = recall;
    }
    Some(area)
}

/// Counts at a decision threshold.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Confusion {
    /// Predicted positive and actually positive.
    pub true_positive: usize,
    /// Predicted positive but actually negative.
    pub false_positive: usize,
    /// Predicted negative and actually negative.
    pub true_negative: usize,
    /// Predicted negative but actually positive.
    pub false_negative: usize,
}

impl Confusion {
    /// Classify `scores` at `threshold` and count outcomes.
    pub fn at(scores: &[f64], labels: &[bool], threshold: f64) -> Self {
        assert_eq!(scores.len(), labels.len(), "confusion: length mismatch");
        let mut c = Self { true_positive: 0, false_positive: 0, true_negative: 0, false_negative: 0 };
        for (score, &label) in scores.iter().zip(labels) {
            match (*score >= threshold, label) {
                (true, true) => c.true_positive += 1,
                (true, false) => c.false_positive += 1,
                (false, false) => c.true_negative += 1,
                (false, true) => c.false_negative += 1,
            }
        }
        c
    }

    /// Matthews correlation coefficient. `None` when a denominator term is
    /// zero, which happens when everything lands in one row or column.
    pub fn mcc(&self) -> Option<f64> {
        let (tp, fp) = (self.true_positive as f64, self.false_positive as f64);
        let (tn, fln) = (self.true_negative as f64, self.false_negative as f64);
        let denom = ((tp + fp) * (tp + fln) * (tn + fp) * (tn + fln)).sqrt();
        if denom == 0.0 {
            return None;
        }
        Some((tp * tn - fp * fln) / denom)
    }

    /// Harmonic mean of precision and recall. `None` if both are zero.
    pub fn f1(&self) -> Option<f64> {
        let (tp, fp) = (self.true_positive as f64, self.false_positive as f64);
        let fln = self.false_negative as f64;
        if tp == 0.0 {
            return None;
        }
        let precision = tp / (tp + fp);
        let recall = tp / (tp + fln);
        Some(2.0 * precision * recall / (precision + recall))
    }

    /// Mean of sensitivity and specificity, which is not distorted by a
    /// skewed class balance the way plain accuracy is.
    pub fn balanced_accuracy(&self) -> Option<f64> {
        let positives = self.true_positive + self.false_negative;
        let negatives = self.true_negative + self.false_positive;
        if positives == 0 || negatives == 0 {
            return None;
        }
        let sensitivity = self.true_positive as f64 / positives as f64;
        let specificity = self.true_negative as f64 / negatives as f64;
        Some(f64::midpoint(sensitivity, specificity))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn pearson_matches_a_worked_example() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [2.0, 4.0, 6.0, 8.0, 10.0];
        assert!(close(pearson(&x, &y).unwrap(), 1.0));
        let z = [10.0, 8.0, 6.0, 4.0, 2.0];
        assert!(close(pearson(&x, &z).unwrap(), -1.0));
        assert_eq!(pearson(&x, &[1.0; 5]), None, "zero variance is undefined");
    }

    #[test]
    fn spearman_is_rank_based_not_value_based() {
        // Monotone but strongly non-linear: Pearson < 1, Spearman == 1.
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 4.0, 9.0, 16.0, 2500.0];
        assert!(close(spearman(&x, &y).unwrap(), 1.0));
        assert!(pearson(&x, &y).unwrap() < 0.8);
    }

    #[test]
    fn ranks_average_ties() {
        assert_eq!(ranks(&[10.0, 20.0, 20.0, 30.0]), vec![1.0, 2.5, 2.5, 4.0]);
    }

    #[test]
    fn errors_match_hand_computation() {
        let p = [1.0, 2.0, 3.0];
        let a = [1.0, 4.0, 3.0];
        assert!(close(mae(&p, &a), 2.0 / 3.0));
        assert!(close(rmse(&p, &a), (4.0f64 / 3.0).sqrt()));
    }

    #[test]
    fn auroc_handles_perfect_random_and_inverted_rankings() {
        let labels = [false, false, true, true];
        assert!(close(auroc(&[0.1, 0.2, 0.3, 0.4], &labels).unwrap(), 1.0));
        assert!(close(auroc(&[0.4, 0.3, 0.2, 0.1], &labels).unwrap(), 0.0));
        // All scores tied: every pair is a coin flip.
        assert!(close(auroc(&[0.5; 4], &labels).unwrap(), 0.5));
        assert_eq!(auroc(&[0.1, 0.2], &[true, true]), None, "needs both classes");
    }

    #[test]
    fn auprc_of_a_perfect_ranking_is_one() {
        let labels = [true, true, false, false];
        assert!(close(auprc(&[0.9, 0.8, 0.2, 0.1], &labels).unwrap(), 1.0));
        // Always-positive baseline equals the positive rate.
        let skewed = [true, false, false, false];
        assert!(close(auprc(&[0.5; 4], &skewed).unwrap(), 0.25));
    }

    #[test]
    fn confusion_counts_and_derived_scores() {
        let scores = [0.9, 0.8, 0.3, 0.1];
        let labels = [true, false, true, false];
        let c = Confusion::at(&scores, &labels, 0.5);
        assert_eq!(c, Confusion { true_positive: 1, false_positive: 1, true_negative: 1, false_negative: 1 });
        assert!(close(c.mcc().unwrap(), 0.0));
        assert!(close(c.f1().unwrap(), 0.5));
        assert!(close(c.balanced_accuracy().unwrap(), 0.5));
    }

    #[test]
    fn perfect_classifier_has_unit_scores() {
        let c = Confusion::at(&[0.9, 0.1], &[true, false], 0.5);
        assert!(close(c.mcc().unwrap(), 1.0));
        assert!(close(c.f1().unwrap(), 1.0));
        assert!(close(c.balanced_accuracy().unwrap(), 1.0));
    }
}
