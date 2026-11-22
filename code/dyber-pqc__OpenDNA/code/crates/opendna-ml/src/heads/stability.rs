use super::Mlp;
use candle_core::{Module, Tensor};
use candle_nn::VarBuilder;

/// Number of standard amino acids, the width of each one-hot block.
const AMINO_ACIDS: usize = 20;
/// Wild-type one-hot, mutant one-hot, and the language model's log-odds.
const EXTRA_FEATURES: usize = AMINO_ACIDS * 2 + 1;
/// Predicted relative solvent accessibility at the mutated position, when the
/// head was trained with it.
///
/// This carries no information the embedding does not already contain — it is
/// predicted *from* that embedding by [`super::BurialHead`]. What it adds is
/// supervision from somewhere else: the burial head learned from deposited
/// coordinates, and distilling that into one number the small head can use
/// directly is cheaper than expecting this MLP to rediscover it from 389,068
/// mutations. Whether that helps is an empirical question, which is why it is
/// optional and why both variants get benchmarked.
const BURIAL_FEATURE: usize = 1;

/// Predicts the change in folding free energy, ΔΔG, for a point mutation.
///
/// Features per mutation, concatenated:
/// - the wild-type sequence's embedding **at the mutated position**;
/// - one-hot encodings of the wild-type and mutant residues;
/// - the masked-language-model log-odds `log P(mutant) - log P(wild type)` at
///   that position, which is a strong predictor on its own (Meier et al., 2021)
///   and lets the head start from a sensible zero-shot signal.
///
/// The sign convention is the dataset's: negative ΔΔG is destabilising.
#[derive(Debug)]
pub struct StabilityHead {
    mlp: Mlp,
    /// Whether this head's first layer includes the burial column.
    ///
    /// Held because [`Self::predict`] validates the row width against what
    /// *this* head takes, not against the plain formula. Without it, a head
    /// built with burial rejected the very rows it needs — which is how this
    /// was caught.
    with_burial: bool,
}

/// Assembled features for a batch of mutations.
#[derive(Debug, Clone)]
pub struct StabilityFeatures {
    /// Row-major `(batch, embedding_dim + 41)`.
    pub values: Vec<f32>,
    /// Number of mutations.
    pub batch: usize,
    /// Width of one row.
    pub width: usize,
}

impl StabilityHead {
    /// Width of one input row for a given embedding size.
    pub fn input_width(embedding_dim: usize) -> usize {
        embedding_dim + EXTRA_FEATURES
    }

    /// Width of one input row when burial is included.
    pub fn input_width_with_burial(embedding_dim: usize) -> usize {
        embedding_dim + EXTRA_FEATURES + BURIAL_FEATURE
    }

    /// Column holding the language model's log-odds.
    ///
    /// Exposed because the training code needs it to report the zero-shot
    /// baseline, and it used to find it by taking the last column. That
    /// worked until burial was appended after it, at which point the
    /// "baseline" silently became burial-against-ΔΔG — a plausible number
    /// measuring the wrong thing. The layout is defined here, so the index
    /// is answered here.
    pub fn log_odds_index(embedding_dim: usize) -> usize {
        embedding_dim + AMINO_ACIDS * 2
    }

    /// Hidden widths used by the released model.
    pub fn hidden_widths() -> [usize; 2] {
        [256, 64]
    }

    /// Build the head.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn new(embedding_dim: usize, vb: VarBuilder<'_>) -> candle_core::Result<Self> {
        Self::new_with(embedding_dim, false, vb)
    }

    /// Build the head, with or without the burial feature.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn new_with(
        embedding_dim: usize,
        with_burial: bool,
        vb: VarBuilder<'_>,
    ) -> candle_core::Result<Self> {
        let width = if with_burial {
            Self::input_width_with_burial(embedding_dim)
        } else {
            Self::input_width(embedding_dim)
        };
        let [a, b] = Self::hidden_widths();
        Ok(Self { mlp: Mlp::new(&[width, a, b, 1], vb)?, with_burial })
    }

    /// Predict ΔΔG for each row of `features`, returning `(batch,)`.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn forward(&self, features: &Tensor) -> candle_core::Result<Tensor> {
        self.mlp.forward(features)?.squeeze(1)
    }

    /// Predict ΔΔG for each assembled feature row.
    ///
    /// Takes and returns plain slices so callers do not need candle: the
    /// desktop application depends on this crate, not on a tensor library.
    ///
    /// # Errors
    /// Propagates tensor errors, or reports a row width that does not match
    /// the head.
    pub fn predict(
        &self,
        rows: &[Vec<f32>],
        embedding_dim: usize,
        device: &candle_core::Device,
    ) -> candle_core::Result<Vec<f32>> {
        if rows.is_empty() {
            return Ok(Vec::new());
        }
        let width = if self.with_burial {
            Self::input_width_with_burial(embedding_dim)
        } else {
            Self::input_width(embedding_dim)
        };
        if let Some(bad) = rows.iter().position(|r| r.len() != width) {
            return Err(candle_core::Error::Msg(format!(
                "stability feature row {bad} is {} wide; this head expects {width}",
                rows[bad].len()
            )));
        }
        let flat: Vec<f32> = rows.iter().flatten().copied().collect();
        let tensor = Tensor::from_vec(flat, (rows.len(), width), device)?;
        self.forward(&tensor)?.to_vec1::<f32>()
    }

    /// Assemble one feature row.
    ///
    /// `embedding` is the wild-type embedding at the mutated position.
    /// `wild_type` and `mutant` are indices into the 20 standard amino acids.
    /// `log_odds` is the language model's `log P(mutant) - log P(wild type)`.
    pub fn features(embedding: &[f32], wild_type: usize, mutant: usize, log_odds: f32) -> Vec<f32> {
        Self::features_with(embedding, wild_type, mutant, log_odds, None)
    }

    /// Features, optionally with predicted accessibility at the mutated
    /// position appended.
    ///
    /// `burial` must be `Some` exactly when the head was built with
    /// `with_burial`, or the row will be the wrong width and the matrix
    /// multiply will fail — loudly, which is the intent.
    pub fn features_with(
        embedding: &[f32],
        wild_type: usize,
        mutant: usize,
        log_odds: f32,
        burial: Option<f32>,
    ) -> Vec<f32> {
        let mut row = Vec::with_capacity(embedding.len() + EXTRA_FEATURES + BURIAL_FEATURE);
        row.extend_from_slice(embedding);
        for index in [wild_type, mutant] {
            let start = row.len();
            row.resize(start + AMINO_ACIDS, 0.0);
            if index < AMINO_ACIDS {
                row[start + index] = 1.0;
            }
        }
        row.push(log_odds);
        if let Some(value) = burial {
            row.push(value);
        }
        row
    }
}

#[cfg(test)]
mod burial_variant_tests {
    use super::*;
    use candle_core::{DType, Device};
    use candle_nn::VarMap;

    fn head(with_burial: bool) -> (StabilityHead, Device) {
        let device = Device::Cpu;
        let map = VarMap::new();
        let vb = VarBuilder::from_varmap(&map, DType::F32, &device);
        (StabilityHead::new_with(8, with_burial, vb).expect("head"), device)
    }

    /// The variant head must accept the wider row and the plain head must
    /// not. Before the head remembered which it was, `predict` validated
    /// against the plain formula and a burial head rejected the only rows it
    /// can use — a benchmark that could not run at all.
    #[test]
    fn each_variant_accepts_only_its_own_row_width() {
        let plain_row = StabilityHead::features(&[0.0; 8], 0, 1, 0.5);
        let burial_row = StabilityHead::features_with(&[0.0; 8], 0, 1, 0.5, Some(0.3));
        assert_eq!(plain_row.len(), StabilityHead::input_width(8));
        assert_eq!(burial_row.len(), StabilityHead::input_width_with_burial(8));
        assert_eq!(burial_row.len(), plain_row.len() + 1);

        let (plain, device) = head(false);
        assert!(plain.predict(std::slice::from_ref(&plain_row), 8, &device).is_ok());
        assert!(
            plain.predict(std::slice::from_ref(&burial_row), 8, &device).is_err(),
            "plain took a burial row"
        );

        let (variant, device) = head(true);
        assert!(variant.predict(std::slice::from_ref(&burial_row), 8, &device).is_ok());
        assert!(
            variant.predict(std::slice::from_ref(&plain_row), 8, &device).is_err(),
            "variant took a plain row"
        );
    }

    /// The log-odds column must stay findable once burial is appended after
    /// it, or the zero-shot baseline silently measures burial instead.
    #[test]
    fn the_log_odds_column_is_where_it_says_it_is() {
        let row = StabilityHead::features_with(&[0.0; 8], 0, 1, 0.75, Some(0.3));
        let at = StabilityHead::log_odds_index(8);
        // Compared with a tolerance only to satisfy the float_cmp lint; these
        // values are stored verbatim, so they are in fact bit-identical.
        assert!((row[at] - 0.75).abs() < 1e-6, "log-odds not at its stated index");
        assert!((*row.last().expect("row") - 0.3).abs() < 1e-6, "burial should be last");
        assert_ne!(at, row.len() - 1, "the last column is burial, not log-odds");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::float_cmp, reason = "one-hot entries are exactly 0.0 or 1.0")]
    fn feature_row_has_expected_layout() {
        let emb = [0.5f32; 8];
        let row = StabilityHead::features(&emb, 3, 7, -1.5);
        assert_eq!(row.len(), StabilityHead::input_width(8));
        assert_eq!(&row[..8], &emb);
        assert_eq!(row[8 + 3], 1.0, "wild-type one-hot");
        assert_eq!(row[8 + AMINO_ACIDS + 7], 1.0, "mutant one-hot");
        assert_eq!(row[row.len() - 1], -1.5, "log-odds is the last feature");
        assert_eq!(row[8..8 + AMINO_ACIDS].iter().sum::<f32>(), 1.0);
    }
}
