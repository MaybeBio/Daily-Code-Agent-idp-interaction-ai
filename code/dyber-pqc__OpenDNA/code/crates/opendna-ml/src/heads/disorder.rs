use super::Mlp;
use candle_core::{Module, Tensor};
use candle_nn::{conv1d, Conv1d, Conv1dConfig, VarBuilder};

/// Width of the convolution that smooths predictions along the chain.
/// Disorder comes in runs rather than isolated residues, so a little local
/// context helps; 9 residues is roughly the span of a short linker.
const KERNEL: usize = 9;
const CHANNELS: usize = 128;

/// Predicts, for each residue, whether it is intrinsically disordered.
///
/// A per-residue projection is followed by one convolution along the
/// sequence, then a per-residue output. The output is a **logit**; apply a
/// sigmoid for a probability.
#[derive(Debug)]
pub struct DisorderHead {
    project: Mlp,
    conv: Conv1d,
    output: Mlp,
    /// Fitted map from logit to calibrated probability, when one has been
    /// measured for these weights. `None` leaves the plain sigmoid.
    calibration: Option<crate::heads::Calibration>,
}

impl DisorderHead {
    /// Build the head.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn new(embedding_dim: usize, vb: VarBuilder<'_>) -> candle_core::Result<Self> {
        let cfg = Conv1dConfig { padding: KERNEL / 2, ..Conv1dConfig::default() };
        Ok(Self {
            project: Mlp::new(&[embedding_dim, CHANNELS], vb.pp("project"))?,
            conv: conv1d(CHANNELS, CHANNELS, KERNEL, cfg, vb.pp("conv"))?,
            output: Mlp::new(&[CHANNELS, 64, 1], vb.pp("output"))?,
            calibration: None,
        })
    }

    /// Per-residue logits for one sequence.
    ///
    /// `embeddings` is `(length, embedding_dim)`; the result is `(length,)`.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn forward(&self, embeddings: &Tensor) -> candle_core::Result<Tensor> {
        let hidden = self.project.forward(embeddings)?.gelu_erf()?;
        // Conv1d wants (batch, channels, length).
        let shaped = hidden.transpose(0, 1)?.unsqueeze(0)?.contiguous()?;
        let convolved = self.conv.forward(&shaped)?.gelu_erf()?;
        let per_residue = convolved.squeeze(0)?.transpose(0, 1)?.contiguous()?;
        self.output.forward(&per_residue)?.squeeze(1)
    }

    #[must_use]
    /// Attach a calibration map, which [`Self::predict`] will then apply.
    ///
    /// Taken at load time from the metadata beside the weights, so a head and
    /// its calibration cannot be separated by accident.
    pub fn with_calibration(mut self, calibration: Option<crate::heads::Calibration>) -> Self {
        self.calibration = calibration;
        self
    }

    /// Per-residue disorder **probabilities** for one sequence.
    ///
    /// `per_residue` is [`crate::esm2::Embedding::per_residue`], row-major
    /// `length x embedding_dim`. The sigmoid is applied here so that callers
    /// cannot forget it and read raw logits as probabilities.
    ///
    /// When the head carries a calibration map, it is applied in place of the
    /// plain sigmoid. The map is monotonic, so it changes the probabilities
    /// without changing their order, and every rank statistic -- AUROC
    /// included -- is unaffected.
    ///
    /// # Errors
    /// Propagates tensor errors, or reports a length that does not divide the
    /// embedding evenly.
    pub fn predict(
        &self,
        per_residue: &[f32],
        length: usize,
        device: &candle_core::Device,
    ) -> candle_core::Result<Vec<f32>> {
        if length == 0 {
            return Ok(Vec::new());
        }
        if !per_residue.len().is_multiple_of(length) {
            return Err(candle_core::Error::Msg(format!(
                "{} embedding values do not divide into {length} residues",
                per_residue.len()
            )));
        }
        let width = per_residue.len() / length;
        let tensor = Tensor::from_vec(per_residue.to_vec(), (length, width), device)?;
        let logits = self.forward(&tensor)?.to_vec1::<f32>()?;
        Ok(match &self.calibration {
            Some(map) => logits.into_iter().map(|l| map.apply(l)).collect(),
            None => logits.into_iter().map(|l| 1.0 / (1.0 + (-l).exp())).collect(),
        })
    }
}
