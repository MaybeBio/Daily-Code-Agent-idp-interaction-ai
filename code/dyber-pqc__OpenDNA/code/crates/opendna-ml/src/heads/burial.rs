use super::Mlp;
use candle_core::{Module, Tensor};
use candle_nn::{conv1d, Conv1d, Conv1dConfig, VarBuilder};

/// Width of the convolution that smooths predictions along the chain.
/// Burial is a property of a residue's neighbourhood as much as of the
/// residue, and a turn of helix is about 4 residues; 9 covers two.
const KERNEL: usize = 9;
const CHANNELS: usize = 128;

/// Predicts how solvent-exposed each residue is.
///
/// The output is **relative solvent accessibility**: the residue's accessible
/// area over the theoretical maximum for that amino acid, so roughly 0 for a
/// residue in the core and roughly 1 for one fully exposed on the surface.
///
/// # Why the output is squashed
///
/// A plain linear output would predict negative areas, which are not a thing.
/// Relative accessibility is bounded below by 0 and in practice sits under
/// about 1.2 — the Tien et al. maxima are theoretical, so a residue with an
/// unusually extended side chain can exceed its own maximum. The head
/// therefore emits a sigmoid scaled to `[0, 1.25]`, which cannot produce an
/// impossible value and still reaches every value the labels contain.
const MAX_RSA: f64 = 1.25;

/// Predicts per-residue relative solvent accessibility.
///
/// A per-residue projection, one convolution along the sequence, then a
/// per-residue output. Deliberately the same shape as
/// [`super::DisorderHead`]: both predict a smooth per-residue quantity from
/// the same embeddings, and using one architecture for both means a
/// difference in results is a difference in the task rather than in the
/// model.
#[derive(Debug)]
pub struct BurialHead {
    project: Mlp,
    conv: Conv1d,
    output: Mlp,
}

impl BurialHead {
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
        })
    }

    /// Per-residue **logits** for one sequence.
    ///
    /// `embeddings` is `(length, embedding_dim)`; the result is `(length,)`.
    /// Training uses these directly so the loss sees an unsquashed value;
    /// [`Self::predict`] applies the squashing.
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

    /// Turn logits into relative accessibility in `[0, MAX_RSA]`.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn accessibility(logits: &Tensor) -> candle_core::Result<Tensor> {
        candle_nn::ops::sigmoid(logits)?.affine(MAX_RSA, 0.0)
    }

    /// Per-residue relative solvent accessibility for one sequence.
    ///
    /// `per_residue` is [`crate::esm2::Embedding::per_residue`], row-major
    /// `length x embedding_dim`. The squashing is applied here so a caller
    /// cannot read a raw logit as an area.
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
        let input = Tensor::from_vec(per_residue.to_vec(), (length, width), device)?;
        Self::accessibility(&self.forward(&input)?)?.to_vec1()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use candle_core::{DType, Device};
    use candle_nn::VarMap;

    /// A negative accessible area is not a thing, and neither is one far
    /// above the theoretical maximum. The squashing has to make both
    /// unreachable whatever the weights do.
    #[test]
    fn accessibility_stays_within_its_physical_range() {
        let device = Device::Cpu;
        let extreme = Tensor::from_vec(vec![-1e3f32, -5.0, 0.0, 5.0, 1e3], 5, &device).expect("tensor");
        let got = BurialHead::accessibility(&extreme).expect("squash").to_vec1::<f32>().expect("vec");
        assert!(got.iter().all(|v| (0.0..=1.25).contains(v)), "{got:?}");
        assert!(got[0] < 1e-3, "a large negative logit should approach 0: {got:?}");
        assert!(got[4] > 1.24, "a large positive logit should approach the maximum: {got:?}");
        // Monotonic, so ordering a set of residues by logit orders them by
        // accessibility.
        assert!(got.windows(2).all(|w| w[0] <= w[1]), "{got:?}");
    }

    #[test]
    fn predict_returns_one_value_per_residue() {
        let device = Device::Cpu;
        let map = VarMap::new();
        let vb = VarBuilder::from_varmap(&map, DType::F32, &device);
        let head = BurialHead::new(8, vb).expect("head");
        let got = head.predict(&[0.1f32; 8 * 5], 5, &device).expect("predict");
        assert_eq!(got.len(), 5);
        assert!(got.iter().all(|v| (0.0..=1.25).contains(v)), "{got:?}");
    }

    #[test]
    fn an_embedding_that_does_not_divide_is_refused() {
        let device = Device::Cpu;
        let map = VarMap::new();
        let vb = VarBuilder::from_varmap(&map, DType::F32, &device);
        let head = BurialHead::new(8, vb).expect("head");
        assert!(head.predict(&[0.1f32; 23], 5, &device).is_err());
        assert_eq!(head.predict(&[], 0, &device).expect("empty").len(), 0);
    }
}
