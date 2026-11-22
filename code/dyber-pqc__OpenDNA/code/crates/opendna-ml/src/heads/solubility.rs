use super::Mlp;
use candle_core::{Module, Tensor};
use candle_nn::VarBuilder;

/// Predicts protein solubility from a mean-pooled sequence embedding.
///
/// The training target is the eSOL solubility percentage, so the output is on
/// that scale. Values above 100 occur in the source data and are not clipped.
#[derive(Debug)]
pub struct SolubilityHead {
    mlp: Mlp,
}

impl SolubilityHead {
    /// Hidden widths used by the released model.
    pub fn hidden_widths() -> [usize; 2] {
        [256, 64]
    }

    /// Build the head.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn new(embedding_dim: usize, vb: VarBuilder<'_>) -> candle_core::Result<Self> {
        let [a, b] = Self::hidden_widths();
        Ok(Self { mlp: Mlp::new(&[embedding_dim, a, b, 1], vb)? })
    }

    /// Predict solubility for each pooled embedding, returning `(batch,)`.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn forward(&self, pooled: &Tensor) -> candle_core::Result<Tensor> {
        self.mlp.forward(pooled)?.squeeze(1)
    }

    /// Predict solubility from one mean-pooled embedding.
    ///
    /// Pass [`crate::esm2::Embedding::pooled`] directly: the head was trained
    /// on that exact pooling, and recomputing it elsewhere risks the two
    /// drifting apart.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn predict(&self, pooled: &[f32], device: &candle_core::Device) -> candle_core::Result<f32> {
        let tensor = Tensor::from_vec(pooled.to_vec(), (1, pooled.len()), device)?;
        let out = self.forward(&tensor)?.to_vec1::<f32>()?;
        out.first()
            .copied()
            .ok_or_else(|| candle_core::Error::Msg("the solubility head returned no value".to_string()))
    }
}
