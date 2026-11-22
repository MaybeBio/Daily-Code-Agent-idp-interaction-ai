//! Model hyperparameters, read from a Hugging Face `config.json`.

use crate::MlError;
use serde::Deserialize;

/// ESM-2 hyperparameters.
///
/// Field names match the keys in `config.json` so the file deserialises
/// directly.
#[derive(Debug, Clone, Deserialize)]
pub struct Config {
    /// Number of tokens in the alphabet (33 for every published ESM-2).
    pub vocab_size: usize,
    /// Embedding width.
    pub hidden_size: usize,
    /// Number of transformer blocks.
    pub num_hidden_layers: usize,
    /// Attention heads per block.
    pub num_attention_heads: usize,
    /// Feed-forward width.
    pub intermediate_size: usize,
    /// Epsilon for every layer normalisation.
    pub layer_norm_eps: f64,
    /// Longest token sequence the position encoding supports.
    #[serde(default = "default_max_positions")]
    pub max_position_embeddings: usize,
    /// Must be `rotary`; ESM-2 has no learned position embeddings.
    pub position_embedding_type: String,
    /// Whether to apply ESM's inference-time embedding rescaling.
    #[serde(default)]
    pub token_dropout: bool,
    /// Whether a layer normalisation precedes the first block.
    #[serde(default)]
    pub emb_layer_norm_before: bool,
}

fn default_max_positions() -> usize {
    1026
}

impl Config {
    /// Width of one attention head.
    pub fn head_dim(&self) -> usize {
        self.hidden_size / self.num_attention_heads
    }

    /// Longest residue sequence that fits, allowing for `<cls>` and `<eos>`.
    pub fn max_residues(&self) -> usize {
        self.max_position_embeddings.saturating_sub(2)
    }

    /// Reject configurations this implementation does not handle, so a
    /// mismatch surfaces at load time instead of as wrong numbers.
    ///
    /// # Errors
    /// [`MlError::Config`] describing the unsupported setting.
    pub fn validate(&self) -> Result<(), MlError> {
        if self.position_embedding_type != "rotary" {
            return Err(MlError::Config(format!(
                "position_embedding_type is {:?}; only \"rotary\" is implemented",
                self.position_embedding_type
            )));
        }
        if !self.hidden_size.is_multiple_of(self.num_attention_heads) {
            return Err(MlError::Config(format!(
                "hidden_size {} is not divisible by num_attention_heads {}",
                self.hidden_size, self.num_attention_heads
            )));
        }
        if !self.head_dim().is_multiple_of(2) {
            return Err(MlError::Config(format!(
                "head dimension {} is odd; rotary embeddings need an even width",
                self.head_dim()
            )));
        }
        if self.vocab_size != super::alphabet::TOKENS.len() {
            return Err(MlError::Config(format!(
                "vocab_size {} does not match the built-in ESM-2 alphabet ({} tokens)",
                self.vocab_size,
                super::alphabet::TOKENS.len()
            )));
        }
        Ok(())
    }
}
