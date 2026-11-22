//! The ESM-2 transformer.
//!
//! Follows the reference implementation (Lin et al., 2023) as published in
//! `transformers`' `modeling_esm.py`. Three details are easy to get wrong and
//! are called out where they appear:
//!
//! 1. **Token dropout.** At inference ESM-2 rescales the token embeddings by
//!    `1 - 0.12`. Omitting it changes every downstream number.
//! 2. **Pre-normalisation.** Each block normalises its input, not its output.
//! 3. **Exact GELU.** The feed-forward uses the erf form, not the tanh
//!    approximation; the two differ by enough to fail a tight comparison.

use super::config::Config;
use candle_core::{DType, Device, IndexOp, Tensor, D};
use candle_nn::{embedding, layer_norm, linear, Embedding, LayerNorm, Linear, Module, VarBuilder};

/// Fraction of tokens masked during ESM-2 pre-training (15% selected, 80% of
/// those replaced by `<mask>`). Used by the inference-time rescaling.
const TRAIN_MASK_RATIO: f64 = 0.15 * 0.8;

/// Rotary position embeddings, applied to queries and keys in every block.
///
/// The inverse frequencies are taken from the checkpoint rather than
/// recomputed. Published ESM-2 checkpoints ship them rounded to bfloat16
/// precision (0.31616 where the formula gives 0.316228), and the reference
/// implementation uses the stored values; recomputing them shifts every
/// downstream activation by ~1e-3.
struct Rotary {
    cos: Tensor,
    sin: Tensor,
}

impl Rotary {
    /// Inverse frequencies from the formula, for checkpoints that omit the
    /// buffer. Matches the reference only up to bfloat16 rounding.
    fn default_inv_freq(head_dim: usize, device: &Device) -> candle_core::Result<Tensor> {
        let half = head_dim / 2;
        let values: Vec<f32> =
            (0..half).map(|i| 1f32 / 10_000f32.powf(2.0 * i as f32 / head_dim as f32)).collect();
        Tensor::from_vec(values, (1, half), device)
    }

    /// Build `cos`/`sin` tables of shape `(1, 1, seq_len, head_dim)` from
    /// `inv_freq`, which has `head_dim / 2` elements.
    fn new(inv_freq: &Tensor, seq_len: usize, dtype: DType, device: &Device) -> candle_core::Result<Self> {
        let half = inv_freq.elem_count();
        let inv_freq = inv_freq.reshape((1, half))?.to_dtype(DType::F32)?;
        let positions: Vec<f32> = (0..seq_len).map(|p| p as f32).collect();
        let positions = Tensor::from_vec(positions, (seq_len, 1), device)?;
        // (seq_len, half), then duplicated to (seq_len, head_dim) to match
        // the reference's `cat((freqs, freqs), -1)`.
        let freqs = positions.matmul(&inv_freq)?;
        let emb = Tensor::cat(&[&freqs, &freqs], D::Minus1)?;
        let shape = (1, 1, seq_len, half * 2);
        Ok(Self {
            cos: emb.cos()?.reshape(shape)?.to_dtype(dtype)?,
            sin: emb.sin()?.reshape(shape)?.to_dtype(dtype)?,
        })
    }

    /// `(-x2, x1)` over the two halves of the last axis.
    fn rotate_half(x: &Tensor) -> candle_core::Result<Tensor> {
        let half = x.dim(D::Minus1)? / 2;
        let x1 = x.narrow(D::Minus1, 0, half)?;
        let x2 = x.narrow(D::Minus1, half, half)?;
        Tensor::cat(&[&x2.neg()?, &x1], D::Minus1)
    }

    /// `x * cos + rotate_half(x) * sin`, broadcast over batch and heads.
    fn apply(&self, x: &Tensor) -> candle_core::Result<Tensor> {
        let rotated = Self::rotate_half(x)?;
        x.broadcast_mul(&self.cos)?.add(&rotated.broadcast_mul(&self.sin)?)
    }
}

struct SelfAttention {
    query: Linear,
    key: Linear,
    value: Linear,
    output: Linear,
    layer_norm: LayerNorm,
    inv_freq: Tensor,
    num_heads: usize,
    head_dim: usize,
}

impl SelfAttention {
    fn load(cfg: &Config, vb: VarBuilder<'_>) -> candle_core::Result<Self> {
        let h = cfg.hidden_size;
        let self_vb = vb.pp("self");
        let head_dim = cfg.head_dim();
        let inv_freq = self_vb
            .pp("rotary_embeddings")
            .get(head_dim / 2, "inv_freq")
            .or_else(|_| Rotary::default_inv_freq(head_dim, vb.device()))?;
        Ok(Self {
            query: linear(h, h, self_vb.pp("query"))?,
            key: linear(h, h, self_vb.pp("key"))?,
            value: linear(h, h, self_vb.pp("value"))?,
            output: linear(h, h, vb.pp("output").pp("dense"))?,
            layer_norm: layer_norm(h, cfg.layer_norm_eps, vb.pp("LayerNorm"))?,
            inv_freq,
            num_heads: cfg.num_attention_heads,
            head_dim,
        })
    }

    /// `(batch, seq, hidden)` -> `(batch, heads, seq, head_dim)`.
    fn split_heads(&self, x: &Tensor) -> candle_core::Result<Tensor> {
        let (b, s, _) = x.dims3()?;
        x.reshape((b, s, self.num_heads, self.head_dim))?.transpose(1, 2)?.contiguous()
    }

    /// `mask` is additive: 0 for real tokens, a large negative value for padding.
    fn forward(&self, hidden: &Tensor, mask: Option<&Tensor>) -> candle_core::Result<Tensor> {
        let (b, s, _) = hidden.dims3()?;
        let rotary = Rotary::new(&self.inv_freq, s, hidden.dtype(), hidden.device())?;
        // Pre-normalisation: the block normalises its input and keeps the
        // unnormalised value for the residual.
        let normed = self.layer_norm.forward(hidden)?;

        // The reference scales the query before the rotation. Rotation is
        // orthogonal so the order does not matter, but matching it keeps the
        // arithmetic identical.
        let scale = 1.0 / (self.head_dim as f64).sqrt();
        let q = (self.split_heads(&self.query.forward(&normed)?)? * scale)?;
        let k = self.split_heads(&self.key.forward(&normed)?)?;
        let v = self.split_heads(&self.value.forward(&normed)?)?;

        let q = rotary.apply(&q)?;
        let k = rotary.apply(&k)?;

        let mut scores = q.matmul(&k.transpose(D::Minus1, D::Minus2)?.contiguous()?)?;
        if let Some(mask) = mask {
            scores = scores.broadcast_add(mask)?;
        }
        let probs = candle_nn::ops::softmax_last_dim(&scores)?;
        let context = probs.matmul(&v)?.transpose(1, 2)?.reshape((b, s, self.num_heads * self.head_dim))?;
        // Residual added to the *unnormalised* input.
        self.output.forward(&context)?.add(hidden)
    }
}

struct Block {
    attention: SelfAttention,
    intermediate: Linear,
    output: Linear,
    layer_norm: LayerNorm,
}

impl Block {
    fn load(cfg: &Config, vb: VarBuilder<'_>) -> candle_core::Result<Self> {
        Ok(Self {
            attention: SelfAttention::load(cfg, vb.pp("attention"))?,
            intermediate: linear(cfg.hidden_size, cfg.intermediate_size, vb.pp("intermediate").pp("dense"))?,
            output: linear(cfg.intermediate_size, cfg.hidden_size, vb.pp("output").pp("dense"))?,
            layer_norm: layer_norm(cfg.hidden_size, cfg.layer_norm_eps, vb.pp("LayerNorm"))?,
        })
    }

    fn forward(&self, hidden: &Tensor, mask: Option<&Tensor>) -> candle_core::Result<Tensor> {
        let attended = self.attention.forward(hidden, mask)?;
        let normed = self.layer_norm.forward(&attended)?;
        // Exact erf GELU, not the tanh approximation.
        let ffn = self.output.forward(&self.intermediate.forward(&normed)?.gelu_erf()?)?;
        ffn.add(&attended)
    }
}

/// Masked-language-model head. Only needed for per-residue token scoring;
/// plain embedding does not touch it.
struct LmHead {
    dense: Linear,
    layer_norm: LayerNorm,
    decoder: Tensor,
    bias: Tensor,
}

impl LmHead {
    fn load(cfg: &Config, vb: VarBuilder<'_>, word_embeddings: &Tensor) -> candle_core::Result<Self> {
        let h = cfg.hidden_size;
        // The decoder is weight-tied to the input embedding and is usually
        // absent from the checkpoint; fall back to the tied tensor.
        let decoder =
            vb.pp("decoder").get((cfg.vocab_size, h), "weight").unwrap_or_else(|_| word_embeddings.clone());
        Ok(Self {
            dense: linear(h, h, vb.pp("dense"))?,
            layer_norm: layer_norm(h, cfg.layer_norm_eps, vb.pp("layer_norm"))?,
            decoder,
            bias: vb.get(cfg.vocab_size, "bias")?,
        })
    }

    fn forward(&self, hidden: &Tensor) -> candle_core::Result<Tensor> {
        let x = self.layer_norm.forward(&self.dense.forward(hidden)?.gelu_erf()?)?;
        x.broadcast_matmul(&self.decoder.t()?)?.broadcast_add(&self.bias)
    }
}

/// A loaded ESM-2 model.
pub struct Model {
    embeddings: Embedding,
    emb_layer_norm_before: Option<LayerNorm>,
    blocks: Vec<Block>,
    emb_layer_norm_after: LayerNorm,
    lm_head: LmHead,
    config: Config,
    device: Device,
    dtype: DType,
}

impl Model {
    /// Load from a [`VarBuilder`] rooted at the checkpoint's top level.
    ///
    /// # Errors
    /// Propagates tensor errors when a weight is missing or mis-shaped.
    pub fn load(config: Config, vb: VarBuilder<'_>) -> candle_core::Result<Self> {
        let device = vb.device().clone();
        let dtype = vb.dtype();
        let esm = vb.pp("esm");
        let emb_vb = esm.pp("embeddings");
        let embeddings = embedding(config.vocab_size, config.hidden_size, emb_vb.pp("word_embeddings"))?;
        let encoder = esm.pp("encoder");
        let emb_layer_norm_before = if config.emb_layer_norm_before {
            Some(layer_norm(config.hidden_size, config.layer_norm_eps, encoder.pp("emb_layer_norm_before"))?)
        } else {
            None
        };
        let blocks = (0..config.num_hidden_layers)
            .map(|i| Block::load(&config, encoder.pp("layer").pp(i)))
            .collect::<candle_core::Result<Vec<_>>>()?;
        let emb_layer_norm_after =
            layer_norm(config.hidden_size, config.layer_norm_eps, encoder.pp("emb_layer_norm_after"))?;
        let lm_head = LmHead::load(&config, vb.pp("lm_head"), embeddings.embeddings())?;
        Ok(Self {
            embeddings,
            emb_layer_norm_before,
            blocks,
            emb_layer_norm_after,
            lm_head,
            config,
            device,
            dtype,
        })
    }

    /// Model hyperparameters.
    pub fn config(&self) -> &Config {
        &self.config
    }

    /// Device the weights live on.
    pub fn device(&self) -> &Device {
        &self.device
    }

    /// Run the encoder.
    ///
    /// `tokens` is `(batch, seq)` of token ids and `attention` is `(batch, seq)`
    /// with 1 for real tokens and 0 for padding. Returns the final hidden
    /// states, `(batch, seq, hidden)`, after the trailing layer normalisation.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn encode(&self, tokens: &Tensor, attention: &Tensor) -> candle_core::Result<Tensor> {
        let mut hidden = self.embeddings.forward(tokens)?;

        if self.config.token_dropout {
            hidden = self.apply_token_dropout(&hidden, tokens, attention)?;
        }
        if let Some(ln) = &self.emb_layer_norm_before {
            hidden = ln.forward(&hidden)?;
        }
        // Zero the padded positions so they cannot leak through the residual stream.
        let keep = attention.to_dtype(self.dtype)?.unsqueeze(D::Minus1)?;
        hidden = hidden.broadcast_mul(&keep)?;

        let mask = self.additive_mask(attention)?;
        for block in &self.blocks {
            hidden = block.forward(&hidden, mask.as_ref())?;
        }
        self.emb_layer_norm_after.forward(&hidden)
    }

    /// Token logits, `(batch, seq, vocab)`.
    ///
    /// # Errors
    /// Propagates tensor errors.
    pub fn logits(&self, tokens: &Tensor, attention: &Tensor) -> candle_core::Result<Tensor> {
        self.lm_head.forward(&self.encode(tokens, attention)?)
    }

    /// ESM-2 rescales token embeddings at inference to match the statistics
    /// seen during masked pre-training: masked positions are zeroed, then
    /// everything is scaled by `(1 - 0.12) / (1 - observed_mask_ratio)`.
    fn apply_token_dropout(
        &self,
        hidden: &Tensor,
        tokens: &Tensor,
        attention: &Tensor,
    ) -> candle_core::Result<Tensor> {
        let is_mask = tokens.eq(super::alphabet::MASK)?.to_dtype(self.dtype)?;
        // Zero the embedding at masked positions.
        let keep = (1.0 - &is_mask)?.unsqueeze(D::Minus1)?;
        let hidden = hidden.broadcast_mul(&keep)?;

        let lengths = attention.to_dtype(self.dtype)?.sum(D::Minus1)?;
        let masked = is_mask.sum(D::Minus1)?;
        let observed = masked.div(&lengths)?;
        let scale = ((1.0 - TRAIN_MASK_RATIO) / (1.0 - observed)?)?;
        hidden.broadcast_mul(&scale.reshape((scale.elem_count(), 1, 1))?)
    }

    /// Additive attention mask, or `None` when nothing is padded.
    fn additive_mask(&self, attention: &Tensor) -> candle_core::Result<Option<Tensor>> {
        let total = attention.sum_all()?.to_scalar::<u32>()? as usize;
        if total == attention.elem_count() {
            return Ok(None);
        }
        let (b, s) = attention.dims2()?;
        // Large finite penalty rather than -inf: a fully padded row would
        // otherwise produce NaN in the softmax.
        let blocked = (1.0 - attention.to_dtype(self.dtype)?)?;
        Ok(Some((blocked * -1e30)?.reshape((b, 1, 1, s))?))
    }

    /// Mean of the per-residue hidden states, excluding `<cls>`, `<eos>` and
    /// padding. This is the vector the prediction heads are trained on.
    ///
    /// # Errors
    /// Propagates tensor errors.
    #[allow(clippy::unused_self, reason = "method for symmetry with the rest of the model API")]
    pub fn pool(&self, hidden: &Tensor, lengths: &[usize]) -> candle_core::Result<Tensor> {
        let mut rows = Vec::with_capacity(lengths.len());
        for (i, &len) in lengths.iter().enumerate() {
            // Residues occupy positions 1..=len.
            rows.push(hidden.i(i)?.narrow(0, 1, len)?.mean(0)?);
        }
        Tensor::stack(&rows, 0)
    }
}
