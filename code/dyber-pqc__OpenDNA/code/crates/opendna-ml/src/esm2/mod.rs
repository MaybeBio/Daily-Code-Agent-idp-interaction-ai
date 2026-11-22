//! ESM-2, a protein language model (Lin et al., 2023), ported to candle.
//!
//! ```no_run
//! # fn main() -> Result<(), opendna_ml::MlError> {
//! use opendna_ml::esm2::{Esm2, ModelSize};
//! let model = Esm2::from_pretrained(ModelSize::Small8M, &opendna_ml::best_device()?)?;
//! let out = model.embed(&["MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLE"])?;
//! assert_eq!(out[0].pooled.len(), model.hidden_size());
//! # Ok(()) }
//! ```

pub mod alphabet;
mod config;
mod model;

pub use config::Config;

use crate::{weights, MlError};
use candle_core::{DType, Device, Tensor};
use candle_nn::VarBuilder;
use std::path::Path;

/// Published ESM-2 checkpoints.
///
/// Larger models embed better but cost proportionally more time and memory.
/// The 8M model is used by the test suite because its activations are small
/// enough to commit as fixtures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModelSize {
    /// 6 layers, 320 dimensions, ~8M parameters.
    Small8M,
    /// 12 layers, 480 dimensions, ~35M parameters.
    Medium35M,
    /// 30 layers, 640 dimensions, ~150M parameters.
    Large150M,
    /// 33 layers, 1280 dimensions, ~650M parameters.
    XLarge650M,
}

impl ModelSize {
    /// Hugging Face repository identifier.
    pub fn repo(self) -> &'static str {
        match self {
            Self::Small8M => "facebook/esm2_t6_8M_UR50D",
            Self::Medium35M => "facebook/esm2_t12_35M_UR50D",
            Self::Large150M => "facebook/esm2_t30_150M_UR50D",
            Self::XLarge650M => "facebook/esm2_t33_650M_UR50D",
        }
    }

    /// Embedding width, for sizing downstream heads without loading weights.
    pub fn hidden_size(self) -> usize {
        match self {
            Self::Small8M => 320,
            Self::Medium35M => 480,
            Self::Large150M => 640,
            Self::XLarge650M => 1280,
        }
    }
}

/// Embeddings for one sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct Embedding {
    /// Mean over residues, length `hidden_size`. Input to whole-sequence heads.
    pub pooled: Vec<f32>,
    /// One vector per residue, `length × hidden_size`, row-major. Input to
    /// per-residue heads such as disorder.
    pub per_residue: Vec<f32>,
    /// Residue count, so `per_residue` can be reshaped.
    pub length: usize,
}

impl Embedding {
    /// Borrow the embedding of residue `index` (0-based).
    pub fn residue(&self, index: usize) -> Option<&[f32]> {
        let width = self.per_residue.len() / self.length.max(1);
        self.per_residue.get(index * width..(index + 1) * width)
    }
}

/// Group sequence indices into batches that stay inside a token budget.
///
/// A batch costs `members x longest member` padded tokens, so the grouping has
/// to consider length and not just count: a 60-residue sequence batched with a
/// 1000-residue one is computed as though it were 1000 long. Sorting by length
/// first puts similar lengths together, which is what keeps the waste small.
///
/// Ties are broken by index so the plan is deterministic, which matters
/// because padding changes floating-point summation order and therefore the
/// last few bits of the result. Given the same lengths and budget this returns
/// the same plan on every machine.
///
/// A sequence larger than the whole budget gets a batch to itself rather than
/// being refused: the budget is advice about memory, not a limit on what can
/// be embedded.
fn plan_batches(lengths: &[usize], token_budget: usize) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..lengths.len()).collect();
    order.sort_by_key(|&i| (lengths[i], i));

    let mut plan: Vec<Vec<usize>> = Vec::new();
    let mut batch: Vec<usize> = Vec::new();
    let mut widest = 0usize;
    for &index in &order {
        let widened = widest.max(lengths[index]);
        if !batch.is_empty() && (batch.len() + 1) * widened > token_budget {
            plan.push(std::mem::take(&mut batch));
            widest = 0;
        }
        widest = widest.max(lengths[index]);
        batch.push(index);
    }
    if !batch.is_empty() {
        plan.push(batch);
    }
    plan
}

/// A loaded ESM-2 model.
///
/// # Concurrency
///
/// Calls are independent and this type is `Sync`, but driving one instance
/// from several threads at once on a CUDA device changes results in the last
/// bits of an f32 — around 2e-7 relative, from cuBLAS choosing differently
/// shaped reductions under concurrent load. Sequential calls are
/// bit-identical. That is far below any meaningful signal, so it does not
/// make a prediction wrong, but it does mean results are only
/// bit-reproducible if inference is serialised. Where reproducibility
/// matters — building an embedding cache, comparing two runs — keep calls on
/// one thread or behind a lock, which is what the embedding pipeline and the
/// desktop application both do.
pub struct Esm2 {
    model: model::Model,
}

impl Esm2 {
    /// Download (or reuse the cached copy of) a published checkpoint.
    ///
    /// # Errors
    /// [`MlError::Download`] if the weights cannot be fetched, or
    /// [`MlError::Config`] if the checkpoint is not a supported ESM-2.
    pub fn from_pretrained(size: ModelSize, device: &Device) -> Result<Self, MlError> {
        let files = weights::fetch(size.repo())?;
        Self::from_files(&files.config, &files.weights, device)
    }

    /// Load from a local `config.json` and `model.safetensors`.
    ///
    /// # Errors
    /// [`MlError::File`] if a file is unreadable, [`MlError::Config`] if the
    /// configuration is unsupported, [`MlError::Tensor`] if a weight is missing.
    pub fn from_files(config: &Path, weights: &Path, device: &Device) -> Result<Self, MlError> {
        let text = std::fs::read_to_string(config)
            .map_err(|source| MlError::File { path: config.display().to_string(), source })?;
        let config: Config =
            serde_json::from_str(&text).map_err(|e| MlError::Config(format!("{}: {e}", config.display())))?;
        config.validate()?;
        // Read the tensors into memory rather than memory-mapping them: the
        // mmap constructor is `unsafe` (the file must not change underneath),
        // and this crate forbids unsafe code. The cost is holding one copy of
        // the weights in RAM while loading.
        let tensors = candle_core::safetensors::load(weights, device)?;
        let vb = VarBuilder::from_tensors(tensors, DType::F32, device);
        Ok(Self { model: model::Model::load(config, vb)? })
    }

    /// Embedding width.
    pub fn hidden_size(&self) -> usize {
        self.model.config().hidden_size
    }

    /// Longest sequence this checkpoint accepts.
    pub fn max_residues(&self) -> usize {
        self.model.config().max_residues()
    }

    /// Device the weights are on.
    pub fn device(&self) -> &Device {
        self.model.device()
    }

    /// Tokenise a batch, padding to the longest member.
    fn tokenize(&self, sequences: &[&str]) -> Result<(Tensor, Tensor, Vec<usize>), MlError> {
        let limit = self.max_residues();
        let mut encoded = Vec::with_capacity(sequences.len());
        for (index, seq) in sequences.iter().enumerate() {
            let residues = seq.chars().filter(|c| !c.is_whitespace()).count();
            if residues == 0 {
                return Err(MlError::EmptySequence { index });
            }
            if residues > limit {
                return Err(MlError::SequenceTooLong { index, length: residues, limit });
            }
            encoded.push((alphabet::encode(seq), residues));
        }
        let width = encoded.iter().map(|(t, _)| t.len()).max().unwrap_or(0);
        let mut tokens = Vec::with_capacity(sequences.len() * width);
        let mut mask = Vec::with_capacity(sequences.len() * width);
        let mut lengths = Vec::with_capacity(sequences.len());
        for (ids, residues) in encoded {
            let pad = width - ids.len();
            mask.extend(std::iter::repeat_n(1u32, ids.len()));
            mask.extend(std::iter::repeat_n(0u32, pad));
            tokens.extend(ids);
            tokens.extend(std::iter::repeat_n(alphabet::PAD, pad));
            lengths.push(residues);
        }
        let device = self.device();
        let shape = (sequences.len(), width);
        Ok((Tensor::from_vec(tokens, shape, device)?, Tensor::from_vec(mask, shape, device)?, lengths))
    }

    /// Embed a batch of sequences.
    ///
    /// Every sequence in one call is padded to the longest, so batching
    /// similar lengths together is much cheaper than mixing them.
    ///
    /// # Errors
    /// [`MlError::EmptySequence`], [`MlError::SequenceTooLong`], or a tensor error.
    pub fn embed(&self, sequences: &[&str]) -> Result<Vec<Embedding>, MlError> {
        let (tokens, mask, lengths) = self.tokenize(sequences)?;
        let hidden = self.model.encode(&tokens, &mask)?;
        let pooled = self.model.pool(&hidden, &lengths)?.to_vec2::<f32>()?;
        let mut out = Vec::with_capacity(sequences.len());
        for (i, (&length, pooled)) in lengths.iter().zip(pooled).enumerate() {
            let residues = hidden.i(i)?.narrow(0, 1, length)?.flatten_all()?.to_vec1::<f32>()?;
            out.push(Embedding { pooled, per_residue: residues, length });
        }
        Ok(out)
    }

    /// Embed many sequences, batching them to bound the padded work.
    ///
    /// Every member of a batch is padded to the longest, so a batch chosen by
    /// count alone can waste almost all of its work: a 60-residue sequence
    /// sharing a batch with a 1000-residue one is computed as though it were
    /// 1000 long. So sequences are sorted by length and each batch is filled
    /// until `token_budget` padded tokens, which keeps the padding small and
    /// bounds the activation memory a batch needs — the quantity that has to
    /// be bounded, since it is what runs a device out of memory.
    ///
    /// Results come back in the order given, not in length order.
    ///
    /// A budget smaller than one sequence still embeds that sequence alone:
    /// refusing would be choosing not to do work the model can do.
    ///
    /// # Errors
    /// [`MlError::EmptySequence`], [`MlError::SequenceTooLong`] for a sequence
    /// past the checkpoint's limit, or a tensor error. Use
    /// [`Self::embed_windowed`] for longer sequences.
    pub fn embed_many(&self, sequences: &[&str], token_budget: usize) -> Result<Vec<Embedding>, MlError> {
        let lengths: Vec<usize> =
            sequences.iter().map(|s| s.chars().filter(|c| !c.is_whitespace()).count()).collect();
        let mut out: Vec<Option<Embedding>> = vec![None; sequences.len()];
        for batch in plan_batches(&lengths, token_budget) {
            let group: Vec<&str> = batch.iter().map(|&i| sequences[i]).collect();
            for (&index, embedding) in batch.iter().zip(self.embed(&group)?) {
                out[index] = Some(embedding);
            }
        }
        out.into_iter().enumerate().map(|(index, e)| e.ok_or(MlError::EmptySequence { index })).collect()
    }

    /// Per-residue embeddings for a sequence of any length.
    ///
    /// Sequences within the model's limit are embedded directly. Longer ones
    /// are covered by overlapping windows, and each residue takes its
    /// embedding from the window in which it sits most centrally, so no
    /// residue is represented by a vector computed at a window edge where
    /// half its context is missing.
    ///
    /// Returns `length × hidden_size` values, row-major.
    ///
    /// # Errors
    /// [`MlError::EmptySequence`], or a tensor error.
    pub fn embed_windowed(&self, sequence: &str) -> Result<Embedding, MlError> {
        let residues: Vec<char> = sequence.chars().filter(|c| !c.is_whitespace()).collect();
        if residues.is_empty() {
            return Err(MlError::EmptySequence { index: 0 });
        }
        let window = self.max_residues();
        if residues.len() <= window {
            return Ok(self.embed(&[sequence])?.remove(0));
        }

        let width = self.hidden_size();
        let length = residues.len();
        // A quarter-window of overlap on each side gives every interior
        // residue at least `window / 4` residues of context.
        let stride = window - window / 2;
        let mut per_residue = vec![0f32; length * width];
        // How central the residue's current embedding was in its window;
        // a later window only overwrites a residue if it covers it better.
        let mut best = vec![usize::MAX; length];

        let mut start = 0;
        loop {
            let end = (start + window).min(length);
            let chunk: String = residues[start..end].iter().collect();
            let embedded = self.embed(&[&chunk])?.remove(0);
            let centre = start + (end - start) / 2;
            for (offset, position) in (start..end).enumerate() {
                let distance = position.abs_diff(centre);
                if distance < best[position] {
                    best[position] = distance;
                    let src = offset * width;
                    per_residue[position * width..(position + 1) * width]
                        .copy_from_slice(&embedded.per_residue[src..src + width]);
                }
            }
            if end == length {
                break;
            }
            start += stride;
        }

        let mut pooled = vec![0f32; width];
        for row in per_residue.chunks_exact(width) {
            for (acc, v) in pooled.iter_mut().zip(row) {
                *acc += v;
            }
        }
        for v in &mut pooled {
            *v /= length as f32;
        }
        Ok(Embedding { pooled, per_residue, length })
    }

    /// Masked log-probabilities over the alphabet at the given positions.
    ///
    /// For each 0-based residue `position`, the residue is replaced by
    /// `<mask>` and the model is run, giving the distribution it expects
    /// there from the rest of the sequence. The difference between two
    /// entries is the log-odds used by the stability head (Meier et al., 2021).
    ///
    /// Returns one `vocab_size`-long row per requested position, in order.
    ///
    /// # Errors
    /// [`MlError::SequenceTooLong`] if the sequence exceeds the model limit,
    /// or a tensor error.
    pub fn masked_log_probs(&self, sequence: &str, positions: &[usize]) -> Result<Vec<Vec<f32>>, MlError> {
        let residues: Vec<char> = sequence.chars().filter(|c| !c.is_whitespace()).collect();
        if residues.is_empty() {
            return Err(MlError::EmptySequence { index: 0 });
        }
        if residues.len() > self.max_residues() {
            return Err(MlError::SequenceTooLong {
                index: 0,
                length: residues.len(),
                limit: self.max_residues(),
            });
        }
        if positions.is_empty() {
            return Ok(Vec::new());
        }

        let base: Vec<u32> = alphabet::encode(sequence);
        let width = base.len();
        let mut tokens = Vec::with_capacity(positions.len() * width);
        for &position in positions {
            let mut row = base.clone();
            // Token 0 is <cls>, so residue `position` is at index `position + 1`.
            let index =
                position.checked_add(1).filter(|i| *i < width - 1).ok_or(MlError::SequenceTooLong {
                    index: position,
                    length: residues.len(),
                    limit: residues.len(),
                })?;
            row[index] = alphabet::MASK;
            tokens.extend(row);
        }
        let shape = (positions.len(), width);
        let device = self.device();
        let tokens = Tensor::from_vec(tokens, shape, device)?;
        let mask = Tensor::ones(shape, candle_core::DType::U32, device)?;
        let logits = self.model.logits(&tokens, &mask)?;
        let log_probs = candle_nn::ops::log_softmax(&logits, candle_core::D::Minus1)?;

        let mut out = Vec::with_capacity(positions.len());
        for (row, &position) in positions.iter().enumerate() {
            out.push(log_probs.i((row, position + 1))?.to_vec1::<f32>()?);
        }
        Ok(out)
    }

    /// Final hidden states including the `<cls>` and `<eos>` positions,
    /// `(batch, seq, hidden)`. Mainly useful for comparing against reference
    /// implementations.
    ///
    /// # Errors
    /// As [`Esm2::embed`].
    pub fn hidden_states(&self, sequences: &[&str]) -> Result<Tensor, MlError> {
        let (tokens, mask, _) = self.tokenize(sequences)?;
        Ok(self.model.encode(&tokens, &mask)?)
    }

    /// Token logits, `(batch, seq, vocab)`.
    ///
    /// # Errors
    /// As [`Esm2::embed`].
    pub fn logits(&self, sequences: &[&str]) -> Result<Tensor, MlError> {
        let (tokens, mask, _) = self.tokenize(sequences)?;
        Ok(self.model.logits(&tokens, &mask)?)
    }
}

use candle_core::IndexOp;

#[cfg(test)]
mod tests {
    use super::*;

    /// The budget exists to bound memory, and a batch's cost is members times
    /// the longest member. If the plan ignores that, a long sequence joining a
    /// batch of short ones silently multiplies the activation size.
    #[test]
    fn no_batch_exceeds_the_token_budget_unless_one_sequence_already_does() {
        let lengths = [60, 300, 75, 1000, 80, 62, 290];
        for budget in [0, 1, 100, 512, 8192, 100_000] {
            let plan = plan_batches(&lengths, budget);
            for batch in &plan {
                let widest = batch.iter().map(|&i| lengths[i]).max().unwrap_or(0);
                let cost = batch.len() * widest;
                assert!(
                    cost <= budget || batch.len() == 1,
                    "budget {budget}: a batch of {} costs {cost}",
                    batch.len()
                );
            }
        }
    }

    /// Every sequence must be embedded exactly once. A planner that drops one
    /// produces a table quietly missing rows, and one that duplicates it
    /// produces rows nobody asked for.
    #[test]
    fn every_sequence_appears_exactly_once() {
        let lengths = [60, 300, 75, 1000, 80, 62, 290];
        for budget in [0, 1, 100, 512, 8192, 100_000] {
            let mut seen: Vec<usize> = plan_batches(&lengths, budget).concat();
            seen.sort_unstable();
            assert_eq!(seen, (0..lengths.len()).collect::<Vec<_>>(), "budget {budget}");
        }
    }

    /// Padding changes floating-point summation order, so a plan that varied
    /// between runs would make results irreproducible for no reason.
    #[test]
    fn the_plan_is_deterministic_and_groups_similar_lengths() {
        let lengths = [300, 60, 301, 61, 302, 62];
        let first = plan_batches(&lengths, 700);
        assert_eq!(first, plan_batches(&lengths, 700), "the same input must give the same plan");

        // The three short sequences should not be sharing a batch with the
        // three long ones, which is the whole point of sorting by length.
        let shorts: Vec<&Vec<usize>> = first.iter().filter(|b| b.iter().any(|&i| lengths[i] < 100)).collect();
        assert!(
            shorts.iter().all(|b| b.iter().all(|&i| lengths[i] < 100)),
            "a batch mixes short and long: {first:?}"
        );
    }

    /// A budget of zero means one sequence per batch, which pads nothing and
    /// therefore reproduces the unbatched numbers exactly. That is the default
    /// the CLI uses, so it has to hold.
    #[test]
    fn a_zero_budget_puts_every_sequence_in_its_own_batch() {
        let plan = plan_batches(&[60, 300, 75], 0);
        assert_eq!(plan.len(), 3, "{plan:?}");
        assert!(plan.iter().all(|b| b.len() == 1), "{plan:?}");
    }

    #[test]
    fn an_empty_input_plans_nothing() {
        assert_eq!(plan_batches(&[], 8192), [] as [Vec<usize>; 0]);
    }
}
