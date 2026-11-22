//! Builds the embedding caches the training runs read.
//!
//! This is the expensive step: every unique sequence in a dataset is passed
//! through ESM-2 once. Results are cached, so re-running only covers
//! sequences that are not there yet.

use crate::cache::MatrixCache;
use crate::data::{self, Split};
use crate::Task;
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use opendna_ml::esm2::{Esm2, ModelSize};
use std::collections::BTreeSet;
use std::path::Path;

/// Sequences are embedded in batches of roughly this many residues, so long
/// and short sequences both make reasonable use of the GPU.
const RESIDUES_PER_BATCH: usize = 8_192;

/// Cache file name for a task's embeddings.
pub fn embedding_cache_name(task: Task, model: ModelSize) -> String {
    format!("{}_{}_embeddings", task.name(), model_slug(model))
}

/// Cache file name for the stability task's masked log-probabilities.
pub fn log_prob_cache_name(model: ModelSize) -> String {
    format!("stability_{}_logprobs", model_slug(model))
}

fn model_slug(model: ModelSize) -> String {
    model.repo().rsplit('/').next().unwrap_or("model").to_string()
}

/// Every unique sequence a task needs, across all three splits.
fn sequences_for(task: Task, root: &Path) -> Result<BTreeSet<String>> {
    let mut out = BTreeSet::new();
    for split in Split::all() {
        match task {
            Task::Stability => {
                out.extend(data::stability(root, split)?.into_iter().map(|r| r.wt_sequence));
            }
            Task::Solubility => {
                out.extend(data::solubility(root, split)?.into_iter().map(|r| r.sequence));
            }
            Task::Disorder => {
                out.extend(data::disorder(root, split)?.into_iter().map(|r| r.sequence));
            }
            Task::Burial => {
                out.extend(data::burial(root, split)?.into_iter().map(|r| r.sequence));
            }
        }
    }
    Ok(out)
}

fn progress(len: u64, label: &str) -> ProgressBar {
    let bar = ProgressBar::new(len);
    bar.set_style(
        ProgressStyle::with_template("{msg} [{bar:40}] {pos}/{len} ({eta} left)")
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars("=> "),
    );
    bar.set_message(label.to_string());
    bar
}

/// Embed every sequence a task needs and write the caches.
///
/// # Errors
/// If the datasets cannot be read, the model cannot be loaded, or the cache
/// cannot be written.
pub fn run(task: Task, model_size: ModelSize, datasets_root: &Path, cache_dir: &Path) -> Result<()> {
    let device = opendna_ml::best_device().context("selecting a compute device")?;
    println!("Embedding with {} on {}", model_size.repo(), opendna_ml::device_description(&device));
    let model = Esm2::from_pretrained(model_size, &device)
        .with_context(|| format!("loading {}", model_size.repo()))?;

    let sequences = sequences_for(task, datasets_root)?;
    println!("{} unique sequences for {}", sequences.len(), task.name());

    let mut cache = MatrixCache::open(
        cache_dir,
        &embedding_cache_name(task, model_size),
        model_size.repo(),
        model.hidden_size(),
    )?;

    let pending: Vec<&String> = sequences.iter().filter(|s| !cache.contains(s)).collect();
    if pending.is_empty() {
        println!("All {} sequences are already cached.", sequences.len());
    } else {
        println!("{} to embed, {} already cached", pending.len(), sequences.len() - pending.len());
        embed_all(&model, &pending, &mut cache, task)?;
    }

    if task == Task::Stability {
        build_log_probs(&model, model_size, &sequences, datasets_root, cache_dir)?;
    }
    Ok(())
}

/// Embed `pending` into `cache`, batching short sequences together.
fn embed_all(model: &Esm2, pending: &[&String], cache: &mut MatrixCache, task: Task) -> Result<()> {
    let bar = progress(pending.len() as u64, "embedding");
    let limit = model.max_residues();
    let mut batch: Vec<&str> = Vec::new();
    let mut residues = 0;

    // Group by length so a batch is not padded to a much longer member.
    let mut ordered: Vec<&&String> = pending.iter().collect();
    ordered.sort_by_key(|s| s.len());

    for sequence in ordered {
        // Sequences past the model limit need windowing, which handles one
        // sequence at a time.
        if sequence.len() > limit {
            flush(model, &mut batch, &mut residues, cache, task, &bar)?;
            let embedded = model.embed_windowed(sequence)?;
            store(cache, sequence, &embedded, task)?;
            bar.inc(1);
            continue;
        }
        if !batch.is_empty() && residues + sequence.len() > RESIDUES_PER_BATCH {
            flush(model, &mut batch, &mut residues, cache, task, &bar)?;
        }
        residues += sequence.len();
        batch.push(sequence);
    }
    flush(model, &mut batch, &mut residues, cache, task, &bar)?;
    cache.save()?;
    bar.finish_with_message("embedding complete");
    Ok(())
}

fn flush(
    model: &Esm2,
    batch: &mut Vec<&str>,
    residues: &mut usize,
    cache: &mut MatrixCache,
    task: Task,
    bar: &ProgressBar,
) -> Result<()> {
    if batch.is_empty() {
        return Ok(());
    }
    let embedded = model.embed(batch)?;
    for (sequence, embedding) in batch.iter().zip(&embedded) {
        store(cache, sequence, embedding, task)?;
    }
    bar.inc(batch.len() as u64);
    // Save periodically so an interrupted run keeps its progress.
    if cache.len().is_multiple_of(512) {
        cache.save()?;
    }
    batch.clear();
    *residues = 0;
    Ok(())
}

/// Solubility only needs the pooled vector; the other tasks need every residue.
fn store(
    cache: &mut MatrixCache,
    sequence: &str,
    embedding: &opendna_ml::esm2::Embedding,
    task: Task,
) -> Result<()> {
    if task == Task::Solubility {
        cache.append(sequence, 1, &embedding.pooled)
    } else {
        cache.append(sequence, embedding.length, &embedding.per_residue)
    }
}

/// Masked log-probabilities at every position that the stability dataset
/// mutates. These give the head the language model's own opinion about a
/// substitution, which is a strong predictor on its own.
fn build_log_probs(
    model: &Esm2,
    model_size: ModelSize,
    sequences: &BTreeSet<String>,
    datasets_root: &Path,
    cache_dir: &Path,
) -> Result<()> {
    let vocab = opendna_ml::esm2::alphabet::TOKENS.len();
    let mut cache = MatrixCache::open(cache_dir, &log_prob_cache_name(model_size), model_size.repo(), vocab)?;
    let pending: Vec<&String> = sequences.iter().filter(|s| !cache.contains(s)).collect();
    if pending.is_empty() {
        println!("Masked log-probabilities already cached.");
        return Ok(());
    }

    // Only positions that are actually mutated are needed, but the rows are
    // small and a complete matrix keeps indexing simple.
    let _ = datasets_root;
    println!("Computing masked log-probabilities for {} sequences", pending.len());
    let bar = progress(pending.len() as u64, "masking");
    for sequence in pending {
        let positions: Vec<usize> = (0..sequence.len()).collect();
        let rows = model.masked_log_probs(sequence, &positions)?;
        let flat: Vec<f32> = rows.into_iter().flatten().collect();
        cache.append(sequence, sequence.len(), &flat)?;
        bar.inc(1);
        if cache.len().is_multiple_of(64) {
            cache.save()?;
        }
    }
    cache.save()?;
    bar.finish_with_message("masking complete");
    Ok(())
}
