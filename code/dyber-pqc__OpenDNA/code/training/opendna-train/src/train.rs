//! Training loops for the three prediction heads.
//!
//! Every head is trained on frozen ESM-2 embeddings read from the cache, so
//! a run takes minutes rather than hours and is reproducible from a seed.
//! The validation split chooses when to stop; the test split is read once,
//! at the end, for the numbers that go in the model card.

// Each training function is a complete, linear recipe: load, train,
// evaluate, save. Splitting them would scatter a procedure that reads better
// in one piece.
#![allow(clippy::too_many_lines)]

use crate::cache::MatrixCache;
use crate::data::{self, Split};
use crate::metrics;
use crate::Task;
use anyhow::{bail, Context, Result};
use candle_core::{DType, Tensor};
use candle_nn::{AdamW, Optimizer, ParamsAdamW, VarBuilder, VarMap};
use opendna_core::AminoAcid;
use opendna_ml::esm2::ModelSize;
use opendna_ml::heads::{BurialHead, DisorderHead, HeadMetadata, SolubilityHead, StabilityHead};
use rand::prelude::*;
use rand::rngs::StdRng;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// Knobs shared by the three training runs.
#[derive(Debug, Clone)]
pub struct TrainConfig {
    /// Maximum passes over the training split.
    pub epochs: usize,
    /// Rows per optimisation step.
    pub batch_size: usize,
    /// `AdamW` learning rate.
    pub learning_rate: f64,
    /// Stop after this many epochs without a better validation score.
    pub patience: usize,
    /// Seed for shuffling and weight initialisation.
    pub seed: u64,
    /// Give the stability head predicted solvent accessibility at the mutated
    /// position. Ignored by the other tasks.
    pub with_burial: bool,
}

impl Default for TrainConfig {
    fn default() -> Self {
        Self {
            epochs: 100,
            batch_size: 256,
            learning_rate: 1e-3,
            patience: 10,
            seed: 20_260_101,
            with_burial: false,
        }
    }
}

/// Metrics for one split.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplitMetrics {
    /// Rows evaluated.
    pub n: usize,
    /// Metric name to value. Absent metrics are omitted rather than reported
    /// as zero.
    #[serde(flatten)]
    pub values: std::collections::BTreeMap<String, f64>,
}

/// Everything a training run produces, written next to the weights.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    /// Task name.
    pub task: String,
    /// Embedding model used.
    pub embedding_model: String,
    /// Embedding width.
    pub embedding_dim: usize,
    /// Epoch whose weights were kept.
    pub best_epoch: usize,
    /// Training configuration, flattened for readability.
    pub epochs_run: usize,
    /// Learning rate used.
    pub learning_rate: f64,
    /// Batch size used.
    pub batch_size: usize,
    /// Seed used.
    pub seed: u64,
    /// Whether weights were drawn from the seed. Always true now that
    /// initialisation does not depend on the device's own generator.
    pub seed_applied: bool,
    /// Validation metrics at the kept epoch.
    pub validation: SplitMetrics,
    /// Test metrics at the kept epoch.
    pub test: SplitMetrics,
    /// Metrics for the trivial baseline a useful model must beat.
    pub baseline: SplitMetrics,
}

fn metric_map(pairs: &[(&str, Option<f64>)]) -> std::collections::BTreeMap<String, f64> {
    pairs.iter().filter_map(|(k, v)| v.map(|v| ((*k).to_string(), (v * 1e6).round() / 1e6))).collect()
}

/// Re-initialise every variable from a seeded generator.
///
/// candle's CPU backend has no seedable device generator, so a run would
/// otherwise start from different weights each time — on this task that moved
/// the test score by about 0.03 Spearman, which is enough to make two runs
/// disagree about which model is better. Weights are therefore drawn here,
/// from this crate's own generator, which makes a run reproducible on any
/// device.
///
/// The distribution matches `PyTorch`'s default for linear and convolutional
/// layers: uniform over ±1/sqrt(`fan_in`), where `fan_in` is the product of every
/// dimension after the first.
///
/// # Errors
/// Propagates tensor errors.
fn seed_weights(varmap: &VarMap, seed: u64) -> Result<()> {
    let mut rng = StdRng::seed_from_u64(seed);
    // Sort by name so the draw order does not depend on hash iteration order.
    let data = varmap.data().lock().expect("varmap mutex");
    let mut names: Vec<&String> = data.keys().collect();
    names.sort();

    for name in names {
        let var = &data[name];
        let shape = var.shape().clone();
        let dims = shape.dims();
        let fan_in: usize = if dims.len() > 1 { dims[1..].iter().product() } else { dims[0] };
        let bound = 1.0 / (fan_in as f64).sqrt();
        let count = shape.elem_count();
        let values: Vec<f32> = (0..count).map(|_| rng.random_range(-bound..bound) as f32).collect();
        var.set(&Tensor::from_vec(values, shape, var.device())?)?;
    }
    Ok(())
}

fn amino_index(residue: char) -> Option<usize> {
    let aa = AminoAcid::from_one_letter(residue)?;
    AminoAcid::STANDARD.iter().position(|s| *s == aa)
}

/// Index of a residue in the ESM-2 alphabet, for reading log-probability rows.
fn token_index(residue: char) -> usize {
    opendna_ml::esm2::alphabet::token_id(residue) as usize
}

/// Save weights, metadata and metrics.
/// Published version of each head.
///
/// Semantic and per head, because they are separate models: the crate's own
/// version is the same string for all four and says `dev`, which is no use to
/// anyone pinning what they ran. Bump the head's version when a retrain
/// changes what it predicts; the build id below changes on its own whenever
/// any input does.
// The arms are identical today and must stay separate: these are four
// independent models, and merging them would mean a retrain of one silently
// bumped the others.
#[allow(clippy::match_same_arms)]
fn head_version(task: Task) -> &'static str {
    match task {
        Task::Stability => "1.0.0",
        Task::Solubility => "1.0.0",
        // 1.1.0: the same weights, now carrying a calibration map and a
        // measured decision threshold. The probabilities it reports are
        // different numbers from 1.0.0's, so it is not 1.0.0.
        Task::Disorder => "1.1.0",
        Task::Burial => "0.1.0",
    }
}

/// A digest over everything that determines these weights.
///
/// Two runs agreeing on this produced the same model from the same inputs.
/// It covers the dataset the run read, the hyperparameters, the embedding
/// model and the training code's version, so a retrain that forgets to bump
/// the semantic version is still distinguishable.
fn build_id(task: Task, model_size: ModelSize, report: &Report, datasets_root: &Path) -> String {
    build_id_with(task, model_size, report, datasets_root, None)
}

/// As [`build_id`], folding in a calibration when one is recorded.
///
/// A calibrated head produces different numbers from the same weights, so it
/// is a different model and must not share an identity with the uncalibrated
/// one. Without this, applying a calibration would leave `models --check`
/// passing against a lockfile taken before it.
fn build_id_with(
    task: Task,
    model_size: ModelSize,
    report: &Report,
    datasets_root: &Path,
    calibration: Option<&opendna_ml::heads::Calibration>,
) -> String {
    use sha2::{Digest, Sha256};
    let mut h = Sha256::new();
    h.update(task.name().as_bytes());
    h.update(head_version(task).as_bytes());
    h.update(model_size.repo().as_bytes());
    h.update(env!("CARGO_PKG_VERSION").as_bytes());
    h.update(
        format!(
            "lr={} batch={} seed={} epoch={}",
            report.learning_rate, report.batch_size, report.seed, report.best_epoch
        )
        .as_bytes(),
    );
    // The dataset's own manifest already hashes every Parquet file, so
    // folding it in pins the data without rehashing it here.
    let manifest = datasets_root.join(task.name()).join("manifest.json");
    match std::fs::read(&manifest) {
        Ok(bytes) => h.update(&bytes),
        // A build id that cannot see the data says so rather than silently
        // describing less than it claims to.
        Err(_) => h.update(b"dataset-manifest-unavailable"),
    }
    if let Some(c) = calibration {
        h.update(
            format!(
                "calibration={} slope={:.9} intercept={:.9} threshold={:.9}",
                c.method, c.slope, c.intercept, c.threshold
            )
            .as_bytes(),
        );
    }
    format!("{:x}", h.finalize())[..16].to_string()
}

/// Rewrite an existing head's metadata with its identity fields.
///
/// Reads the head's current metadata for the things only training knows
/// (the embedding width, and whether the stability head takes burial), and
/// recomputes the version, build id and weights hash exactly as [`save`]
/// would. The weights are not touched.
///
/// # Errors
/// If the head's files cannot be read or the metadata cannot be written.
pub fn stamp(
    task: Task,
    out_dir: &Path,
    model_size: ModelSize,
    report: &Report,
    datasets_root: &Path,
) -> Result<()> {
    let files = opendna_ml::heads::HeadFiles::in_dir(out_dir, task.name());
    let existing: HeadMetadata = serde_json::from_str(
        &std::fs::read_to_string(&files.metadata)
            .with_context(|| format!("reading {}", files.metadata.display()))?,
    )
    .with_context(|| format!("parsing {}", files.metadata.display()))?;

    let metadata = HeadMetadata {
        task: task.name().to_string(),
        embedding_model: model_size.repo().to_string(),
        embedding_dim: existing.embedding_dim,
        version: head_version(task).to_string(),
        build: Some(build_id_with(task, model_size, report, datasets_root, existing.calibration.as_ref())),
        weights_sha256: Some(opendna_ml::heads::sha256_file(&files.weights)?),
        with_burial: existing.with_burial,
        // Stamping restates an identity without retraining, so the weights
        // are the ones the calibration was fitted against and it still holds.
        calibration: existing.calibration,
        abstention: existing.abstention,
        intervals: existing.intervals,
    };
    std::fs::write(&files.metadata, serde_json::to_string_pretty(&metadata)?)
        .with_context(|| format!("writing {}", files.metadata.display()))?;
    println!(
        "{}: version {} build {} sha256 {}",
        task.name(),
        metadata.version,
        metadata.build.as_deref().unwrap_or("-"),
        &metadata.weights_sha256.as_deref().unwrap_or("-")[..16]
    );
    Ok(())
}

/// Everything `save` needs that is not the weights themselves.
struct Saved<'a> {
    out_dir: &'a Path,
    task: Task,
    model_size: ModelSize,
    embedding_dim: usize,
    report: &'a Report,
    with_burial: bool,
    datasets_root: &'a Path,
}

fn save(varmap: &VarMap, s: &Saved<'_>) -> Result<()> {
    let (out_dir, task, model_size) = (s.out_dir, s.task, s.model_size);
    let (embedding_dim, report, with_burial) = (s.embedding_dim, s.report, s.with_burial);
    let datasets_root = s.datasets_root;
    std::fs::create_dir_all(out_dir)?;
    let weights = out_dir.join(format!("{}.safetensors", task.name()));
    varmap.save(&weights).with_context(|| format!("writing {}", weights.display()))?;

    let metadata = HeadMetadata {
        task: task.name().to_string(),
        embedding_model: model_size.repo().to_string(),
        embedding_dim,
        version: head_version(task).to_string(),
        build: Some(build_id(task, model_size, report, datasets_root)),
        // Hashed after writing, so the recorded digest is of the bytes that
        // actually landed rather than of what was meant to.
        weights_sha256: Some(opendna_ml::heads::sha256_file(&weights)?),
        // Recorded so loading reconstructs the right input width rather than
        // inferring it from the weight shapes.
        with_burial,
        // Deliberately none. A calibration is fitted against particular
        // weights, and these are new ones -- carrying the old map over would
        // apply a correction measured for a different model. `calibrate`
        // measures it again.
        calibration: None,
        // Void for the same reason: a confidence/error trade-off measured
        // against other weights says nothing about these.
        abstention: None,
        // And an interval fitted to another model's residuals is not this
        // model's interval.
        intervals: None,
    };
    std::fs::write(
        out_dir.join(format!("{}.metadata.json", task.name())),
        serde_json::to_string_pretty(&metadata)?,
    )?;
    std::fs::write(
        out_dir.join(format!("{}.metrics.json", task.name())),
        serde_json::to_string_pretty(report)?,
    )?;
    println!("Wrote {}", weights.display());
    Ok(())
}

/// Mean of a slice, or 0 when empty.
fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Stability
// ---------------------------------------------------------------------------

struct StabilityData {
    features: Vec<f32>,
    targets: Vec<f32>,
    width: usize,
}

impl StabilityData {
    fn len(&self) -> usize {
        self.targets.len()
    }
}

/// Build the feature matrix for one split from the cached embeddings.
fn stability_features(
    rows: &[data::StabilityRow],
    embeddings: &HashMap<String, (Vec<f32>, usize)>,
    log_probs: &HashMap<String, (Vec<f32>, usize)>,
    hidden: usize,
    vocab: usize,
    burial: Option<&HashMap<String, Vec<f32>>>,
) -> Result<StabilityData> {
    let width = if burial.is_some() {
        StabilityHead::input_width_with_burial(hidden)
    } else {
        StabilityHead::input_width(hidden)
    };
    let mut features = Vec::with_capacity(rows.len() * width);
    let mut targets = Vec::with_capacity(rows.len());

    for row in rows {
        let (embedding, _) = embeddings
            .get(&row.wt_sequence)
            .with_context(|| "a wild-type sequence is missing from the embedding cache")?;
        let (log_prob, _) = log_probs
            .get(&row.wt_sequence)
            .with_context(|| "a wild-type sequence is missing from the log-probability cache")?;

        let start = row.position * hidden;
        let site = embedding
            .get(start..start + hidden)
            .with_context(|| format!("position {} is past the cached embedding", row.position))?;

        let base = row.position * vocab;
        let wt_logp = log_prob[base + token_index(row.wild_type)];
        let mut_logp = log_prob[base + token_index(row.mutant)];

        let (Some(wt), Some(mt)) = (amino_index(row.wild_type), amino_index(row.mutant)) else {
            bail!("mutation {}{}{} uses a non-standard residue", row.wild_type, row.position + 1, row.mutant);
        };
        // Accessibility at the mutated position, when the variant uses it.
        // A missing sequence is an error rather than a zero: zero means
        // "fully buried", so substituting it would feed the head a confident
        // wrong answer.
        let exposure =
            match burial {
                Some(map) => {
                    let values = map
                        .get(&row.wt_sequence)
                        .with_context(|| "a wild-type sequence is missing from the burial predictions")?;
                    Some(*values.get(row.position).with_context(|| {
                        format!("position {} is past the burial prediction", row.position)
                    })?)
                }
                None => None,
            };
        features.extend(StabilityHead::features_with(site, wt, mt, mut_logp - wt_logp, exposure));
        targets.push(row.ddg as f32);
    }
    Ok(StabilityData { features, targets, width })
}

/// Train the ΔΔG head.
///
/// # Errors
/// If the datasets or caches cannot be read, or training fails.
pub fn stability(
    datasets_root: &Path,
    cache_dir: &Path,
    out_dir: &Path,
    model_size: ModelSize,
    config: &TrainConfig,
) -> Result<Report> {
    let device = opendna_ml::best_device()?;
    let hidden = model_size.hidden_size();
    let vocab = opendna_ml::esm2::alphabet::TOKENS.len();

    let embeddings_cache = MatrixCache::open(
        cache_dir,
        &crate::embed::embedding_cache_name(Task::Stability, model_size),
        model_size.repo(),
        hidden,
    )?;
    let log_probs_cache = MatrixCache::open(
        cache_dir,
        &crate::embed::log_prob_cache_name(model_size),
        model_size.repo(),
        vocab,
    )?;

    let splits: Vec<Vec<data::StabilityRow>> =
        Split::all().iter().map(|s| data::stability(datasets_root, *s)).collect::<Result<_>>()?;

    // Only a few hundred unique wild-type sequences back hundreds of
    // thousands of mutations, so everything fits comfortably in memory.
    let mut embeddings = HashMap::new();
    let mut log_probs = HashMap::new();
    for split in &splits {
        for row in split {
            if !embeddings.contains_key(&row.wt_sequence) {
                embeddings.insert(row.wt_sequence.clone(), embeddings_cache.get(&row.wt_sequence)?);
                log_probs.insert(row.wt_sequence.clone(), log_probs_cache.get(&row.wt_sequence)?);
            }
        }
    }
    println!("{} unique wild-type sequences in memory", embeddings.len());

    // Predicted accessibility per wild-type sequence, from the trained
    // burial head beside the weights being written. It is loaded rather than
    // retrained so the feature is exactly what inference would supply.
    let burial = if config.with_burial {
        let heads = opendna_ml::heads::Heads::load(out_dir, model_size, &device)?;
        let head = heads.burial.as_ref().with_context(|| {
            format!("--with-burial needs a trained burial head in {}; train one first", out_dir.display())
        })?;
        let mut map = HashMap::new();
        for (sequence, (values, rows)) in &embeddings {
            map.insert(sequence.clone(), head.predict(values, *rows, &device)?);
        }
        println!("burial predicted for {} wild-type sequences", map.len());
        Some(map)
    } else {
        None
    };

    let data: Vec<StabilityData> = splits
        .iter()
        .map(|rows| stability_features(rows, &embeddings, &log_probs, hidden, vocab, burial.as_ref()))
        .collect::<Result<_>>()?;
    let (train, validation, test) = (&data[0], &data[1], &data[2]);
    println!("{} train / {} validation / {} test mutations", train.len(), validation.len(), test.len());

    let mut varmap = VarMap::new();
    let vb = VarBuilder::from_varmap(&varmap, DType::F32, &device);
    let head = StabilityHead::new_with(hidden, config.with_burial, vb)?;
    seed_weights(&varmap, config.seed)?;
    let mut optimizer =
        AdamW::new(varmap.all_vars(), ParamsAdamW { lr: config.learning_rate, ..Default::default() })?;

    let predict = |head: &StabilityHead, d: &StabilityData| -> Result<Vec<f64>> {
        let mut out = Vec::with_capacity(d.len());
        for chunk_start in (0..d.len()).step_by(4096) {
            let n = (d.len() - chunk_start).min(4096);
            let slice = &d.features[chunk_start * d.width..(chunk_start + n) * d.width];
            let x = Tensor::from_slice(slice, (n, d.width), &device)?;
            out.extend(head.forward(&x)?.to_vec1::<f32>()?.into_iter().map(f64::from));
        }
        Ok(out)
    };

    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut order: Vec<usize> = (0..train.len()).collect();
    let mut best = f64::NEG_INFINITY;
    let mut best_epoch = 0;
    let mut epochs_run = 0;
    let mut best_weights: Option<Vec<u8>> = None;

    for epoch in 1..=config.epochs {
        epochs_run = epoch;
        order.shuffle(&mut rng);
        let mut losses = Vec::new();
        for batch in order.chunks(config.batch_size) {
            let mut xs = Vec::with_capacity(batch.len() * train.width);
            let mut ys = Vec::with_capacity(batch.len());
            for &i in batch {
                xs.extend_from_slice(&train.features[i * train.width..(i + 1) * train.width]);
                ys.push(train.targets[i]);
            }
            let x = Tensor::from_vec(xs, (batch.len(), train.width), &device)?;
            let y = Tensor::from_vec(ys, batch.len(), &device)?;
            let loss = candle_nn::loss::mse(&head.forward(&x)?, &y)?;
            optimizer.backward_step(&loss)?;
            losses.push(f64::from(loss.to_scalar::<f32>()?));
        }

        let predicted = predict(&head, validation)?;
        let actual: Vec<f64> = validation.targets.iter().map(|v| f64::from(*v)).collect();
        let score = metrics::spearman(&predicted, &actual).unwrap_or(f64::NEG_INFINITY);
        println!(
            "epoch {epoch:3}  loss {:.4}  validation spearman {score:.4}{}",
            mean(&losses),
            if score > best { "  *" } else { "" }
        );
        if score > best {
            best = score;
            best_epoch = epoch;
            // Keep the best weights in memory; the file is written once, at the end.
            let tmp =
                std::env::temp_dir().join(format!("opendna-stability-{}.safetensors", std::process::id()));
            varmap.save(&tmp)?;
            best_weights = Some(std::fs::read(&tmp)?);
            let _ = std::fs::remove_file(&tmp);
        } else if epoch - best_epoch >= config.patience {
            println!("No improvement for {} epochs; stopping.", config.patience);
            break;
        }
    }

    if let Some(bytes) = best_weights {
        let tmp =
            std::env::temp_dir().join(format!("opendna-stability-best-{}.safetensors", std::process::id()));
        std::fs::write(&tmp, bytes)?;
        varmap.load(&tmp)?;
        let _ = std::fs::remove_file(&tmp);
    }

    let evaluate = |d: &StabilityData| -> Result<SplitMetrics> {
        let predicted = predict(&head, d)?;
        let actual: Vec<f64> = d.targets.iter().map(|v| f64::from(*v)).collect();
        Ok(SplitMetrics {
            n: d.len(),
            values: metric_map(&[
                ("spearman", metrics::spearman(&predicted, &actual)),
                ("pearson", metrics::pearson(&predicted, &actual)),
                ("rmse", Some(metrics::rmse(&predicted, &actual))),
                ("mae", Some(metrics::mae(&predicted, &actual))),
            ]),
        })
    };

    // The language model's log-odds alone, with no trained head: a useful
    // head has to beat this.
    // Not `width - 1`: with burial appended the last column is burial, and
    // taking it would report burial-against-ΔΔG as the zero-shot baseline.
    let log_odds = StabilityHead::log_odds_index(hidden);
    let zero_shot: Vec<f64> =
        (0..test.len()).map(|i| f64::from(test.features[i * test.width + log_odds])).collect();
    let test_actual: Vec<f64> = test.targets.iter().map(|v| f64::from(*v)).collect();
    let baseline = SplitMetrics {
        n: test.len(),
        values: metric_map(&[
            ("spearman", metrics::spearman(&zero_shot, &test_actual)),
            ("pearson", metrics::pearson(&zero_shot, &test_actual)),
        ]),
    };

    let report = Report {
        task: Task::Stability.name().to_string(),
        embedding_model: model_size.repo().to_string(),
        embedding_dim: hidden,
        best_epoch,
        epochs_run,
        learning_rate: config.learning_rate,
        batch_size: config.batch_size,
        seed: config.seed,
        seed_applied: true,
        validation: evaluate(validation)?,
        test: evaluate(test)?,
        baseline,
    };
    save(
        &varmap,
        &Saved {
            out_dir,
            task: Task::Stability,
            model_size,
            embedding_dim: hidden,
            report: &report,
            with_burial: config.with_burial,
            datasets_root,
        },
    )?;
    Ok(report)
}

// ---------------------------------------------------------------------------
// Solubility
// ---------------------------------------------------------------------------

/// Train the solubility head.
///
/// # Errors
/// If the datasets or caches cannot be read, or training fails.
pub fn solubility(
    datasets_root: &Path,
    cache_dir: &Path,
    out_dir: &Path,
    model_size: ModelSize,
    config: &TrainConfig,
) -> Result<Report> {
    let device = opendna_ml::best_device()?;
    let hidden = model_size.hidden_size();
    let cache = MatrixCache::open(
        cache_dir,
        &crate::embed::embedding_cache_name(Task::Solubility, model_size),
        model_size.repo(),
        hidden,
    )?;

    let mut features = Vec::new();
    let mut targets = Vec::new();
    for split in Split::all() {
        let rows = data::solubility(datasets_root, split)?;
        let mut x = Vec::with_capacity(rows.len() * hidden);
        let mut y = Vec::with_capacity(rows.len());
        for row in &rows {
            let (values, _) = cache.get(&row.sequence)?;
            x.extend(values);
            y.push(row.solubility as f32);
        }
        features.push(x);
        targets.push(y);
    }
    println!(
        "{} train / {} validation / {} test proteins",
        targets[0].len(),
        targets[1].len(),
        targets[2].len()
    );

    let mut varmap = VarMap::new();
    let vb = VarBuilder::from_varmap(&varmap, DType::F32, &device);
    let head = SolubilityHead::new(hidden, vb)?;
    seed_weights(&varmap, config.seed)?;
    let mut optimizer =
        AdamW::new(varmap.all_vars(), ParamsAdamW { lr: config.learning_rate, ..Default::default() })?;

    let predict = |head: &SolubilityHead, x: &[f32], n: usize| -> Result<Vec<f64>> {
        let t = Tensor::from_slice(x, (n, hidden), &device)?;
        Ok(head.forward(&t)?.to_vec1::<f32>()?.into_iter().map(f64::from).collect())
    };

    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut order: Vec<usize> = (0..targets[0].len()).collect();
    let mut best = f64::NEG_INFINITY;
    let mut best_epoch = 0;
    let mut epochs_run = 0;
    let mut best_weights: Option<Vec<u8>> = None;

    for epoch in 1..=config.epochs {
        epochs_run = epoch;
        order.shuffle(&mut rng);
        let mut losses = Vec::new();
        for batch in order.chunks(config.batch_size) {
            let mut xs = Vec::with_capacity(batch.len() * hidden);
            let mut ys = Vec::with_capacity(batch.len());
            for &i in batch {
                xs.extend_from_slice(&features[0][i * hidden..(i + 1) * hidden]);
                ys.push(targets[0][i]);
            }
            let x = Tensor::from_vec(xs, (batch.len(), hidden), &device)?;
            let y = Tensor::from_vec(ys, batch.len(), &device)?;
            let loss = candle_nn::loss::mse(&head.forward(&x)?, &y)?;
            optimizer.backward_step(&loss)?;
            losses.push(f64::from(loss.to_scalar::<f32>()?));
        }

        let predicted = predict(&head, &features[1], targets[1].len())?;
        let actual: Vec<f64> = targets[1].iter().map(|v| f64::from(*v)).collect();
        let score = metrics::spearman(&predicted, &actual).unwrap_or(f64::NEG_INFINITY);
        println!(
            "epoch {epoch:3}  loss {:.2}  validation spearman {score:.4}{}",
            mean(&losses),
            if score > best { "  *" } else { "" }
        );
        if score > best {
            best = score;
            best_epoch = epoch;
            let tmp = std::env::temp_dir().join(format!("opendna-sol-{}.safetensors", std::process::id()));
            varmap.save(&tmp)?;
            best_weights = Some(std::fs::read(&tmp)?);
            let _ = std::fs::remove_file(&tmp);
        } else if epoch - best_epoch >= config.patience {
            println!("No improvement for {} epochs; stopping.", config.patience);
            break;
        }
    }

    if let Some(bytes) = best_weights {
        let tmp = std::env::temp_dir().join(format!("opendna-sol-best-{}.safetensors", std::process::id()));
        std::fs::write(&tmp, bytes)?;
        varmap.load(&tmp)?;
        let _ = std::fs::remove_file(&tmp);
    }

    let evaluate = |index: usize| -> Result<SplitMetrics> {
        let predicted = predict(&head, &features[index], targets[index].len())?;
        let actual: Vec<f64> = targets[index].iter().map(|v| f64::from(*v)).collect();
        Ok(SplitMetrics {
            n: actual.len(),
            values: metric_map(&[
                ("spearman", metrics::spearman(&predicted, &actual)),
                ("pearson", metrics::pearson(&predicted, &actual)),
                ("rmse", Some(metrics::rmse(&predicted, &actual))),
                ("mae", Some(metrics::mae(&predicted, &actual))),
            ]),
        })
    };

    // Predicting the training mean for every protein.
    let train_mean = f64::from(targets[0].iter().sum::<f32>()) / targets[0].len() as f64;
    let test_actual: Vec<f64> = targets[2].iter().map(|v| f64::from(*v)).collect();
    let constant = vec![train_mean; test_actual.len()];
    let baseline = SplitMetrics {
        n: test_actual.len(),
        values: metric_map(&[
            ("rmse", Some(metrics::rmse(&constant, &test_actual))),
            ("mae", Some(metrics::mae(&constant, &test_actual))),
        ]),
    };

    let report = Report {
        task: Task::Solubility.name().to_string(),
        embedding_model: model_size.repo().to_string(),
        embedding_dim: hidden,
        best_epoch,
        epochs_run,
        learning_rate: config.learning_rate,
        batch_size: config.batch_size,
        seed: config.seed,
        seed_applied: true,
        validation: evaluate(1)?,
        test: evaluate(2)?,
        baseline,
    };
    save(
        &varmap,
        &Saved {
            out_dir,
            task: Task::Solubility,
            model_size,
            embedding_dim: hidden,
            report: &report,
            with_burial: false,
            datasets_root,
        },
    )?;
    Ok(report)
}

// ---------------------------------------------------------------------------
// Disorder
// ---------------------------------------------------------------------------

/// Train the per-residue disorder head.
///
/// Proteins are streamed from the cache one at a time: the full set of
/// per-residue embeddings is several gigabytes and does not need to be
/// resident.
///
/// # Errors
/// If the datasets or caches cannot be read, or training fails.
/// The disorder head's logits and labels for one split, plus the sequences.
struct Predictions {
    /// Per-protein logits and labels.
    proteins: Vec<crate::calibrate::Protein>,
    /// The sequences, in the same order, for sequence-level signals.
    sequences: Vec<String>,
}

/// Run the disorder head over a split and keep its logits.
///
/// # Errors
/// If the head is missing, its embedding cache is incomplete, or a cached
/// embedding disagrees with its label array.
fn disorder_predictions(
    head: &opendna_ml::heads::DisorderHead,
    cache: &MatrixCache,
    datasets_root: &Path,
    split: Split,
    hidden: usize,
    device: &candle_core::Device,
) -> Result<Predictions> {
    let rows = data::disorder(datasets_root, split)?;
    let mut proteins = Vec::with_capacity(rows.len());
    let mut sequences = Vec::with_capacity(rows.len());
    for row in &rows {
        let (values, residues) = cache.get(&row.sequence)?;
        if residues != row.disorder.len() {
            bail!("cached embedding has {residues} residues but the label array has {}", row.disorder.len());
        }
        let x = Tensor::from_vec(values, (residues, hidden), device)?;
        // Logits, not probabilities: the correction is fitted in logit space,
        // and recovering a logit from a rounded probability loses precision
        // exactly where the model is most confident.
        let logits = head.forward(&x)?.to_vec1::<f32>()?;
        proteins.push(crate::calibrate::Protein {
            logits: logits.into_iter().map(f64::from).collect(),
            labels: row.disorder.clone(),
        });
        sequences.push(row.sequence.clone());
    }
    Ok(Predictions { proteins, sequences })
}

/// Predictions against truth for one regressor, on two splits.
///
/// The validation split calibrates and the test split measures, which is the
/// same division the rest of the training code uses. Grouping matters:
/// mutations in one protein and residues in one chain move together, so they
/// are kept together rather than treated as independent observations.
///
/// # Errors
/// If the head is missing, or a cache does not cover the split.
fn regressor_predictions(
    task: Task,
    datasets_root: &Path,
    cache_dir: &Path,
    out_dir: &Path,
    model_size: ModelSize,
) -> Result<(Vec<crate::conformal::Group>, Vec<crate::conformal::Group>, &'static str)> {
    use crate::conformal::{Group, Prediction};

    let device = opendna_ml::best_device()?;
    let hidden = model_size.hidden_size();
    let heads = opendna_ml::heads::Heads::load(out_dir, model_size, &device)?;

    match task {
        Task::Solubility => {
            let head = heads
                .solubility
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("no solubility head in {}", out_dir.display()))?;
            let cache = MatrixCache::open(
                cache_dir,
                &crate::embed::embedding_cache_name(Task::Solubility, model_size),
                model_size.repo(),
                hidden,
            )?;
            let read = |split: Split| -> Result<Vec<Group>> {
                let mut out = Vec::new();
                for row in &data::solubility(datasets_root, split)? {
                    let (values, _) = cache.get(&row.sequence)?;
                    let predicted = f64::from(head.predict(&values, &device)?);
                    // One protein, one number: nothing to group.
                    out.push(Group(vec![Prediction { predicted, actual: row.solubility }]));
                }
                Ok(out)
            };
            Ok((read(Split::Validation)?, read(Split::Test)?, "percentage points"))
        }
        Task::Burial => {
            let head = heads
                .burial
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("no burial head in {}", out_dir.display()))?;
            let cache = MatrixCache::open(
                cache_dir,
                &crate::embed::embedding_cache_name(Task::Burial, model_size),
                model_size.repo(),
                hidden,
            )?;
            let read = |split: Split| -> Result<Vec<Group>> {
                let mut out = Vec::new();
                for row in &data::burial(datasets_root, split)? {
                    let (values, residues) = cache.get(&row.sequence)?;
                    if residues != row.rsa.len() {
                        bail!(
                            "cached embedding has {residues} residues but the label array has {}",
                            row.rsa.len()
                        );
                    }
                    let predicted = head.predict(&values, residues, &device)?;
                    // One group per chain: neighbouring residues are not
                    // independent of each other.
                    out.push(Group(
                        predicted
                            .iter()
                            .zip(&row.rsa)
                            .map(|(p, a)| Prediction { predicted: f64::from(*p), actual: f64::from(*a) })
                            .collect(),
                    ));
                }
                Ok(out)
            };
            Ok((read(Split::Validation)?, read(Split::Test)?, "RSA"))
        }
        Task::Stability => {
            let head = heads
                .stability
                .as_ref()
                .ok_or_else(|| anyhow::anyhow!("no stability head in {}", out_dir.display()))?;
            let vocab = opendna_ml::esm2::alphabet::TOKENS.len();
            let embeddings_cache = MatrixCache::open(
                cache_dir,
                &crate::embed::embedding_cache_name(Task::Stability, model_size),
                model_size.repo(),
                hidden,
            )?;
            let log_probs_cache = MatrixCache::open(
                cache_dir,
                &crate::embed::log_prob_cache_name(model_size),
                model_size.repo(),
                vocab,
            )?;
            let splits: Vec<Vec<data::StabilityRow>> = [Split::Validation, Split::Test]
                .iter()
                .map(|s| data::stability(datasets_root, *s))
                .collect::<Result<_>>()?;

            let mut embeddings = HashMap::new();
            let mut log_probs = HashMap::new();
            for split in &splits {
                for row in split {
                    if !embeddings.contains_key(&row.wt_sequence) {
                        embeddings.insert(row.wt_sequence.clone(), embeddings_cache.get(&row.wt_sequence)?);
                        log_probs.insert(row.wt_sequence.clone(), log_probs_cache.get(&row.wt_sequence)?);
                    }
                }
            }

            let build = |rows: &[data::StabilityRow]| -> Result<Vec<Group>> {
                let built = stability_features(rows, &embeddings, &log_probs, hidden, vocab, None)?;
                let tensor = Tensor::from_slice(&built.features, (built.len(), built.width), &device)?;
                let predicted: Vec<f64> =
                    head.forward(&tensor)?.to_vec1::<f32>()?.into_iter().map(f64::from).collect();
                // Grouped by wild-type sequence: hundreds of mutations in one
                // small protein are not hundreds of independent tests.
                let mut by_protein: std::collections::BTreeMap<&str, Vec<Prediction>> =
                    std::collections::BTreeMap::new();
                for (index, row) in rows.iter().enumerate() {
                    by_protein.entry(row.wt_sequence.as_str()).or_default().push(Prediction {
                        predicted: predicted[index],
                        actual: f64::from(built.targets[index]),
                    });
                }
                Ok(by_protein.into_values().map(Group).collect())
            };
            Ok((build(&splits[0])?, build(&splits[1])?, "kcal/mol"))
        }
        Task::Disorder => bail!(
            "disorder is a classifier, not a regressor. Its uncertainty is a calibrated \
             probability -- see `calibrate` and `abstain`."
        ),
    }
}

/// Fit conformal prediction intervals and measure the coverage they achieve.
///
/// Calibrates on the validation split and measures on the test split, which
/// is deliberately the hard case: the dataset builder clusters by sequence
/// identity so near-identical proteins cannot straddle a split, which is
/// exactly the exchangeability conformal assumes. Whether the guarantee
/// survives that is the question, and it is measured rather than asserted.
///
/// # Errors
/// If the head is missing or a cache does not cover the split.
pub fn intervals(
    task: Task,
    datasets_root: &Path,
    cache_dir: &Path,
    out_dir: &Path,
    model_size: ModelSize,
    apply: bool,
) -> Result<()> {
    let (calibration, held_out, units) =
        regressor_predictions(task, datasets_root, cache_dir, out_dir, model_size)?;
    println!(
        "{}: {} calibration groups, {} held-out groups\n",
        task.name(),
        calibration.len(),
        held_out.len()
    );

    let report = crate::conformal::run(&calibration, &held_out, units);
    print!("{}", report.summary());

    if !apply {
        println!("\nNothing written. Re-run with --apply to record these intervals beside the head.");
        return Ok(());
    }

    let files = opendna_ml::heads::HeadFiles::in_dir(out_dir, task.name());
    let mut metadata = opendna_ml::heads::read_metadata(&files, task.name(), model_size)?;
    metadata.intervals = Some(opendna_ml::heads::Intervals {
        units: units.to_string(),
        levels: report
            .levels
            .iter()
            .map(|l| opendna_ml::heads::IntervalLevel {
                nominal: l.nominal,
                half_width: l.half_width,
                achieved: l.achieved,
            })
            .collect(),
        held_out_mae: report.mae,
    });
    std::fs::write(&files.metadata, serde_json::to_string_pretty(&metadata)? + "\n")
        .with_context(|| format!("writing {}", files.metadata.display()))?;
    println!("\nRecorded in {}.", files.metadata.display());
    println!("No prediction changes: an interval is reported beside the number, not instead of it.");
    Ok(())
}

/// Measure whether abstaining on low-confidence predictions would help.
///
/// Reports a risk-coverage curve for each candidate signal against two
/// references: discarding residues at random, which is what no signal looks
/// like, and an oracle that discards exactly the errors, which is what there
/// was to win. Nothing is changed by this command. Its output is the evidence
/// for or against writing an abstention rule at all, and "against" is a
/// result.
///
/// # Errors
/// If the head is missing or its embedding cache is incomplete.
pub fn abstain(
    task: Task,
    datasets_root: &Path,
    cache_dir: &Path,
    out_dir: &Path,
    model_size: ModelSize,
    apply: bool,
) -> Result<()> {
    if task != Task::Disorder {
        bail!(
            "only the disorder head has a calibrated probability to abstain on. A regressor \
             needs a prediction interval, which is a different measurement."
        );
    }
    let device = opendna_ml::best_device()?;
    let hidden = model_size.hidden_size();
    let heads = opendna_ml::heads::Heads::load(out_dir, model_size, &device)?;
    let head = heads
        .disorder
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("no disorder head in {}", out_dir.display()))?;
    let cache = MatrixCache::open(
        cache_dir,
        &crate::embed::embedding_cache_name(Task::Disorder, model_size),
        model_size.repo(),
        hidden,
    )?;

    let files = opendna_ml::heads::HeadFiles::in_dir(out_dir, task.name());
    let metadata = opendna_ml::heads::read_metadata(&files, task.name(), model_size)?;
    let calibration = metadata.calibration.as_ref();
    // The head applies its map inside `predict`, so the measurement has to
    // apply the same one to the raw logits it reads here.
    let mapping = calibration.map_or_else(crate::calibrate::Mapping::default, |c| {
        crate::calibrate::Mapping { slope: c.slope, intercept: c.intercept }
    });
    let threshold = calibration.map_or(0.5, |c| c.threshold);
    println!(
        "Head {} at threshold {threshold:.2}, map p = sigmoid({:.4} * logit + {:.4})\n",
        metadata.version, mapping.slope, mapping.intercept
    );

    let test = disorder_predictions(head, &cache, datasets_root, Split::Test, hidden, &device)?;
    // Novelty is measured against the training split, which is what the head
    // actually saw.
    let training: Vec<String> =
        data::disorder(datasets_root, Split::Train)?.into_iter().map(|r| r.sequence).collect();
    let novel = crate::selective::novelty(&training, &test.sequences, 5);

    let report = crate::selective::run(&test.proteins, mapping, threshold, &novel, 20_260_101);
    print!("{}", report.summary());

    if !apply {
        println!(
            "\nNothing was changed. Re-run with --apply to record the trade-off beside the \
             head, which lets the CLI quote a measured error rate for a chosen cut-off."
        );
        return Ok(());
    }

    let trade_off = crate::selective::Report::trade_off(&test.proteins, mapping, threshold);
    let curve: Vec<opendna_ml::heads::Reliability> = trade_off
        .iter()
        .map(|(min_confidence, coverage, error)| opendna_ml::heads::Reliability {
            min_confidence: *min_confidence,
            coverage: *coverage,
            error: *error,
        })
        .collect();
    let mut metadata = opendna_ml::heads::read_metadata(&files, task.name(), model_size)?;
    metadata.abstention = Some(opendna_ml::heads::Abstention {
        curve,
        base_error: report.base_risk,
        captured: report.captured(),
        rejected: report.rejected(),
    });
    std::fs::write(&files.metadata, serde_json::to_string_pretty(&metadata)? + "\n")
        .with_context(|| format!("writing {}", files.metadata.display()))?;
    println!("\nRecorded in {}.", files.metadata.display());
    println!(
        "No prediction changes: the confidence is derived from the probability the head \
         already reports, so the version and build id stay as they are. What is recorded is \
         what a cut-off was measured to buy."
    );
    Ok(())
}

/// Measure a head's calibration and, with `apply`, record the correction.
///
/// The validation split fits, the test split reports, and the test split is
/// read once. Choosing the form of the correction by its test score would turn
/// the test number into a training score, so the choice is made on validation.
///
/// Only the disorder head is in scope. The other three are regressors: asking
/// whether a predicted ddG of 1.2 kcal/mol "means" anything probabilistic is a
/// different question with a different answer, and conflating the two would be
/// the sort of thing this project exists not to do.
///
/// # Errors
/// If the head is not installed, if its embedding cache is missing, or if the
/// metadata cannot be rewritten.
pub fn calibrate(
    task: Task,
    datasets_root: &Path,
    cache_dir: &Path,
    out_dir: &Path,
    model_size: ModelSize,
    apply: bool,
) -> Result<()> {
    if task != Task::Disorder {
        bail!(
            "only the disorder head emits probabilities, so only it can be calibrated this way. \
             {} is a regressor: what it needs is a prediction interval, which is a different \
             measurement.",
            task.name()
        );
    }

    let device = opendna_ml::best_device()?;
    let hidden = model_size.hidden_size();
    let heads = opendna_ml::heads::Heads::load(out_dir, model_size, &device)?;
    let head = heads
        .disorder
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("no disorder head in {}; train it first", out_dir.display()))?;

    let cache = MatrixCache::open(
        cache_dir,
        &crate::embed::embedding_cache_name(Task::Disorder, model_size),
        model_size.repo(),
        hidden,
    )?;

    // The same splits the training run used, so "validation" and "test" mean
    // here exactly what they meant there.
    let validation =
        disorder_predictions(head, &cache, datasets_root, Split::Validation, hidden, &device)?.proteins;
    let test = disorder_predictions(head, &cache, datasets_root, Split::Test, hidden, &device)?.proteins;
    let residues: usize = validation.iter().chain(&test).map(|p| p.logits.len()).sum();
    println!(
        "{} validation and {} test proteins, {residues} assessed residues\n",
        validation.len(),
        test.len()
    );

    let report = crate::calibrate::run(&validation, &test, 20_260_101);
    print!("{}", report.summary());

    if !apply {
        println!("\nNothing written. Re-run with --apply to record this in the head's metadata.");
        return Ok(());
    }
    // Recorded either way. A correction that does not beat its own noise is
    // stored as the identity map, which changes no number but records that
    // the question was asked and answered. "Measured and unnecessary" and
    // "never measured" are different states, and a model card has to be able
    // to tell a reader which one it is looking at.
    let worth = report.worth_applying();
    let mapping = if worth { report.mapping } else { crate::calibrate::Mapping::default() };

    let files = opendna_ml::heads::HeadFiles::in_dir(out_dir, task.name());
    let mut metadata = opendna_ml::heads::read_metadata(&files, task.name(), model_size)?;
    metadata.calibration = Some(opendna_ml::heads::Calibration {
        slope: mapping.slope,
        intercept: mapping.intercept,
        method: if worth { report.form.name().to_string() } else { "none".to_string() },
        test_ece_before: report.test.0.point,
        test_ece_after: if worth { report.test.1.point } else { report.test.0.point },
        threshold: report.recommended.0.threshold,
        test_mcc: report.recommended.1.mcc,
    });
    // The head now reports different numbers from the same weights, so it is a
    // different model and needs a different build id. Recomputed here rather
    // than left for a later `stamp`, because a calibrated head carrying the
    // uncalibrated one's identity would pass a lockfile check it should fail.
    let metrics_path = out_dir.join(format!("{}.metrics.json", task.name()));
    if let Ok(text) = std::fs::read_to_string(&metrics_path) {
        if let Ok(metrics) = serde_json::from_str::<Report>(&text) {
            metadata.version = head_version(task).to_string();
            metadata.build =
                Some(build_id_with(task, model_size, &metrics, datasets_root, metadata.calibration.as_ref()));
        }
    }
    let json = serde_json::to_string_pretty(&metadata)?;
    std::fs::write(&files.metadata, json + "\n")
        .with_context(|| format!("writing {}", files.metadata.display()))?;
    println!("\nRecorded in {}.", files.metadata.display());
    if worth {
        println!("The head's predictions will now pass through this mapping.");
    } else {
        println!(
            "Recorded as the identity, so no prediction changes. What is recorded is that \
             calibration was measured on these weights and found not worth correcting."
        );
    }
    Ok(())
}

pub fn disorder(
    datasets_root: &Path,
    cache_dir: &Path,
    out_dir: &Path,
    model_size: ModelSize,
    config: &TrainConfig,
) -> Result<Report> {
    let device = opendna_ml::best_device()?;
    let hidden = model_size.hidden_size();
    let cache = MatrixCache::open(
        cache_dir,
        &crate::embed::embedding_cache_name(Task::Disorder, model_size),
        model_size.repo(),
        hidden,
    )?;

    let splits: Vec<Vec<data::DisorderRow>> =
        Split::all().iter().map(|s| data::disorder(datasets_root, *s)).collect::<Result<_>>()?;
    println!(
        "{} train / {} validation / {} test proteins",
        splits[0].len(),
        splits[1].len(),
        splits[2].len()
    );

    // Disordered residues are the minority; weight them up so the head does
    // not simply predict "ordered" everywhere.
    let positives: usize = splits[0].iter().map(|r| r.disorder.iter().filter(|d| **d).count()).sum();
    let total: usize = splits[0].iter().map(|r| r.disorder.len()).sum();
    let positive_weight = (total - positives) as f64 / positives.max(1) as f64;
    println!("positive rate {:.3}, positive weight {positive_weight:.2}", positives as f64 / total as f64);

    let mut varmap = VarMap::new();
    let vb = VarBuilder::from_varmap(&varmap, DType::F32, &device);
    let head = DisorderHead::new(hidden, vb)?;
    seed_weights(&varmap, config.seed)?;
    let mut optimizer =
        AdamW::new(varmap.all_vars(), ParamsAdamW { lr: config.learning_rate, ..Default::default() })?;

    let weight = positive_weight as f32;
    let forward = |head: &DisorderHead, row: &data::DisorderRow| -> Result<Tensor> {
        let (values, rows) = cache.get(&row.sequence)?;
        if rows != row.disorder.len() {
            bail!("cached embedding has {rows} residues but the label array has {}", row.disorder.len());
        }
        let x = Tensor::from_vec(values, (rows, hidden), &device)?;
        Ok(head.forward(&x)?)
    };

    let evaluate = |head: &DisorderHead, rows: &[data::DisorderRow]| -> Result<(Vec<f64>, Vec<bool>)> {
        let mut scores = Vec::new();
        let mut labels = Vec::new();
        for row in rows {
            let logits = forward(head, row)?.to_vec1::<f32>()?;
            scores.extend(logits.into_iter().map(|v| f64::from(1.0 / (1.0 + (-v).exp()))));
            labels.extend(row.disorder.iter().copied());
        }
        Ok((scores, labels))
    };

    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut order: Vec<usize> = (0..splits[0].len()).collect();
    let mut best = f64::NEG_INFINITY;
    let mut best_epoch = 0;
    let mut epochs_run = 0;
    let mut best_weights: Option<Vec<u8>> = None;

    for epoch in 1..=config.epochs {
        epochs_run = epoch;
        order.shuffle(&mut rng);
        let mut losses = Vec::new();
        for &i in &order {
            let row = &splits[0][i];
            let logits = forward(&head, row)?;
            let targets: Vec<f32> = row.disorder.iter().map(|d| f32::from(u8::from(*d))).collect();
            let y = Tensor::from_vec(targets, row.disorder.len(), &device)?;
            // Weighted binary cross-entropy, written out because candle's
            // helper has no class-weight parameter.
            let max = logits.relu()?;
            let stable = logits.abs()?.neg()?.exp()?.affine(1.0, 1.0)?.log()?;
            let base = (max + stable)?;
            let weights = (y.affine(f64::from(weight) - 1.0, 1.0))?;
            let loss = (((base - (&logits * &y)?)? * weights)?).mean_all()?;
            optimizer.backward_step(&loss)?;
            losses.push(f64::from(loss.to_scalar::<f32>()?));
        }

        let (scores, labels) = evaluate(&head, &splits[1])?;
        let score = metrics::auroc(&scores, &labels).unwrap_or(f64::NEG_INFINITY);
        println!(
            "epoch {epoch:3}  loss {:.4}  validation auroc {score:.4}{}",
            mean(&losses),
            if score > best { "  *" } else { "" }
        );
        if score > best {
            best = score;
            best_epoch = epoch;
            let tmp = std::env::temp_dir().join(format!("opendna-dis-{}.safetensors", std::process::id()));
            varmap.save(&tmp)?;
            best_weights = Some(std::fs::read(&tmp)?);
            let _ = std::fs::remove_file(&tmp);
        } else if epoch - best_epoch >= config.patience {
            println!("No improvement for {} epochs; stopping.", config.patience);
            break;
        }
    }

    if let Some(bytes) = best_weights {
        let tmp = std::env::temp_dir().join(format!("opendna-dis-best-{}.safetensors", std::process::id()));
        std::fs::write(&tmp, bytes)?;
        varmap.load(&tmp)?;
        let _ = std::fs::remove_file(&tmp);
    }

    let summarise = |rows: &[data::DisorderRow]| -> Result<SplitMetrics> {
        let (scores, labels) = evaluate(&head, rows)?;
        let confusion = metrics::Confusion::at(&scores, &labels, 0.5);
        Ok(SplitMetrics {
            n: labels.len(),
            values: metric_map(&[
                ("auroc", metrics::auroc(&scores, &labels)),
                ("auprc", metrics::auprc(&scores, &labels)),
                ("mcc", confusion.mcc()),
                ("f1", confusion.f1()),
                ("balanced_accuracy", confusion.balanced_accuracy()),
            ]),
        })
    };

    // Always predicting "disordered" — the floor any real model must clear.
    let (_, test_labels) = evaluate(&head, &splits[2])?;
    let ones = vec![1.0; test_labels.len()];
    let confusion = metrics::Confusion::at(&ones, &test_labels, 0.5);
    let baseline = SplitMetrics {
        n: test_labels.len(),
        values: metric_map(&[
            ("auroc", Some(0.5)),
            ("auprc", Some(test_labels.iter().filter(|l| **l).count() as f64 / test_labels.len() as f64)),
            ("f1", confusion.f1()),
            ("balanced_accuracy", confusion.balanced_accuracy()),
        ]),
    };

    let report = Report {
        task: Task::Disorder.name().to_string(),
        embedding_model: model_size.repo().to_string(),
        embedding_dim: hidden,
        best_epoch,
        epochs_run,
        learning_rate: config.learning_rate,
        batch_size: config.batch_size,
        seed: config.seed,
        seed_applied: true,
        validation: summarise(&splits[1])?,
        test: summarise(&splits[2])?,
        baseline,
    };
    save(
        &varmap,
        &Saved {
            out_dir,
            task: Task::Disorder,
            model_size,
            embedding_dim: hidden,
            report: &report,
            with_burial: false,
            datasets_root,
        },
    )?;
    Ok(report)
}

/// Train the per-residue solvent accessibility head.
///
/// Deliberately the disorder loop with the loss swapped, because the two
/// tasks differ only in what a residue's label is. The differences that do
/// matter:
///
/// * The target is continuous, so this is squared error on the squashed
///   output rather than weighted cross-entropy. There is no class imbalance
///   to correct — roughly half of all residues are buried.
/// * Spearman selects the epoch. Accessibility is used to rank positions
///   ("which residues are exposed enough to mutate"), and a model that
///   ranks well with a small systematic offset is more useful than one with
///   lower error and worse ordering.
/// * The baseline is predicting the training mean everywhere, which is the
///   floor for any regression and the one the solubility head uses.
///
/// # Errors
/// If the dataset or the embedding cache cannot be read, or training fails.
pub fn burial(
    datasets_root: &Path,
    cache_dir: &Path,
    out_dir: &Path,
    model_size: ModelSize,
    config: &TrainConfig,
) -> Result<Report> {
    let device = opendna_ml::best_device()?;
    let hidden = model_size.hidden_size();
    let cache = MatrixCache::open(
        cache_dir,
        &crate::embed::embedding_cache_name(Task::Burial, model_size),
        model_size.repo(),
        hidden,
    )?;

    let splits: Vec<Vec<data::BurialRow>> =
        Split::all().iter().map(|s| data::burial(datasets_root, *s)).collect::<Result<_>>()?;
    println!("{} train / {} validation / {} test domains", splits[0].len(), splits[1].len(), splits[2].len());
    let residues: usize = splits[0].iter().map(|r| r.rsa.len()).sum();
    let training_mean =
        splits[0].iter().flat_map(|r| &r.rsa).map(|v| f64::from(*v)).sum::<f64>() / residues as f64;
    println!("{residues} training residues, mean accessibility {training_mean:.4}");

    let mut varmap = VarMap::new();
    let vb = VarBuilder::from_varmap(&varmap, DType::F32, &device);
    let head = BurialHead::new(hidden, vb)?;
    seed_weights(&varmap, config.seed)?;
    let mut optimizer =
        AdamW::new(varmap.all_vars(), ParamsAdamW { lr: config.learning_rate, ..Default::default() })?;

    let forward = |head: &BurialHead, row: &data::BurialRow| -> Result<Tensor> {
        let (values, rows) = cache.get(&row.sequence)?;
        if rows != row.rsa.len() {
            bail!("cached embedding has {rows} residues but the label array has {}", row.rsa.len());
        }
        let x = Tensor::from_vec(values, (rows, hidden), &device)?;
        Ok(head.forward(&x)?)
    };

    let evaluate = |head: &BurialHead, rows: &[data::BurialRow]| -> Result<(Vec<f64>, Vec<f64>)> {
        let mut predicted = Vec::new();
        let mut actual = Vec::new();
        for row in rows {
            let values = BurialHead::accessibility(&forward(head, row)?)?.to_vec1::<f32>()?;
            predicted.extend(values.into_iter().map(f64::from));
            actual.extend(row.rsa.iter().map(|v| f64::from(*v)));
        }
        Ok((predicted, actual))
    };

    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut order: Vec<usize> = (0..splits[0].len()).collect();
    let mut best = f64::NEG_INFINITY;
    let mut best_epoch = 0;
    let mut epochs_run = 0;
    let mut best_weights: Option<Vec<u8>> = None;

    for epoch in 1..=config.epochs {
        epochs_run = epoch;
        order.shuffle(&mut rng);
        let mut losses = Vec::new();
        for &i in &order {
            let row = &splits[0][i];
            let predicted = BurialHead::accessibility(&forward(&head, row)?)?;
            let y = Tensor::from_vec(row.rsa.clone(), row.rsa.len(), &device)?;
            let loss = (predicted - y)?.sqr()?.mean_all()?;
            optimizer.backward_step(&loss)?;
            losses.push(f64::from(loss.to_scalar::<f32>()?));
        }

        let (predicted, actual) = evaluate(&head, &splits[1])?;
        let score = metrics::spearman(&predicted, &actual).unwrap_or(f64::NEG_INFINITY);
        println!(
            "epoch {epoch:3}  loss {:.5}  validation spearman {score:.4}{}",
            mean(&losses),
            if score > best { "  *" } else { "" }
        );
        if score > best {
            best = score;
            best_epoch = epoch;
            let tmp = std::env::temp_dir().join(format!("opendna-bur-{}.safetensors", std::process::id()));
            varmap.save(&tmp)?;
            best_weights = Some(std::fs::read(&tmp)?);
            let _ = std::fs::remove_file(&tmp);
        } else if epoch - best_epoch >= config.patience {
            println!("No improvement for {} epochs; stopping.", config.patience);
            break;
        }
    }

    if let Some(bytes) = best_weights {
        let tmp = std::env::temp_dir().join(format!("opendna-bur-best-{}.safetensors", std::process::id()));
        std::fs::write(&tmp, bytes)?;
        varmap.load(&tmp)?;
        let _ = std::fs::remove_file(&tmp);
    }

    let summarise = |rows: &[data::BurialRow]| -> Result<SplitMetrics> {
        let (predicted, actual) = evaluate(&head, rows)?;
        Ok(SplitMetrics {
            n: actual.len(),
            values: metric_map(&[
                ("spearman", metrics::spearman(&predicted, &actual)),
                ("pearson", metrics::pearson(&predicted, &actual)),
                ("rmse", Some(metrics::rmse(&predicted, &actual))),
                ("mae", Some(metrics::mae(&predicted, &actual))),
            ]),
        })
    };

    // Predicting the training mean for every residue: the floor.
    let (_, test_actual) = evaluate(&head, &splits[2])?;
    let flat = vec![training_mean; test_actual.len()];
    let baseline = SplitMetrics {
        n: test_actual.len(),
        values: metric_map(&[
            ("rmse", Some(metrics::rmse(&flat, &test_actual))),
            ("mae", Some(metrics::mae(&flat, &test_actual))),
        ]),
    };

    let report = Report {
        task: Task::Burial.name().to_string(),
        embedding_model: model_size.repo().to_string(),
        embedding_dim: hidden,
        best_epoch,
        epochs_run,
        learning_rate: config.learning_rate,
        batch_size: config.batch_size,
        seed: config.seed,
        seed_applied: true,
        validation: summarise(&splits[1])?,
        test: summarise(&splits[2])?,
        baseline,
    };
    save(
        &varmap,
        &Saved {
            out_dir,
            task: Task::Burial,
            model_size,
            embedding_dim: hidden,
            report: &report,
            with_burial: false,
            datasets_root,
        },
    )?;
    Ok(report)
}
