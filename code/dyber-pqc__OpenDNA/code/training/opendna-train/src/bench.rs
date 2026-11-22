//! External benchmarks: running the heads on reference sets that other
//! published methods also report on.
//!
//! Our own test splits measure how well a head does on data drawn the way its
//! training data was drawn. They cannot say how the head compares to the
//! field, because no other method reports a number on them. That makes every
//! figure in our model cards unsituated: 0.883 AUROC is only meaningful next
//! to what else scores on the same proteins, under the same protocol.
//!
//! # No benchmark data is redistributed
//!
//! Reference sets are read from a path you supply; nothing here downloads
//! them, and no built dataset contains them. Some of these sets carry no
//! licence at all — `protddg-bench` has no `LICENSE` file — which permits
//! evaluation but not republication, the same reasoning that kept `FireProtDB`
//! out of the stability dataset. The CAID3 reference is fetched by the
//! `datasets` workspace, which records its licence and checksum.
//!
//! # The protocol is the benchmark
//!
//! A reference set is not just a list of proteins; it is a list of proteins
//! plus rules about which residues count. CAID masks residues it did not
//! assess, and scoring those residues anyway would produce a number that
//! looks like an AUROC and is not comparable to any published one. Where this
//! module departs from a published protocol, it says so in the report rather
//! than in a comment nobody reads.

use crate::metrics::{self, Confusion};
use anyhow::{anyhow, bail, Context, Result};
use opendna_ml::esm2::{Esm2, ModelSize};
use opendna_ml::heads::{read_metadata, HeadFiles, Heads};
use serde::{Deserialize, Serialize};
use std::fmt::Write as _;
use std::path::Path;

/// One protein of a CAID reference set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    /// `DisProt` identifier, as the reference file gives it.
    pub id: String,
    /// Amino-acid sequence.
    pub sequence: String,
    /// Per-residue truth. `None` where the reference did not assess the
    /// residue, which is a third of CAID3 Disorder-PDB and must be excluded
    /// from scoring rather than treated as ordered.
    pub labels: Vec<Option<bool>>,
}

/// Parse a CAID reference set.
///
/// The format is three lines per protein: `>id`, the sequence, and a string of
/// `1` (disordered), `0` (ordered) and `-` (not assessed), one character per
/// residue.
///
/// # Errors
/// If a record is truncated, a label character is unrecognised, or the label
/// string is not the same length as the sequence — any of which would silently
/// shift every label along the chain.
pub fn read_caid_reference(text: &str) -> Result<Vec<Reference>> {
    let lines: Vec<&str> = text.lines().map(str::trim_end).filter(|l| !l.is_empty()).collect();
    if !lines.len().is_multiple_of(3) {
        bail!(
            "{} non-empty lines is not a whole number of 3-line records; the file looks truncated",
            lines.len()
        );
    }
    let mut out = Vec::with_capacity(lines.len() / 3);
    for (record, chunk) in lines.as_chunks::<3>().0.iter().enumerate() {
        let id = chunk[0]
            .strip_prefix('>')
            .ok_or_else(|| anyhow!("record {record}: expected a '>' header, found {:?}", chunk[0]))?;
        let (sequence, labels) = (chunk[1], chunk[2]);
        if sequence.chars().count() != labels.chars().count() {
            bail!(
                "{id}: {} residues but {} labels; one of the two is wrong and every label after \
                 the difference would be attributed to the wrong residue",
                sequence.chars().count(),
                labels.chars().count()
            );
        }
        let labels = labels
            .chars()
            .map(|c| match c {
                '1' => Ok(Some(true)),
                '0' => Ok(Some(false)),
                '-' => Ok(None),
                other => Err(anyhow!("{id}: unrecognised label character {other:?}")),
            })
            .collect::<Result<Vec<_>>>()?;
        out.push(Reference { id: id.to_string(), sequence: sequence.to_string(), labels });
    }
    if out.is_empty() {
        bail!("the reference file contains no records");
    }
    Ok(out)
}

/// Assays excluded from the stability suite because they are the head's own
/// training data.
///
/// `ProteinGym`'s stability assays are dominated by Tsuboyama 2023 — 64 of its
/// 217 assays use `cDNA display proteolysis`, which is that paper's method —
/// and that is exactly the dataset the stability head was trained on. Scoring
/// them would measure memorisation and report it as generalisation.
const TRAINING_DATA_ASSAY: &str = "cDNA display proteolysis";

/// Assay types in the `ProteinGym` reference file that measure something
/// stability-like: variant abundance, thermostability, or explicit protein
/// stability. None of them is ΔΔG, which is stated in the report's caveats.
const STABILITY_LIKE: [&str; 4] =
    ["Protein stability", "thermostability", "protein abundance", "VAMP-seq, drug sensitivity"];

/// Load the `ProteinGym` stability-like assays that are not training data.
///
/// `reference` is `DMS_substitutions.csv`, which carries each assay's
/// wild-type sequence; `assay_dir` holds the per-assay score files. An assay
/// named in the reference but missing from `assay_dir` is skipped with a note
/// rather than failing, so the suite can be run with a subset downloaded.
///
/// # Errors
/// If the reference file cannot be read or parsed, if no usable assay is
/// found, or if an assay's mutant numbering disagrees with its sequence.
pub fn read_proteingym(reference: &Path, assay_dir: &Path) -> Result<(Vec<Assay>, Vec<String>)> {
    let mut index =
        csv::Reader::from_path(reference).with_context(|| format!("reading {}", reference.display()))?;
    let mut assays = Vec::new();
    let mut notes = Vec::new();
    let mut excluded_training = 0usize;
    for record in index.deserialize::<std::collections::HashMap<String, String>>() {
        let row = record.context("parsing a row of the ProteinGym reference file")?;
        let get = |key: &str| row.get(key).map(String::as_str).unwrap_or_default();
        let selection = get("selection_type").to_string();
        if selection == TRAINING_DATA_ASSAY {
            excluded_training += 1;
            continue;
        }
        if !STABILITY_LIKE.contains(&selection.as_str()) {
            continue;
        }
        let id = get("DMS_id").to_string();
        let path = assay_dir.join(get("DMS_filename"));
        if !path.is_file() {
            notes.push(format!("{id}: not downloaded ({})", path.display()));
            continue;
        }
        let sequence = get("target_seq").to_string();
        let residues: Vec<char> = sequence.chars().collect();
        let mut scores =
            csv::Reader::from_path(&path).with_context(|| format!("reading {}", path.display()))?;
        let mut mutants = Vec::new();
        for row in scores.deserialize::<std::collections::HashMap<String, String>>() {
            let row = row.with_context(|| format!("parsing {}", path.display()))?;
            let label = row.get("mutant").map(String::as_str).unwrap_or_default();
            // Multi-substitution entries are joined with ':' and the head
            // only models single substitutions.
            if label.contains(':') {
                continue;
            }
            let score: f64 = match row.get("DMS_score").and_then(|s| s.parse().ok()) {
                Some(s) => s,
                None => continue,
            };
            let (position, mutant) =
                parse_mutant(label, &residues).with_context(|| format!("{id}: in {}", path.display()))?;
            mutants.push((position, mutant, score));
        }
        if mutants.is_empty() {
            notes.push(format!("{id}: no single substitutions"));
            continue;
        }
        assays.push(Assay { id, selection, sequence, mutants });
    }
    notes.push(format!("{excluded_training} assays excluded as training data ({TRAINING_DATA_ASSAY})"));
    if assays.is_empty() {
        bail!("no usable assay found in {}; notes: {}", assay_dir.display(), notes.join("; "));
    }
    Ok((assays, notes))
}

/// What a benchmark run measured.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchReport {
    /// Which reference set.
    pub suite: String,
    /// Embedding model the head was trained against.
    pub embedding_model: String,
    /// Head version from its metadata.
    pub head_version: String,
    /// Whether the head that was scored takes burial features.
    ///
    /// Two heads can be trained on the same embedder and scored on the same
    /// suite and differ only in this. Without it a report cannot say which of
    /// them it measured, and a card showing both presents a variant that was
    /// rejected beside the one that ships. Defaults to false so reports
    /// written before the field existed still parse; the one of those that was
    /// a with-burial run has been corrected by hand.
    #[serde(default)]
    pub with_burial: bool,
    /// Proteins scored.
    pub proteins: usize,
    /// Residues the reference assessed, and which therefore count.
    pub residues_scored: usize,
    /// Residues the reference did not assess, excluded from every metric.
    pub residues_masked: usize,
    /// Fraction of scored residues that are positive, for a classification
    /// suite. `None` for a regression suite, where it has no meaning —
    /// serialising a `NaN` here would write a JSON `null` that reads as a
    /// measurement someone forgot to take.
    pub positive_rate: Option<f64>,
    /// Proteins longer than the embedding model's limit, covered by
    /// overlapping windows rather than a single pass.
    pub windowed: usize,
    /// Metrics over all scored residues pooled, which is what CAID ranks on.
    pub pooled: Vec<(String, f64)>,
    /// Mean over per-protein AUROCs, and the count it averages over. Reported
    /// because the pooled figure treats 100k correlated residues as 100k
    /// independent observations, so it has a far tighter apparent confidence
    /// than the ~300 proteins behind it justify.
    pub per_protein_auroc: Option<f64>,
    /// Proteins that could not get a per-protein AUROC because every assessed
    /// residue had the same label.
    pub per_protein_undefined: usize,
    /// Anything about this run that a reader needs in order to not
    /// over-read it.
    pub caveats: Vec<String>,
}

/// One protein's predictions, reduced to the residues the reference assessed.
struct Scored {
    scores: Vec<f64>,
    labels: Vec<bool>,
    masked: usize,
    windowed: bool,
}

/// Predict disorder for one protein and keep only the assessed residues.
///
/// Dropping the masked residues here, rather than at the call site, is what
/// keeps the mask from being forgotten: there is no path through this
/// function that returns a score without its matching label.
fn score_one_protein(
    entry: &Reference,
    head: &opendna_ml::heads::DisorderHead,
    model: &Esm2,
    device: &candle_core::Device,
    limit: usize,
) -> Result<Scored> {
    let windowed = entry.sequence.chars().count() > limit;
    let embedding = if windowed {
        model.embed_windowed(&entry.sequence)?
    } else {
        model.embed(&[entry.sequence.as_str()])?.remove(0)
    };
    let probabilities = head
        .predict(&embedding.per_residue, embedding.length, device)
        .with_context(|| format!("predicting disorder for {} ({} residues)", entry.id, embedding.length))?;
    if probabilities.len() != entry.labels.len() {
        bail!("{}: {} predictions for {} labels", entry.id, probabilities.len(), entry.labels.len());
    }

    // Only assessed residues contribute. This is the whole reason the
    // reference ships a mask, and ignoring it would inflate every metric by
    // scoring predictions against labels that do not exist.
    let mut out = Scored { scores: Vec::new(), labels: Vec::new(), masked: 0, windowed };
    for (probability, label) in probabilities.iter().zip(&entry.labels) {
        match label {
            Some(truth) => {
                out.scores.push(f64::from(*probability));
                out.labels.push(*truth);
            }
            None => out.masked += 1,
        }
    }
    Ok(out)
}

/// Score the disorder head against a CAID reference set.
///
/// # Errors
/// If the head or embedding model cannot be loaded, if the head was trained
/// against a different embedding model, or if inference fails.
pub fn disorder(
    reference: &[Reference],
    models_dir: &Path,
    embedder: ModelSize,
    threshold: f64,
) -> Result<BenchReport> {
    let device = opendna_ml::best_device()?;
    println!("Loading on {}", opendna_ml::device_description(&device));
    let heads = Heads::load(models_dir, embedder, &device)
        .with_context(|| format!("loading heads from {}", models_dir.display()))?;
    let head =
        heads.disorder.as_ref().ok_or_else(|| anyhow!("no disorder head in {}", models_dir.display()))?;
    // The version identifies which trained head produced these numbers, so a
    // benchmark report cannot be mistaken for one from a different run.
    let metadata = read_metadata(&HeadFiles::in_dir(models_dir, "disorder"), "disorder", embedder)?;
    let model = Esm2::from_pretrained(embedder, &device)?;
    let limit = model.max_residues();

    let mut scores = Vec::new();
    let mut labels = Vec::new();
    let mut per_protein = Vec::new();
    let mut undefined = 0usize;
    let mut masked = 0usize;
    let mut windowed = 0usize;

    let bar = indicatif::ProgressBar::new(reference.len() as u64);
    for entry in reference {
        let one = score_one_protein(entry, head, &model, &device, limit)?;
        if one.windowed {
            windowed += 1;
        }
        masked += one.masked;
        match metrics::auroc(&one.scores, &one.labels) {
            Some(a) => per_protein.push(a),
            None => undefined += 1,
        }
        scores.extend_from_slice(&one.scores);
        labels.extend_from_slice(&one.labels);
        bar.inc(1);
    }
    bar.finish_and_clear();

    if scores.is_empty() {
        bail!("the reference set assessed no residues, so there is nothing to score");
    }
    let positives = labels.iter().filter(|l| **l).count();
    let confusion = Confusion::at(&scores, &labels, threshold);
    let mut pooled = Vec::new();
    for (name, value) in [
        ("auroc", metrics::auroc(&scores, &labels)),
        ("auprc", metrics::auprc(&scores, &labels)),
        ("mcc", confusion.mcc()),
        ("f1", confusion.f1()),
        ("balanced_accuracy", confusion.balanced_accuracy()),
    ] {
        if let Some(v) = value {
            pooled.push((name.to_string(), v));
        }
    }

    let per_protein_auroc = if per_protein.is_empty() {
        None
    } else {
        Some(per_protein.iter().sum::<f64>() / per_protein.len() as f64)
    };

    Ok(BenchReport {
        suite: "caid3-disorder-pdb".into(),
        embedding_model: embedder.repo().into(),
        head_version: metadata.version,
        with_burial: false,
        proteins: reference.len(),
        residues_scored: scores.len(),
        residues_masked: masked,
        positive_rate: Some(positives as f64 / scores.len() as f64),
        windowed,
        pooled,
        per_protein_auroc,
        per_protein_undefined: undefined,
        caveats: vec![
            format!(
                "Residues the reference did not assess ('-') are excluded, which is the CAID \
                 protocol. {masked} of {} residues were excluded.",
                masked + scores.len()
            ),
            format!(
                "{windowed} of {} proteins are longer than the embedding model's {limit}-residue \
                 limit and were covered by overlapping windows. No published CAID3 entrant uses \
                 this windowing, so those proteins are not scored the way another method would \
                 score them.",
                reference.len()
            ),
            format!(
                "The threshold for MCC, F1 and balanced accuracy is {threshold}, taken from \
                 the head rather than tuned on this set: tuning it here would make the \
                 number incomparable to a published one and would be fitting to the test \
                 set. It therefore comes from a different distribution than this one -- a \
                 calibrated head's threshold is fitted on its own validation split, whose \
                 class balance and definition of disorder are not this reference set's -- so \
                 a threshold chosen here would score better. The rank metrics above are \
                 unaffected by any of this, which is why AUROC and AUPRC are the figures to \
                 compare."
            ),
            "CAID3 rankings are computed by the organisers' own code on a server-side \
             ground truth. This is our reimplementation of the protocol against the published \
             reference file, so a small difference from an official number is expected and this \
             figure should not be presented as an official CAID3 result."
                .into(),
        ],
    })
}

/// One deep-mutational-scanning assay from `ProteinGym`.
#[derive(Debug, Clone)]
pub struct Assay {
    /// `ProteinGym` assay identifier.
    pub id: String,
    /// What the assay measured, from the reference file. Worth carrying
    /// because "thermostability" and "protein abundance" are different
    /// quantities and neither is ΔΔG.
    pub selection: String,
    /// Wild-type sequence the mutants are numbered against.
    pub sequence: String,
    /// `(0-based position, mutant residue, measured score)`. Higher scores are
    /// fitter in `ProteinGym`'s convention, so for a stability-like assay a
    /// correctly ordered ΔΔG prediction correlates **positively**.
    pub mutants: Vec<(usize, char, f64)>,
}

/// Parse a `ProteinGym` mutant label such as `A32N` against its wild-type
/// sequence.
///
/// Returns the 0-based position and the mutant residue.
///
/// # Errors
/// If the label is malformed, the position is outside the sequence, or the
/// wild-type letter in the label disagrees with the sequence at that
/// position. The last of those is the valuable one: it is what catches an
/// off-by-one in the numbering convention, which would otherwise produce a
/// benchmark number that looks reasonable and scores the wrong residues.
pub fn parse_mutant(label: &str, sequence: &[char]) -> Result<(usize, char)> {
    let mut chars = label.chars();
    let wild = chars.next().ok_or_else(|| anyhow!("empty mutant label"))?;
    let mutant = label.chars().next_back().ok_or_else(|| anyhow!("empty mutant label"))?;
    let digits: String = chars.clone().take(label.chars().count() - 2).collect();
    let position: usize = digits.parse().with_context(|| format!("{label}: {digits:?} is not a position"))?;
    if position == 0 || position > sequence.len() {
        bail!("{label}: position {position} is outside a {}-residue sequence", sequence.len());
    }
    let index = position - 1;
    if sequence[index] != wild {
        bail!(
            "{label}: the label says the wild-type residue at position {position} is {wild}, but \
             the sequence has {}. The numbering conventions disagree, so every score from this \
             assay would be attributed to the wrong residue.",
            sequence[index]
        );
    }
    Ok((index, mutant))
}

/// Predict ΔΔG for every single substitution in one assay.
///
/// Returns the predictions and the measurements they line up with, in the
/// same order, so a correlation cannot be computed across mismatched vectors.
/// Substitutions involving a non-standard residue are dropped from both.
fn score_one_assay(
    assay: &Assay,
    head: &opendna_ml::heads::StabilityHead,
    burial: Option<&opendna_ml::heads::BurialHead>,
    model: &Esm2,
    device: &candle_core::Device,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let residues: Vec<char> = assay.sequence.chars().collect();
    let embedding = model.embed(&[assay.sequence.as_str()])?.remove(0);

    // One masked forward pass per mutated position, reused across every
    // substitution at that position — otherwise this would be one pass per
    // mutant, roughly nineteen times the work.
    let mut positions: Vec<usize> = assay.mutants.iter().map(|(i, _, _)| *i).collect();
    positions.sort_unstable();
    positions.dedup();
    println!("  {}: {} mutants at {} positions", assay.id, assay.mutants.len(), positions.len());
    let log_probs = model.masked_log_probs(&assay.sequence, &positions)?;
    let index_of: std::collections::HashMap<usize, usize> =
        positions.iter().enumerate().map(|(n, p)| (*p, n)).collect();

    // Accessibility for the whole sequence, from the same burial head
    // inference would use, when this variant takes it.
    let exposure = match burial {
        Some(b) => Some(b.predict(&embedding.per_residue, embedding.length, device)?),
        None => None,
    };

    let standard = opendna_core::AminoAcid::STANDARD;
    let mut rows = Vec::new();
    let mut measured = Vec::new();
    for (position, mutant, score) in &assay.mutants {
        let wild = residues[*position];
        let (Some(wt_index), Some(mut_index)) = (
            standard.iter().position(|a| a.one_letter() == wild),
            standard.iter().position(|a| a.one_letter() == *mutant),
        ) else {
            // A non-standard residue on either side has no one-hot encoding
            // the head was trained with.
            continue;
        };
        let lp = &log_probs[index_of[position]];
        let logp = |letter: char| -> f32 { lp[opendna_ml::esm2::alphabet::token_id(letter) as usize] };
        let site = embedding
            .residue(*position)
            .ok_or_else(|| anyhow!("{}: no embedding at {position}", assay.id))?;
        let at_position = match &exposure {
            Some(values) => Some(
                *values
                    .get(*position)
                    .ok_or_else(|| anyhow!("{}: no burial prediction at {position}", assay.id))?,
            ),
            None => None,
        };
        rows.push(opendna_ml::heads::StabilityHead::features_with(
            site,
            wt_index,
            mut_index,
            logp(*mutant) - logp(wild),
            at_position,
        ));
        measured.push(*score);
    }
    let predicted = head.predict(&rows, model.hidden_size(), device)?;
    Ok((predicted.into_iter().map(f64::from).collect(), measured))
}

/// Score the stability head against `ProteinGym` assays.
///
/// Each assay's Spearman correlation is reported separately, because the
/// assays measure different quantities in different proteins and a single
/// pooled number across them would not mean anything.
///
/// # Errors
/// If a head or the embedding model cannot be loaded, if an assay's numbering
/// does not match its sequence, or if inference fails.
pub fn stability(assays: &[Assay], models_dir: &Path, embedder: ModelSize) -> Result<BenchReport> {
    let device = opendna_ml::best_device()?;
    println!("Loading on {}", opendna_ml::device_description(&device));
    let heads = Heads::load(models_dir, embedder, &device)
        .with_context(|| format!("loading heads from {}", models_dir.display()))?;
    let head =
        heads.stability.as_ref().ok_or_else(|| anyhow!("no stability head in {}", models_dir.display()))?;
    let metadata = read_metadata(&HeadFiles::in_dir(models_dir, "stability"), "stability", embedder)?;
    // The head's metadata says whether it takes predicted burial. Supplying
    // it when it was not trained with it, or omitting it when it was, is a
    // width mismatch — which is how this was caught rather than quietly
    // scoring a different model than the one on disk.
    let with_burial = metadata.with_burial;
    if with_burial && heads.burial.is_none() {
        bail!(
            "the stability head in {} was trained with burial, but there is no burial head there \
             to supply it",
            models_dir.display()
        );
    }
    let model = Esm2::from_pretrained(embedder, &device)?;

    let mut pooled = Vec::new();
    let mut scored = 0usize;
    let mut skipped_assays = Vec::new();
    for assay in assays {
        if assay.sequence.chars().count() > model.max_residues() {
            // Windowing changes which context each residue sees, and the
            // stability head was trained on single-pass embeddings. Rather
            // than quietly score it differently, drop it and say so.
            skipped_assays.push(format!(
                "{} ({} residues, over the model's {} limit)",
                assay.id,
                assay.sequence.chars().count(),
                model.max_residues()
            ));
            continue;
        }
        let (predicted, measured) = score_one_assay(
            assay,
            head,
            if with_burial { heads.burial.as_ref() } else { None },
            &model,
            &device,
        )?;
        scored += predicted.len();
        if let Some(rho) = metrics::spearman(&predicted, &measured) {
            pooled.push((format!("spearman::{}", assay.id), rho));
        }
    }

    if pooled.is_empty() {
        bail!("no assay produced a correlation; nothing was scored");
    }
    let mean = pooled.iter().map(|(_, v)| *v).sum::<f64>() / pooled.len() as f64;
    pooled.push(("spearman::mean_over_assays".into(), mean));

    let mut caveats = vec![
        "These assays measure variant abundance or thermostability in whole proteins, not ΔΔG \
         by calorimetry. They are correlates of folding stability, so a positive Spearman means \
         the head orders destabilising substitutions correctly; the magnitude is not comparable \
         to a ΔΔG-against-ΔΔG correlation."
            .into(),
        "Every ProteinGym assay derived from Tsuboyama 2023 is excluded. Those 64 assays are the \
         head's own training data, and scoring them would measure memorisation."
            .into(),
        "ProteinGym scores are higher-is-fitter, and the head predicts higher-is-more-stable, so \
         the expected sign is positive. A negative correlation here would mean the head is \
         ordering mutations backwards on this protein, not that the convention is flipped."
            .into(),
        format!(
            "{} assays, {scored} mutants scored. This is a handful of proteins, so the mean \
             across assays has no useful confidence interval; read the per-assay numbers.",
            pooled.len() - 1
        ),
    ];
    if !skipped_assays.is_empty() {
        caveats.push(format!("Skipped: {}.", skipped_assays.join(", ")));
    }
    // What each assay actually measured. Without this the per-assay numbers
    // read as five goes at one quantity, when they are different quantities
    // that happen to correlate with folding stability.
    for assay in assays {
        caveats.push(format!(
            "{}: {} in a {}-residue protein, {} single substitutions.",
            assay.id,
            assay.selection,
            assay.sequence.chars().count(),
            assay.mutants.len()
        ));
    }

    Ok(BenchReport {
        suite: "proteingym-stability-out-of-domain".into(),
        embedding_model: embedder.repo().into(),
        head_version: metadata.version,
        with_burial,
        proteins: assays.len() - skipped_assays.len(),
        residues_scored: scored,
        residues_masked: 0,
        positive_rate: None,
        windowed: 0,
        pooled,
        per_protein_auroc: None,
        per_protein_undefined: 0,
        caveats,
    })
}

/// Which head a suite scores.
///
/// A report records the suite it ran, not the head, because a suite is the
/// fixed thing: the head is whatever was in `--models`. The association is
/// here, spelled out, rather than inferred from a filename.
#[must_use]
pub fn task_of(suite: &str) -> Option<&'static str> {
    match suite {
        "caid3-disorder-pdb" => Some("disorder"),
        "proteingym-stability-out-of-domain" => Some("stability"),
        _ => None,
    }
}

/// Load the benchmark reports in `dir` that describe `task` scored with
/// `embedding_model`.
///
/// The embedding model has to match, or a card would quote a number measured
/// on a different model than the one it documents. A directory that does not
/// exist yields no reports rather than an error: benchmarks are optional, and
/// a card without them is incomplete rather than wrong.
///
/// # Errors
/// If a file in `dir` ends in `.json` but is not a readable report. A
/// malformed report is a mistake worth stopping for, since the alternative is
/// silently publishing a card with a benchmark missing.
pub fn load_reports(
    dir: &Path,
    task: &str,
    embedding_model: &str,
    with_burial: bool,
) -> Result<Vec<BenchReport>> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Ok(Vec::new());
    };
    let mut out = Vec::new();
    for entry in entries {
        let path = entry.context("listing the benchmarks directory")?.path();
        if path.extension().is_none_or(|e| e != "json") {
            continue;
        }
        let text = std::fs::read_to_string(&path).with_context(|| format!("reading {}", path.display()))?;
        let report: BenchReport = serde_json::from_str(&text)
            .with_context(|| format!("parsing {} as a benchmark report", path.display()))?;
        // All three must match. A report from a variant that was measured
        // and not published is a measurement of a different model, and putting
        // it on this model's card would be describing something nobody can
        // download.
        if task_of(&report.suite) == Some(task)
            && report.embedding_model == embedding_model
            && report.with_burial == with_burial
        {
            out.push(report);
        }
    }
    out.sort_by(|a, b| a.suite.cmp(&b.suite));
    Ok(out)
}

/// Render a report as a short markdown section.
pub fn render(report: &BenchReport) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "## {}\n", report.suite);
    let balance = match report.positive_rate {
        Some(rate) => format!(", {:.1}% positive", rate * 100.0),
        None => String::new(),
    };
    let _ = writeln!(
        out,
        "{} proteins, {} measurements scored ({} masked){balance}.\n",
        report.proteins, report.residues_scored, report.residues_masked
    );
    let _ = writeln!(out, "| Metric | Value |\n|---|---:|");
    for (name, value) in &report.pooled {
        let _ = writeln!(out, "| {name} | {value:.4} |");
    }
    if let Some(mean) = report.per_protein_auroc {
        let _ = writeln!(
            out,
            "| auroc (mean over {} proteins) | {mean:.4} |",
            report.proteins - report.per_protein_undefined
        );
    }
    let _ = writeln!(out, "\nHead `{}` on `{}`.\n", report.head_version, report.embedding_model);
    for caveat in &report.caveats {
        let _ = writeln!(out, "- {caveat}");
    }
    out
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    const THREE: &str = ">A\nMKV\n10-\n>B\nQQ\n00\n";

    #[test]
    fn reads_sequences_labels_and_the_mask() {
        let r = read_caid_reference(THREE).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].id, "A");
        assert_eq!(r[0].sequence, "MKV");
        assert_eq!(r[0].labels, vec![Some(true), Some(false), None]);
        assert_eq!(r[1].labels, vec![Some(false), Some(false)]);
    }

    /// A label string of the wrong length would shift every label after the
    /// difference onto a neighbouring residue, producing a plausible-looking
    /// and completely wrong benchmark number.
    #[test]
    fn a_label_string_of_the_wrong_length_is_rejected() {
        let err = read_caid_reference(">A\nMKV\n10\n").unwrap_err().to_string();
        assert!(err.contains("3 residues but 2 labels"), "{err}");
    }

    #[test]
    fn an_unknown_label_character_is_rejected() {
        let err = read_caid_reference(">A\nMKV\n1X0\n").unwrap_err().to_string();
        assert!(err.contains("unrecognised label character"), "{err}");
    }

    #[test]
    fn a_truncated_record_is_rejected() {
        let err = read_caid_reference(">A\nMKV\n").unwrap_err().to_string();
        assert!(err.contains("truncated"), "{err}");
    }

    #[test]
    fn a_missing_header_is_rejected() {
        let err = read_caid_reference("A\nMKV\n10-\n").unwrap_err().to_string();
        assert!(err.contains("expected a '>' header"), "{err}");
    }

    #[test]
    fn an_empty_file_is_rejected() {
        assert!(read_caid_reference("").unwrap_err().to_string().contains("no records"));
    }

    /// Blank lines between records are tolerated, because a reference file
    /// written by someone else is not ours to reformat.
    #[test]
    fn blank_lines_between_records_are_ignored() {
        let r = read_caid_reference(">A\nMKV\n10-\n\n\n>B\nQQ\n00\n").unwrap();
        assert_eq!(r.len(), 2);
    }

    #[test]
    fn a_mutant_label_resolves_to_a_zero_based_position() {
        let seq: Vec<char> = "MKVLA".chars().collect();
        assert_eq!(parse_mutant("M1A", &seq).unwrap(), (0, 'A'));
        assert_eq!(parse_mutant("A5P", &seq).unwrap(), (4, 'P'));
        // Two-digit positions must not be read one character at a time.
        let long: Vec<char> = "ACDEFGHIKLMNP".chars().collect();
        assert_eq!(parse_mutant("P13W", &long).unwrap(), (12, 'W'));
    }

    /// The check that earns its keep: if the benchmark's numbering convention
    /// ever disagrees with ours, every score is attributed to the wrong
    /// residue and the resulting correlation looks plausible. A 0-based read
    /// of `K2V` would land on V, and this is what refuses it.
    #[test]
    fn a_wild_type_mismatch_is_refused_rather_than_scored() {
        let seq: Vec<char> = "MKVLA".chars().collect();
        let err = parse_mutant("V2L", &seq).unwrap_err().to_string();
        assert!(err.contains("the sequence has K"), "{err}");
        assert!(err.contains("wrong residue"), "{err}");
    }

    #[test]
    fn a_position_outside_the_sequence_is_refused() {
        let seq: Vec<char> = "MKV".chars().collect();
        assert!(parse_mutant("M9A", &seq).unwrap_err().to_string().contains("outside"));
        assert!(parse_mutant("M0A", &seq).unwrap_err().to_string().contains("outside"));
    }
}
