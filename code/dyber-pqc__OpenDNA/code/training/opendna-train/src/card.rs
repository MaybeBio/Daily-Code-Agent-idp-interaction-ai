//! Generates a Hugging Face model card from a training report.
//!
//! The card is assembled from the metrics the run actually produced, so the
//! numbers in it cannot drift from the numbers that were measured. The parts
//! that need judgement — what the model is for, where it fails — are written
//! here as prose, per task, and reviewed by hand before publishing.

use crate::train::Report;
use crate::Task;
use anyhow::{Context, Result};
use std::fmt::Write;
use std::path::Path;

/// Text that cannot be derived from the metrics.
struct Prose {
    summary: &'static str,
    dataset: &'static str,
    dataset_url: &'static str,
    license: &'static str,
    /// What the *weights* are released under, as a Hugging Face licence
    /// identifier, with `weights_license_note` explaining the choice.
    ///
    /// This is not always `apache-2.0`. Whether model weights are a derivative
    /// work of their training data is unsettled, and a head trained on
    /// share-alike data is not something this project can relicense on the
    /// strength of an assumption. Where the upstream licence is share-alike,
    /// the weights carry it forward.
    weights_license: &'static str,
    weights_license_note: &'static str,
    inputs: &'static str,
    outputs: &'static str,
    baseline: &'static str,
    limitations: &'static [&'static str],
    citation: &'static str,
}

// One `match` of per-task prose. Splitting it per task would hide the fact
// that the three cards are the same document with different words, which is
// the property that keeps them comparable.
#[allow(clippy::too_many_lines)]
fn prose(task: Task) -> Prose {
    match task {
        Task::Stability => Prose {
            summary: "Predicts the change in folding free energy (ΔΔG) caused by a single \
                      amino-acid substitution.",
            dataset: "dyber-pqc/opendna-stability",
            dataset_url: "https://huggingface.co/datasets/dyber-pqc/opendna-stability",
            license: "CC BY 4.0 (Tsuboyama et al. 2023)",
            weights_license: "apache-2.0",
            weights_license_note: "The training data is CC BY 4.0, which asks for attribution \
                                   and places no condition on the licence of derived work, so \
                                   the weights are Apache-2.0. Cite the dataset.",
            inputs: "The wild-type sequence's ESM-2 embedding at the mutated position, one-hot \
                     encodings of the wild-type and mutant residues, and the masked \
                     log-odds `log P(mutant) - log P(wild type)` at that position.",
            outputs: "A single ΔΔG value in kcal/mol, following the dataset's sign convention: \
                      negative is destabilising.",
            baseline: "The language model's masked log-odds used directly, with no trained head.",
            limitations: &[
                "Trained almost entirely on small, single-domain proteins of 32-72 residues. \
                 Behaviour on large multi-domain proteins, membrane proteins and complexes is \
                 untested.",
                "The source assay infers ΔΔG from protease resistance over a range of roughly \
                 -1 to 5 kcal/mol for ΔG. Mutations outside that range are absent from training, \
                 so strongly destabilising substitutions are under-represented.",
                "Only single substitutions. The model says nothing about combined mutations, \
                 insertions or deletions.",
                "Trained on frozen embeddings, so accuracy is bounded by what ESM-2 already \
                 encodes about the position.",
                "Trained with a squared-error loss, so predictions are compressed toward the \
                 training mean and span a much narrower range than the measurements do. The rank \
                 order is the usable signal, which is why Spearman is the headline metric; treat \
                 an individual predicted value as a weak estimate.",
            ],
            citation: "Tsuboyama, K. et al. Mega-scale experimental analysis of protein folding \
                       stability in biology and design. Nature 620, 434-444 (2023).",
        },
        Task::Solubility => Prose {
            summary: "Predicts protein solubility from sequence.",
            dataset: "dyber-pqc/opendna-solubility",
            dataset_url: "https://huggingface.co/datasets/dyber-pqc/opendna-solubility",
            license: "CC BY-SA 2.1 JP (eSOL, Niwa et al. 2009); the UniProtKB sequences it is \
                      joined against are CC BY 4.0",
            weights_license: "other\nlicense_name: cc-by-sa-2.1-jp\nlicense_link: \
                              https://creativecommons.org/licenses/by-sa/2.1/jp/deed.en",
            weights_license_note: "The eSOL measurements are licensed CC BY-SA 2.1 JP, which is \
                                   share-alike. Whether trained weights are a derivative work of \
                                   their training data is an open legal question, and this \
                                   project is not going to answer it in its own favour, so these \
                                   weights are released under the same share-alike terms rather \
                                   than the Apache-2.0 used for the other heads.",
            inputs: "The mean-pooled ESM-2 embedding of the sequence.",
            outputs: "A solubility percentage on the source assay's scale. Values above 100 occur \
                      in the source data and are not clipped.",
            baseline: "Predicting the training-set mean for every protein.",
            limitations: &[
                "Trained on E. coli cell-free expression measurements. Solubility in other hosts, \
                 or under different expression conditions, is a different quantity.",
                "Sequences are UniProtKB canonical sequences, not the exact expressed constructs; \
                 tags and truncations in the real experiments are not represented.",
                "A few thousand proteins is a small training set. Treat the output as a ranking \
                 signal for prioritising constructs, not as a predicted assay value.",
            ],
            citation: "Niwa, T. et al. Bimodal protein solubility distribution revealed by an \
                       aggregation analysis of the entire ensemble of Escherichia coli proteins. \
                       PNAS 106, 4201-4206 (2009).",
        },
        Task::Burial => Prose {
            summary: "Predicts, for each residue, how solvent-exposed it is.",
            dataset: "dyber-pqc/opendna-burial",
            dataset_url: "https://huggingface.co/datasets/dyber-pqc/opendna-burial",
            license: "CC BY 4.0 (CATH v4.4.0); the coordinates it is computed from are CC0 1.0 \
                      (wwPDB)",
            weights_license: "apache-2.0",
            weights_license_note: "CATH is CC BY 4.0 and wwPDB coordinates are CC0, neither of \
                                   which places a condition on the licence of derived work, so \
                                   the weights are Apache-2.0. Cite the dataset.",
            inputs: "Per-residue ESM-2 embeddings. Sequences longer than the embedding model's \
                     limit are covered by overlapping windows.",
            outputs: "Relative solvent accessibility per residue: the residue's accessible area \
                      over the theoretical maximum for that amino acid, so roughly 0 in the core \
                      and roughly 1 fully exposed. The output is squashed into [0, 1.25], \
                      because a negative area is not a thing and the theoretical maxima are \
                      occasionally exceeded.",
            baseline: "Predicting the training-set mean accessibility for every residue.",
            limitations: &[
                "Labels come from a single chain of a crystal or NMR model, so a residue that \
                 forms an inter-chain interface in an obligate oligomer is labelled exposed. The \
                 head learns that convention and will do the same.",
                "Trained on single-segment CATH domains whose sequence could be located in \
                 current coordinates, which excludes discontinuous domains — about one in eight \
                 of CATH — and entries substantially revised since CATH v4.4.0. The training \
                 distribution is biased toward simple, stable structures.",
                "Accessibility is a property of a fold. For a sequence with no stable structure \
                 the number is still produced and still means nothing; read it alongside the \
                 disorder head rather than alone.",
                "The metrics are pooled over residues, and residues within a domain are strongly \
                 correlated because neighbours share a local environment. The effective sample \
                 size is nearer the number of domains than the number of residues, so small \
                 differences between splits are noise.",
                "Trained on frozen embeddings, so accuracy is bounded by what ESM-2 already \
                 encodes about the position.",
            ],
            // Two references, and the card has to say which is which: the
            // labels are computed from CATH domains and wwPDB coordinates,
            // and Tien supplies the per-residue maxima they are divided by.
            // Crediting Tien as the source of the data was an attribution
            // error.
            citation: "Sillitoe, I. et al. CATH: increased structural coverage of functional \
                       space. Nucleic Acids Research 49, D266-D273 (2021).\n>\n> The labels are \
                       computed from wwPDB coordinates and normalised by the theoretical maxima \
                       of: Tien, M. Z. et al. Maximum allowed solvent accessibilites of residues \
                       in proteins. PLOS ONE 8, e80635 (2013).",
        },
        Task::Disorder => Prose {
            summary: "Predicts, for each residue, whether it lies in an intrinsically disordered \
                      region.",
            dataset: "dyber-pqc/opendna-disorder",
            dataset_url: "https://huggingface.co/datasets/dyber-pqc/opendna-disorder",
            license: "CC BY 4.0 (DisProt 2025_12); the CAID3 reference set used only to hold \
                      proteins out is CC BY 3.0",
            weights_license: "apache-2.0",
            weights_license_note: "Both sources ask for attribution and place no condition on \
                                   the licence of derived work, so the weights are Apache-2.0. \
                                   Cite the dataset.",
            inputs: "Per-residue ESM-2 embeddings. Sequences longer than the embedding model's \
                     limit are covered by overlapping windows.",
            outputs: "One calibrated probability per residue. The threshold for a disorder \
                      call is recorded in the head's metadata and stated under Calibration \
                      below, because it is a property of the trained head rather than a \
                      convention.",
            baseline: "Labelling every residue disordered.",
            limitations: &[
                "DisProt records disorder that has been experimentally demonstrated. Residues \
                 that are simply unannotated are treated as ordered during training, so the \
                 negative class contains an unknown number of undiscovered disordered residues.",
                "Annotation is biased toward proteins people chose to study, and toward longer \
                 disordered regions. Short linkers are under-represented.",
                "The CAID3 reference proteins are held out of training, and their whole \
                 identity clusters with them, which is what makes the CAID3 figure above a \
                 measurement rather than a leak. The results table beside it is this dataset's \
                 own test split and is not comparable to anybody's published number: our \
                 negative class contains undiscovered disorder and CAID masks what it did not \
                 assess, so the two are scoring different things.",
                "The metrics are pooled over residues, but residues within one protein are \
                 strongly correlated — neighbouring positions almost always share a label. The \
                 effective sample size is therefore closer to the few hundred proteins in each \
                 split than to the hundred thousand or more residues the counts show, so small \
                 differences between splits should not be read as real. On this run the test \
                 split scored higher than validation on every metric, which is what that noise \
                 looks like rather than evidence the model generalises unusually well.",
            ],
            citation: "Aspromonte, M. C. et al. DisProt in 2024: improving function annotation of \
                       intrinsically disordered proteins. Nucleic Acids Research 52, D434-D441 (2024).",
        },
    }
}

/// Display name for a metric key, so the card reads as prose rather than as
/// field names.
fn label(key: &str) -> &str {
    match key {
        "spearman" => "Spearman ρ",
        "pearson" => "Pearson r",
        "rmse" => "RMSE",
        "mae" => "MAE",
        "auroc" => "AUROC",
        "auprc" => "AUPRC",
        "mcc" => "MCC",
        "f1" => "F1",
        "balanced_accuracy" => "Balanced accuracy",
        other => other,
    }
}

fn metrics_table(report: &Report) -> String {
    let mut keys: Vec<&String> = report.test.values.keys().collect();
    for k in report.validation.values.keys() {
        if !keys.contains(&k) {
            keys.push(k);
        }
    }
    keys.sort();

    let mut out = String::from("| Metric | Validation | Test | Baseline (test) |\n|---|---|---|---|\n");
    // Which metrics the baseline has no value for, because the reason differs
    // and the two reasons are near opposites.
    let mut missing_correlation = false;
    let mut missing_error = false;
    for key in keys {
        let cell = |m: &crate::train::SplitMetrics| {
            m.values.get(key).map_or_else(|| "—".to_string(), |v| format!("{v:.3}"))
        };
        if !report.baseline.values.contains_key(key) {
            if matches!(key.as_str(), "pearson" | "spearman") {
                missing_correlation = true;
            } else {
                missing_error = true;
            }
        }
        let _ = writeln!(
            out,
            "| {} | {} | {} | {} |",
            label(key),
            cell(&report.validation),
            cell(&report.test),
            cell(&report.baseline)
        );
    }
    if missing_correlation || missing_error {
        out.push_str(
            "\nA — in the baseline column means the metric is undefined for that baseline, not \
             unmeasured. ",
        );
    }
    // A constant baseline — the training mean, or the majority class — has no
    // correlation to report, because correlation with a constant is undefined
    // rather than zero.
    if missing_correlation {
        out.push_str(
            "This baseline predicts the same value for every input, so it has no correlation \
             with the truth. ",
        );
    }
    // A baseline that is a real predictor on another scale is the opposite
    // case: it correlates, and an absolute error against the measurements is
    // not a thing it has.
    if missing_error {
        out.push_str(
            "This baseline is a score rather than a prediction in the same units as the \
             measurements, so an absolute or squared error against them is not defined for it \
             — only how well it orders them. ",
        );
    }
    if missing_correlation || missing_error {
        // Each sentence above ends in a space so they join on one line; the
        // last one must not leave trailing whitespace in a published file.
        while out.ends_with(' ') {
            out.pop();
        }
        out.push('\n');
    }
    out
}

/// The external-benchmark section, or an explanation of why there is none.
///
/// The table above it is measured on our own test split, which says nothing
/// about how the head compares to published work, because nobody else reports
/// a number on our splits. A card that gives only that number invites the
/// reader to treat it as a standing rather than a measurement, so the absence
/// of an external figure is stated rather than left as a gap.
fn external_section(reports: &[crate::bench::BenchReport]) -> String {
    let mut out = String::from("### Against published benchmarks\n\n");
    if reports.is_empty() {
        out.push_str(
            "None. The external benchmarks for this task come from sources whose terms do not \
             permit redistribution, or from services that have been retired, so there is no \
             set this head can be scored on alongside other published methods. The figures \
             above are therefore not evidence of where it stands relative to the field — only \
             of how it does on data drawn the way its training data was drawn.\n\n",
        );
        return out;
    }
    for report in reports {
        let _ = writeln!(out, "**{}** — {} items scored.\n", report.suite, report.residues_scored);
        let _ = writeln!(out, "| Metric | Value |\n|---|---:|");
        for (name, value) in &report.pooled {
            let _ = writeln!(out, "| {name} | {value:.4} |");
        }
        if let Some(mean) = report.per_protein_auroc {
            let _ = writeln!(out, "| auroc (mean over proteins) | {mean:.4} |");
        }
        out.push('\n');
        for caveat in &report.caveats {
            let _ = writeln!(out, "- {caveat}");
        }
        out.push('\n');
    }
    out
}

/// Say outright whether the intervals cover less than they promise.
///
/// A reader can compare two columns, and a reader skimming will not. Burial
/// falls short at every level and that is the fact most worth not missing on
/// its card, so it is stated in words rather than left in the table.
fn coverage_direction(intervals: &opendna_ml::heads::Intervals) -> String {
    let short = |l: &opendna_ml::heads::IntervalLevel| l.achieved < l.nominal;
    if intervals.levels.is_empty() {
        return String::new();
    }
    if intervals.levels.iter().all(short) {
        return String::from(
            "On this head it points the same way at every level: the intervals cover less than \
             they say, so read each one as the lower bound of a range that ought to be wider.",
        );
    }
    if !intervals.levels.iter().any(short) {
        return String::from(
            "On this head it points the other way at every level: the intervals cover more than \
             they say, so they are conservative rather than optimistic.",
        );
    }
    String::from("On this head it falls short at some levels and not others; read the column.")
}

/// What the interval's quantile actually rests on, which is not the residual
/// count.
///
/// Different for each regressor, and the difference is not cosmetic. Burial
/// and stability produce many residuals per protein and those residuals are
/// correlated, so the effective count is nearer the number of proteins.
/// Solubility produces exactly one per protein, so its residuals are
/// independent and the only caveat is how few there are — stating a
/// correlation there would be inventing a problem, and stating it in
/// stability's vocabulary on the burial card, which has no wild types and no
/// substitutions, was how this paragraph read when it was written once for all
/// three.
fn independence_note(task: Task) -> String {
    match task {
        Task::Stability => String::from(
            "It rests on fewer observations than the residual count suggests. Residuals from \
             one wild-type protein are correlated — substitutions in the same fold share \
             whatever the head has got wrong about it — so the effective number of independent \
             observations behind each quantile is nearer the number of distinct wild types in \
             the split than the number of substitutions. The dataset card records both.",
        ),
        Task::Burial => String::from(
            "It rests on fewer observations than the residual count suggests. Residuals from \
             one domain are correlated — neighbouring residues share a local environment and \
             whatever the head has got wrong about that fold — so the effective number of \
             independent observations behind each quantile is nearer the number of domains in \
             the split than the number of residues. The dataset card records both.",
        ),
        Task::Solubility => String::from(
            "Each protein contributes exactly one residual, so unlike the per-residue heads \
             these observations are independent. There are simply few of them: the split the \
             quantile is fitted on holds a few hundred proteins, which is why the interval is \
             wide and why it should be read as the honest width rather than a tight one.",
        ),
        // No intervals on a classifier, so this is never asked for.
        Task::Disorder => String::new(),
    }
}

/// Which threshold the table's threshold-dependent metrics used.
///
/// Only for the classifier, and only because it has two different Matthews
/// correlations on one card. The training report is written before any
/// calibration exists, so its F1, balanced accuracy and MCC are computed at
/// the conventional 0.5 on raw scores. The figure under Calibration is at the
/// head's own threshold on calibrated ones. Both are honest and they are not
/// the same measurement, and a reader comparing either to a published number
/// has to know which they are holding.
fn threshold_note(task: Task) -> String {
    match task {
        Task::Disorder => String::from(
            "F1, balanced accuracy and MCC above are computed at a probability of 0.5 on the \
             head's raw output, which is how the training report records them. The head's own \
             decision threshold is lower and is stated under Calibration, along with the \
             Matthews correlation it achieves — a different and better number, measured a \
             different way. AUROC and AUPRC are rank statistics and do not depend on either \
             threshold.\n\n",
        ),
        // The regressors' metrics are MAE, RMSE and correlations. No threshold
        // enters them, so a note about one would be noise.
        Task::Burial | Task::Solubility | Task::Stability => String::new(),
    }
}

/// The Hub metadata that decides whether anybody finds the repository.
///
/// Returns the `pipeline_tag` line, empty when no tag is honest, followed by
/// the task's own tags.
///
/// Hugging Face recognises a fixed list of pipeline tags and puts a repository
/// under a task filter only if it carries one. One of these four heads fits a
/// tag on that list: disorder assigns a label to each residue, which is token
/// classification. The others are regressions — accessibility per residue,
/// solubility per protein, stability per substitution — and the list has no
/// entry for a continuous output over a sequence. `tabular-regression`
/// describes feature tables, and `text-classification` would advertise classes
/// these heads do not have. A wrong tag also drives a widget that would submit
/// a sequence and show nothing useful.
///
/// Free-form tags carry the rest. They are what text search matches on, and
/// they can say "regression" without claiming a pipeline.
fn hub_metadata(task: Task) -> String {
    let (pipeline, tags) = match task {
        Task::Disorder => (
            Some("token-classification"),
            &["intrinsic-disorder", "idp", "per-residue", "token-classification"][..],
        ),
        Task::Burial => (None, &["solvent-accessibility", "burial", "rsa", "per-residue", "regression"][..]),
        Task::Solubility => (None, &["solubility", "expression", "regression"][..]),
        Task::Stability => {
            (None, &["stability", "ddg", "mutation-effect", "variant-effect", "regression"][..])
        }
    };
    let mut out = String::new();
    if let Some(pipeline) = pipeline {
        let _ = writeln!(out, "pipeline_tag: {pipeline}");
    }
    out.push_str("tags:\n  - protein\n  - biology\n  - protein-language-model\n");
    for tag in tags {
        let _ = writeln!(out, "  - {tag}");
    }
    out
}

/// How to actually run these weights, in the card itself.
///
/// Hugging Face renders a "Use this model" snippet for libraries it knows how
/// to load. `library_name: opendna` is not one, which is accurate -- these are
/// not transformers weights -- and leaves a reader with metrics and no way to
/// act on them. So the card carries the commands.
///
/// The prediction line differs by task, and the difference is not cosmetic:
/// stability scores named substitutions and cannot do anything without being
/// told which ones.
fn usage_section(task: Task) -> String {
    let name = task.name();
    let predict = match task {
        // A mutation list is required, so a card that showed the bare command
        // would be showing one that exits 2.
        Task::Stability => format!(
            "# a CSV of `name,mutation` rows, for example `ubiquitin,L8P`\n\
             opendna predict --task {name} -i candidates.fasta \\\n\
             \u{20}   --mutations muts.csv -o {name}.csv"
        ),
        Task::Disorder | Task::Burial | Task::Solubility => {
            format!("opendna predict --task {name} -i candidates.fasta -o {name}.csv")
        }
    };
    let caid =
        if matches!(task, Task::Disorder) { "\nFor a CAID submission, add `--format caid`.\n" } else { "" };
    format!(
        "## Using it\n\n\
         These are not transformers weights and no `AutoModel` will load them. They are read by \
         [OpenDNA](https://github.com/dyber-pqc/OpenDNA), which is a single binary with no \
         runtime to install — take it from the \
         [releases page](https://github.com/dyber-pqc/OpenDNA/releases) or build it with \
         `cargo build --release -p opendna-cli`.\n\n\
         ```bash\n\
         hf download dyber-pqc/opendna-{name} --local-dir models\n\
         export OPENDNA_MODELS=$PWD/models\n\n\
         # confirm the download is intact: the weights are hashed against the\n\
         # digest recorded in their metadata\n\
         opendna models\n\n\
         {predict}\n\
         ```\n\n\
         The output carries the head version, the embedding model, the device and the \
         limitations below as `#` comments, and a `{name}.manifest.json` beside it records the \
         digests of what went in and what came out. The first run downloads the embedding model; \
         after that `opendna stage` makes a directory that works with no network at all.\n\n\
         `opendna models` checks that the file is the one it claims to be. `opendna selftest` \
         goes further and re-runs the heads over a known sequence to check they reproduce \
         recorded values — worth doing, and it needs all four heads installed.\n\
         {caid}\n"
    )
}

/// Render a model card.
pub fn render(task: Task, report: &Report, external: &[crate::bench::BenchReport]) -> String {
    let p = prose(task);
    let mut out = String::new();
    let _ = write!(
        out,
        "---\nlicense: {weights_license}\nlibrary_name: opendna\n{hub}datasets:\n  - \
         {dataset}\nbase_model: {base}\n---\n\n\
         # OpenDNA {task} head\n\n{summary}\n\n\
         This is a small head trained on **frozen** [{base}](https://huggingface.co/{base}) \
         embeddings. It is not a fine-tuned language model: the embedding model is unchanged, and \
         only the head below was trained.\n\n\
         {usage}\
         ## Inputs and outputs\n\n**Input.** {inputs}\n\n**Output.** {outputs}\n\n\
         ## Results\n\n{table}\n\
         Baseline: {baseline}\n\n\
         {threshold_note}\
         Validation chose the epoch that was kept (epoch {best_epoch} of {epochs_run} run). The \
         test split was evaluated once, afterwards. Splits come from the dataset builder, which \
         clusters by sequence identity so near-identical proteins cannot appear on both sides.\n\n\
         {external}\
         ## Training\n\n\
         | | |\n|---|---|\n\
         | Embedding model | `{base}` ({dim} dimensions, frozen) |\n\
         | Optimiser | AdamW, learning rate {lr} |\n\
         | Batch size | {batch} |\n\
         | Seed | {seed}{seed_note} |\n\
         | Dataset | [`{dataset}`]({dataset_url}) |\n\
         | Data licence | {license} |\n\n\
         ### Licence of these weights\n\n{weights_license_note}\n\n\
         ### Reproducing this run\n\n\
         Download the dataset from [`{dataset}`]({dataset_url}) into `datasets/out/{task}/`, or \
         rebuild it from its upstream sources with the `datasets/` workspace of the \
         [OpenDNA repository](https://github.com/dyber-pqc/OpenDNA) — the dataset card lists \
         every source file with its SHA-256. Then, from `training/`:\n\n\
         ```bash\ncargo run --release -p opendna-train -- embed --task {task}\n\
         cargo run --release -p opendna-train -- train --task {task} \\\n\
         \u{20}   --batch-size {batch} --learning-rate {lr} --seed {seed}\n```\n\n\
         The hyperparameters are spelled out above even where they match the current defaults, \
         so the command reproduces *this* run rather than whatever the defaults happen to be \
         later.\n\n\
         Both commands take `--datasets <dir>` if the dataset is not at the default \
         `../datasets/out`. Embedding is cached, so the second run of either command is cheap.\n\n\
         ## Limitations\n\n",
        dataset = p.dataset,
        base = report.embedding_model,
        task = task.name(),
        summary = p.summary,
        inputs = p.inputs,
        outputs = p.outputs,
        table = metrics_table(report),
        baseline = p.baseline,
        best_epoch = report.best_epoch,
        epochs_run = report.epochs_run,
        dim = report.embedding_dim,
        lr = report.learning_rate,
        batch = report.batch_size,
        seed = report.seed,
        seed_note = if report.seed_applied {
            ""
        } else {
            " (data order only; the training device could not seed weight initialisation)"
        },
        hub = hub_metadata(task),
        threshold_note = threshold_note(task),
        usage = usage_section(task),
        dataset_url = p.dataset_url,
        license = p.license,
        weights_license = p.weights_license,
        weights_license_note = p.weights_license_note,
        external = external_section(external),
    );
    for limitation in p.limitations {
        let _ = writeln!(out, "- {limitation}");
    }
    let _ = write!(
        out,
        "\n## Intended use\n\nRanking and triage during protein engineering: comparing candidate \
         mutations or constructs against each other. It is not a substitute for measurement, and \
         it should not be used to make claims about a specific protein without experimental \
         confirmation.\n\n\
         ## Citation\n\nThe training data comes from:\n\n> {citation}\n\n\
         Please cite that work as well as this model.\n",
        citation = p.citation
    );
    out
}

/// Write a model card next to the weights.
///
/// # Errors
/// If the file cannot be written.
/// What the card says about calibration.
///
/// Three states, and they are different claims: a correction was fitted and
/// applied, the question was measured and no correction was warranted, or
/// nobody has looked. An absent section would collapse the last two.
fn calibration_section(task: Task, calibration: Option<&opendna_ml::heads::Calibration>) -> String {
    if task != Task::Disorder {
        // The regressors emit quantities, not probabilities. Saying nothing
        // about calibration is correct for them; claiming it was checked
        // would not be.
        return String::new();
    }
    let Some(c) = calibration else {
        return "\n## Calibration\n\nNot measured. The probabilities this head reports have not \
                been checked against observed frequencies, so a reported 0.7 should be read as a \
                score for ranking rather than as a seven-in-ten chance.\n"
            .to_string();
    };
    if c.method == "none" {
        return format!(
            "\n## Calibration\n\nMeasured and left alone. Expected calibration error on the \
             held-out split is {:.4}, and no fitted correction improved on that by more than the \
             width of its own bootstrap interval, so none is applied. Residues are called \
             disordered at {:.2}.\n",
            c.test_ece_before, c.threshold
        );
    }
    format!(
        "\n## Calibration\n\nThe probabilities are corrected before they are reported. This head \
         is trained with weighted cross-entropy so that it does not predict \"ordered\" \
         everywhere, which helps the ranking and inflates the probabilities; a {} map fitted on \
         the validation split removes the inflation.\n\n\
         | | |\n|---|---|\n\
         | Map | `p = sigmoid({:.4} * logit + {:.4})` |\n\
         | Expected calibration error, held-out split | {:.4} before, {:.4} after |\n\
         | Threshold for a disorder call | {:.2} |\n\
         | Matthews correlation at that threshold, held-out split | {:.4} |\n\n\
         The map is monotonic, so it changes no ranking and no rank statistic: AUROC is \
         identical before and after. What it changes is the meaning of the number, and therefore \
         where the decision threshold belongs -- which is why the threshold is recorded here \
         rather than left at the conventional 0.5. The threshold was chosen on the validation \
         split, so the correlation above is an honest estimate rather than the best figure \
         findable with the held-out labels in hand.\n\n\
         Calibration is against the labels. DisProt records disorder that has been demonstrated \
         experimentally and unannotated residues are treated as ordered, so these probabilities \
         are calibrated against the annotation rather than against nature. If the annotation \
         under-reports disorder, and it does by an unknown amount, the calibrated probabilities \
         inherit that.\n\n\
         The threshold is a property of a distribution, not of the head. It was chosen on the \
         validation split of this dataset, and it does not suit every set: on the CAID3 \
         reference above, MCC, F1 and balanced accuracy are lower at this threshold than they \
         were at the conventional 0.5, because {:.2} calibrated is about 0.78 raw and that set \
         prefers a lower cut. The rank metrics are identical either way. If you threshold \
         these probabilities for a decision of your own, the right cut depends on your class \
         balance and on what a false positive costs you; the calibration is what makes it \
         possible to reason about that at all.\n",
        c.method,
        c.slope,
        c.intercept,
        c.test_ece_before,
        c.test_ece_after,
        c.threshold,
        c.test_mcc,
        c.threshold
    )
}

/// What the card says about confidence and abstention.
///
/// Includes the signals that were measured and rejected. A card that lists
/// only what worked invites the next person to spend an afternoon re-testing
/// what did not.
fn confidence_section(abstention: Option<&opendna_ml::heads::Abstention>) -> String {
    let Some(a) = abstention else {
        return String::new();
    };
    let mut out = String::from(
        "\n## Confidence, and whether to abstain\n\n\
         How far the calibrated probability sits from the decision threshold orders residues by \
         how likely the call is to be wrong, so a caller can trade coverage for accuracy. \
         Measured on the held-out split:\n\n\
         | Minimum confidence | Residues kept | Error among them |\n|---|---|---|\n",
    );
    for point in &a.curve {
        let _ = writeln!(
            out,
            "| {:.2} | {:.0}% | {:.1}% |",
            point.min_confidence,
            100.0 * point.coverage,
            100.0 * point.error
        );
    }
    let _ = write!(
        out,
        "\nNothing abstains by default. The confidence is reported for every residue and the \
         cut-off is the caller's, because the right one depends on what a wrong call costs them. \
         `opendna predict --min-confidence` marks the rows that do not clear a cut-off rather than \
         removing them: a table quietly missing rows looks exactly like a complete one.\n\n\
         This signal captures {:.0}% of the error reduction an oracle shows is available -- an \
         oracle being a rule that discards exactly the wrong calls. Half of what is there, which \
         is worth having and is not a guarantee.\n",
        100.0 * a.captured
    );
    if !a.rejected.is_empty() {
        let _ = write!(
            out,
            "\nMeasured and rejected: {}. Ordering residues by this was no better than \
             discarding them at random, so no rule is built on it. It is recorded because the \
             obvious next idea is to try it.\n",
            a.rejected.join(", ")
        );
    }
    out
}

/// What the card says about prediction intervals.
fn interval_section(task: Task, intervals: Option<&opendna_ml::heads::Intervals>) -> String {
    let Some(i) = intervals else {
        return String::new();
    };
    let mut out = String::from(
        "\n## How wrong the number is likely to be\n\n\
         A point estimate invites being read as a measurement. These intervals say how wide the \
         range around it has to be, fitted by split conformal prediction on the validation split: \
         the half-width is whatever the residuals require, with no assumption about the shape of \
         the error.\n\n\
         | Confidence | Interval | Coverage achieved on held-out data |\n|---|---|---|\n",
    );
    for level in &i.levels {
        let _ = writeln!(
            out,
            "| {:.0}% | ±{:.3} {} | {:.1}% |",
            100.0 * level.nominal,
            level.half_width,
            i.units,
            100.0 * level.achieved
        );
    }
    let _ = write!(
        out,
        "\nMean absolute error on the same split is {:.3} {}.\n\n\
         The achieved column is there because it differs from the nominal one, and it can \
         differ either way. A conformal guarantee assumes the calibration and held-out sets are \
         exchangeable, and this project's splits are deliberately not: the dataset builder \
         clusters by sequence identity so that near-identical proteins cannot straddle a split, \
         which makes the held-out set systematically more distant than a random sample. That is \
         the right way to build the splits, and the cost is that the guarantee is approximate. \
         The measured column is the one to believe, in whichever direction it points. \
         {direction}\n\n\
         {independence}\n\n\
         The half-width is constant, which is the method's main limitation: the interval is \
         wider than it needs to be where the head is reliable and narrower than it should be \
         where it is not.\n",
        i.held_out_mae,
        i.units,
        direction = coverage_direction(i),
        independence = independence_note(task)
    );
    out
}

pub fn write(task: Task, report: &Report, out_dir: &Path, benchmarks: &Path) -> Result<()> {
    let path = out_dir.join(format!("{}.README.md", task.name()));
    // Read from the metadata beside the weights, so the card cannot describe
    // a calibration the head is not carrying.
    let files = opendna_ml::heads::HeadFiles::in_dir(out_dir, task.name());
    // Benchmark reports are matched on the embedding model and on whether the
    // head takes burial features. Both are needed: the with-burial stability
    // variant scores the same suite on the same embedder, so matching on the
    // embedder alone put a rejected variant's figures on the published head's
    // card, listed first.
    let with_burial = std::fs::read_to_string(&files.metadata)
        .ok()
        .and_then(|t| serde_json::from_str::<opendna_ml::heads::HeadMetadata>(&t).ok())
        .is_some_and(|m| m.with_burial);
    let external = crate::bench::load_reports(benchmarks, task.name(), &report.embedding_model, with_burial)?;
    let calibration = std::fs::read_to_string(&files.metadata)
        .ok()
        .and_then(|t| serde_json::from_str::<opendna_ml::heads::HeadMetadata>(&t).ok())
        .and_then(|m| m.calibration);
    let abstention = std::fs::read_to_string(&files.metadata)
        .ok()
        .and_then(|t| serde_json::from_str::<opendna_ml::heads::HeadMetadata>(&t).ok())
        .and_then(|m| m.abstention);
    let intervals = std::fs::read_to_string(&files.metadata)
        .ok()
        .and_then(|t| serde_json::from_str::<opendna_ml::heads::HeadMetadata>(&t).ok())
        .and_then(|m| m.intervals);
    let body = render(task, report, &external)
        + &calibration_section(task, calibration.as_ref())
        + &confidence_section(abstention.as_ref())
        + &interval_section(task, intervals.as_ref());
    std::fs::write(&path, body).with_context(|| format!("writing {}", path.display()))?;
    println!(
        "Wrote {} ({} external benchmark{})",
        path.display(),
        external.len(),
        if external.len() == 1 { "" } else { "s" }
    );
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::train::SplitMetrics;
    use std::collections::BTreeMap;

    fn report() -> Report {
        let values = |pairs: &[(&str, f64)]| SplitMetrics {
            n: 100,
            values: pairs.iter().map(|(k, v)| ((*k).to_string(), *v)).collect::<BTreeMap<_, _>>(),
        };
        Report {
            task: "stability".into(),
            embedding_model: "facebook/esm2_t30_150M_UR50D".into(),
            embedding_dim: 640,
            best_epoch: 7,
            epochs_run: 17,
            learning_rate: 1e-3,
            batch_size: 256,
            seed: 1,
            seed_applied: true,
            validation: values(&[("spearman", 0.61)]),
            test: values(&[("spearman", 0.6)]),
            baseline: values(&[("spearman", 0.28)]),
        }
    }

    #[test]
    fn card_reports_the_measured_numbers_and_the_baseline() {
        let card = render(Task::Stability, &report(), &[]);
        assert!(card.contains("| Spearman ρ | 0.610 | 0.600 | 0.280 |"), "{card}");
        assert!(card.contains("epoch 7 of 17"));
        // The run used a non-default batch size once, and the card's command
        // omitted it — so following the card produced a different model than
        // the one it described. The command now carries every value.
        assert!(card.contains("--batch-size 256"), "{card}");
        assert!(card.contains("--seed 1"), "{card}");
        assert!(card.contains("facebook/esm2_t30_150M_UR50D"));
        assert!(card.starts_with("---\nlicense: apache-2.0"), "needs YAML front matter");
        assert!(card.contains("## Limitations"));
    }

    fn levels(pairs: &[(f64, f64)]) -> opendna_ml::heads::Intervals {
        opendna_ml::heads::Intervals {
            units: "RSA".into(),
            levels: pairs
                .iter()
                .map(|(nominal, achieved)| opendna_ml::heads::IntervalLevel {
                    nominal: *nominal,
                    half_width: 0.1,
                    achieved: *achieved,
                })
                .collect(),
            held_out_mae: 0.1,
        }
    }

    /// Whether an interval covers what it promises is the fact most worth not
    /// missing, and it is the one a reader skimming two columns misses.
    #[test]
    fn the_card_says_which_way_the_coverage_misses() {
        let under = levels(&[(0.5, 0.433), (0.8, 0.760), (0.9, 0.875)]);
        assert!(
            coverage_direction(&under).contains("cover less than they say"),
            "burial's measured coverage is short at every level"
        );
        let over = levels(&[(0.5, 0.546), (0.8, 0.852), (0.9, 0.937)]);
        assert!(
            coverage_direction(&over).contains("cover more than they say"),
            "stability's is long at every level"
        );
        let mixed = levels(&[(0.5, 0.502), (0.8, 0.770), (0.9, 0.911)]);
        assert!(
            coverage_direction(&mixed).contains("some levels and not others"),
            "solubility's goes both ways"
        );
    }

    /// The per-residue heads produce correlated residuals; solubility produces
    /// one per protein, so claiming a correlation there invents a problem.
    #[test]
    fn each_head_gets_its_own_account_of_what_the_quantile_rests_on() {
        let burial = independence_note(Task::Burial);
        assert!(burial.contains("one domain are correlated"), "{burial}");
        assert!(!burial.contains("wild-type"), "burial has no wild types: {burial}");
        assert!(!burial.contains("substitution"), "burial has no substitutions: {burial}");

        let stability = independence_note(Task::Stability);
        assert!(stability.contains("wild-type protein are correlated"), "{stability}");

        let solubility = independence_note(Task::Solubility);
        assert!(solubility.contains("exactly one residual"), "{solubility}");
        assert!(!solubility.contains("correlated"), "one residual per protein is independent: {solubility}");
    }

    /// The burial card credited Tien et al. as the source of its training
    /// data. Tien supplies the maxima the labels are divided by; the data is
    /// CATH domains and wwPDB coordinates.
    #[test]
    fn the_burial_card_credits_the_data_and_not_only_the_normalisation() {
        let card = render(Task::Burial, &report(), &[]);
        assert!(card.contains("Sillitoe"), "CATH is the source of the domains:\n{card}");
        assert!(card.contains("normalised by the theoretical maxima"), "{card}");
        assert!(card.contains("Tien"), "Tien still has to be cited, for what it is:\n{card}");
    }

    /// The dash in the baseline column has two opposite causes.
    ///
    /// Burial and solubility are compared against the training mean, which
    /// correlates with nothing, so the correlations are the dashes. Stability
    /// is compared against the language model's masked log-odds, a real
    /// predictor on another scale, so the correlations have values and the
    /// error metrics are the dashes. One sentence cannot explain both, and the
    /// one that was there explained the wrong one on the stability card.
    #[test]
    fn the_baseline_dash_is_explained_by_what_the_baseline_is() {
        // The real stability report: a baseline that correlates but has no
        // error in kcal/mol, because it is the log-odds.
        let mut ranking_baseline = report();
        for split in [&mut ranking_baseline.validation, &mut ranking_baseline.test] {
            split.values.insert("mae".into(), 0.5);
            split.values.insert("rmse".into(), 0.7);
        }
        let card = render(Task::Stability, &ranking_baseline, &[]);
        assert!(
            card.contains("a score rather than a prediction in the same units"),
            "a log-odds baseline has no error against the measurements:
{card}"
        );
        assert!(
            !card.contains("predicts the same value for every input"),
            "the log-odds baseline is not a constant:
{card}"
        );

        // And the other way round: a constant baseline has the error metrics
        // and no correlation, which is burial's and solubility's case.
        let mut constant_baseline = report();
        constant_baseline.baseline.values.remove("spearman");
        constant_baseline.baseline.values.insert("mae".into(), 0.19);
        let card = render(Task::Burial, &constant_baseline, &[]);
        assert!(card.contains("predicts the same value for every input"), "{card}");
        assert!(!card.contains("a score rather than a prediction"), "{card}");
    }

    /// Two Matthews correlations on one card need telling apart.
    ///
    /// The results table is at 0.5 on raw scores, because that is when the
    /// training report is written; the calibration table is at the head's own
    /// threshold on calibrated ones. Without the note a reader sees 0.493 and
    /// 0.5213 for the same-sounding metric and cannot reconcile them, and
    /// would not know which to compare against a published figure.
    #[test]
    fn the_classifier_card_says_which_threshold_its_metrics_used() {
        let card = render(Task::Disorder, &report(), &[]);
        assert!(card.contains("at a probability of 0.5 on the head's raw output"), "{card}");
        assert!(card.contains("rank statistics and do not depend on either threshold"), "{card}");
        // The regressors report MAE, RMSE and correlations, which no threshold
        // enters, so the note would be noise on their cards.
        for task in [Task::Burial, Task::Solubility, Task::Stability] {
            let card = render(task, &report(), &[]);
            assert!(
                !card.contains("raw output"),
                "{} has no decision threshold and should carry no note about one",
                task.name()
            );
        }
    }

    /// Only the head that classifies claims a classification pipeline.
    ///
    /// A `pipeline_tag` puts the repository under a task filter and drives an
    /// inference widget. Disorder labels each residue, so token classification
    /// is what it does. The other three return a number, and the Hub's list has
    /// no tag for that over a sequence, so they carry none — a tag that
    /// advertises classes a regressor does not have would be a false claim
    /// about the output, for a bit of discoverability.
    #[test]
    fn the_pipeline_tag_goes_only_where_it_is_true() {
        let classifier = render(Task::Disorder, &report(), &[]);
        assert!(
            classifier.contains("pipeline_tag: token-classification"),
            "the disorder card should claim the pipeline it fits:\n{classifier}"
        );
        for task in [Task::Burial, Task::Solubility, Task::Stability] {
            let card = render(task, &report(), &[]);
            assert!(
                !card.contains("pipeline_tag:"),
                "{} is a regressor and should claim no pipeline tag",
                task.name()
            );
            assert!(
                card.contains("  - regression"),
                "{} should still be findable as a regression",
                task.name()
            );
        }
    }

    /// Front matter the Hub cannot parse is front matter the Hub ignores, and
    /// then none of the metadata counts — licence included.
    #[test]
    fn the_front_matter_is_one_well_formed_block() {
        for task in [Task::Stability, Task::Solubility, Task::Disorder, Task::Burial] {
            let card = render(task, &report(), &[]);
            let rest = card.strip_prefix("---\n").expect("front matter opens the card");
            let (front, body) = rest.split_once("\n---\n").expect("front matter closes");
            assert!(body.starts_with("\n# OpenDNA"), "{} card: body follows the block", task.name());
            for line in front.lines() {
                // A list item, `key: value`, or `key:` opening a block —
                // which is what `tags:` and `datasets:` are.
                let shaped = line.starts_with("  - ")
                    || line.strip_suffix(':').is_some_and(|key| !key.contains(' '))
                    || line.split_once(": ").is_some_and(|(key, _)| !key.contains(' '));
                assert!(shaped, "{} card: `{line}` is neither a key nor a list item", task.name());
            }
            assert!(front.contains("library_name: opendna"), "{} card", task.name());
            assert!(front.contains("base_model: "), "{} card", task.name());
        }
    }

    /// A reader who lands on the card has to be able to run the thing.
    ///
    /// Hugging Face renders no "Use this model" snippet for this library, so
    /// the card is the only place the commands can come from. Each one has to
    /// name its own repository, and the verification step has to be one that
    /// works for somebody who downloaded a single head: `selftest` needs all
    /// four and exits 3 with fewer, so a card that told them to run it would
    /// be telling them to run a command that fails.
    #[test]
    fn every_card_says_how_to_run_the_weights() {
        for task in [Task::Stability, Task::Solubility, Task::Disorder, Task::Burial] {
            let card = render(task, &report(), &[]);
            let name = task.name();
            assert!(card.contains("## Using it"), "{name} card: no usage section");
            assert!(
                card.contains(&format!("hf download dyber-pqc/opendna-{name}")),
                "{name} card: does not say how to download itself"
            );
            assert!(
                card.contains(&format!("opendna predict --task {name}")),
                "{name} card: does not show a prediction command"
            );
            assert!(card.contains("opendna models"), "{name} card: no way to verify the download");
            assert!(
                !card.contains("\nopendna selftest\n"),
                "{name} card: tells a reader with one head to run a command needing four"
            );
        }
    }

    /// Stability scores named substitutions, so the bare command exits 2.
    #[test]
    fn the_stability_card_shows_the_mutation_file() {
        let card = render(Task::Stability, &report(), &[]);
        assert!(card.contains("--mutations muts.csv"), "{card}");
        assert!(card.contains("`name,mutation`"), "{card}");
    }

    /// The dataset link has to resolve for someone who only has the card.
    ///
    /// It used to point at `datasets/` inside the GitHub repository, a
    /// directory that is not on the published branch — so every card shipped
    /// a dead link to the one thing a reader most needs. It now points at the
    /// dataset's own Hugging Face repository, which is published alongside.
    #[test]
    fn the_dataset_link_points_at_a_published_dataset() {
        for task in [Task::Stability, Task::Solubility, Task::Disorder, Task::Burial] {
            let card = render(task, &report(), &[]);
            let expected = format!("https://huggingface.co/datasets/dyber-pqc/opendna-{}", task.name());
            assert!(card.contains(&expected), "{} card: missing {expected}", task.name());
            assert!(
                !card.contains("OpenDNA/tree/main/datasets"),
                "{} card still links to a directory that is not on the published branch",
                task.name()
            );
        }
    }

    /// Share-alike data must not be laundered into an Apache-2.0 model.
    ///
    /// eSOL is CC BY-SA 2.1 JP. Declaring the head Apache-2.0 would assert
    /// that weights are not a derivative work of their training data, which
    /// is an open question and not one to settle in our own favour.
    #[test]
    fn a_head_trained_on_share_alike_data_does_not_claim_apache() {
        let solubility = render(Task::Solubility, &report(), &[]);
        assert!(solubility.starts_with("---\nlicense: other\n"), "{solubility}");
        assert!(solubility.contains("license_name: cc-by-sa-2.1-jp"), "{solubility}");
        assert!(solubility.contains("### Licence of these weights"), "{solubility}");
        // The other two are CC BY, which carries no such condition.
        for task in [Task::Stability, Task::Disorder] {
            let card = render(task, &report(), &[]);
            assert!(card.starts_with("---\nlicense: apache-2.0"), "{} card", task.name());
        }
    }

    fn bench(suite: &str) -> crate::bench::BenchReport {
        crate::bench::BenchReport {
            suite: suite.into(),
            embedding_model: "facebook/esm2_t30_150M_UR50D".into(),
            head_version: "0.6.0-dev".into(),
            with_burial: false,
            proteins: 319,
            residues_scored: 99_239,
            residues_masked: 61_183,
            positive_rate: Some(0.3164),
            windowed: 34,
            pooled: vec![("auroc".into(), 0.941_02)],
            per_protein_auroc: Some(0.925_2),
            per_protein_undefined: 86,
            caveats: vec!["Not an official CAID3 result.".into()],
        }
    }

    /// An in-domain number on its own invites a reader to treat it as a
    /// standing rather than a measurement.
    #[test]
    fn external_benchmarks_appear_on_the_card() {
        let card = render(Task::Disorder, &report(), &[bench("caid3-disorder-pdb")]);
        assert!(card.contains("### Against published benchmarks"), "{card}");
        assert!(card.contains("caid3-disorder-pdb"), "{card}");
        assert!(card.contains("| auroc | 0.9410 |"), "{card}");
        assert!(card.contains("Not an official CAID3 result."), "{card}");
    }

    /// Having no external benchmark is a fact about the task, not a gap to
    /// pass over in silence: solubility has none, and the card says why.
    #[test]
    fn the_absence_of_an_external_benchmark_is_stated() {
        let card = render(Task::Solubility, &report(), &[]);
        assert!(card.contains("### Against published benchmarks"), "{card}");
        assert!(card.contains("None."), "{card}");
        assert!(card.contains("not evidence of where it stands"), "{card}");
    }

    /// A report measured on another embedding model must not be quoted on a
    /// card documenting this one, and a suite belonging to another head must
    /// not be picked up at all.
    #[test]
    fn a_report_from_a_different_embedder_is_not_loaded() {
        let dir = std::env::temp_dir().join(format!("opendna-card-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut other = bench("caid3-disorder-pdb");
        other.embedding_model = "facebook/esm2_t33_650M_UR50D".into();
        std::fs::write(dir.join("caid3.json"), serde_json::to_string(&other).unwrap()).unwrap();

        let loaded = |task, model| crate::bench::load_reports(&dir, task, model, false).unwrap().len();
        assert_eq!(loaded("disorder", "facebook/esm2_t30_150M_UR50D"), 0, "650M quoted on a 150M card");
        assert_eq!(loaded("disorder", "facebook/esm2_t33_650M_UR50D"), 1);
        assert_eq!(loaded("stability", "facebook/esm2_t33_650M_UR50D"), 0, "wrong head's suite");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// The stability card carried the `ProteinGym` section twice.
    ///
    /// Same suite, same embedder, same item count, same caveats, different
    /// numbers: the with-burial variant was measured, found not to help, and
    /// never published, but it scored the same suite on the same embedder as
    /// the head that ships. Matching on those two alone put it on the card,
    /// listed first, so the figure a reader saw before any other was a model
    /// nobody can download.
    #[test]
    fn a_variant_that_was_not_published_stays_off_the_card() {
        let dir = std::env::temp_dir().join(format!("opendna-variant-{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).unwrap();

        let shipped = bench("proteingym-stability-out-of-domain");
        let mut variant = shipped.clone();
        variant.with_burial = true;
        std::fs::write(dir.join("plain.json"), serde_json::to_string(&shipped).unwrap()).unwrap();
        std::fs::write(dir.join("burial.json"), serde_json::to_string(&variant).unwrap()).unwrap();

        let model = "facebook/esm2_t30_150M_UR50D";
        let for_head =
            |with_burial| crate::bench::load_reports(&dir, "stability", model, with_burial).unwrap();
        let plain = for_head(false);
        assert_eq!(plain.len(), 1, "a head without burial features gets one report, not both");
        assert!(!plain[0].with_burial);
        let with = for_head(true);
        assert_eq!(with.len(), 1, "and the variant's card gets the variant's report");
        assert!(with[0].with_burial);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_missing_benchmarks_directory_is_not_an_error() {
        let got = crate::bench::load_reports(
            Path::new("definitely/not/here"),
            "disorder",
            "facebook/esm2_t30_150M_UR50D",
            false,
        )
        .unwrap();
        assert!(got.is_empty());
    }

    /// `FireProtDB` was dropped from the stability dataset for licensing
    /// reasons, but the card went on crediting it — attributing data that is
    /// not in the model.
    #[test]
    fn the_stability_card_credits_only_the_source_that_is_in_the_data() {
        let card = render(Task::Stability, &report(), &[]);
        assert!(!card.contains("FireProtDB"), "{card}");
        assert!(card.contains("Tsuboyama"), "{card}");
    }

    #[test]
    fn a_missing_baseline_metric_is_explained_rather_than_left_as_a_dash() {
        let mut r = report();
        r.test.values.insert("rmse".into(), 23.4);
        // The constant baseline has rmse but no spearman.
        r.baseline.values.remove("spearman");
        r.baseline.values.insert("rmse".into(), 32.3);
        let card = render(Task::Solubility, &r, &[]);
        assert!(card.contains("| Spearman ρ | 0.610 | 0.600 | — |"), "{card}");
        assert!(card.contains("undefined for that baseline"), "{card}");
    }

    #[test]
    fn a_baseline_with_every_metric_gets_no_footnote() {
        assert!(!render(Task::Stability, &report(), &[]).contains("undefined for that baseline"));
    }

    #[test]
    fn an_unseeded_run_says_so() {
        let mut r = report();
        r.seed_applied = false;
        assert!(render(Task::Stability, &r, &[]).contains("could not seed"));
    }

    #[test]
    fn every_task_renders() {
        for task in [Task::Stability, Task::Solubility, Task::Disorder] {
            let mut r = report();
            r.task = task.name().to_string();
            let card = render(task, &r, &[]);
            assert!(card.contains("## Limitations"), "{} has no limitations", task.name());
            assert!(card.contains("## Citation"));
        }
    }
}
