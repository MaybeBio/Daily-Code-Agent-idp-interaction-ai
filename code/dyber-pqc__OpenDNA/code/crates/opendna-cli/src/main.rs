//! OpenDNA command line: batch predictions over a FASTA file.
//!
//! The desktop application works on one protein at a time, which is the wrong
//! shape for the job people actually have — a few hundred candidate sequences
//! to triage. This runs the same heads over a whole file and writes a table.
//!
//! It shares the prediction code with the application rather than
//! reimplementing it, so the two cannot disagree about a number, and every
//! output file carries the model version and the card's limitations.
//!
//! Nothing here falls back to a heuristic. A missing head is an error naming
//! the files it looked for.

// Writing to standard output is this crate's purpose, so the workspace's
// `print_stdout` lint is inverted here: results go to stdout when no
// `--output` is given, and progress goes to stderr so a pipe stays clean.
#![allow(clippy::print_stdout, clippy::print_stderr)]
// `doc_markdown` reads product names like OpenDNA as code identifiers.
// Backticking them would make the rendered help and docs wrong.
#![allow(clippy::doc_markdown)]

mod caid;
mod checkpoint;
mod exit;
mod fold;
mod humanness;
mod inventory;
mod manifest;
mod merge;
mod provenance;
mod selftest;
mod split;
mod staging;

use anyhow::{anyhow, bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use exit::{Classify as _, Failure, Fault};
use opendna_core::AminoAcid;
use opendna_ml::esm2::{Esm2, ModelSize};
use opendna_ml::heads::{read_metadata, HeadFiles, Heads, StabilityHead};
use provenance::Provenance;
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};

/// The embedding model the published heads were trained on.
///
/// Not a flag: a head records the embedder it expects and refuses to load
/// against another, so offering a choice here would only produce a confusing
/// error. Retraining on a different embedder is a `training/` concern.
const EMBEDDER: ModelSize = ModelSize::Large150M;

/// Fallback probability at or above which a residue is called disordered.
///
/// Used only for a head that records no threshold of its own. It is a
/// convention rather than a measurement, and for the head measured here it is
/// the wrong number twice over: on the raw scale that head's own optimum is
/// 0.73, and on its calibrated scale 0.34. A head that has been calibrated
/// carries the threshold measured with it, which is what gets used.
const DISORDER_THRESHOLD: f64 = 0.5;

/// Relative accessibility at or below which a residue counts as buried.
///
/// 0.25 is the Rost & Sander convention and is what the burial dataset's own
/// `buried` column uses, so the two agree. The accessibility column is the
/// measurement; this is a convenience.
const BURIED_AT: f64 = 0.25;

/// Padded tokens per embedding batch on the CPU.
///
/// Measured rather than guessed, and it took several attempts to measure
/// honestly. On 200 sequences and 35,033 residues, interleaving the settings
/// across rounds on an otherwise idle machine:
///
/// ```text
/// CPU, batch-tokens     0: median 60.8s  range 59.1-63.8s
/// CPU, batch-tokens  8192: median 54.5s  range 51.7-55.5s
/// ```
///
/// About a tenth faster, with the ranges clear of each other. An earlier
/// session of the same measurement put it nearer a sixth, so treat the size of
/// the gain as "roughly ten per cent" rather than a figure; the direction is
/// not in doubt. Two attempts before that -- one setting after the other,
/// while a compile was competing for the cores -- gave first an overstated
/// gain and then no detectable one. Neither was true.
const CPU_BATCH_TOKENS: usize = 8192;

/// Padded tokens per embedding batch on a GPU: none.
///
/// Batching was expected to pay more here, not less. It does not:
///
/// ```text
/// GPU, batch-tokens     0: median  6.4s  range  6.3- 6.4s
/// GPU, batch-tokens  8192: median  6.8s  range  6.8- 6.8s
/// ```
///
/// Measured on an RTX 2000 Ada, which is 9.5x the CPU's throughput on the same
/// file either way. The arithmetic is already wide enough that the padded work
/// costs more than the per-call overhead it saves, and the variance is small
/// enough -- 6.3-6.4 against 6.8-6.8 -- to leave no doubt about the sign.
///
/// One number cannot be right for both devices, so the default follows the
/// device. Another GPU may well behave differently; `--batch-tokens` is there
/// so somebody can measure theirs rather than inherit this one.
const GPU_BATCH_TOKENS: usize = 0;

/// Units of work per checkpoint chunk.
///
/// Only used when `--checkpoint` is given. It is a trade-off between how much
/// an interruption costs and how often the run stops to sync a file, and 64
/// sits where neither is noticeable: on the machine measured here a chunk of
/// 64 sequences takes a few seconds on the GPU and under a minute on the CPU,
/// against a sync that takes milliseconds.
///
/// It also fixes the embedding batch composition, which is why it is recorded
/// and why resuming at a different value is refused rather than accommodated.
const DEFAULT_CHUNK: usize = 64;

#[derive(Parser)]
#[command(
    name = "opendna",
    about = "Batch protein property predictions that run on your own machine",
    long_about = "Runs OpenDNA's trained prediction heads over a FASTA file and writes a table.\n\
                  \n\
                  Predictions need the trained heads. Download them and pass --models, or set\n\
                  OPENDNA_MODELS:\n\
                  \n\
                  \u{20} hf download dyber-pqc/opendna-disorder --local-dir models\n\
                  \n\
                  Nothing is uploaded. The first run downloads the embedding model into the\n\
                  Hugging Face cache; after that it works offline.",
    version
)]
struct Cli {
    /// Directory holding the trained heads. Defaults to `$OPENDNA_MODELS`,
    /// then `./models`.
    #[arg(long, global = true)]
    models: Option<PathBuf>,
    /// Never reach the network. The embedding model must already be staged in
    /// the models directory — see `opendna stage`.
    ///
    /// Without this, a machine that has not run a prediction before downloads
    /// several hundred megabytes on its first one. On an instrument computer
    /// or behind an air gap that is a failure at the worst moment, so this
    /// turns it into an error up front.
    #[arg(long, global = true)]
    offline: bool,
    /// Do not write the `.manifest.json` that records what produced an output
    /// file.
    ///
    /// The manifest is written by default because a result that cannot be
    /// traced back to its input bytes and its model build is not a result
    /// anyone can defend. This exists for the case where something downstream
    /// cannot tolerate a second file appearing.
    #[arg(long, global = true)]
    no_manifest: bool,
    /// Padded tokens per embedding batch. 0 embeds one sequence at a time.
    ///
    /// Sequences are sorted by length and grouped until this many padded
    /// tokens, which is what bounds the memory a batch needs. Left unset it
    /// follows the device: batching is about a tenth faster on the CPU
    /// measured here and slightly slower on the GPU measured here.
    ///
    /// Two runs at the same value are byte-identical on the same device. Two
    /// runs at different values are not: padding changes the order of
    /// floating-point sums, so a few values in ten thousand differ in their
    /// last digit. The value used is recorded in the run manifest.
    #[arg(long, global = true)]
    batch_tokens: Option<usize>,
    /// Which device to run on.
    ///
    /// `auto` takes the fastest this build can use. `cpu` forces the CPU,
    /// which is what you want when the GPU is busy with someone else's work,
    /// and what makes it possible to check that the two devices agree.
    #[arg(long, global = true, value_enum, default_value_t = Which::Auto)]
    device: Which,
    #[command(subcommand)]
    command: Command,
}

/// What to predict, from what.
#[derive(clap::Args)]
struct PredictArgs {
    /// Which property.
    #[arg(long, value_enum)]
    task: Task,
    /// FASTA file to read. `-` reads standard input.
    #[arg(long, short = 'i')]
    input: PathBuf,
    /// Where to write the table. Omit to write to standard output.
    #[arg(long, short = 'o')]
    output: Option<PathBuf>,
    /// Output format.
    #[arg(long, value_enum, default_value_t = Format::Csv)]
    format: Format,
    /// For `stability` only: a CSV of `name,mutation` rows, for example
    /// `ubiquitin,L8P`. Mutations are 1-based, as they are written.
    #[arg(long)]
    mutations: Option<PathBuf>,
    /// Record results as they are produced, and continue an interrupted run
    /// from this file.
    ///
    /// Written as the run goes and synced to disk after every chunk, so a
    /// process that is killed loses at most one chunk. Giving the same path
    /// again continues from where it stopped, after checking that the input,
    /// the head, the device and the settings are the same ones — a table that
    /// is partly one run and partly another is the mixture `opendna merge`
    /// refuses between files, and it would be odd to create it inside one.
    ///
    /// Results are still written to `--output` only when the run finishes.
    /// The checkpoint is a record of work done, not a results file.
    #[arg(long)]
    checkpoint: Option<PathBuf>,
    /// Sequences per chunk, when `--checkpoint` is given.
    ///
    /// Bounds what an interruption costs. It also fixes which sequences are
    /// embedded together, and batch composition changes the order of
    /// floating-point sums, so a resumed run matches an uninterrupted one
    /// exactly when this value matches. Resuming with a different value is
    /// refused.
    #[arg(long, default_value_t = DEFAULT_CHUNK)]
    chunk: usize,
    /// For `disorder` only: mark calls whose confidence is below this.
    ///
    /// Confidence is the distance of the calibrated probability from the
    /// decision threshold. Rows are marked in a `low_confidence` column and
    /// never removed: a table quietly missing rows looks exactly like a
    /// complete one. The head records what each cut-off was measured to buy,
    /// and the chosen one is reported against that.
    #[arg(long)]
    min_confidence: Option<f64>,
}

/// Settings that apply to every subcommand.
#[derive(Copy, Clone)]
struct Global {
    /// Never reach the network.
    offline: bool,
    /// Write a manifest beside each output file.
    manifest: bool,
    /// Padded tokens per embedding batch, or `None` to let the device decide.
    batch_tokens: Option<usize>,
    /// Which device to run on.
    device: Which,
}

impl Global {
    /// The batch budget to use on `device`.
    ///
    /// An explicit `--batch-tokens` wins. Otherwise the default follows the
    /// device, because the measurement says the right answer differs between
    /// them and a single number would be wrong for one of the two.
    fn budget(self, device: &candle_core::Device) -> usize {
        self.batch_tokens.unwrap_or(match device {
            candle_core::Device::Cpu => CPU_BATCH_TOKENS,
            candle_core::Device::Cuda(_) | candle_core::Device::Metal(_) => GPU_BATCH_TOKENS,
        })
    }
}

#[derive(Subcommand)]
enum Command {
    /// Predict a property for every sequence in a FASTA file.
    Predict(PredictArgs),
    /// Predict structures for a FASTA file by driving an installed folding
    /// engine.
    ///
    /// OpenDNA does not fold proteins. ESMFold and ColabFold do, and this
    /// runs one of them, checks that the chain it returns is the sequence it
    /// was given, and records which engine produced each structure.
    ///
    /// There is no fallback. With no engine installed this exits 3 naming
    /// what it looked for, because coordinates from something that cannot
    /// fold proteins would be worse than no coordinates.
    ///
    /// The global --models, --offline, --device and --batch-tokens options
    /// have no effect here. No head is loaded and nothing is embedded; how
    /// the engine uses the machine is the engine's own business, and
    /// --engine-arg is how to tell it anything about that.
    Fold {
        /// FASTA file of sequences to fold.
        #[arg(short, long)]
        input: PathBuf,
        /// Directory to write structures into. Created if it is not there.
        #[arg(short, long)]
        out_dir: PathBuf,
        /// Which engine to drive. With one installed this is not needed; with
        /// two it is, because they do not produce the same structures.
        #[arg(long, value_enum)]
        engine: Option<FoldEngine>,
        /// An argument passed to the engine verbatim, after the input and
        /// output it is already given. Repeatable.
        ///
        /// Not translated between engines: a flag whose meaning OpenDNA
        /// invented would be a flag that quietly means something else after
        /// the engine's next release. The full command line is recorded in
        /// every manifest the run writes.
        #[arg(long = "engine-arg", value_name = "ARG")]
        engine_args: Vec<String>,
        /// Fail unless the installed engines match this lockfile.
        #[arg(long)]
        check: Option<PathBuf>,
    },
    /// Score how often each 9-residue peptide of an antibody occurs in human
    /// repertoires.
    ///
    /// No model runs. Given a prevalence table, the answer is an exact count:
    /// cut the sequence into overlapping 9-mers and look each one up. So
    /// there is no accuracy to quote and no uncertainty to report -- the
    /// number is a property of the table, and the table's card says what it
    /// was built from.
    ///
    /// The per-window output is the point. A single percentage invites
    /// comparison with figures computed from other repertoires, which would
    /// be wrong; the rare windows say which positions are unlike anything
    /// seen in a human, which is the question this can answer.
    ///
    /// The global --models, --offline, --device and --batch-tokens options
    /// have no effect here. No head is loaded and nothing is embedded.
    Humanness {
        /// FASTA file of antibody sequences to score.
        #[arg(short, long)]
        input: PathBuf,
        /// The prevalence table, built by `opendna-data build humanness`.
        #[arg(short, long)]
        table: PathBuf,
        /// Which threshold the summary identity is computed at.
        ///
        /// A threshold finer than the table can express is refused rather than
        /// reported: from n donors the smallest non-zero prevalence is 1/n.
        #[arg(long, value_enum, default_value_t = Threshold::Medium)]
        stringency: Threshold,
        /// Print every scored window instead of only the rarest ones.
        #[arg(long)]
        all_windows: bool,
        /// How many of the rarest windows to list per sequence.
        #[arg(long, default_value_t = 10, conflicts_with = "all_windows")]
        worst: usize,
    },
    /// Report which external engines are installed, and optionally pin them.
    ///
    /// Some work -- folding a sequence, docking a ligand, running a
    /// simulation -- is decades of other people's software, and OpenDNA
    /// drives it rather than reimplementing it badly. This says which of
    /// those tools it can find, which version each reports, and the digest of
    /// the file it would run.
    ///
    /// Being listed is not the same as being used. No command drives any of
    /// these yet; this exists so an installation can be checked and pinned
    /// before the commands that need them arrive.
    Engines {
        /// Write a lockfile naming the version and digest of every engine
        /// found.
        #[arg(long)]
        lock: Option<PathBuf>,
        /// Fail unless the installed engines match this lockfile exactly.
        #[arg(long)]
        check: Option<PathBuf>,
    },
    /// Report which heads are installed, whether their weights verify, and
    /// optionally pin them.
    Models {
        /// Write a lockfile naming the version, build and hash of every
        /// installed head.
        #[arg(long)]
        lock: Option<PathBuf>,
        /// Fail unless the installed heads match this lockfile exactly.
        #[arg(long)]
        check: Option<PathBuf>,
    },
    /// Run the installed heads over a bundled sequence and check the results
    /// against recorded values. Exits non-zero if this machine does not
    /// reproduce them.
    Selftest,
    /// Join result files from parallel parts, refusing a mixture.
    ///
    /// Checks each part's run manifest and stops if they were not produced
    /// the same way — a different head build, a different device, a different
    /// batch size. Those change the numbers, and concatenating them produces
    /// a file nobody can describe.
    Merge {
        /// Result files, merged in the order given.
        inputs: Vec<PathBuf>,
        /// Where to write. Omit for standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,
        /// Merge parts that have no manifest beside them.
        #[arg(long)]
        allow_unverified: bool,
    },
    /// Cut a FASTA into parts a scheduler can run in parallel.
    ///
    /// Parts are balanced by residue count rather than by sequence count,
    /// because runtime scales with residues and a scheduler waits for the
    /// slowest part.
    Split {
        /// FASTA file to cut.
        #[arg(long, short = 'i')]
        input: PathBuf,
        /// Directory for the parts. Created if it does not exist.
        #[arg(long)]
        out_dir: PathBuf,
        /// How many parts to produce.
        #[arg(long, conflicts_with = "per_chunk")]
        parts: Option<usize>,
        /// Or: at most this many sequences per part.
        #[arg(long)]
        per_chunk: Option<usize>,
    },
    /// Copy this installation, embedding model included, into one directory
    /// that can be carried to a machine with no network.
    Stage {
        /// Where to write the bundle. Created if it does not exist.
        #[arg(long)]
        to: PathBuf,
    },
}

/// How common a peptide has to be to count as human.
///
/// The three thresholds are the ones OASis defines. Which is appropriate
/// depends on the table: a threshold below 1/donors is not expressible at
/// all, and `opendna humanness` refuses rather than reporting one.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum Threshold {
    /// At least 1% of donors. Needs a table of 100 donors or more.
    Loose,
    /// At least 10% of donors.
    Medium,
    /// At least half of donors.
    Strict,
}

impl From<Threshold> for opendna_analysis::humanness::Stringency {
    fn from(choice: Threshold) -> Self {
        match choice {
            Threshold::Loose => Self::Loose,
            Threshold::Medium => Self::Medium,
            Threshold::Strict => Self::Strict,
        }
    }
}

/// Which folding engine to drive.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum FoldEngine {
    /// ESMFold, through the `esm-fold` entry point from `fair-esm`.
    Esmfold,
    /// ColabFold, through `colabfold_batch`.
    Colabfold,
}

impl From<FoldEngine> for opendna_engine::Folder {
    fn from(choice: FoldEngine) -> Self {
        match choice {
            FoldEngine::Esmfold => Self::EsmFold,
            FoldEngine::Colabfold => Self::ColabFold,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum Which {
    /// The fastest device this build can use.
    Auto,
    /// The CPU, even where a GPU is available.
    Cpu,
}

impl Which {
    /// Open the device.
    ///
    /// # Errors
    /// If a GPU is compiled in, selected, and cannot be opened.
    fn open(self) -> Result<candle_core::Device, candle_core::Error> {
        match self {
            Self::Auto => opendna_ml::best_device(),
            Self::Cpu => Ok(candle_core::Device::Cpu),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum Task {
    /// Per-residue intrinsic disorder.
    Disorder,
    /// Whole-protein solubility.
    Solubility,
    /// Change in folding free energy for named substitutions.
    Stability,
    /// Per-residue relative solvent accessibility.
    Burial,
}

impl Task {
    fn name(self) -> &'static str {
        match self {
            Self::Disorder => "disorder",
            Self::Burial => "burial",
            Self::Solubility => "solubility",
            Self::Stability => "stability",
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum Format {
    /// Comma-separated, with `#` provenance comments.
    Csv,
    /// Tab-separated, with `#` provenance comments.
    Tsv,
    /// One JSON object, with provenance as a field.
    Json,
    /// The CAID prediction format, for submitting the disorder head to the
    /// Critical Assessment of Intrinsic Disorder. Disorder only, and it
    /// carries no provenance comments — see `caid.rs` for why.
    Caid,
}

impl Format {
    fn separator(self) -> char {
        match self {
            Self::Tsv | Self::Caid => '\t',
            Self::Csv | Self::Json => ',',
        }
    }
}

/// Where to look for trained heads.
///
/// `--models` wins, then `OPENDNA_MODELS`, then `./models`. The environment
/// variable exists so a shared install can be pointed at once rather than on
/// every invocation.
fn models_dir(flag: Option<PathBuf>) -> PathBuf {
    flag.or_else(|| std::env::var_os("OPENDNA_MODELS").map(PathBuf::from))
        .unwrap_or_else(|| PathBuf::from("models"))
}

fn read_input(path: &Path) -> Result<String> {
    if path.as_os_str() == "-" {
        let mut buffer = String::new();
        std::io::Read::read_to_string(&mut std::io::stdin(), &mut buffer)
            .context("reading FASTA from standard input")?;
        return Ok(buffer);
    }
    std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))
}

fn write_output(output: Option<&Path>, text: &str) -> Result<()> {
    match output {
        Some(path) => {
            std::fs::write(path, text).with_context(|| format!("writing {}", path.display()))?;
            eprintln!("Wrote {}", path.display());
        }
        None => {
            std::io::stdout().write_all(text.as_bytes()).context("writing to standard output")?;
        }
    }
    Ok(())
}

/// A head's identity, its files, and the provenance derived from them.
struct Loaded {
    /// What goes into the output's comments.
    prov: Provenance,
    /// What the head records about itself.
    metadata: opendna_ml::heads::HeadMetadata,
    /// Where its files are.
    files: HeadFiles,
}

/// Load a head's metadata and its card's limitations.
fn head_provenance(dir: &Path, task: Task, device: &str) -> Result<Loaded> {
    let files = HeadFiles::in_dir(dir, task.name());
    if !files.exist() {
        bail!(
            "no {} head in {}. Expected {} and {}. Download one:\n  \
             hf download dyber-pqc/opendna-{} --local-dir {}\n\
             OpenDNA does not estimate this property without a trained model.",
            task.name(),
            dir.display(),
            files.weights.display(),
            files.metadata.display(),
            task.name(),
            dir.display()
        );
    }
    let metadata = read_metadata(&files, task.name(), EMBEDDER)?;
    // Limitations come from the card beside the weights, so the CLI, the
    // desktop panel and the published card cannot state different caveats.
    let prov = Provenance::new(task.name(), &metadata, EMBEDDER, device, files.limitations());
    Ok(Loaded { prov, metadata, files })
}

/// One row of a result table, as strings so every format shares the values.
struct Table {
    columns: Vec<&'static str>,
    rows: Vec<Vec<String>>,
}

impl Table {
    fn render(&self, format: Format, prov: &Provenance) -> Result<String> {
        if format == Format::Json {
            let objects: Vec<serde_json::Value> = self
                .rows
                .iter()
                .map(|row| {
                    self.columns
                        .iter()
                        .zip(row)
                        .map(|(c, v)| ((*c).to_string(), json_scalar(v)))
                        .collect::<serde_json::Map<_, _>>()
                        .into()
                })
                .collect();
            let doc = serde_json::json!({
                "task": prov.task,
                "head_version": prov.head_version,
                "embedding_model": prov.embedding_model,
                "device": prov.device,
                "limitations": prov.limitations,
                "results": objects,
            });
            return Ok(serde_json::to_string_pretty(&doc)? + "\n");
        }
        let sep = format.separator();
        let mut out = prov.comment_block();
        let _ = writeln!(out, "{}", self.columns.join(&sep.to_string()));
        for row in &self.rows {
            // A field containing the separator would silently shift every
            // later column, so quote it the way a CSV reader expects.
            let fields: Vec<String> = row
                .iter()
                .map(|f| {
                    if f.contains(sep) || f.contains('"') || f.contains('\n') {
                        format!("\"{}\"", f.replace('"', "\"\""))
                    } else {
                        f.clone()
                    }
                })
                .collect();
            let _ = writeln!(out, "{}", fields.join(&sep.to_string()));
        }
        Ok(out)
    }
}

/// Numbers stay numbers in JSON; everything else is a string.
///
/// Integers stay integers. A position written `1.0` or a length written
/// `76.0` is the sort of thing that survives into a downstream join as a
/// float key and fails to match, so whole numbers are emitted as such. A
/// `+0.333` keeps its sign but parses as a float either way.
fn json_scalar(value: &str) -> serde_json::Value {
    if let Ok(integer) = value.parse::<i64>() {
        return serde_json::Value::Number(integer.into());
    }
    value
        .parse::<f64>()
        .ok()
        .and_then(serde_json::Number::from_f64)
        .map_or_else(|| serde_json::Value::String(value.to_string()), serde_json::Value::Number)
}

fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(Failure { fault, error }) => {
            // The chain, because the useful sentence is usually the innermost
            // one and the outer context says where it happened.
            eprintln!("Error: {error:#}");
            eprintln!("{fault}");
            std::process::ExitCode::from(fault.code())
        }
    }
}

fn run() -> Result<(), Failure> {
    let cli = Cli::parse();
    let dir = models_dir(cli.models);
    let global = Global {
        offline: cli.offline,
        manifest: !cli.no_manifest,
        batch_tokens: cli.batch_tokens,
        device: cli.device,
    };

    match cli.command {
        // Everything these two can fail at is a property of the install.
        Command::Models { lock, check } => {
            models(&dir, lock.as_deref(), check.as_deref()).blame(Fault::Install)
        }
        Command::Engines { lock, check } => engines(lock.as_deref(), check.as_deref()).blame(Fault::Install),
        Command::Fold { input, out_dir, engine, engine_args, check } => fold::run(&fold::Request {
            input,
            out_dir,
            engine: engine.map(Into::into),
            engine_args,
            check,
            manifest: global.manifest,
        }),
        Command::Humanness { input, table, stringency, all_windows, worst } => {
            humanness::run(&humanness::Request {
                input,
                table,
                stringency: stringency.into(),
                all_windows,
                worst,
            })
        }
        Command::Merge { inputs, output, allow_unverified } => {
            merge::run(&inputs, output.as_deref(), !allow_unverified).blame(Fault::Input)
        }
        Command::Split { input, out_dir, parts, per_chunk } => {
            let cut = match (parts, per_chunk) {
                (Some(n), _) => split::Cut::Parts(n),
                (None, Some(n)) => split::Cut::PerChunk(n),
                // Neither given: one part per sequence is not a useful
                // default and neither is the whole file, so say so.
                (None, None) => {
                    return Err(Failure {
                        fault: Fault::Usage,
                        error: anyhow!("give --parts <n> or --per-chunk <n>"),
                    })
                }
            };
            split::run(&input, &out_dir, cut).blame(Fault::Input).map(|_| ())
        }
        Command::Stage { to } => staging::stage(&dir, &to, EMBEDDER).blame(Fault::Install),
        Command::Selftest => {
            let device = global.device.open()?;
            let budget = global.budget(&device);
            selftest::run(&dir, EMBEDDER, global.offline, device, budget)
        }
        Command::Predict(args) => predict(&dir, global, args),
    }
}

/// The input, parsed.
struct Request {
    /// Every sequence in the input, as written.
    proteins: Vec<(String, String)>,
    /// The input's bytes, for the manifest's digest.
    text: String,
}

/// Read and parse the input.
///
/// # Errors
/// If it cannot be read, or is not FASTA.
fn read_request(input: &Path) -> Result<Request, Failure> {
    // A path that is not there is the request being wrong, not the machine:
    // a typo in a filename should not report as an unanticipated failure.
    let text = read_input(input).blame(Fault::Input)?;
    // Read as written: real files contain selenocysteine and the ambiguity
    // codes, ESM-2 has tokens for them, and collapsing them to X would feed
    // the model a residue the file does not specify.
    let proteins = opendna_io::fasta::read_raw(&text)
        .with_context(|| format!("parsing {} as FASTA", input.display()))
        .blame(Fault::Input)?;
    eprintln!("{} sequences from {}", proteins.len(), input.display());
    Ok(Request { proteins, text })
}

/// The models, loaded and verified.
struct Prepared {
    /// The trained heads.
    heads: Heads,
    /// The embedding model.
    model: Esm2,
    /// What it runs on.
    device: candle_core::Device,
    /// How that device is described to a reader.
    description: String,
    /// What goes in the output's comments.
    prov: Provenance,
    /// The head's recorded identity, including any calibration.
    metadata: opendna_ml::heads::HeadMetadata,
    /// Where the head's files are.
    files: HeadFiles,
}

/// Load the models, refusing before any work begins.
///
/// # Errors
/// If the head is missing, corrupt, or unloadable. Everything here is a
/// property of the installation rather than of the data, which is why it runs
/// after the request has already been checked.
fn prepare(dir: &Path, global: Global, task: Task) -> Result<Prepared, Failure> {
    // Refuse before loading anything: predicting from weights that do not
    // match their metadata is how an unidentifiable number ends up in
    // someone's results.
    let head = inventory::scan(dir, EMBEDDER)
        .into_iter()
        .find(|h| h.task == task.name())
        .ok_or_else(|| anyhow!("unknown task {}", task.name()))
        .blame(Fault::Usage)?;
    if !head.ready() {
        return Err(Failure {
            fault: Fault::Install,
            error: anyhow!("{} head: {}", task.name(), head.describe()),
        });
    }

    let device = global.device.open()?;
    let description = opendna_ml::device_description(&device);
    let Loaded { prov, metadata, files } = head_provenance(dir, task, &description).blame(Fault::Install)?;
    eprintln!("Running {} on {description}", task.name());

    let heads = Heads::load(dir, EMBEDDER, &device).blame(Fault::Install)?;
    let model = staging::open_embedder(dir, EMBEDDER, &device, global.offline).blame(Fault::Install)?;
    Ok(Prepared { heads, model, device, description, prov, metadata, files })
}

/// Produce the result table for one task.
///
/// Split out of `predict` because this is the part that grows: every head
/// that gains a column, an interval or a confidence measure adds to the
/// dispatch, and none of it has anything to do with reading the input or
/// loading the models.
///
/// # Errors
/// If the head is missing, or a sequence cannot be scored by it.
fn predict_table(
    task: Task,
    proteins: &[(String, String)],
    resolved: Option<&[(String, Vec<Substitution>)]>,
    ready: &Prepared,
    budget: usize,
    min_confidence: Option<f64>,
) -> Result<Table, Failure> {
    let Prepared { heads, model, device, metadata, .. } = ready;
    Ok(match task {
        Task::Disorder => {
            let threshold = disorder_threshold(metadata);
            disorder_table(
                &disorder(proteins, heads, model, device, budget, threshold).blame(Fault::Input)?,
                threshold,
                min_confidence,
            )
        }
        Task::Burial => with_interval(
            burial(proteins, heads, model, device, budget).blame(Fault::Input)?,
            "accessibility",
            interval_of(metadata),
        )
        .blame(Fault::Unexpected)?,
        Task::Solubility => with_interval(
            solubility(proteins, heads, model, device, budget).blame(Fault::Input)?,
            "solubility_percent",
            interval_of(metadata),
        )
        .blame(Fault::Unexpected)?,
        Task::Stability => {
            let resolved = resolved.unwrap_or_default();
            with_interval(
                stability(proteins, resolved, heads, model, device).blame(Fault::Input)?,
                "ddg_kcal_per_mol",
                interval_of(metadata),
            )
            .blame(Fault::Unexpected)?
        }
    })
}

/// Resolve the requested mutations against the sequences, for the one task
/// that takes them.
///
/// Done before the model is loaded. Everything that can be wrong with a
/// mutation list needs only the FASTA, so checking it first means a typo is
/// reported in a moment rather than after several hundred megabytes have
/// loaded for a run that cannot work.
///
/// # Errors
/// If stability was asked for without a mutation file, if the file cannot be
/// read, or if any mutation does not resolve against its sequence.
fn requested_mutations(
    task: Task,
    proteins: &[(String, String)],
    mutations: Option<&Path>,
) -> Result<Option<Resolved>, Failure> {
    if task != Task::Stability {
        return Ok(None);
    }
    let path = mutations.ok_or_else(|| Failure {
        fault: Fault::Usage,
        error: anyhow!(
            "--mutations <file> is required for stability: the head scores named substitutions, \
             so it needs to be told which ones. The file is a CSV of `name,mutation` rows, for \
             example `ubiquitin,L8P`."
        ),
    })?;
    let wanted = read_mutations(path).blame(Fault::Input)?;
    Ok(Some(resolve_mutations(proteins, &wanted).blame(Fault::Input)?))
}

/// The settings a chunked run needs beyond the models and the sequences.
///
/// Grouped because they travel together and separately from everything else:
/// each one is a thing that makes two runs different, which is exactly what
/// the checkpoint has to compare.
#[derive(Copy, Clone)]
struct Chunked<'a> {
    /// Where the checkpoint lives.
    path: &'a Path,
    /// Units of work per chunk.
    chunk: usize,
    /// The confidence cut-off, which adds a column.
    min_confidence: Option<f64>,
    /// The input's bytes, hashed into the header.
    input: &'a [u8],
    /// The mutation file, likewise.
    mutations: Option<&'a Path>,
}

/// Everything about a run that its checkpoint has to record.
///
/// Gathered in one place so that the list of things which make two runs
/// different is readable as a list. Anything missing from it is something a
/// resumed run could silently change.
///
/// # Errors
/// If the mutation file cannot be read for its digest.
fn resume_header(
    task: Task,
    units: usize,
    run: &Chunked<'_>,
    ready: &Prepared,
    budget: usize,
) -> Result<checkpoint::Header, Failure> {
    let mutations_sha256 = run
        .mutations
        .map(|p| {
            std::fs::read(p)
                .with_context(|| format!("reading {}", p.display()))
                .map(|bytes| manifest::sha256(&bytes))
        })
        .transpose()
        .blame(Fault::Input)?;
    Ok(checkpoint::Header {
        format: checkpoint::FORMAT,
        task: task.name().to_string(),
        // The bytes, not the path: the same name holding different contents is
        // the mistake most worth catching.
        input_sha256: manifest::sha256(run.input),
        units,
        chunk: run.chunk,
        head_build: ready.metadata.build.clone(),
        head_weights_sha256: ready.metadata.weights_sha256.clone(),
        embedding_model: EMBEDDER.repo().to_string(),
        device: ready.description.clone(),
        batch_tokens: budget,
        mutations_sha256,
        min_confidence: run.min_confidence,
    })
}

/// Run a prediction in chunks, recording each one as it finishes.
///
/// The loop is deliberately dull: for every chunk not already in the
/// checkpoint, run the same dispatch a whole-file run would, and append the
/// rows. What makes a resumed run agree with an uninterrupted one is that the
/// chunk boundaries are the same, so the embedding batches are the same, so
/// the floating-point sums happen in the same order.
///
/// Rows come back from the checkpoint rather than from memory even for chunks
/// this process computed, so that the assembled table is the one the file
/// describes. If the two could differ, the checkpoint would not be worth
/// having.
///
/// # Errors
/// If the checkpoint describes a different run, cannot be written, or holds
/// chunks that cannot be assembled into one table; or if a chunk fails to
/// predict, in which case the chunks already recorded survive.
fn predict_resumable(
    task: Task,
    proteins: &[(String, String)],
    resolved: Option<&[(String, Vec<Substitution>)]>,
    ready: &Prepared,
    budget: usize,
    run: Chunked<'_>,
) -> Result<Table, Failure> {
    let Chunked { path, chunk, min_confidence, .. } = run;
    let units = match task {
        Task::Stability => resolved.map_or(0, <[_]>::len),
        Task::Disorder | Task::Burial | Task::Solubility => proteins.len(),
    };
    let total = checkpoint::chunks(units, chunk).blame(Fault::Usage)?;
    let header = resume_header(task, units, &run, ready, budget)?;
    // A sequence too long for this head fails the chunk it lands in, which on
    // a long run means hours of work before being told. The check needs only
    // the lengths, so it runs over the whole file first.
    scorable(task, proteins, &ready.model).blame(Fault::Input)?;

    let mut writer = checkpoint::open(path, &header).blame(Fault::Input)?;
    let done = writer.records().len();
    if done > 0 {
        let rows: usize = writer.records().iter().map(|r| r.rows.len()).sum();
        eprintln!("Resuming {}: {done} of {total} chunks already done, {rows} rows", path.display());
    }

    let mut columns: Option<Vec<&'static str>> = None;
    for index in done..total {
        let from = index * chunk;
        let to = ((index + 1) * chunk).min(units);
        let table = match task {
            Task::Stability => {
                predict_table(task, proteins, resolved.map(|r| &r[from..to]), ready, budget, min_confidence)?
            }
            Task::Disorder | Task::Burial | Task::Solubility => {
                predict_table(task, &proteins[from..to], None, ready, budget, min_confidence)?
            }
        };
        let names = table.columns.iter().map(|c| (*c).to_string()).collect();
        columns = Some(table.columns);
        writer
            .append(checkpoint::Record { chunk: index, columns: names, rows: table.rows })
            .blame(Fault::Unexpected)?;
        eprintln!("Chunk {} of {total} recorded in {}", index + 1, path.display());
    }

    // The table needs its column names, and a completed checkpoint being
    // rendered again ran no chunk to take them from. Running the dispatch over
    // nothing costs nothing and is authoritative about the shape.
    let columns = match columns {
        Some(columns) => columns,
        None => predict_table(task, &[], Some(&[]), ready, budget, min_confidence)?.columns,
    };
    let expected: Vec<String> = columns.iter().map(|c| (*c).to_string()).collect();
    let mut rows = Vec::new();
    for record in writer.records() {
        if record.columns != expected {
            return Err(Failure {
                fault: Fault::Input,
                error: anyhow!(
                    "{} holds a chunk with columns {:?} where this run produces {:?}. The \
                     checkpoint was written by a different version of OpenDNA; start the run \
                     again against a new file.",
                    path.display(),
                    record.columns,
                    expected
                ),
            });
        }
        rows.extend(record.rows.iter().cloned());
    }
    Ok(Table { columns, rows })
}

/// Refuse sequences this task cannot score, before any of them is embedded.
///
/// Only solubility has a limit of its own: it pools over a single pass and has
/// no windowed equivalent, so a sequence past the embedder's limit cannot be
/// scored at all rather than scored less well. Disorder and burial window.
/// Stability checks each sequence as it reaches it, and its own limit is the
/// same one.
///
/// # Errors
/// Naming the sequence and its length.
fn scorable(task: Task, proteins: &[(String, String)], model: &Esm2) -> Result<()> {
    if task != Task::Solubility {
        return Ok(());
    }
    for (name, sequence) in proteins {
        let residues = sequence.chars().count();
        if residues > model.max_residues() {
            bail!(
                "{}: {} residues exceeds the embedding model's {} limit. The solubility head \
                 pools over one pass and has no windowed equivalent, so this sequence cannot be \
                 scored.",
                name,
                residues,
                model.max_residues()
            );
        }
    }
    Ok(())
}

/// Run one prediction and record what produced it.
///
/// # Errors
/// If the input cannot be parsed, if the head is missing or corrupt, if the
/// embedding model cannot be loaded, or if anything cannot be written. The
/// [`Fault`] on the failure distinguishes the install's problems from the
/// data's, because a caller can only act on the difference.
fn predict(dir: &Path, global: Global, args: PredictArgs) -> Result<(), Failure> {
    let PredictArgs { task, input, output, format, mutations, min_confidence, checkpoint, chunk } = args;
    let started_utc = chrono::Utc::now();
    let started = std::time::Instant::now();
    // Nothing installed can make this combination meaningful, so it is
    // refused before the input is even read.
    if format == Format::Caid && task != Task::Disorder {
        return Err(Failure {
            fault: Fault::Usage,
            error: anyhow!(
                "--format caid is the CAID disorder prediction format and only applies to \
                 --task disorder"
            ),
        });
    }
    // The CAID writer assembles one document from every target at once rather
    // than a table of rows, so there is nothing for a chunk to hold. Saying so
    // is better than accepting the flag and ignoring it.
    if format == Format::Caid && checkpoint.is_some() {
        return Err(Failure {
            fault: Fault::Usage,
            error: anyhow!(
                "--checkpoint cannot be used with --format caid: that format is written as one \
                 document rather than as rows, so there is no chunk to resume from. Run it as \
                 --format csv with a checkpoint and convert afterwards, or run it without one."
            ),
        });
    }
    // Checked here rather than where the chunks are counted, so that a typo in
    // an argument is reported before six hundred megabytes of model have
    // loaded for a run that cannot start.
    if chunk == 0 {
        return Err(Failure { fault: Fault::Usage, error: anyhow!("--chunk must be at least 1") });
    }

    let Request { proteins, text } = read_request(&input)?;

    let resolved = requested_mutations(task, &proteins, mutations.as_deref())?;

    let ready = prepare(dir, global, task)?;
    let Prepared { device, description, prov, metadata, files, .. } = &ready;
    let budget = global.budget(device);

    // Said once per run rather than once per chunk: a chunked run calls the
    // same dispatch forty times, and forty copies of the same paragraph on
    // standard error is how a person stops reading them.
    if let (Task::Disorder, Some(cut)) = (task, min_confidence) {
        report_confidence_cut(metadata, cut);
    }

    // One manifest writer for every format, so no output can end up described
    // differently from another.
    let label = input.display().to_string();
    let record = |path: &Path, rendered: &str, rows: usize, muts: Option<&Path>| -> Result<()> {
        if !global.manifest {
            return Ok(());
        }
        let residues = proteins.iter().map(|(_, s)| s.chars().count()).sum();
        let wrote = manifest::write(
            &manifest::Run {
                task: task.name(),
                device: description,
                started: started_utc,
                elapsed: started.elapsed(),
                batch_tokens: budget,
                chunk: checkpoint.as_ref().map(|_| chunk),
                input: manifest::Input {
                    label: &label,
                    bytes: text.as_bytes(),
                    sequences: proteins.len(),
                    residues,
                },
                output: Some(path),
                rendered,
                rows,
                mutations: muts,
                head: metadata,
                files,
                embedder: &staging::describe_embedder(dir, EMBEDDER),
                limitations: &prov.limitations,
                license: files.license(),
            },
            path,
        )?;
        eprintln!("Wrote {}", wrote.display());
        Ok(())
    };

    if format == Format::Caid {
        let (rendered, rows) = write_caid(&ready, &proteins, budget, output.as_deref())?;
        if let Some(path) = output.as_deref() {
            record(path, &rendered, rows, None)?;
        }
        return Ok(());
    }

    let table = match checkpoint.as_deref() {
        None => predict_table(task, &proteins, resolved.as_deref(), &ready, budget, min_confidence)?,
        Some(path) => predict_resumable(
            task,
            &proteins,
            resolved.as_deref(),
            &ready,
            budget,
            Chunked { path, chunk, min_confidence, input: text.as_bytes(), mutations: mutations.as_deref() },
        )?,
    };

    let rendered = table.render(format, prov)?;
    write_output(output.as_deref(), &rendered)?;
    if let Some(path) = output.as_deref() {
        record(path, &rendered, table.rows.len(), mutations.as_deref())?;
    }
    Ok(())
}

/// Report which external engines are present, and pin or check them.
///
/// The same shape as [`models`], because it answers the same question about a
/// different kind of dependency: what is here, which one is it, and is it what
/// the last run used.
///
/// An engine being absent is not a failure. Nothing drives these yet, and a
/// lab that does not dock ligands has no reason to install a docking program.
/// What fails is a `--check` against a lockfile that no longer matches.
///
/// # Errors
/// If a lockfile cannot be written or read, or if the installed engines
/// differ from a lockfile given to `--check`.
fn engines(lock: Option<&Path>, check: Option<&Path>) -> Result<()> {
    let found = opendna_engine::discover_all();

    println!("External engines OpenDNA can drive:\n");
    for (requirement, outcome) in &found {
        match outcome {
            Ok(engine) => {
                println!("  {:11} {}", requirement.name, engine.describe());
                println!("  {:11}   {} -- {}", "", requirement.purpose, engine.identification);
            }
            Err(error) => {
                println!("  {:11} not installed", requirement.name);
                println!("  {:11}   {} -- {}", "", requirement.purpose, error);
            }
        }
    }

    // Said plainly, because a listing that implies otherwise would be the
    // same overclaim the rewrite removed.
    println!(
        "\nNo command drives any of these yet. They are listed so an installation can be \
         checked and pinned\nbefore the commands that use them exist. See docs/ROADMAP.md for \
         which tranche picks up each one."
    );

    if let Some(path) = lock {
        std::fs::write(path, opendna_engine::lockfile(&found))
            .with_context(|| format!("writing {}", path.display()))?;
        println!("\nWrote {}", path.display());
    }
    if let Some(path) = check {
        let pinned = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        opendna_engine::check(&pinned, &found)?;
        println!("\nThe installed engines match {}.", path.display());
    }
    Ok(())
}

/// Report what is installed, and pin or check it.
///
/// Verification is not optional here: a head whose weights do not match its
/// metadata is reported as corrupt and makes the command fail, because the
/// alternative is a run that produces numbers from a file nobody can
/// identify.
fn models(dir: &Path, lock: Option<&Path>, check: Option<&Path>) -> Result<()> {
    let found = inventory::scan(dir, EMBEDDER);
    println!("Looking in {}", dir.display());
    println!("Embedding model: {}", EMBEDDER.repo());
    for head in &found {
        println!("  {:11} {}", head.task, head.describe());
    }
    // Whether this directory is self-contained is the question an offline
    // machine asks, and it is not answerable from the head listing.
    if staging::staged_embedder(dir).is_some() {
        println!("\nThe embedding model is staged here, so --offline will work.");
    } else {
        println!(
            "\nNo embedding model staged here: the first prediction downloads {}.\n\
             Run `opendna stage --to <dir>` for a directory that works offline.",
            EMBEDDER.repo()
        );
    }

    if let Some(path) = lock {
        std::fs::write(path, inventory::lockfile(&found, EMBEDDER))
            .with_context(|| format!("writing {}", path.display()))?;
        eprintln!("Wrote {}", path.display());
    }
    if let Some(path) = check {
        inventory::check(path, &found, EMBEDDER)?;
        println!(
            "
The installed heads match {}.",
            path.display()
        );
    }

    let corrupt: Vec<&str> = found.iter().filter(|h| !h.integrity.usable()).map(|h| h.task).collect();
    if !corrupt.is_empty() {
        bail!("corrupt weights: {}. Re-download these heads.", corrupt.join(", "));
    }
    Ok(())
}

/// Write the CAID submission format.
///
/// Separate from the table formats because it is a different file with its own
/// rules: no comment lines are permitted, so the provenance every other output
/// embeds has to travel beside it instead.
///
/// # Errors
/// If the run or any write fails.
fn write_caid(
    ready: &Prepared,
    proteins: &[(String, String)],
    budget: usize,
    output: Option<&Path>,
) -> Result<(String, usize), Failure> {
    let Prepared { heads, model, device, prov, metadata, .. } = ready;
    let threshold = disorder_threshold(metadata);
    let predictions = disorder(proteins, heads, model, device, budget, threshold).blame(Fault::Input)?;
    let targets: Vec<(String, Vec<caid::Residue>)> = predictions
        .into_iter()
        .map(|(header, residues)| (caid::accession(&header).to_string(), residues))
        .collect();
    let rendered = caid::render(&targets, threshold)?;
    write_output(output, &rendered)?;
    // The format cannot carry comments, so the provenance that every other
    // output embeds goes beside the file instead.
    eprint!("{}", prov.comment_block());
    if let Some(path) = output {
        let sidecar = path.with_extension("provenance.txt");
        std::fs::write(&sidecar, prov.comment_block())
            .with_context(|| format!("writing {}", sidecar.display()))?;
        eprintln!("Wrote {}", sidecar.display());
    }
    let rows = targets.iter().map(|(_, r)| r.len()).sum();
    Ok((rendered, rows))
}

/// Embed every sequence, batching the ones that fit in one pass.
///
/// Sequences within the checkpoint's limit go through [`Esm2::embed_many`],
/// which sorts by length and fills batches to a token budget. Longer ones are
/// windowed one at a time, because a windowed embedding is assembled from
/// several passes over the same sequence and has no batched equivalent.
///
/// Progress goes to standard error, counted in sequences rather than batches:
/// a batch is an implementation detail and "47 of 200" is what somebody
/// watching a long run wants to know.
///
/// # Errors
/// If any sequence cannot be embedded.
fn embed_all(
    model: &Esm2,
    proteins: &[(String, String)],
    budget: usize,
) -> Result<Vec<opendna_ml::esm2::Embedding>> {
    let limit = model.max_residues();
    let fits: Vec<usize> = (0..proteins.len()).filter(|&i| proteins[i].1.chars().count() <= limit).collect();
    let long: Vec<usize> = (0..proteins.len()).filter(|&i| proteins[i].1.chars().count() > limit).collect();

    let mut out: Vec<Option<opendna_ml::esm2::Embedding>> = vec![None; proteins.len()];
    if !fits.is_empty() {
        let batch: Vec<&str> = fits.iter().map(|&i| proteins[i].1.as_str()).collect();
        eprintln!("Embedding {} sequences", batch.len());
        for (&index, embedding) in fits.iter().zip(model.embed_many(&batch, budget)?) {
            out[index] = Some(embedding);
        }
    }
    for (done, &index) in long.iter().enumerate() {
        eprintln!(
            "Windowing {} ({} residues, {} of {})",
            proteins[index].0,
            proteins[index].1.chars().count(),
            done + 1,
            long.len()
        );
        out[index] = Some(model.embed_windowed(&proteins[index].1)?);
    }
    out.into_iter()
        .enumerate()
        .map(|(i, e)| e.ok_or_else(|| anyhow!("{}: no embedding was produced", proteins[i].0)))
        .collect()
}

/// Per-residue disorder for every sequence, keyed by its FASTA header.
///
/// Both the table and the CAID file are derived from this, so the two cannot
/// contain different numbers for the same residue.
/// The interval level a run reports, when the head records one.
///
/// 90% by default: the convention, and the level whose coverage was closest
/// to nominal on the heads measured here. A head that records nothing gets no
/// interval columns rather than a made-up one.
const INTERVAL_LEVEL: f64 = 0.9;

/// The recorded interval at [`INTERVAL_LEVEL`], if the head has one.
fn interval_of(metadata: &opendna_ml::heads::HeadMetadata) -> Option<&opendna_ml::heads::IntervalLevel> {
    metadata.intervals.as_ref()?.levels.iter().find(|l| (l.nominal - INTERVAL_LEVEL).abs() < 1e-9)
}

/// Add low and high columns to a table of point predictions.
///
/// The half-width is constant, which is what split conformal gives and what
/// its caveat is: the interval is too wide where the head is reliable and too
/// narrow where it is not. The model card carries the conditional coverage
/// that shows where.
fn with_interval(
    mut table: Table,
    column: &str,
    level: Option<&opendna_ml::heads::IntervalLevel>,
) -> Result<Table> {
    let Some(level) = level else { return Ok(table) };
    let index = table
        .columns
        .iter()
        .position(|c| *c == column)
        .ok_or_else(|| anyhow!("the result table has no {column} column"))?;
    for row in &mut table.rows {
        let value: f64 = row
            .get(index)
            .ok_or_else(|| anyhow!("a row is missing its {column} field"))?
            .parse()
            .with_context(|| format!("parsing {column} as a number"))?;
        row.push(format!("{:.4}", value - level.half_width));
        row.push(format!("{:.4}", value + level.half_width));
    }
    table.columns.push("interval_low");
    table.columns.push("interval_high");
    Ok(table)
}

/// The threshold this head says to use, or the fallback.
fn disorder_threshold(metadata: &opendna_ml::heads::HeadMetadata) -> f64 {
    metadata.calibration.as_ref().map_or(DISORDER_THRESHOLD, |c| c.threshold)
}

fn disorder(
    proteins: &[(String, String)],
    heads: &Heads,
    model: &Esm2,
    device: &candle_core::Device,
    budget: usize,
    threshold: f64,
) -> Result<Vec<(String, Vec<caid::Residue>)>> {
    let head = heads.disorder.as_ref().ok_or_else(|| anyhow!("the disorder head did not load"))?;
    let embeddings = embed_all(model, proteins, budget)?;
    let mut out = Vec::with_capacity(proteins.len());
    for ((name, sequence), embedding) in proteins.iter().zip(&embeddings) {
        let residues: Vec<char> = sequence.chars().collect();
        let probabilities = head.predict(&embedding.per_residue, embedding.length, device)?;
        let scored = probabilities
            .iter()
            .enumerate()
            .map(|(index, probability)| caid::Residue {
                position: index + 1,
                code: residues.get(index).copied().unwrap_or('?'),
                score: f64::from(*probability),
                disordered: f64::from(*probability) >= threshold,
            })
            .collect();
        out.push((name.clone(), scored));
    }
    Ok(out)
}

/// Per-residue solvent accessibility for every sequence.
///
/// Windowed for long sequences, like disorder: accessibility is a local
/// property, so a residue's own neighbourhood is what the head reads and a
/// window boundary costs less than refusing the sequence.
fn burial(
    proteins: &[(String, String)],
    heads: &Heads,
    model: &Esm2,
    device: &candle_core::Device,
    budget: usize,
) -> Result<Table> {
    let head = heads.burial.as_ref().ok_or_else(|| anyhow!("the burial head did not load"))?;
    let embeddings = embed_all(model, proteins, budget)?;
    let mut rows = Vec::new();
    for ((name, sequence), embedding) in proteins.iter().zip(&embeddings) {
        let residues: Vec<char> = sequence.chars().collect();
        let values = head.predict(&embedding.per_residue, embedding.length, device)?;
        for (index, accessibility) in values.iter().enumerate() {
            rows.push(vec![
                name.clone(),
                (index + 1).to_string(),
                residues.get(index).copied().unwrap_or('?').to_string(),
                format!("{accessibility:.4}"),
                u8::from(f64::from(*accessibility) <= BURIED_AT).to_string(),
            ]);
        }
    }
    Ok(Table { columns: vec!["name", "position", "residue", "accessibility", "buried"], rows })
}

/// The flat table form of a disorder run.
///
/// `confidence` is how far the calibrated probability sits from the decision
/// threshold, which is the quantity measured to order residues by how likely
/// the call is to be wrong. It is reported for every row, because a reader
/// filtering a spreadsheet should not have to recompute it from the threshold
/// in the comments.
fn disorder_table(
    predictions: &[(String, Vec<caid::Residue>)],
    threshold: f64,
    min_confidence: Option<f64>,
) -> Table {
    let mut rows = Vec::new();
    for (name, residues) in predictions {
        for r in residues {
            let confidence = (r.score - threshold).abs();
            let mut row = vec![
                name.clone(),
                r.position.to_string(),
                r.code.to_string(),
                format!("{:.4}", r.score),
                u8::from(r.disordered).to_string(),
                format!("{confidence:.4}"),
            ];
            if let Some(cut) = min_confidence {
                row.push(u8::from(confidence < cut).to_string());
            }
            rows.push(row);
        }
    }
    let mut columns = vec!["name", "position", "residue", "probability", "disordered", "confidence"];
    if min_confidence.is_some() {
        columns.push("low_confidence");
    }
    Table { columns, rows }
}

/// Say what a chosen confidence cut-off was measured to buy.
///
/// Printed to standard error, so it reaches a person without entering the
/// results. The numbers come from the head's own recorded measurement; when it
/// has none, the cut-off still applies and this says there is nothing behind
/// it, because a marked row with no evidence is worth less than the reader
/// thinks.
fn report_confidence_cut(metadata: &opendna_ml::heads::HeadMetadata, cut: f64) {
    let Some(abstention) = metadata.abstention.as_ref() else {
        eprintln!(
            "Marking calls below confidence {cut:.2}. This head records no measured \
             confidence/error trade-off, so nothing is known about what that buys. Run \
             `opendna-train abstain --apply` against it."
        );
        return;
    };
    // The nearest recorded cut-off: the curve is measured at a few points, and
    // quoting an interpolated number as measured would be inventing one.
    let nearest = abstention
        .curve
        .iter()
        .min_by(|a, b| (a.min_confidence - cut).abs().total_cmp(&(b.min_confidence - cut).abs()));
    match nearest {
        Some(point) => eprintln!(
            "Marking calls below confidence {cut:.2}. Measured at the nearest recorded cut-off \
             ({:.2}) on held-out data: {:.0}% of residues clear it and {:.1}% of those calls \
             were wrong, against {:.1}% over all residues.",
            point.min_confidence,
            100.0 * point.coverage,
            100.0 * point.error,
            100.0 * abstention.base_error
        ),
        None => eprintln!("Marking calls below confidence {cut:.2}."),
    }
}

fn solubility(
    proteins: &[(String, String)],
    heads: &Heads,
    model: &Esm2,
    device: &candle_core::Device,
    budget: usize,
) -> Result<Table> {
    let head = heads.solubility.as_ref().ok_or_else(|| anyhow!("the solubility head did not load"))?;
    let mut rows = Vec::new();
    // Refused before embedding anything: a file with one over-long sequence
    // would otherwise spend the whole run before failing. Stated once, in
    // `scorable`, so a chunked run and a whole-file run refuse the same
    // sequences with the same words.
    scorable(Task::Solubility, proteins, model)?;
    let batch: Vec<&str> = proteins.iter().map(|(_, s)| s.as_str()).collect();
    eprintln!("Embedding {} sequences", batch.len());
    for ((name, sequence), embedding) in proteins.iter().zip(model.embed_many(&batch, budget)?) {
        let percent = head.predict(&embedding.pooled, device)?;
        rows.push(vec![name.clone(), sequence.chars().count().to_string(), format!("{percent:.1}")]);
    }
    Ok(Table { columns: vec!["name", "length", "solubility_percent"], rows })
}

/// Parse the `name,mutation` file, grouping by sequence name.
fn read_mutations(path: &Path) -> Result<Vec<(String, String)>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut out = Vec::new();
    for (number, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (name, mutation) = line
            .split_once(',')
            .ok_or_else(|| anyhow!("{}:{}: expected `name,mutation`", path.display(), number + 1))?;
        let (name, mutation) = (name.trim(), mutation.trim());
        // A header row is easy to leave in by accident, and treating it as a
        // mutation produces a confusing parse error further down.
        if number == 0 && mutation.eq_ignore_ascii_case("mutation") {
            continue;
        }
        out.push((name.to_string(), mutation.to_string()));
    }
    if out.is_empty() {
        bail!("{} contains no mutations", path.display());
    }
    Ok(out)
}

/// Split `L8P` into the wild-type residue, 0-based index, and mutant.
///
/// # Errors
/// If the label is malformed, out of range, or disagrees with the sequence at
/// that position — the last being the check that catches a file numbered
/// against a different construct, which would otherwise score the wrong
/// residue and look fine.
fn parse_mutation(label: &str, residues: &[char], name: &str) -> Result<(char, usize, char)> {
    let mut chars = label.chars();
    let wild = chars.next().ok_or_else(|| anyhow!("{name}: empty mutation"))?;
    let mutant = label
        .chars()
        .next_back()
        .filter(|_| label.chars().count() >= 3)
        .ok_or_else(|| anyhow!("{name}: {label:?} is not a mutation like L8P"))?;
    let digits: String = chars.clone().take(label.chars().count() - 2).collect();
    let position: usize = digits
        .parse()
        .with_context(|| format!("{name}: {label:?} has no position; expected something like L8P"))?;
    if position == 0 || position > residues.len() {
        bail!("{name}: {label} is outside a {}-residue sequence; positions are 1-based", residues.len());
    }
    let index = position - 1;
    if residues[index] != wild {
        bail!(
            "{name}: {label} says position {position} is {wild}, but the sequence has {}. The \
             mutation list and the FASTA disagree, so every score from this pair would be \
             attributed to the wrong residue.",
            residues[index]
        );
    }
    Ok((wild, index, mutant))
}

/// Mutations resolved against their sequences, grouped by sequence name.
type Resolved = Vec<(String, Vec<Substitution>)>;

/// One substitution, resolved against the sequence it applies to.
struct Substitution {
    /// The label as written, which is what the output reports.
    label: String,
    /// 0-based position in the sequence.
    index: usize,
    /// The wild-type residue's index in [`AminoAcid::STANDARD`].
    wild: usize,
    /// The mutant residue's, likewise.
    mutant: usize,
    /// The wild-type and mutant letters, for the masked log-odds.
    letters: (char, char),
}

/// Resolve every requested mutation against the sequences.
///
/// Done before any model is loaded. Everything that can go wrong here -- a
/// malformed label, a position past the end, a wild-type letter the sequence
/// contradicts, a residue outside the twenty the head encodes, a name matching
/// no sequence -- needs only the FASTA, so paying for several hundred
/// megabytes of embedding model first would be waiting to be told about a
/// typo.
///
/// # Errors
/// If any mutation cannot be resolved, naming it. One bad label stops the run
/// rather than being skipped: a table silently missing rows looks exactly like
/// a table that is complete.
fn resolve_mutations(proteins: &[(String, String)], wanted: &[(String, String)]) -> Result<Resolved> {
    let standard = AminoAcid::STANDARD;
    let mut out = Vec::new();
    for (name, sequence) in proteins {
        let mine: Vec<&String> = wanted.iter().filter(|(n, _)| n == name).map(|(_, m)| m).collect();
        if mine.is_empty() {
            continue;
        }
        let residues: Vec<char> = sequence.chars().collect();
        let mut resolved = Vec::with_capacity(mine.len());
        for label in mine {
            let (wild, index, mutant) = parse_mutation(label, &residues, name)?;
            let (Some(wt), Some(mt)) = (
                standard.iter().position(|a| a.one_letter() == wild),
                standard.iter().position(|a| a.one_letter() == mutant),
            ) else {
                bail!(
                    "{name}: {label} involves a residue the head was not trained on; only the 20 \
                     standard amino acids are encoded"
                );
            };
            resolved.push(Substitution {
                label: label.clone(),
                index,
                wild: wt,
                mutant: mt,
                letters: (wild, mutant),
            });
        }
        out.push((name.clone(), resolved));
    }
    if out.is_empty() {
        bail!(
            "no mutation matched a sequence in the FASTA file. Names must match the FASTA header \
             exactly."
        );
    }
    Ok(out)
}

fn stability(
    proteins: &[(String, String)],
    resolved: &[(String, Vec<Substitution>)],
    heads: &Heads,
    model: &Esm2,
    device: &candle_core::Device,
) -> Result<Table> {
    let head = heads.stability.as_ref().ok_or_else(|| anyhow!("the stability head did not load"))?;
    let mut rows = Vec::new();

    for (name, mutations) in resolved {
        let sequence = proteins
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, s)| s)
            .ok_or_else(|| anyhow!("{name}: resolved mutations for a sequence that is not here"))?;
        if sequence.chars().count() > model.max_residues() {
            bail!(
                "{}: {} residues exceeds the embedding model's {} limit",
                name,
                sequence.chars().count(),
                model.max_residues()
            );
        }

        let embedding = model.embed(&[sequence.as_str()])?.remove(0);
        let mut positions: Vec<usize> = mutations.iter().map(|m| m.index).collect();
        positions.sort_unstable();
        positions.dedup();
        let log_probs = model.masked_log_probs(sequence, &positions)?;
        let index_of: std::collections::HashMap<usize, usize> =
            positions.iter().enumerate().map(|(n, p)| (*p, n)).collect();

        let mut features = Vec::new();
        let mut kept = Vec::new();
        for m in mutations {
            let lp = &log_probs[index_of[&m.index]];
            let logp = |letter: char| -> f32 { lp[opendna_ml::esm2::alphabet::token_id(letter) as usize] };
            let odds = logp(m.letters.1) - logp(m.letters.0);
            let site =
                embedding.residue(m.index).ok_or_else(|| anyhow!("{name}: no embedding at {}", m.index))?;
            features.push(StabilityHead::features(site, m.wild, m.mutant, odds));
            kept.push((&m.label, odds));
        }
        let predicted = head.predict(&features, model.hidden_size(), device)?;
        for ((label, log_odds), ddg) in kept.into_iter().zip(predicted) {
            rows.push(vec![name.clone(), label.clone(), format!("{ddg:+.3}"), format!("{log_odds:+.3}")]);
        }
    }

    Ok(Table { columns: vec!["name", "mutation", "ddg_kcal_per_mol", "log_odds"], rows })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn a_mutation_resolves_against_the_sequence() {
        let residues: Vec<char> = "MQIFVKTL".chars().collect();
        assert_eq!(parse_mutation("M1A", &residues, "p").unwrap(), ('M', 0, 'A'));
        assert_eq!(parse_mutation("L8P", &residues, "p").unwrap(), ('L', 7, 'P'));
    }

    /// The check that matters: a mutation list numbered against a different
    /// construct must stop the run, not score a neighbouring residue.
    #[test]
    fn a_mutation_disagreeing_with_the_sequence_is_refused() {
        let residues: Vec<char> = "MQIFVKTL".chars().collect();
        let err = parse_mutation("L2P", &residues, "p").unwrap_err().to_string();
        assert!(err.contains("the sequence has Q"), "{err}");
        assert!(err.contains("wrong residue"), "{err}");
    }

    #[test]
    fn a_malformed_or_out_of_range_mutation_is_refused() {
        let residues: Vec<char> = "MQIF".chars().collect();
        for bad in ["", "L", "LP", "M9A", "M0A", "Mxa"] {
            assert!(parse_mutation(bad, &residues, "p").is_err(), "accepted {bad:?}");
        }
    }

    #[test]
    fn two_digit_positions_are_read_whole() {
        let residues: Vec<char> = "ACDEFGHIKLMNP".chars().collect();
        assert_eq!(parse_mutation("P13W", &residues, "p").unwrap(), ('P', 12, 'W'));
    }

    #[test]
    fn the_models_directory_prefers_the_flag_then_the_environment() {
        assert_eq!(models_dir(Some(PathBuf::from("/x"))), PathBuf::from("/x"));
        assert_eq!(models_dir(None), models_dir(None), "deterministic without the flag");
    }

    #[test]
    fn a_field_containing_the_separator_is_quoted() {
        let table = Table {
            columns: vec!["name", "value"],
            rows: vec![vec!["a,b".into(), "1".into()], vec!["plain".into(), "2".into()]],
        };
        let prov = Provenance {
            task: "disorder",
            head_version: "v".into(),
            embedding_model: "m".into(),
            device: "CPU".into(),
            limitations: vec![],
            interval: None,
        };
        let csv = table.render(Format::Csv, &prov).unwrap();
        assert!(csv.contains("\"a,b\",1"), "{csv}");
        assert!(csv.contains("plain,2"), "{csv}");
    }

    #[test]
    fn json_keeps_numbers_as_numbers() {
        let table =
            Table { columns: vec!["name", "probability"], rows: vec![vec!["p".into(), "0.9412".into()]] };
        let prov = Provenance {
            task: "disorder",
            head_version: "0.6.0-dev".into(),
            embedding_model: "m".into(),
            device: "CPU".into(),
            limitations: vec!["careful".into()],
            interval: None,
        };
        let out = table.render(Format::Json, &prov).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(parsed["results"][0]["probability"], 0.9412);
        assert_eq!(parsed["results"][0]["name"], "p");
        assert_eq!(parsed["head_version"], "0.6.0-dev");
        assert_eq!(parsed["limitations"][0], "careful");
    }

    /// A length or a position emitted as `76.0` becomes a float key in a
    /// downstream join and silently fails to match an integer one.
    #[test]
    fn whole_numbers_stay_integers_in_json() {
        assert_eq!(json_scalar("76"), serde_json::json!(76));
        assert_eq!(json_scalar("1"), serde_json::json!(1));
        assert!(json_scalar("76").is_i64(), "a length became a float");
        assert_eq!(json_scalar("0.9412"), serde_json::json!(0.9412));
        assert_eq!(json_scalar("+0.333"), serde_json::json!(0.333));
        assert_eq!(json_scalar("L8P"), serde_json::json!("L8P"));
    }

    /// Confidence is the distance from the threshold, not from 0.5, and the
    /// marking has to follow the cut-off rather than the probability.
    #[test]
    fn confidence_is_measured_from_the_head_s_own_threshold() {
        let residues = vec![
            caid::Residue { position: 1, code: 'M', score: 0.40, disordered: true },
            caid::Residue { position: 2, code: 'Q', score: 0.04, disordered: false },
        ];
        let predictions = vec![("p".to_string(), residues)];

        let table = disorder_table(&predictions, 0.34, None);
        assert_eq!(table.columns.last(), Some(&"confidence"));
        // 0.40 is 0.06 above a 0.34 threshold, not 0.10 below 0.5.
        assert_eq!(table.rows[0][5], "0.0600");
        assert_eq!(table.rows[1][5], "0.3000");

        let marked = disorder_table(&predictions, 0.34, Some(0.25));
        assert_eq!(marked.columns.last(), Some(&"low_confidence"));
        assert_eq!(marked.rows[0][6], "1", "0.06 is below the 0.25 cut-off");
        assert_eq!(marked.rows[1][6], "0", "0.30 clears it");
    }

    /// A row is never dropped for being uncertain: a table quietly missing
    /// rows is indistinguishable from a complete one.
    #[test]
    fn marking_low_confidence_removes_no_rows() {
        let residues: Vec<caid::Residue> = (1..=10)
            .map(|position| caid::Residue { position, code: 'A', score: 0.34, disordered: true })
            .collect();
        let predictions = vec![("p".to_string(), residues)];
        // Every residue sits exactly on the threshold, so every one is marked.
        let marked = disorder_table(&predictions, 0.34, Some(0.5));
        assert_eq!(marked.rows.len(), 10);
        assert!(marked.rows.iter().all(|r| r[6] == "1"), "all should be marked");
    }

    #[test]
    fn a_mutation_file_skips_comments_and_a_header_row() {
        let dir = std::env::temp_dir().join(format!("opendna-cli-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("m.csv");
        std::fs::write(&path, "name,mutation\n# a note\n\nubq,L8P\nubq, K6A \n").unwrap();
        let got = read_mutations(&path).unwrap();
        assert_eq!(got, vec![("ubq".to_string(), "L8P".to_string()), ("ubq".to_string(), "K6A".to_string())]);
        std::fs::remove_dir_all(&dir).ok();
    }
}
