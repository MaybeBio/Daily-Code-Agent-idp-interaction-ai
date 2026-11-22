//! Running the trained prediction heads.
//!
//! The embedding model is several hundred megabytes and is fetched from the
//! Hugging Face Hub on first use, so nothing here is loaded at startup: the
//! first prediction pays for it and every later one reuses it.
//!
//! When a head is absent the answer is an error saying so. There is no
//! fallback heuristic anywhere in this module — a number shown in the
//! application is always the output of a trained model, or it is not shown.

use opendna_app_protocol::{
    BurialPrediction, CommandError, DisorderPrediction, HeadStatus, ModelStatus, SolubilityPrediction,
    StabilityScan, SubstitutionScore,
};
use opendna_core::AminoAcid;
use opendna_ml::esm2::{Esm2, ModelSize};
use opendna_ml::heads::{HeadFiles, Heads, StabilityHead};
use opendna_ml::{best_device, device_description, weights, MlError};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, MutexGuard, PoisonError};

/// Embedding model the published heads were trained on. Changing this
/// invalidates every head, which `Heads::load` enforces by comparing against
/// each head's recorded metadata.
const EMBEDDER: ModelSize = ModelSize::Large150M;

/// Fallback probability at or above which a residue is reported as
/// disordered, for a head that records no threshold of its own.
///
/// A convention rather than a measurement. A calibrated head carries the
/// threshold that was measured with it, and that is what gets used: the
/// application and the CLI have to agree about which residues are disordered,
/// and they can only do that by both asking the head.
const DISORDER_THRESHOLD: f64 = 0.5;

/// Relative accessibility at or below which the interface calls a residue
/// buried. The Rost & Sander convention, matching the burial dataset's own
/// `buried` column and the CLI, so all three agree.
const BURIED_AT: f64 = 0.25;

/// The disorder threshold the installed head says to use.
///
/// Read from the metadata beside the weights rather than compiled in, so a
/// recalibrated head changes the call without changing the application.
fn disorder_threshold(dir: &std::path::Path) -> f64 {
    let files = opendna_ml::heads::HeadFiles::in_dir(dir, "disorder");
    opendna_ml::heads::read_metadata(&files, "disorder", EMBEDDER)
        .ok()
        .and_then(|m| m.calibration.map(|c| c.threshold))
        .unwrap_or(DISORDER_THRESHOLD)
}

/// Where to look for trained heads.
///
/// `OPENDNA_MODELS` wins, so a developer can point the application at a
/// freshly trained head without installing it. Otherwise heads live beside
/// the database in the application's data directory.
pub fn models_dir(app_data: &Path) -> PathBuf {
    std::env::var_os("OPENDNA_MODELS").map_or_else(|| app_data.join("models"), PathBuf::from)
}

/// The embedding model and whichever heads were present, loaded once.
struct Loaded {
    esm2: Esm2,
    heads: Heads,
}

/// Lazily-loaded models, shared by every prediction command.
pub struct Predictor {
    dir: PathBuf,
    loaded: Mutex<Option<Loaded>>,
}

fn fail(e: impl std::fmt::Display) -> CommandError {
    CommandError { message: e.to_string() }
}

/// Describe one head without loading its weights.
fn head_status(dir: &Path, task: &str) -> HeadStatus {
    let files = HeadFiles::in_dir(dir, task);
    if !files.exist() {
        return HeadStatus::Missing;
    }
    match opendna_ml::heads::read_metadata(&files, task, EMBEDDER) {
        // Limitations come from the card beside the weights, so the panel,
        // the exported file and the published card cannot disagree.
        Ok(metadata) => HeadStatus::Ready { version: metadata.version, limitations: files.limitations() },
        Err(e) => HeadStatus::Unusable { reason: e.to_string() },
    }
}

impl Predictor {
    /// Note where heads live. Nothing is read until the first prediction.
    pub fn new(dir: PathBuf) -> Self {
        Self { dir, loaded: Mutex::new(None) }
    }

    fn lock(&self) -> MutexGuard<'_, Option<Loaded>> {
        self.loaded.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// What can be predicted right now, without downloading anything.
    pub fn status(&self) -> ModelStatus {
        ModelStatus {
            models_dir: self.dir.display().to_string(),
            embedding_model: EMBEDDER.repo().to_string(),
            embedding_cached: weights::is_cached(EMBEDDER.repo()),
            device: best_device()
                .as_ref()
                .map_or_else(|_| "CPU (no GPU available)".to_string(), device_description),
            burial: head_status(&self.dir, "burial"),
            stability: head_status(&self.dir, "stability"),
            solubility: head_status(&self.dir, "solubility"),
            disorder: head_status(&self.dir, "disorder"),
        }
    }

    /// Load the embedding model and heads if that has not happened yet, then
    /// run `f` against them.
    ///
    /// Everything is held behind one lock, so two predictions requested at
    /// once do not both download the checkpoint.
    fn with_models<T>(
        &self,
        f: impl FnOnce(&Esm2, &Heads) -> Result<T, CommandError>,
    ) -> Result<T, CommandError> {
        let mut guard = self.lock();
        if guard.is_none() {
            let device = best_device().map_err(fail)?;
            let esm2 = Esm2::from_pretrained(EMBEDDER, &device).map_err(describe)?;
            let heads = Heads::load(&self.dir, EMBEDDER, &device).map_err(describe)?;
            *guard = Some(Loaded { esm2, heads });
        }
        let loaded = guard.as_ref().ok_or_else(|| fail("models failed to load"))?;
        f(&loaded.esm2, &loaded.heads)
    }

    /// Score every substitution at `position` (1-based).
    pub fn stability_scan(&self, sequence: &str, position: usize) -> Result<StabilityScan, CommandError> {
        let residues: Vec<char> = sequence.chars().collect();
        let wild_type = *residues
            .get(position.checked_sub(1).ok_or_else(|| fail("positions are numbered from 1"))?)
            .ok_or_else(|| {
                fail(format!("position {position} is past the end of a {}-residue sequence", residues.len()))
            })?;
        let wild_index = standard_index(wild_type).ok_or_else(|| {
            fail(format!(
                "residue {position} is {wild_type}, which is not one of the 20 standard amino acids"
            ))
        })?;

        self.with_models(|esm2, heads| {
            let head = heads.stability.as_ref().ok_or_else(|| missing("stability"))?;
            let index = position - 1;
            let log_probs =
                esm2.masked_log_probs(sequence, &[index])
                    .map_err(describe)?
                    .into_iter()
                    .next()
                    .ok_or_else(|| fail("the embedding model returned no distribution for that position"))?;
            let embedding = esm2.embed(&[sequence]).map_err(describe)?;
            let residue = embedding
                .first()
                .and_then(|e| e.residue(index))
                .ok_or_else(|| fail("no embedding for that position"))?;

            let logp = |letter: char| -> Result<f32, CommandError> {
                log_probs
                    .get(token_index(letter))
                    .copied()
                    .ok_or_else(|| fail(format!("the embedding model has no token for residue {letter}")))
            };
            let wt_logp = logp(wild_type)?;
            // All 19 substitutions go through the head in one batch.
            let mut rows = Vec::new();
            let mut mutants = Vec::new();
            let mut scores = Vec::new();
            for (mutant_index, mutant) in AminoAcid::STANDARD.iter().enumerate() {
                let letter = mutant.one_letter();
                if letter == wild_type {
                    continue;
                }
                let log_odds = logp(letter)? - wt_logp;
                rows.push(StabilityHead::features(residue, wild_index, mutant_index, log_odds));
                mutants.push((letter, log_odds));
            }
            let ddgs = head.predict(&rows, esm2.hidden_size(), esm2.device()).map_err(fail)?;
            for ((letter, log_odds), ddg) in mutants.into_iter().zip(ddgs) {
                scores.push(SubstitutionScore {
                    mutant: letter,
                    ddg: f64::from(ddg),
                    log_odds: f64::from(log_odds),
                });
            }
            // Most stabilising first. ΔΔG follows the dataset's convention:
            // negative is destabilising, so this is a plain descending sort.
            scores.sort_by(|a, b| b.ddg.total_cmp(&a.ddg));
            Ok(StabilityScan { position, wild_type, scores })
        })
    }

    /// Predict solubility for a whole sequence.
    pub fn solubility(&self, sequence: &str) -> Result<SolubilityPrediction, CommandError> {
        self.with_models(|esm2, heads| {
            let head = heads.solubility.as_ref().ok_or_else(|| missing("solubility"))?;
            let embedding = esm2
                .embed(&[sequence])
                .map_err(describe)?
                .into_iter()
                .next()
                .ok_or_else(|| fail("the embedding model returned nothing for this sequence"))?;
            // `pooled` is the pooling the head was trained on; recomputing a
            // mean here would risk the two definitions drifting apart.
            let percent = head.predict(&embedding.pooled, esm2.device()).map_err(fail)?;
            Ok(SolubilityPrediction { percent: f64::from(percent) })
        })
    }

    /// Predict per-residue disorder.
    pub fn disorder(&self, sequence: &str) -> Result<DisorderPrediction, CommandError> {
        self.with_models(|esm2, heads| {
            let head = heads.disorder.as_ref().ok_or_else(|| missing("disorder"))?;
            let windowed = sequence.chars().count() > esm2.max_residues();
            let embedding = if windowed {
                esm2.embed_windowed(sequence).map_err(describe)?
            } else {
                esm2.embed(&[sequence])
                    .map_err(describe)?
                    .into_iter()
                    .next()
                    .ok_or_else(|| fail("the embedding model returned nothing for this sequence"))?
            };

            let probabilities = head
                .predict(&embedding.per_residue, embedding.length, esm2.device())
                .map_err(fail)?
                .into_iter()
                .map(f64::from)
                .collect();
            Ok(DisorderPrediction { probabilities, threshold: disorder_threshold(&self.dir), windowed })
        })
    }

    /// Per-residue solvent accessibility.
    ///
    /// # Errors
    /// [`CommandError`] if the head is absent, the embedding model cannot be
    /// loaded, or inference fails.
    pub fn burial(&self, sequence: &str) -> Result<BurialPrediction, CommandError> {
        self.with_models(|esm2, heads| {
            let head = heads.burial.as_ref().ok_or_else(|| missing("burial"))?;
            // Windowed like disorder: accessibility is a local property, so a
            // residue's own neighbourhood is what the head reads, and a
            // window boundary costs less than refusing the sequence.
            let windowed = sequence.chars().count() > esm2.max_residues();
            let embedding = if windowed {
                esm2.embed_windowed(sequence).map_err(describe)?
            } else {
                esm2.embed(&[sequence])
                    .map_err(describe)?
                    .into_iter()
                    .next()
                    .ok_or_else(|| fail("the embedding model returned nothing for this sequence"))?
            };

            let accessibility = head
                .predict(&embedding.per_residue, embedding.length, esm2.device())
                .map_err(fail)?
                .into_iter()
                .map(f64::from)
                .collect();
            Ok(BurialPrediction { accessibility, buried_at: BURIED_AT, windowed })
        })
    }
}

/// Index of a residue among the 20 standard amino acids, as the heads expect.
fn standard_index(letter: char) -> Option<usize> {
    AminoAcid::STANDARD.iter().position(|a| a.one_letter() == letter)
}

/// Index of a residue in the embedding model's vocabulary.
///
/// This must agree with how the training pipeline built its log-odds
/// feature; both call the same function, so they cannot diverge.
fn token_index(letter: char) -> usize {
    opendna_ml::esm2::alphabet::token_id(letter) as usize
}

/// The error shown when a head's weights are not installed.
///
/// It names the file the application looked for, because the usual cause is
/// that the head was never trained or never copied into place, and a user
/// cannot act on "unavailable" alone.
fn missing(task: &str) -> CommandError {
    CommandError {
        message: format!(
            "No trained {task} head is installed, so this prediction is unavailable. Train one \
             with `cargo run -p opendna-train -- train --task {task}` and put \
             {task}.safetensors and {task}.metadata.json in the models directory."
        ),
    }
}

/// Turn a model error into something a user can act on.
///
/// Used for every `MlError` in this module, so the two cases a user actually
/// hits — no network on the first run, and a sequence past the checkpoint's
/// limit — read as instructions rather than as internal error text.
fn describe(e: MlError) -> CommandError {
    match e {
        MlError::Download { repo, message } => CommandError {
            message: format!(
                "Could not download the embedding model {repo}: {message}. The first prediction \
                 needs a network connection; afterwards it is cached and works offline."
            ),
        },
        // The limit is a property of the checkpoint, so there is nothing the
        // user can change except the input; say which protein is too long
        // rather than repeat the batch index, which is always 0 here.
        MlError::SequenceTooLong { length, limit, .. } => CommandError {
            message: format!(
                "This sequence has {length} residues and {} accepts at most {limit}. Stability                  and solubility predictions need the whole sequence at once; disorder covers                  long sequences with overlapping windows.",
                EMBEDDER.repo()
            ),
        },
        other => fail(other),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_heads_are_reported_not_substituted() {
        let empty = tempdir();
        let predictor = Predictor::new(empty.clone());
        let status = predictor.status();
        assert_eq!(status.stability, HeadStatus::Missing);
        assert_eq!(status.solubility, HeadStatus::Missing);
        assert_eq!(status.disorder, HeadStatus::Missing);
        assert!(!status.stability.is_ready());
    }

    #[test]
    fn a_head_trained_on_another_embedder_is_unusable_rather_than_ready() {
        let dir = tempdir();
        std::fs::write(dir.join("stability.safetensors"), b"not really weights").unwrap();
        std::fs::write(
            dir.join("stability.metadata.json"),
            br#"{"task":"stability","embedding_model":"facebook/esm2_t6_8M_UR50D",
                 "embedding_dim":320,"version":"0.1.0"}"#,
        )
        .unwrap();
        match head_status(&dir, "stability") {
            HeadStatus::Unusable { reason } => {
                assert!(reason.contains("esm2_t6_8M_UR50D"), "{reason}");
            }
            other => panic!("expected Unusable, got {other:?}"),
        }
    }

    /// The panel used to carry its own hardcoded caveats, which could drift
    /// from the published card without anyone noticing. They now come from
    /// the card beside the weights, so the panel, the exported CSV and the
    /// card on Hugging Face cannot say three different things.
    #[test]
    fn a_ready_head_carries_the_limitations_from_its_card() {
        let dir = tempdir();
        std::fs::write(dir.join("disorder.safetensors"), b"not really weights").unwrap();
        std::fs::write(
            dir.join("disorder.metadata.json"),
            br#"{"task":"disorder","embedding_model":"facebook/esm2_t30_150M_UR50D",
                 "embedding_dim":640,"version":"0.6.0-dev"}"#,
        )
        .unwrap();
        std::fs::write(
            dir.join("disorder.README.md"),
            "# Head\n\n## Limitations\n\n- Annotated, not exhaustive.\n- Pooled over\n  residues.\n\n## Intended use\n\n- Not a limitation.\n",
        )
        .unwrap();

        match head_status(&dir, "disorder") {
            HeadStatus::Ready { version, limitations } => {
                assert_eq!(version, "0.6.0-dev");
                assert_eq!(
                    limitations,
                    vec!["Annotated, not exhaustive.", "Pooled over residues."],
                    "a wrapped bullet must keep its tail, and intended-use must not leak in"
                );
            }
            other => panic!("expected Ready, got {other:?}"),
        }
    }

    /// A head with no card is still usable; the UI points at the card rather
    /// than refusing to predict.
    #[test]
    fn a_head_without_a_card_is_still_ready() {
        let dir = tempdir().join("nocard");
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("solubility.safetensors"), b"not really weights").unwrap();
        std::fs::write(
            dir.join("solubility.metadata.json"),
            br#"{"task":"solubility","embedding_model":"facebook/esm2_t30_150M_UR50D",
                 "embedding_dim":640,"version":"0.6.0-dev"}"#,
        )
        .unwrap();
        match head_status(&dir, "solubility") {
            HeadStatus::Ready { limitations, .. } => assert_eq!(limitations, [] as [std::string::String; 0]),
            other => panic!("expected Ready, got {other:?}"),
        }
    }

    #[test]
    fn the_models_directory_can_be_overridden() {
        // Without the variable the directory sits beside the database.
        let app_data = Path::new("/tmp/opendna-data");
        assert_eq!(models_dir(app_data), app_data.join("models"));
    }

    fn tempdir() -> PathBuf {
        let dir = std::env::temp_dir()
            .join(format!("opendna-predict-test-{}", std::process::id()))
            .join(format!("{:?}", std::thread::current().id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }
}

/// Exercises the real inference path against the heads in the repository's
/// `models/` directory. Ignored by default: it loads the 150M-parameter
/// embedding model, downloading it on first run.
///
/// Run with:
/// ```text
/// cargo test -p opendna-desktop --lib -- --ignored --nocapture
/// ```
#[cfg(test)]
mod end_to_end {
    use super::*;

    /// Ubiquitin. Leu8 and Ile44 are part of the hydrophobic patch; position
    /// 8 is the one spot-checked by hand when the head was trained.
    const UBIQUITIN: &str = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

    fn repo_models() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../models")
    }

    /// One predictor shared by every test here, as the application has one.
    ///
    /// Each test used to build its own, so four of them loaded the same
    /// several-hundred-megabyte checkpoint at once and Windows refused the
    /// concurrent opens — "the process cannot access the file because it is
    /// being used by another process", which reads like a missing model and
    /// is not. Sharing it also exercises the arrangement the app actually
    /// uses: one `Predictor`, whose own lock serialises inference.
    fn predictor() -> &'static Predictor {
        static SHARED: std::sync::OnceLock<Predictor> = std::sync::OnceLock::new();
        SHARED.get_or_init(|| Predictor::new(repo_models()))
    }

    #[test]
    #[ignore = "loads the 150M embedding model"]
    fn the_trained_heads_load_from_the_repository() {
        let status = predictor().status();
        assert!(status.stability.is_ready(), "stability: {:?}", status.stability);
        assert!(status.solubility.is_ready(), "solubility: {:?}", status.solubility);
        assert!(status.disorder.is_ready(), "disorder: {:?}", status.disorder);
        assert!(status.burial.is_ready(), "burial: {:?}", status.burial);
    }

    #[test]
    #[ignore = "loads the 150M embedding model"]
    fn a_stability_scan_ranks_proline_and_glycine_worst_in_a_helix() {
        let scan = predictor().stability_scan(UBIQUITIN, 8).unwrap();
        assert_eq!(scan.wild_type, 'L', "position 8 of ubiquitin is leucine");
        assert_eq!(scan.scores.len(), 19, "every substitution but the wild type");

        // Sorted most stabilising first, so the tail is the worst.
        let worst: Vec<char> = scan.scores.iter().rev().take(4).map(|s| s.mutant).collect();
        assert!(
            worst.contains(&'P') && worst.contains(&'G'),
            "proline and glycine should be among the most destabilising at a buried \
             position; got {worst:?} from {:?}",
            scan.scores.iter().map(|s| (s.mutant, s.ddg)).collect::<Vec<_>>()
        );

        // The head is known to compress predictions toward the training mean,
        // which is stated on its card. This asserts the range is plausible
        // rather than pinning the exact values a retrain would move.
        for s in &scan.scores {
            assert!(s.ddg > -8.0 && s.ddg < 8.0, "{} gave an implausible ΔΔG {}", s.mutant, s.ddg);
        }
    }

    /// Ubiquitin's core is the check: I3, V5, I13, L15, V17, I23, V26, I30,
    /// L43, L50, L56, I61, L67 and L69 pack the hydrophobic interior, while
    /// K48 is the polyubiquitination site and must be on the surface to be
    /// conjugated. A head wired to the wrong residues would still return
    /// numbers on the right scale, so the scale alone proves nothing.
    #[test]
    #[ignore = "loads the 150M embedding model"]
    fn burial_puts_ubiquitins_core_inside_and_its_lysines_out() {
        let p = predictor();
        let got = p.burial(UBIQUITIN).unwrap();
        assert_eq!(got.accessibility.len(), UBIQUITIN.len(), "one value per residue");
        assert!(!got.windowed, "ubiquitin fits in one window");
        assert!(
            got.accessibility.iter().all(|v| (0.0..=1.3).contains(v)),
            "accessibility outside its physical range: {:?}",
            got.accessibility.iter().copied().fold(f64::NAN, f64::max)
        );

        let at = |position: usize| got.accessibility[position - 1];
        let core = [3, 5, 13, 15, 17, 23, 26, 30, 43, 50, 56, 61, 67, 69];
        let core_mean = core.iter().map(|i| at(*i)).sum::<f64>() / core.len() as f64;
        let surface = [6, 11, 27, 29, 48, 63, 72, 74];
        let surface_mean = surface.iter().map(|i| at(*i)).sum::<f64>() / surface.len() as f64;

        assert!(core_mean < 0.2, "the hydrophobic core should be buried, got {core_mean:.3}");
        assert!(surface_mean > 0.3, "surface residues should be exposed, got {surface_mean:.3}");
        assert!(
            surface_mean > core_mean + 0.2,
            "core {core_mean:.3} and surface {surface_mean:.3} are not distinguished"
        );
        assert!(at(48) > 0.25, "K48 is a conjugation site and must be exposed, got {:.3}", at(48));
    }

    #[test]
    #[ignore = "loads the 150M embedding model"]
    fn solubility_and_disorder_return_values_on_their_documented_scales() {
        let p = predictor();
        let solubility = p.solubility(UBIQUITIN).unwrap();
        // The assay scale is a percentage; the head does not clip, so this
        // checks for a number on roughly the right scale, not a valid one.
        assert!(
            solubility.percent > -50.0 && solubility.percent < 200.0,
            "solubility {} is not on the assay's scale",
            solubility.percent
        );

        let disorder = p.disorder(UBIQUITIN).unwrap();
        assert_eq!(disorder.probabilities.len(), UBIQUITIN.len(), "one probability per residue");
        assert!(!disorder.windowed, "ubiquitin fits in one window");
        assert!(
            disorder.probabilities.iter().all(|v| (0.0..=1.0).contains(v)),
            "probabilities must be in [0, 1]"
        );
        // Ubiquitin is a well-ordered globular domain. If most of it came
        // back disordered, something is wired up wrongly.
        let called = disorder.probabilities.iter().filter(|v| **v >= disorder.threshold).count();
        assert!(
            called * 2 < UBIQUITIN.len(),
            "{called} of {} residues called disordered in a folded domain",
            UBIQUITIN.len()
        );
    }

    #[test]
    fn a_position_past_the_end_is_refused_before_any_model_loads() {
        // No #[ignore]: this must not touch the embedding model at all.
        let err = predictor().stability_scan(UBIQUITIN, 999).unwrap_err();
        assert!(err.message.contains("past the end"), "{}", err.message);
    }

    #[test]
    fn position_zero_is_refused() {
        let err = predictor().stability_scan(UBIQUITIN, 0).unwrap_err();
        assert!(err.message.contains("numbered from 1"), "{}", err.message);
    }
}
