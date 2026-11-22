//! Driving a folding engine, and checking that what came back was asked for.
//!
//! OpenDNA does not fold anything. ESMFold and ColabFold do, and this runs
//! one of them and then does the part that is easy to skip: confirms that the
//! chain in the file the engine wrote has the sequence that was sent to it.
//!
//! That check is not ceremony. A folding run writes its outputs into a
//! directory under names of its own choosing, several per sequence in
//! ColabFold's case, and the ordinary way a batch goes wrong is that a
//! structure ends up filed against the wrong sequence — a stale file from a
//! previous run, a name mangled on the way through, a failure that produced
//! fewer outputs than inputs without saying so. Matching on the residues
//! rather than on the filename makes that class of mistake impossible to
//! reach, and costs a parse of a file that was going to be parsed anyway.
//!
//! What the engine wrote is kept byte for byte. Nothing here re-serialises a
//! structure: the file copied into place is the file the engine produced, so
//! its digest can be compared against the engine's own output directory,
//! which is also kept.

use crate::{Engine, Requirement, KNOWN};
use opendna_core::Structure;
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::Command;

/// A folding engine OpenDNA knows how to drive.
///
/// Separate from [`Requirement`], which only says how to find a tool. Being
/// findable and being drivable are different facts, and only these two are
/// both.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Folder {
    /// ESMFold, through the `esm-fold` entry point shipped with `fair-esm`.
    EsmFold,
    /// ColabFold, through `colabfold_batch`.
    ColabFold,
}

impl Folder {
    /// Every engine that can be driven, in the order they are tried.
    pub const ALL: &'static [Self] = &[Self::EsmFold, Self::ColabFold];

    /// The name used on the command line, in a lockfile and in a manifest.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::EsmFold => "esmfold",
            Self::ColabFold => "colabfold",
        }
    }

    /// How to find it.
    ///
    /// # Panics
    /// Never in a built binary: both names are in [`KNOWN`], and a test here
    /// holds that true.
    #[must_use]
    pub fn requirement(self) -> &'static Requirement {
        KNOWN
            .iter()
            .find(|r| r.name == self.name())
            .unwrap_or_else(|| unreachable!("{} is not in KNOWN", self.name()))
    }

    /// The arguments that ask it to fold `input` into `out`.
    ///
    /// Both engines take the same two pieces of information in different
    /// shapes and neither has a stable flag for anything else, so anything
    /// further is passed through verbatim by the caller rather than being
    /// translated here. A flag this code invented a meaning for would be a
    /// flag that quietly means something else after the next release.
    fn arguments(self, input: &Path, out: &Path) -> Vec<OsString> {
        match self {
            Self::EsmFold => vec!["-i".into(), input.into(), "-o".into(), out.into()],
            Self::ColabFold => vec![input.into(), out.into()],
        }
    }
}

/// Why a fold did not produce a usable structure.
#[derive(Debug, thiserror::Error)]
pub enum FoldError {
    /// The engine ran and failed.
    #[error("{engine} exited {status}. Its last output was:\n{tail}")]
    Failed {
        /// Which engine.
        engine: &'static str,
        /// The exit status, as a short string.
        status: String,
        /// The tail of what it wrote, which is where the reason usually is.
        tail: String,
    },
    /// The engine succeeded but returned nothing that is this sequence.
    #[error("{engine} did not return a structure matching {name}. {looked}")]
    Missing {
        /// Which engine.
        engine: &'static str,
        /// The sequence that went unanswered.
        name: String,
        /// What was searched and what was found instead.
        looked: String,
    },
    /// Several outputs match and nothing distinguishes them.
    #[error("{engine} produced {} structures for {name} and nothing distinguishes them:\n{}", candidates.len(), candidates.join("\n"))]
    Ambiguous {
        /// Which engine.
        engine: &'static str,
        /// The sequence in question.
        name: String,
        /// The files that all match it.
        candidates: Vec<String>,
    },
    /// An output was unreadable or was not a single protein chain.
    #[error("{} is not a structure this can use: {reason}", path.display())]
    Unusable {
        /// The file.
        path: PathBuf,
        /// What was wrong with it.
        reason: String,
    },
    /// The engine could not be started, or a file could not be read.
    #[error("{what}: {reason}")]
    Io {
        /// What was being attempted.
        what: String,
        /// The underlying message.
        reason: String,
    },
}

/// One sequence, folded and checked.
#[derive(Debug, Clone)]
pub struct Prediction {
    /// The name from the FASTA header.
    pub name: String,
    /// The sequence that was sent, as written.
    pub sequence: String,
    /// The file the engine wrote, where it wrote it.
    pub source: PathBuf,
    /// Its bytes, unchanged. What gets written out is this.
    pub bytes: Vec<u8>,
    /// How many amino-acid residues the chain has.
    pub residues: usize,
    /// Mean B-factor over the chain's Cα atoms.
    ///
    /// Both engines put their per-residue confidence (pLDDT, 0–100) in this
    /// column, which is why it is worth reporting. It is not renamed to
    /// pLDDT here: the column is the B-factor column, this is what was in it,
    /// and nothing in this crate checked what the engine meant by it.
    pub mean_ca_b_factor: f64,
}

/// Build the command that runs the engine, without running it.
///
/// Separate from [`run`] so the caller can record the exact argv in a
/// manifest and so a test can substitute a stand-in for the engine.
///
/// `input` and `out` should be absolute: ColabFold writes some of its
/// intermediates relative to the working directory, and an engine that
/// resolves a relative path against a directory the caller did not choose is
/// an engine whose outputs land somewhere nobody looks.
#[must_use]
pub fn command(folder: Folder, engine: &Engine, input: &Path, out: &Path, extra: &[String]) -> Command {
    let mut command = Command::new(&engine.path);
    command.args(folder.arguments(input, out));
    command.args(extra);
    command
}

/// The command as a list of strings, for a manifest and for an error message.
///
/// Lossy on a path that is not valid Unicode, which is the right trade for a
/// field whose only job is to be read.
#[must_use]
pub fn command_line(command: &Command) -> Vec<String> {
    std::iter::once(command.get_program())
        .chain(command.get_args())
        .map(|piece| piece.to_string_lossy().into_owned())
        .collect()
}

/// Run the engine and match its outputs to the sequences that were requested.
///
/// `requested` is `(name, sequence)` in input order; the result is in the same
/// order, one entry per request, or an error naming the first request that
/// could not be satisfied.
///
/// # Errors
/// [`FoldError::Failed`] if the engine exits non-zero, and the other variants
/// if it exits cleanly without having answered the question.
pub fn run(
    folder: Folder,
    mut command: Command,
    requested: &[(String, String)],
    out: &Path,
) -> Result<Vec<Prediction>, FoldError> {
    let output = command.output().map_err(|e| FoldError::Io {
        what: format!("running {}", command_line(&command).join(" ")),
        reason: e.to_string(),
    })?;

    if !output.status.success() {
        return Err(FoldError::Failed {
            engine: folder.name(),
            status: crate::status_of(output.status),
            tail: tail_of(&output.stderr, &output.stdout),
        });
    }

    let mut candidates = collect(out)?;
    let mut predictions = Vec::with_capacity(requested.len());
    for (name, sequence) in requested {
        let chosen = select(folder, name, sequence, &candidates)?;
        let candidate = candidates.remove(chosen);
        predictions.push(Prediction {
            name: name.clone(),
            sequence: sequence.clone(),
            residues: candidate.sequence.len(),
            mean_ca_b_factor: candidate.mean_ca_b_factor,
            source: candidate.path,
            bytes: candidate.bytes,
        });
    }
    Ok(predictions)
}

/// A structure found in the output directory, parsed far enough to match.
#[derive(Debug)]
struct Candidate {
    path: PathBuf,
    /// The one-letter sequence of its single protein chain.
    sequence: String,
    mean_ca_b_factor: f64,
    bytes: Vec<u8>,
}

/// Every PDB file under `out`, parsed.
///
/// Files that are not single-chain protein structures are skipped rather than
/// refused: a folding run drops other things in its output directory, and a
/// scores file or a template is not a failure. A requested sequence with no
/// match is what gets reported, and it reports what was there.
fn collect(out: &Path) -> Result<Vec<Candidate>, FoldError> {
    let mut found = Vec::new();
    let mut files = Vec::new();
    walk(out, 0, &mut files)?;
    files.sort();
    for path in files {
        let Ok(bytes) = std::fs::read(&path) else { continue };
        let Ok(text) = std::str::from_utf8(&bytes) else { continue };
        let Ok(structure) = opendna_io::pdb::read(text) else { continue };
        let Some((sequence, mean)) = single_chain(&structure) else { continue };
        found.push(Candidate { path, sequence, mean_ca_b_factor: mean, bytes });
    }
    Ok(found)
}

/// Files ending in `.pdb`, to a bounded depth.
///
/// Bounded because the directory being walked was written by a program this
/// one does not control, and an unbounded walk of somebody else's output is
/// how a symlink loop becomes a hang.
fn walk(dir: &Path, depth: usize, into: &mut Vec<PathBuf>) -> Result<(), FoldError> {
    if depth > 3 {
        return Ok(());
    }
    let entries = std::fs::read_dir(dir)
        .map_err(|e| FoldError::Io { what: format!("reading {}", dir.display()), reason: e.to_string() })?;
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            walk(&path, depth + 1, into)?;
        } else if path.extension().is_some_and(|e| e.eq_ignore_ascii_case("pdb")) {
            into.push(path);
        }
    }
    Ok(())
}

/// The one-letter sequence and mean Cα B-factor of a lone protein chain.
///
/// `None` for anything with no protein chain or more than one: a monomer was
/// requested, and a file holding something else is not an answer to it.
fn single_chain(structure: &Structure) -> Option<(String, f64)> {
    let mut chains = structure.chains.iter().filter(|c| c.amino_acids().next().is_some());
    let chain = chains.next()?;
    if chains.next().is_some() {
        return None;
    }
    let b_factors: Vec<f64> =
        chain.amino_acids().filter_map(|r| r.atom("CA")).map(|a| f64::from(a.b_factor)).collect();
    let mean =
        if b_factors.is_empty() { 0.0 } else { b_factors.iter().sum::<f64>() / b_factors.len() as f64 };
    Some((chain.sequence().to_string(), mean))
}

/// Which candidate answers this request.
///
/// Matching is on the residues first, which is the part that cannot be faked
/// by a filename. Only when that leaves more than one — ColabFold ranks, or
/// two inputs that are the same sequence under different names — do the
/// filenames break the tie, and if they cannot, this refuses rather than
/// taking the first.
fn select(folder: Folder, name: &str, sequence: &str, candidates: &[Candidate]) -> Result<usize, FoldError> {
    let matching: Vec<usize> = candidates
        .iter()
        .enumerate()
        .filter(|(_, c)| same_sequence(sequence, &c.sequence))
        .map(|(i, _)| i)
        .collect();

    if matching.is_empty() {
        let looked = if candidates.is_empty() {
            "No readable single-chain structure was left in its output directory.".to_string()
        } else {
            format!(
                "{} structure(s) were there, none with this sequence: {}",
                candidates.len(),
                candidates
                    .iter()
                    .map(|c| format!("{} ({} residues)", c.path.display(), c.sequence.len()))
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        };
        return Err(FoldError::Missing { engine: folder.name(), name: name.to_string(), looked });
    }
    if let [only] = matching[..] {
        return Ok(only);
    }

    // The sequence did not separate them. Try the name the caller gave, then
    // ColabFold's own ranking, in that order.
    let narrowed = narrow(&matching, candidates, |stem| stem.contains(&sanitise(name)));
    if let [only] = narrowed[..] {
        return Ok(only);
    }
    let ranked = narrow(&narrowed, candidates, |stem| stem.contains("rank_001"));
    if let [only] = ranked[..] {
        return Ok(only);
    }

    Err(FoldError::Ambiguous {
        engine: folder.name(),
        name: name.to_string(),
        candidates: ranked.iter().map(|&i| candidates[i].path.display().to_string()).collect(),
    })
}

/// Keep the entries whose filename satisfies `keep`, or all of them if none
/// does. Narrowing to nothing would turn a tie into a false absence.
fn narrow(indices: &[usize], candidates: &[Candidate], keep: impl Fn(&str) -> bool) -> Vec<usize> {
    let kept: Vec<usize> = indices
        .iter()
        .copied()
        .filter(|&i| {
            let stem = candidates[i].path.file_name().unwrap_or_default().to_string_lossy();
            keep(&stem)
        })
        .collect();
    if kept.is_empty() {
        indices.to_vec()
    } else {
        kept
    }
}

/// A name reduced to characters that are a filename everywhere.
#[must_use]
pub fn sanitise(name: &str) -> String {
    name.chars().map(|c| if c.is_ascii_alphanumeric() || "._-".contains(c) { c } else { '_' }).collect()
}

/// Whether a returned chain is the sequence that was asked for.
///
/// Equal length and equal residues, except where the request was ambiguous.
/// `X`, `B`, `Z` and `J` do not name a residue, so an engine has to choose
/// one to build, and holding it to a letter that does not specify an answer
/// would reject a correct structure.
fn same_sequence(requested: &str, got: &str) -> bool {
    if requested.len() != got.len() {
        return false;
    }
    requested.bytes().zip(got.bytes()).all(|(want, have)| {
        want.eq_ignore_ascii_case(&have) || matches!(want.to_ascii_uppercase(), b'X' | b'B' | b'Z' | b'J')
    })
}

/// The last few lines of whatever the engine wrote, for an error message.
fn tail_of(stderr: &[u8], stdout: &[u8]) -> String {
    let err = String::from_utf8_lossy(stderr);
    let out = String::from_utf8_lossy(stdout);
    let text = if err.trim().is_empty() { out } else { err };
    let lines: Vec<&str> = text.lines().filter(|l| !l.trim().is_empty()).collect();
    let start = lines.len().saturating_sub(8);
    let tail = lines[start..].join("\n");
    if tail.is_empty() {
        "  (nothing)".to_string()
    } else {
        tail
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// A Cα-only PDB line in fixed columns, which is all the matcher reads.
    fn atom(serial: usize, residue: &str, number: usize, b: f64) -> String {
        format!(
            "ATOM  {serial:>5}  CA  {residue} A{number:>4}    {:8.3}{:8.3}{:8.3}  1.00{b:6.2}           C",
            f64::from(u16::try_from(serial).unwrap_or(0)) * 3.8,
            0.0,
            0.0
        )
    }

    fn pdb(residues: &[(&str, f64)]) -> String {
        let mut text = String::new();
        for (i, (residue, b)) in residues.iter().enumerate() {
            text.push_str(&atom(i + 1, residue, i + 1, *b));
            text.push('\n');
        }
        text.push_str("END\n");
        text
    }

    /// MQ, with the two pLDDT values that should average to 83.
    fn mq(b: [f64; 2]) -> String {
        pdb(&[("MET", b[0]), ("GLN", b[1])])
    }

    fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, text).unwrap();
        path
    }

    fn scratch(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "opendna-fold-{tag}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn both_drivable_engines_are_findable_ones() {
        for folder in Folder::ALL {
            assert_eq!(folder.requirement().name, folder.name());
            assert_eq!(folder.requirement().purpose, "structure prediction from sequence");
        }
    }

    /// The two engines take the same two facts in different shapes, and
    /// getting either shape wrong is a run that fails or writes somewhere
    /// nobody looks.
    #[test]
    fn each_engine_is_asked_in_its_own_shape() {
        let input = Path::new("/w/in.fasta");
        let out = Path::new("/w/out");
        let esm: Vec<String> =
            Folder::EsmFold.arguments(input, out).iter().map(|a| a.to_string_lossy().into_owned()).collect();
        assert_eq!(esm, ["-i", "/w/in.fasta", "-o", "/w/out"]);
        let colab: Vec<String> = Folder::ColabFold
            .arguments(input, out)
            .iter()
            .map(|a| a.to_string_lossy().into_owned())
            .collect();
        assert_eq!(colab, ["/w/in.fasta", "/w/out"]);
    }

    #[test]
    fn a_structure_is_matched_by_its_residues_and_its_confidence_read() {
        let dir = scratch("match");
        write(&dir, "anything-at-all.pdb", &mq([88.0, 78.0]));
        let found = collect(&dir).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].sequence, "MQ");
        assert!((found[0].mean_ca_b_factor - 83.0).abs() < 1e-9, "{}", found[0].mean_ca_b_factor);
        assert_eq!(select(Folder::EsmFold, "p1", "MQ", &found).unwrap(), 0);
    }

    /// The filename is deliberately not what identifies a structure. A file
    /// named after the sequence it is not is the mistake this guards.
    #[test]
    fn a_structure_filed_under_the_wrong_name_is_refused_not_used() {
        let dir = scratch("wrongname");
        write(&dir, "p1.pdb", &pdb(&[("ALA", 50.0), ("ALA", 50.0)]));
        let found = collect(&dir).unwrap();
        let error = select(Folder::EsmFold, "p1", "MQ", &found).unwrap_err().to_string();
        assert!(error.contains("none with this sequence"), "{error}");
        assert!(error.contains("p1.pdb"), "{error}");
    }

    #[test]
    fn a_sequence_with_no_output_at_all_says_so() {
        let dir = scratch("empty");
        let found = collect(&dir).unwrap();
        let error = select(Folder::ColabFold, "p1", "MQ", &found).unwrap_err().to_string();
        assert!(error.contains("did not return a structure matching p1"), "{error}");
        assert!(error.contains("No readable single-chain structure"), "{error}");
    }

    /// ColabFold writes one file per ranked model. The top-ranked one is the
    /// answer, and picking it by name is safe only after the residues have
    /// already narrowed the field to that sequence's models.
    #[test]
    fn colabfold_ranks_are_resolved_to_rank_001() {
        let dir = scratch("ranks");
        write(&dir, "p1_unrelaxed_rank_002_model_3.pdb", &mq([60.0, 60.0]));
        write(&dir, "p1_unrelaxed_rank_001_model_1.pdb", &mq([90.0, 90.0]));
        let found = collect(&dir).unwrap();
        let chosen = select(Folder::ColabFold, "p1", "MQ", &found).unwrap();
        assert!(found[chosen].path.to_string_lossy().contains("rank_001"), "{:?}", found[chosen].path);
        assert!((found[chosen].mean_ca_b_factor - 90.0).abs() < 1e-9);
    }

    /// Two different proteins that happen to have the same sequence are
    /// separated by name; nothing is guessed.
    #[test]
    fn identical_sequences_under_different_names_are_separated_by_name() {
        let dir = scratch("twins");
        write(&dir, "alpha.pdb", &mq([70.0, 70.0]));
        write(&dir, "beta.pdb", &mq([40.0, 40.0]));
        let found = collect(&dir).unwrap();
        let alpha = select(Folder::EsmFold, "alpha", "MQ", &found).unwrap();
        assert!((found[alpha].mean_ca_b_factor - 70.0).abs() < 1e-9);
    }

    /// When nothing separates two matches, refusing is the only honest
    /// answer: taking the first would file a structure against a sequence on
    /// the strength of directory order.
    #[test]
    fn an_unresolvable_tie_refuses_rather_than_taking_the_first() {
        let dir = scratch("tie");
        write(&dir, "one.pdb", &mq([70.0, 70.0]));
        write(&dir, "two.pdb", &mq([40.0, 40.0]));
        let found = collect(&dir).unwrap();
        let error = select(Folder::EsmFold, "p1", "MQ", &found).unwrap_err().to_string();
        assert!(error.contains("nothing distinguishes them"), "{error}");
        assert!(error.contains("one.pdb") && error.contains("two.pdb"), "{error}");
    }

    /// An ambiguity code does not name a residue, so the engine had to pick
    /// one and holding it to the letter would reject a correct structure.
    #[test]
    fn an_ambiguity_code_in_the_request_matches_whatever_was_built() {
        assert!(same_sequence("MXQ", "MAQ"));
        assert!(same_sequence("MBQ", "MNQ"));
        // Everything else still has to be the residue that was asked for.
        assert!(!same_sequence("MAQ", "MGQ"));
        assert!(!same_sequence("MQ", "MQA"));
        assert!(same_sequence("mq", "MQ"));
    }

    /// A scores file, a log, a template: an output directory holds more than
    /// structures, and none of that is a failure.
    #[test]
    fn non_structures_in_the_output_directory_are_ignored() {
        let dir = scratch("clutter");
        write(&dir, "scores.json", "{\"plddt\": [90, 90]}");
        write(&dir, "log.txt", "done\n");
        write(&dir, "broken.pdb", "not a pdb file at all\n");
        write(&dir, "p1.pdb", &mq([90.0, 90.0]));
        let found = collect(&dir).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].sequence, "MQ");
    }

    /// A monomer was requested. A complex is not an answer to it, and
    /// silently taking one of its chains would be.
    #[test]
    fn a_multi_chain_file_is_not_an_answer_to_a_monomer() {
        let dir = scratch("complex");
        let mut text = mq([90.0, 90.0]);
        text = text.replace("END\n", "");
        text.push_str("ATOM      3  CA  MET B   1       0.000   0.000   0.000  1.00 90.00           C\n");
        write(&dir, "complex.pdb", &text);
        assert!(collect(&dir).unwrap().is_empty());
    }

    #[test]
    fn outputs_in_subdirectories_are_found() {
        let dir = scratch("nested");
        let inner = dir.join("predictions");
        std::fs::create_dir_all(&inner).unwrap();
        write(&inner, "p1.pdb", &mq([90.0, 90.0]));
        assert_eq!(collect(&dir).unwrap().len(), 1);
    }

    /// The whole run, with a stand-in for the engine: a command that exits
    /// cleanly over an output directory that already holds a structure, which
    /// is the state a real engine leaves behind.
    #[test]
    fn a_successful_run_returns_the_bytes_the_engine_wrote() {
        let out = scratch("run");
        let structure = mq([88.0, 78.0]);
        std::fs::write(out.join("p1.pdb"), &structure).unwrap();

        let predictions =
            run(Folder::EsmFold, succeeding_command(), &[("p1".into(), "MQ".into())], &out).unwrap();
        assert_eq!(predictions.len(), 1);
        assert_eq!(predictions[0].name, "p1");
        assert_eq!(predictions[0].residues, 2);
        assert_eq!(predictions[0].bytes, structure.as_bytes(), "the engine's bytes, unchanged");
        assert!((predictions[0].mean_ca_b_factor - 83.0).abs() < 1e-9);
    }

    /// An engine that fails must not look like one that produced nothing:
    /// the two have different causes and different fixes.
    #[test]
    fn an_engine_that_fails_reports_what_it_said() {
        let out = scratch("failure");
        let predictions = run(Folder::EsmFold, failing_command(), &[("p1".into(), "MQ".into())], &out);
        let error = predictions.unwrap_err().to_string();
        assert!(error.contains("esmfold exited"), "{error}");
        assert!(error.contains("CUDA out of memory"), "{error}");
    }

    #[test]
    fn the_command_line_is_recoverable_for_a_manifest() {
        let engine = Engine {
            name: "esmfold",
            path: PathBuf::from("/usr/bin/esm-fold"),
            version: "1.0".into(),
            sha256: "ab".repeat(32),
            identification: crate::Identification::EntryPoint,
        };
        let command = command(
            Folder::EsmFold,
            &engine,
            Path::new("/w/in.fasta"),
            Path::new("/w/out"),
            &["--num-recycles".to_string(), "4".to_string()],
        );
        let line = command_line(&command);
        assert_eq!(line[0], "/usr/bin/esm-fold");
        assert_eq!(&line[1..], ["-i", "/w/in.fasta", "-o", "/w/out", "--num-recycles", "4"]);
    }

    #[test]
    fn a_name_becomes_a_filename_without_surprises() {
        assert_eq!(sanitise("sp|P0DTC2|SPIKE"), "sp_P0DTC2_SPIKE");
        assert_eq!(sanitise("../../etc/passwd"), ".._.._etc_passwd");
        assert_eq!(sanitise("p1"), "p1");
    }

    #[cfg(windows)]
    fn succeeding_command() -> Command {
        let mut command = Command::new("cmd");
        command.args(["/c", "exit 0"]);
        command
    }

    #[cfg(not(windows))]
    fn succeeding_command() -> Command {
        let mut command = Command::new("sh");
        command.args(["-c", "exit 0"]);
        command
    }

    #[cfg(windows)]
    fn failing_command() -> Command {
        let mut command = Command::new("cmd");
        command.args(["/c", "echo CUDA out of memory 1>&2 & exit 1"]);
        command
    }

    #[cfg(not(windows))]
    fn failing_command() -> Command {
        let mut command = Command::new("sh");
        command.args(["-c", "echo 'CUDA out of memory' >&2; exit 1"]);
        command
    }
}
