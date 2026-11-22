//! `opendna fold`: structure prediction through an installed engine.
//!
//! This is the first command that drives external software rather than a head
//! trained here, and it holds that software to the contract the heads are
//! held to. The engine is discovered, identified and recorded; what it
//! returns is checked against what it was asked for; and every structure
//! written gets a manifest saying which engine produced it, from which bytes.
//!
//! What it does not do is fold anything. There is no fallback, no bundled
//! mini-model, nothing that produces coordinates when no engine is installed.
//! A missing engine is a missing dependency and exits saying so, because a
//! plausible-looking structure from a program that cannot fold proteins is
//! the exact failure this project was rewritten to remove.

use crate::exit::{Classify as _, Failure, Fault};
use crate::manifest;
use anyhow::{anyhow, bail, Context, Result};
use opendna_engine::{Engine, FoldError, Folder};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// What `fold` was asked to do.
pub struct Request {
    /// FASTA file of sequences to fold.
    pub input: PathBuf,
    /// Directory to write structures into.
    pub out_dir: PathBuf,
    /// Which engine, or `None` to use whichever one is installed.
    pub engine: Option<Folder>,
    /// Arguments passed to the engine verbatim, after the ones OpenDNA sets.
    pub engine_args: Vec<String>,
    /// Fail unless the installed engines match this lockfile.
    pub check: Option<PathBuf>,
    /// Write a manifest beside each structure.
    pub manifest: bool,
}

/// Fold every sequence in the input.
///
/// # Errors
/// Classified by whose problem it is: a missing or mismatched engine is the
/// installation's, unreadable or unusable FASTA is the input's, an engine
/// that returned the wrong sequence is a check that did not hold, and an
/// engine that simply failed is none of those — see [`run`]'s body for why
/// that one is deliberately left unclassified.
pub fn run(args: &Request) -> Result<(), Failure> {
    let text = std::fs::read_to_string(&args.input)
        .with_context(|| format!("reading {}", args.input.display()))
        .blame(Fault::Input)?;
    let sequences = opendna_io::fasta::read_raw(&text)
        .with_context(|| format!("parsing {} as FASTA", args.input.display()))
        .blame(Fault::Input)?;
    let requested = names(&sequences).blame(Fault::Input)?;

    let (folder, engine) = resolve(args.engine)?;
    if let Some(path) = &args.check {
        verify(path).blame(Fault::Install)?;
    }

    eprintln!(
        "{} sequences from {} through {} ({})",
        requested.len(),
        args.input.display(),
        folder.name(),
        engine.version
    );

    let work = prepare_directories(&args.out_dir).blame(Fault::Unexpected)?;
    let fasta = work.raw.join("input.fasta");
    std::fs::write(&fasta, &text)
        .with_context(|| format!("writing {}", fasta.display()))
        .blame(Fault::Unexpected)?;

    let command = opendna_engine::command(folder, &engine, &fasta, &work.engine_out, &args.engine_args);
    let line = opendna_engine::command_line(&command);
    eprintln!("Running: {}", line.join(" "));

    let started = chrono::Utc::now();
    let clock = std::time::Instant::now();
    let folded = opendna_engine::run(folder, command, &requested, &work.engine_out)
        .map_err(|e| Failure { fault: fault_for(&e), error: anyhow!(e) })?;
    let elapsed = clock.elapsed();

    let residues: usize = sequences.iter().map(|(_, s)| s.len()).sum();
    let label = args.input.display().to_string();
    let mut written = Vec::new();
    for prediction in &folded {
        let output = args.out_dir.join(format!("{}.pdb", opendna_engine::sanitise(&prediction.name)));
        std::fs::write(&output, &prediction.bytes)
            .with_context(|| format!("writing {}", output.display()))
            .blame(Fault::Unexpected)?;

        if args.manifest {
            manifest::write_fold(
                &manifest::FoldRun {
                    folder: folder.name(),
                    engine: &engine,
                    command: &line,
                    started,
                    elapsed,
                    input: manifest::Input {
                        label: &label,
                        bytes: text.as_bytes(),
                        sequences: sequences.len(),
                        residues,
                    },
                    name: &prediction.name,
                    sequence: &prediction.sequence,
                    output: &output,
                    bytes: &prediction.bytes,
                    source: &prediction.source,
                    mean_ca_b_factor: prediction.mean_ca_b_factor,
                    residues: prediction.residues,
                },
                &output,
            )
            .blame(Fault::Unexpected)?;
        }
        written.push((prediction, output));
    }

    report(&written, &work, args.manifest);
    Ok(())
}

/// Where the run's files go.
struct Workspace {
    /// Everything the engine wrote, kept.
    raw: PathBuf,
    /// The directory handed to the engine.
    engine_out: PathBuf,
}

/// Create the output directory and the one the engine writes into.
///
/// The engine's own output is kept rather than deleted. It holds things
/// OpenDNA does not read -- ColabFold's per-model scores, the alignments it
/// used, its log -- and deleting them would leave a structure whose
/// provenance stops at a digest.
fn prepare_directories(out_dir: &Path) -> Result<Workspace> {
    let raw = out_dir.join("engine");
    let engine_out = raw.join("out");
    std::fs::create_dir_all(&engine_out).with_context(|| format!("creating {}", engine_out.display()))?;
    // Absolute, because the engine resolves what it is given against a
    // working directory this process does not set.
    //
    // `absolute` rather than `canonicalize`: on Windows the latter returns an
    // extended-length path beginning `\\?\`, which is correct and which a
    // great deal of software -- including the Python that both these engines
    // are written in -- mishandles. Handing a third-party tool a path shape
    // it may not understand is a failure invented here rather than found.
    let absolute = |p: &Path| std::path::absolute(p).unwrap_or_else(|_| p.to_path_buf());
    Ok(Workspace { raw: absolute(&raw), engine_out: absolute(&engine_out) })
}

/// Pair each sequence with the name its structure will be filed under.
///
/// Refuses duplicates before anything runs. Two sequences that sanitise to
/// one filename would overwrite each other, and discovering that after an
/// hour of folding is discovering it too late.
fn names(sequences: &[(String, String)]) -> Result<Vec<(String, String)>> {
    if sequences.is_empty() {
        bail!("no sequences in the input");
    }
    let mut seen: BTreeMap<String, &str> = BTreeMap::new();
    for (name, _) in sequences {
        let file = opendna_engine::sanitise(name);
        if file.is_empty() {
            bail!("the sequence named {name:?} has no characters usable in a filename");
        }
        if let Some(first) = seen.insert(file.clone(), name) {
            bail!(
                "{first:?} and {name:?} would both be written to {file}.pdb; \
                 give them names that differ in more than punctuation"
            );
        }
    }
    Ok(sequences.to_vec())
}

/// Which engine to drive, and where it is.
///
/// With no choice given, an installed engine is used and two installed
/// engines is a question rather than a guess: they do not produce the same
/// structures, and picking one by list order would make a run's result depend
/// on something nobody wrote down.
fn resolve(wanted: Option<Folder>) -> Result<(Folder, Engine), Failure> {
    if let Some(folder) = wanted {
        let engine = opendna_engine::discover(folder.requirement())
            .with_context(|| format!("{} was asked for with --engine", folder.name()))
            .blame(Fault::Install)?;
        return Ok((folder, engine));
    }

    let mut found = Vec::new();
    let mut missing = Vec::new();
    for folder in Folder::ALL {
        match opendna_engine::discover(folder.requirement()) {
            Ok(engine) => found.push((*folder, engine)),
            Err(error) => missing.push(format!("  {}: {error}", folder.name())),
        }
    }

    match found.len() {
        1 => Ok(found.remove(0)),
        0 => Err(Failure {
            fault: Fault::Install,
            error: anyhow!(
                "no folding engine is installed.\n{}\n\nOpenDNA does not fold proteins itself and \
                 has nothing to fall back on.",
                missing.join("\n")
            ),
        }),
        _ => Err(Failure {
            fault: Fault::Usage,
            error: anyhow!(
                "{} are both installed. Say which with --engine, because they do not produce the \
                 same structures and choosing for you would leave the result depending on \
                 something nobody wrote down.",
                found.iter().map(|(f, _)| f.name()).collect::<Vec<_>>().join(" and ")
            ),
        }),
    }
}

/// Check the installed engines against a lockfile.
fn verify(path: &Path) -> Result<()> {
    let pinned = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    opendna_engine::check(&pinned, &opendna_engine::discover_all())?;
    eprintln!("The installed engines match {}.", path.display());
    Ok(())
}

/// Whose problem a fold failure is.
///
/// An engine that exited non-zero is deliberately [`Fault::Unexpected`]. The
/// cause could be the installation (a missing CUDA library), the input (a
/// sequence longer than the engine's memory allows) or the engine itself, and
/// nothing here can tell those apart. Filing it under a category a caller
/// would act on would be guessing at the action, which is worse than saying
/// the failure was not anticipated and printing what the engine said.
fn fault_for(error: &FoldError) -> Fault {
    match error {
        // The engine ran, exited cleanly, and did not return the sequence it
        // was given. A check ran and did not hold.
        FoldError::Missing { .. } | FoldError::Ambiguous { .. } | FoldError::Unusable { .. } => Fault::Checks,
        FoldError::Failed { .. } | FoldError::Io { .. } => Fault::Unexpected,
    }
}

/// What the run produced, for someone reading the terminal.
fn report(written: &[(&opendna_engine::Prediction, PathBuf)], work: &Workspace, manifests: bool) {
    println!("{:<24} {:>9} {:>12}  structure", "name", "residues", "mean B");
    for (prediction, output) in written {
        println!(
            "{:<24} {:>9} {:>12.1}  {}",
            prediction.name,
            prediction.residues,
            prediction.mean_ca_b_factor,
            output.display()
        );
    }
    println!(
        "\n{} structure(s). The mean B column is the B-factor column averaged over each chain's \
         Ca atoms,\nwhich is where both engines write their per-residue confidence. It is the \
         engine's estimate of\nits own output, not a measurement.",
        written.len()
    );
    println!("What the engine wrote is kept in {}.", work.raw.display());
    if manifests {
        println!("Each structure has a .manifest.json beside it naming the engine that produced it.");
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// Two names that differ only in punctuation become one filename, and
    /// the second would overwrite the first. Caught before the engine runs,
    /// because the engine is the part that takes an hour.
    #[test]
    fn names_that_collide_are_refused_before_anything_runs() {
        let sequences = vec![("p 1".to_string(), "MQ".to_string()), ("p|1".to_string(), "MA".to_string())];
        let error = names(&sequences).unwrap_err().to_string();
        assert!(error.contains("p_1.pdb"), "{error}");
        assert!(error.contains("\"p 1\"") && error.contains("\"p|1\""), "{error}");
    }

    #[test]
    fn distinct_names_are_kept_in_input_order() {
        let sequences = vec![("beta".to_string(), "MA".to_string()), ("alpha".to_string(), "MQ".to_string())];
        let kept = names(&sequences).unwrap();
        assert_eq!(kept[0].0, "beta");
        assert_eq!(kept[1].0, "alpha");
    }

    #[test]
    fn an_empty_input_is_an_input_error_not_an_empty_run() {
        assert!(names(&[]).unwrap_err().to_string().contains("no sequences"));
    }

    #[test]
    fn a_name_with_nothing_usable_in_it_is_refused() {
        let sequences = vec![(String::new(), "MQ".to_string())];
        assert!(names(&sequences).unwrap_err().to_string().contains("no characters usable"));
    }

    /// An engine that exits non-zero could be any of three different
    /// problems, and the exit code must not claim to know which.
    #[test]
    fn a_failing_engine_is_not_filed_as_an_install_or_input_problem() {
        let failed =
            FoldError::Failed { engine: "esmfold", status: "1".into(), tail: "CUDA out of memory".into() };
        assert_eq!(fault_for(&failed), Fault::Unexpected);
        assert_eq!(fault_for(&failed).code(), 1);
    }

    /// A structure that is not the sequence that was asked for is a check
    /// that ran and did not hold, which is its own exit code.
    #[test]
    fn a_wrong_or_missing_structure_exits_as_a_failed_check() {
        let missing =
            FoldError::Missing { engine: "esmfold", name: "p1".into(), looked: "nothing there".into() };
        assert_eq!(fault_for(&missing), Fault::Checks);
        assert_eq!(fault_for(&missing).code(), 5);
    }
}
