//! External tools OpenDNA drives, held to the same standard as its own heads.
//!
//! Some of what a protein engineer wants is decades of other people's work:
//! folding a sequence, docking a ligand, running a simulation. Writing worse
//! versions of those here would be the mistake this project was rewritten to
//! undo. Driving the real ones is the alternative, and it is only honest if a
//! result that came out of an external tool can be traced as precisely as one
//! that came out of a trained head.
//!
//! So an engine gets the same treatment the heads get:
//!
//! - **Found, or an error naming what was looked for.** No search of the
//!   filesystem, no guessing at a sibling install, and never a download. An
//!   absent engine is a missing dependency, which is the caller's to fix.
//! - **Identified by what it is, not by what it is called.** The version
//!   string it prints and the SHA-256 of the file invoked, both recorded.
//! - **Pinnable.** `engines.lock` is to an engine what `models.lock` is to a
//!   head: write down what was validated, and fail later if it changed.
//!
//! What a digest here does and does not cover is worth being plain about. It
//! identifies the entry point — the executable or script actually run — and
//! not the environment behind it. For a self-contained binary that is nearly
//! everything; for a Python entry point whose behaviour lives in the packages
//! it imports, it is a weak identifier, and the version string is doing most
//! of the work. [`Engine::identification`] says which case an engine is in
//! rather than leaving a reader to assume the stronger one.
//!
//! Nothing in this crate runs an engine for its real purpose. It answers
//! "is it here, which one is it, and is it the one we pinned", and leaves
//! doing science with it to whatever command drives it.

#![forbid(unsafe_code)]
// `doc_markdown` reads the product name as a code identifier. Backticking it
// would make the rendered documentation wrong, which is how opendna-cli
// already resolves this.
#![allow(clippy::doc_markdown)]

use std::fmt;
use std::path::{Path, PathBuf};
use std::process::Command;

mod fold;
mod lock;

pub use fold::{command, command_line, run, sanitise, FoldError, Folder, Prediction};
pub use lock::{check, lockfile, LockError};

/// Why an engine could not be used.
#[derive(Debug, thiserror::Error)]
pub enum EngineError {
    /// Nothing matching was found.
    #[error("{name} is not installed. Looked for {looked_for} on PATH{hint}")]
    NotFound {
        /// The engine's name.
        name: &'static str,
        /// What was looked for, as a readable list.
        looked_for: String,
        /// Where else to look, if the engine declares an override variable.
        hint: String,
    },
    /// Found, but it would not say what it is.
    #[error("{name} at {} would not report its version: {reason}", path.display())]
    NoVersion {
        /// The engine's name.
        name: &'static str,
        /// What was invoked.
        path: PathBuf,
        /// What went wrong.
        reason: String,
    },
    /// Found, but unreadable.
    #[error("{name} at {}: {reason}", path.display())]
    Unreadable {
        /// The engine's name.
        name: &'static str,
        /// What was invoked.
        path: PathBuf,
        /// What went wrong.
        reason: String,
    },
}

/// How well the file we invoke identifies the software behind it.
///
/// Recorded rather than assumed, because the two cases deserve different
/// confidence and a manifest that did not distinguish them would overstate
/// one of them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Identification {
    /// A self-contained executable. Its digest identifies what runs.
    Executable,
    /// An entry point whose behaviour lives elsewhere — a script, a launcher,
    /// an interpreter shim. The digest identifies the entry point only, and
    /// the version string carries most of the identity.
    EntryPoint,
}

impl fmt::Display for Identification {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Executable => write!(f, "executable"),
            Self::EntryPoint => write!(f, "entry point; the digest covers this file, not its environment"),
        }
    }
}

/// An engine OpenDNA knows how to look for.
///
/// Data, not behaviour: naming one here does not mean any command uses it.
#[derive(Debug, Clone, Copy)]
pub struct Requirement {
    /// How it is referred to, and the key in a lockfile.
    pub name: &'static str,
    /// Executable names to try, in order.
    pub binaries: &'static [&'static str],
    /// Arguments that make it print its version.
    pub version_args: &'static [&'static str],
    /// Environment variable naming the executable directly, which wins over
    /// the path search.
    pub override_var: &'static str,
    /// What it would be used for, for a listing a person reads.
    pub purpose: &'static str,
    /// How much the digest of the invoked file identifies.
    pub identification: Identification,
}

/// The engines named so far.
///
/// Being listed means OpenDNA can find and pin it, not that anything drives
/// it yet. The roadmap says which tranche picks each one up.
pub const KNOWN: &[Requirement] = &[
    Requirement {
        name: "esmfold",
        binaries: &["esmfold", "esm-fold"],
        version_args: &["--version"],
        override_var: "OPENDNA_ESMFOLD",
        purpose: "structure prediction from sequence",
        identification: Identification::EntryPoint,
    },
    Requirement {
        name: "colabfold",
        binaries: &["colabfold_batch"],
        version_args: &["--version"],
        override_var: "OPENDNA_COLABFOLD",
        purpose: "structure prediction from sequence",
        identification: Identification::EntryPoint,
    },
    Requirement {
        name: "vina",
        binaries: &["vina", "vina_split"],
        version_args: &["--version"],
        override_var: "OPENDNA_VINA",
        purpose: "ligand docking",
        identification: Identification::Executable,
    },
    Requirement {
        name: "openmm",
        // OpenMM is a Python library rather than a program, so what gets
        // pinned is the interpreter that can import it. `check` below
        // confirms the import rather than trusting the name.
        binaries: &["python3", "python"],
        version_args: &["-c", "import openmm; print('openmm', openmm.__version__)"],
        override_var: "OPENDNA_OPENMM_PYTHON",
        purpose: "molecular dynamics",
        identification: Identification::EntryPoint,
    },
];

/// One engine, found and identified.
#[derive(Debug, Clone)]
pub struct Engine {
    /// Which requirement it satisfies.
    pub name: &'static str,
    /// The file that would be invoked.
    pub path: PathBuf,
    /// The first line it printed when asked for its version, verbatim.
    ///
    /// Not parsed into a version number. Third-party tools word this however
    /// they like and change it between releases; a manifest recording what
    /// the tool actually said can be checked by eye, while a parsed number
    /// can be wrong without anyone noticing.
    pub version: String,
    /// SHA-256 of the file at `path`.
    pub sha256: String,
    /// How much that digest identifies.
    pub identification: Identification,
}

impl Engine {
    /// One line for a person reading a listing.
    pub fn describe(&self) -> String {
        format!("{}  {}  {}", self.version, &self.sha256[..16], self.path.display())
    }
}

/// Look for one engine.
///
/// The override variable wins, then the first of `binaries` found on `PATH`.
/// Nothing else is searched: a filesystem crawl would be slow, would surprise
/// whoever ran it, and could find a second copy of a tool that is not the one
/// they meant.
///
/// # Errors
/// [`EngineError::NotFound`] when nothing matches, and the other variants when
/// something matches but cannot be identified. Being unable to identify it is
/// an error rather than a shrug, because an engine that will not say what it
/// is cannot be recorded in a manifest or pinned.
pub fn discover(requirement: &Requirement) -> Result<Engine, EngineError> {
    let path = locate(requirement).ok_or_else(|| EngineError::NotFound {
        name: requirement.name,
        looked_for: requirement.binaries.join(" or "),
        hint: format!(", or set {} to it", requirement.override_var),
    })?;

    let version = ask_version(requirement, &path)?;
    let bytes = std::fs::read(&path).map_err(|e| EngineError::Unreadable {
        name: requirement.name,
        path: path.clone(),
        reason: e.to_string(),
    })?;

    Ok(Engine {
        name: requirement.name,
        path,
        version,
        sha256: digest(&bytes),
        identification: requirement.identification,
    })
}

/// Look for every known engine, reporting each one's outcome.
pub fn discover_all() -> Vec<(&'static Requirement, Result<Engine, EngineError>)> {
    KNOWN.iter().map(|r| (r, discover(r))).collect()
}

/// Where the engine is, if it is anywhere this looks.
fn locate(requirement: &Requirement) -> Option<PathBuf> {
    if let Some(explicit) = std::env::var_os(requirement.override_var) {
        let path = PathBuf::from(explicit);
        // An override that points at nothing is a mistake worth surfacing as
        // "not found" rather than silently falling through to the path search
        // and running a different copy than the one that was asked for.
        return path.is_file().then_some(path);
    }
    requirement.binaries.iter().find_map(|name| which(name))
}

/// The first match for `name` on `PATH`, with the platform's extensions.
fn which(name: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    let extensions: Vec<String> = if cfg!(windows) {
        std::env::var("PATHEXT")
            .unwrap_or_else(|_| ".EXE;.CMD;.BAT;.COM".into())
            .split(';')
            .map(str::to_ascii_lowercase)
            .collect()
    } else {
        vec![String::new()]
    };
    for directory in std::env::split_paths(&path) {
        for extension in &extensions {
            let candidate = directory.join(format!("{name}{extension}"));
            if candidate.is_file() {
                return Some(candidate);
            }
        }
    }
    None
}

/// Ask the engine what it is.
fn ask_version(requirement: &Requirement, path: &Path) -> Result<String, EngineError> {
    let output = Command::new(path).args(requirement.version_args).output().map_err(|e| {
        EngineError::NoVersion { name: requirement.name, path: path.to_path_buf(), reason: e.to_string() }
    })?;

    // A tool that exited non-zero did not answer the question, whatever it
    // wrote while failing. Reporting openmm as installed at version
    // "Traceback (most recent call last):" is what happens without this: the
    // interpreter is on PATH, the import fails, and the first line of the
    // traceback reads like a version string to anything that only looks at
    // the output. All four declared engines exit 0 when they succeed.
    if !output.status.success() {
        // Either stream: a tool that fails while printing to stdout has still
        // told whoever reads this why, and dropping it would leave them with
        // an exit code and nothing else.
        let err = String::from_utf8_lossy(&output.stderr).into_owned();
        let out = String::from_utf8_lossy(&output.stdout).into_owned();
        let first = [err.as_str(), out.as_str()]
            .into_iter()
            .flat_map(|stream| stream.lines())
            .map(str::trim)
            .find(|l| !l.is_empty())
            .unwrap_or("nothing");
        return Err(EngineError::NoVersion {
            name: requirement.name,
            path: path.to_path_buf(),
            reason: format!("it exited {} and said: {first}", status_of(output.status)),
        });
    }

    // Some tools print their version to standard error while succeeding.
    let mut text = String::from_utf8_lossy(&output.stdout).into_owned();
    if text.trim().is_empty() {
        text = String::from_utf8_lossy(&output.stderr).into_owned();
    }
    let line = text.lines().map(str::trim).find(|l| !l.is_empty()).unwrap_or("").to_string();
    if line.is_empty() {
        return Err(EngineError::NoVersion {
            name: requirement.name,
            path: path.to_path_buf(),
            reason: format!("it printed nothing when asked with {:?}", requirement.version_args),
        });
    }
    Ok(line)
}

/// An exit status as a short string, since `ExitStatus` renders differently
/// across platforms and this ends up in a message somebody reads.
pub(crate) fn status_of(status: std::process::ExitStatus) -> String {
    status.code().map_or_else(|| "abnormally".to_string(), |code| code.to_string())
}

/// SHA-256, hex.
fn digest(bytes: &[u8]) -> String {
    use sha2::{Digest as _, Sha256};
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every engine's declaration has to be usable, which is easy to get
    /// wrong by pasting one and editing half of it.
    #[test]
    fn every_known_engine_is_declared_completely() {
        for requirement in KNOWN {
            assert!(!requirement.name.is_empty(), "an engine with no name");
            assert!(!requirement.binaries.is_empty(), "{}: no binaries to look for", requirement.name);
            assert!(!requirement.version_args.is_empty(), "{}: nothing to ask it with", requirement.name);
            assert!(
                requirement.override_var.starts_with("OPENDNA_"),
                "{}: an override variable should be namespaced, got {}",
                requirement.name,
                requirement.override_var
            );
            assert!(!requirement.purpose.is_empty(), "{}: no stated purpose", requirement.name);
        }
    }

    #[test]
    fn engine_names_are_unique() {
        let mut names: Vec<&str> = KNOWN.iter().map(|r| r.name).collect();
        names.sort_unstable();
        let before = names.len();
        names.dedup();
        assert_eq!(names.len(), before, "two engines share a name");
    }

    /// An engine that is not there is an error that says what was looked for,
    /// not a silent absence.
    #[test]
    fn a_missing_engine_names_what_it_looked_for() {
        let requirement = Requirement {
            name: "nothing",
            binaries: &["opendna-definitely-not-a-real-program"],
            version_args: &["--version"],
            override_var: "OPENDNA_NOTHING",
            purpose: "testing",
            identification: Identification::Executable,
        };
        let error = discover(&requirement).unwrap_err();
        let text = error.to_string();
        assert!(text.contains("opendna-definitely-not-a-real-program"), "{text}");
        assert!(text.contains("OPENDNA_NOTHING"), "{text}");
    }

    /// An override pointing at nothing must not fall through to the path
    /// search: running a different copy than the one asked for is worse than
    /// reporting the engine as absent.
    #[test]
    fn an_override_that_points_at_nothing_does_not_fall_back() {
        let requirement = Requirement {
            name: "ghost",
            // Something that does exist on every platform this builds for,
            // so a fall-through would be visible as a success.
            #[cfg(windows)]
            binaries: &["cmd"],
            #[cfg(not(windows))]
            binaries: &["sh"],
            version_args: &["--version"],
            override_var: "OPENDNA_TEST_GHOST",
            purpose: "testing",
            identification: Identification::Executable,
        };
        std::env::set_var("OPENDNA_TEST_GHOST", "/definitely/not/here/at/all");
        let found = locate(&requirement);
        std::env::remove_var("OPENDNA_TEST_GHOST");
        assert!(found.is_none(), "an override pointing at nothing fell back to PATH: {found:?}");
    }

    /// An engine that cannot say what it is must not be reported as present.
    ///
    /// This is the openmm case: an interpreter is on PATH, the import fails,
    /// and the first line of the traceback looks like a version string to
    /// anything that only reads the output.
    #[test]
    fn a_tool_that_fails_while_talking_is_not_identified() {
        #[cfg(windows)]
        let requirement = Requirement {
            name: "noisy-failure",
            binaries: &["cmd"],
            version_args: &["/c", "echo Traceback (most recent call last): & exit 1"],
            override_var: "OPENDNA_NOISY",
            purpose: "testing",
            identification: Identification::Executable,
        };
        #[cfg(not(windows))]
        let requirement = Requirement {
            name: "noisy-failure",
            binaries: &["sh"],
            version_args: &["-c", "echo 'Traceback (most recent call last):'; exit 1"],
            override_var: "OPENDNA_NOISY",
            purpose: "testing",
            identification: Identification::Executable,
        };

        let error = discover(&requirement).expect_err("a failing tool must not identify");
        let text = error.to_string();
        assert!(text.contains("would not report its version"), "{text}");
        // And the message has to carry what it actually said, or whoever
        // reads it cannot tell a missing dependency from a broken one.
        assert!(text.contains("Traceback"), "{text}");
    }

    #[test]
    fn the_digest_is_sha256_of_the_bytes() {
        // Of the empty input, which is the one value worth knowing by sight.
        assert_eq!(digest(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    }
}
