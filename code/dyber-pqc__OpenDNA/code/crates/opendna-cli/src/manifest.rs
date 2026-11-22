//! A machine-readable record of one run.
//!
//! The `#` comments on an output file say which head produced it, which is
//! enough for a person reading the file. It is not enough to reproduce the
//! run or to answer a question about it a year later: which input, which
//! bytes of input, which build of the head, which embedding model, on what,
//! for how long. Those are the questions a reviewer asks, and reconstructing
//! them from a shell history nobody kept is how a result stops being
//! defensible.
//!
//! So every run that writes an output file writes a manifest beside it. It is
//! JSON because something else has to read it — a pipeline collecting
//! provenance, a lab notebook, a script checking that two runs used the same
//! model. One run, one file, no appending: a manifest that several processes
//! write to is a manifest with a race in it.
//!
//! Inputs and outputs are hashed from the bytes in memory rather than by
//! re-reading the files, so the digests describe what was actually used even
//! when the input came from a pipe and has no path at all.
//!
//! The embedding model is the one thing not hashed by default. It is 568 MB
//! and digesting it takes about as long as a small run, so the manifest
//! records its path and size, uses the digest a staged bundle already wrote
//! down when there is one, and otherwise says plainly that it was not
//! recorded. An absent field with a reason is honest; a missing field someone
//! assumes was checked is not.

use anyhow::{Context, Result};
use opendna_ml::heads::{HeadFiles, HeadMetadata};
use serde_json::json;
use std::path::{Path, PathBuf};

/// Where the manifest for an output file goes.
///
/// `out.csv` becomes `out.manifest.json`, matching how the CAID path already
/// places its provenance sidecar, so the two conventions do not differ.
pub fn path_for(output: &Path) -> PathBuf {
    output.with_extension("manifest.json")
}

/// A path as the manifest states it.
///
/// Forward slashes throughout, even on Windows. A path joined from a
/// forward-slash argument comes out mixed — `C:/work/bundle\\disorder.safetensors`
/// — which is valid for the operating system and awkward for everything that
/// reads the manifest afterwards.
fn shown(path: &Path) -> String {
    path.display().to_string().replace('\\', "/")
}

/// SHA-256 of bytes already in hand.
pub fn sha256(bytes: &[u8]) -> String {
    use sha2::{Digest, Sha256};
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

/// What the run read.
pub struct Input<'a> {
    /// The path, or `-` for standard input.
    pub label: &'a str,
    /// The bytes as read, which is what gets hashed.
    pub bytes: &'a [u8],
    /// How many sequences were parsed out of it.
    pub sequences: usize,
    /// Total residues across them.
    pub residues: usize,
}

/// Which embedding model was used, and how well it can be identified.
pub struct Embedder {
    /// The Hugging Face repository the head names.
    pub repo: &'static str,
    /// Where the files came from, for a reader who has to find them again.
    pub source: &'static str,
    /// The weights file, when it is a path we resolved ourselves.
    pub weights: Option<PathBuf>,
    /// Its size in bytes, which is free to record and catches a truncated
    /// copy even without a digest.
    pub bytes: Option<u64>,
    /// Its SHA-256, when something already knew it.
    pub digest: Option<String>,
    /// Where that digest came from, or why there is none.
    pub digest_source: &'static str,
}

/// Everything one run should be able to say about itself.
pub struct Run<'a> {
    /// Which property was predicted.
    pub task: &'a str,
    /// The device description, as reported to the user.
    pub device: &'a str,
    /// When the run started, in UTC. Local time in a manifest compared across
    /// machines is a trap.
    pub started: chrono::DateTime<chrono::Utc>,
    /// How long it took, wall clock.
    pub elapsed: std::time::Duration,
    /// Padded tokens per embedding batch.
    ///
    /// Recorded because it changes the results: padding alters the order of
    /// floating-point sums, so two runs of the same input at different
    /// settings produce files that differ in the last digit of a few values.
    /// Without this the manifest would claim to identify a run it could not
    /// distinguish from a slightly different one.
    pub batch_tokens: usize,
    /// Sequences per chunk, for a run that used a checkpoint.
    ///
    /// Recorded for the same reason as `batch_tokens`: it decides which
    /// sequences are embedded together, and that changes the last digit of a
    /// few values. A run without a checkpoint has no chunking and records
    /// nothing here, which is a different run from one chunked at any size.
    pub chunk: Option<usize>,
    /// The input.
    pub input: Input<'a>,
    /// Where the table was written, if anywhere.
    pub output: Option<&'a Path>,
    /// The rendered table, hashed as the output.
    pub rendered: &'a str,
    /// How many data rows it holds.
    pub rows: usize,
    /// The mutation file, for stability runs.
    pub mutations: Option<&'a Path>,
    /// The head's recorded identity.
    pub head: &'a HeadMetadata,
    /// Where the head's files are.
    pub files: &'a HeadFiles,
    /// The embedding model.
    pub embedder: &'a Embedder,
    /// The card's limitations, carried so the manifest alone is enough to
    /// state what the numbers do not cover.
    pub limitations: &'a [String],
    /// The head's licence, as its card declares it.
    ///
    /// Part of what a methods section has to say about a number: one of these
    /// heads is trained on share-alike data and released on the same terms,
    /// and a result travelling without that is a result somebody may
    /// redistribute without knowing they have an obligation.
    pub license: Option<opendna_ml::heads::License>,
}

/// Render the manifest.
pub fn render(run: &Run<'_>) -> String {
    let doc = json!({
        "opendna_version": env!("CARGO_PKG_VERSION"),
        "manifest_version": 1,
        "task": run.task,
        "run": {
            "started_utc": run.started.to_rfc3339_opts(chrono::SecondsFormat::Secs, true),
            "seconds": (run.elapsed.as_millis() as f64) / 1000.0,
            "device": run.device,
            "batch_tokens": run.batch_tokens,
            "chunk": run.chunk,
            // The arguments as given, so the run can be repeated verbatim.
            // Not the resolved paths: what was typed is what a reader needs to
            // recognise, and the resolved values are recorded separately.
            "arguments": std::env::args().skip(1).collect::<Vec<_>>(),
        },
        "input": {
            "path": run.input.label.replace('\\', "/"),
            "sha256": sha256(run.input.bytes),
            "bytes": run.input.bytes.len(),
            "sequences": run.input.sequences,
            "residues": run.input.residues,
        },
        "mutations": run.mutations.map(|p| json!({ "path": shown(p) })),
        "output": run.output.map(|p| json!({
            "path": shown(p),
            "sha256": sha256(run.rendered.as_bytes()),
            "bytes": run.rendered.len(),
            "rows": run.rows,
        })),
        "head": {
            "task": run.head.task,
            "version": run.head.version,
            "build": run.head.build,
            "weights_sha256": run.head.weights_sha256,
            "weights_path": shown(&run.files.weights),
            "with_burial": run.head.with_burial,
            "license": run.license.as_ref().map(|l| json!({
                "name": l.label(),
                "link": l.link,
                "share_alike": l.share_alike(),
            })),
        },
        "embedding_model": {
            "repo": run.embedder.repo,
            "source": run.embedder.source,
            "weights_path": run.embedder.weights.as_deref().map(shown),
            "bytes": run.embedder.bytes,
            "sha256": run.embedder.digest,
            "sha256_source": run.embedder.digest_source,
        },
        "limitations": run.limitations,
    });
    // Pretty-printed: a manifest ends up in a diff and in a code review.
    format!("{}\n", serde_json::to_string_pretty(&doc).unwrap_or_else(|e| format!("{{\"error\":\"{e}\"}}")))
}

/// Write the manifest beside `output`.
///
/// # Errors
/// If the file cannot be written. This is an error rather than a warning: a
/// run that produced a table and silently failed to record what made it has
/// produced a number nobody can place, which is the thing this exists to
/// prevent.
pub fn write(run: &Run<'_>, output: &Path) -> Result<PathBuf> {
    let path = path_for(output);
    std::fs::write(&path, render(run)).with_context(|| format!("writing {}", path.display()))?;
    Ok(path)
}

/// What one folded structure should be able to say about itself.
///
/// Nothing here is a head or an embedding model, so none of [`Run`]'s fields
/// fit. What has to be identified instead is the engine — which one, which
/// version it reported, the digest of the file invoked and how much that
/// digest covers — and the sequence the returned chain was checked against.
pub struct FoldRun<'a> {
    /// Which engine was driven, by the name OpenDNA uses for it.
    pub folder: &'a str,
    /// The engine as discovered.
    pub engine: &'a opendna_engine::Engine,
    /// The command as run, program first.
    pub command: &'a [String],
    /// When it started, UTC.
    pub started: chrono::DateTime<chrono::Utc>,
    /// Wall clock for the whole batch, not for this one structure: the
    /// engines fold a file at a time and do not report per-sequence timings,
    /// and dividing by the sequence count would be an invention.
    pub elapsed: std::time::Duration,
    /// The input file, hashed whole, so every structure from one run points
    /// back at the same bytes.
    pub input: Input<'a>,
    /// The FASTA name of this sequence.
    pub name: &'a str,
    /// The sequence that was sent, as written.
    pub sequence: &'a str,
    /// Where the structure was written.
    pub output: &'a Path,
    /// Its bytes, which are the engine's bytes unchanged.
    pub bytes: &'a [u8],
    /// The engine's own file, kept rather than deleted.
    pub source: &'a Path,
    /// Mean B-factor over the chain's Ca atoms.
    pub mean_ca_b_factor: f64,
    /// Residues in the returned chain.
    pub residues: usize,
}

/// Render a fold manifest.
pub fn render_fold(run: &FoldRun<'_>) -> String {
    let doc = json!({
        "opendna_version": env!("CARGO_PKG_VERSION"),
        "manifest_version": 1,
        "task": "fold",
        "run": {
            "started_utc": run.started.to_rfc3339_opts(chrono::SecondsFormat::Secs, true),
            "seconds": (run.elapsed.as_millis() as f64) / 1000.0,
            "seconds_covers": "the whole batch; the engines do not report per-sequence timings",
            "arguments": std::env::args().skip(1).collect::<Vec<_>>(),
            "engine_command": run.command,
        },
        "engine": {
            "name": run.folder,
            "version": run.engine.version,
            "path": shown(&run.engine.path),
            "sha256": run.engine.sha256,
            "identification": match run.engine.identification {
                opendna_engine::Identification::Executable => "executable",
                opendna_engine::Identification::EntryPoint => "entry-point",
            },
            "identification_covers": match run.engine.identification {
                opendna_engine::Identification::Executable =>
                    "the digest identifies the program that ran",
                opendna_engine::Identification::EntryPoint =>
                    "the digest identifies the entry point only, not the environment behind it; \
                     the version string carries the rest of the identity",
            },
        },
        "input": {
            "path": run.input.label.replace('\\', "/"),
            "sha256": sha256(run.input.bytes),
            "bytes": run.input.bytes.len(),
            "sequences": run.input.sequences,
            "residues": run.input.residues,
        },
        "sequence": {
            "name": run.name,
            "residues": run.sequence.len(),
            "one_letter": run.sequence,
        },
        "output": {
            "path": shown(run.output),
            "sha256": sha256(run.bytes),
            "bytes": run.bytes.len(),
            "engine_file": shown(run.source),
            "chain_residues": run.residues,
            "mean_ca_b_factor": (run.mean_ca_b_factor * 1000.0).round() / 1000.0,
        },
        "verification": "the chain in the output was parsed and its sequence compared with the \
                         input residue by residue; the file was matched to this sequence by those \
                         residues rather than by its filename",
        "limitations": [
            "This is a predicted model, not an experimental structure.",
            "OpenDNA did not fold anything. It ran the engine named above and checked that the \
             chain that came back is the sequence that was sent.",
            "mean_ca_b_factor is the mean of the B-factor column over the chain's Ca atoms. \
             ESMFold and ColabFold both write their per-residue confidence (pLDDT) there, but \
             nothing here verified what this engine meant by that column, and the value is the \
             engine's own estimate of its own output.",
        ],
    });
    format!("{}\n", serde_json::to_string_pretty(&doc).unwrap_or_else(|e| format!("{{\"error\":\"{e}\"}}")))
}

/// Write the fold manifest beside `output`.
///
/// # Errors
/// If the file cannot be written, for the same reason as [`fn@write`].
pub fn write_fold(run: &FoldRun<'_>, output: &Path) -> Result<PathBuf> {
    let path = path_for(output);
    std::fs::write(&path, render_fold(run)).with_context(|| format!("writing {}", path.display()))?;
    Ok(path)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use opendna_ml::esm2::ModelSize;

    fn metadata() -> HeadMetadata {
        HeadMetadata {
            task: "disorder".into(),
            embedding_model: ModelSize::Large150M.repo().into(),
            embedding_dim: 640,
            version: "1.0.0".into(),
            build: Some("adbdc038cbd171ff".into()),
            weights_sha256: Some("ff".repeat(32)),
            with_burial: false,
            calibration: None,
            abstention: None,
            intervals: None,
        }
    }

    fn embedder() -> Embedder {
        Embedder {
            repo: ModelSize::Large150M.repo(),
            source: "staged bundle",
            weights: Some(PathBuf::from("bundle/embedder/model.safetensors")),
            bytes: Some(595_257_706),
            digest: Some("ab".repeat(32)),
            digest_source: "sha256sums.txt in the staged bundle",
        }
    }

    fn manifest(output: Option<&Path>) -> serde_json::Value {
        let meta = metadata();
        let files = HeadFiles::in_dir(Path::new("models"), "disorder");
        let emb = embedder();
        let limitations = ["not a measurement".to_string()];
        let text = render(&Run {
            task: "disorder",
            device: "CPU",
            started: chrono::DateTime::from_timestamp(1_763_000_000, 0).unwrap(),
            elapsed: std::time::Duration::from_millis(2500),
            batch_tokens: 0,
            chunk: None,
            input: Input { label: "seqs.fasta", bytes: b">a\nMQIF\n", sequences: 1, residues: 4 },
            output,
            rendered: "name,position\na,1\n",
            rows: 1,
            mutations: None,
            head: &meta,
            files: &files,
            embedder: &emb,
            limitations: &limitations,
            license: Some(opendna_ml::heads::License {
                identifier: "other".into(),
                name: Some("cc-by-sa-2.1-jp".into()),
                link: None,
            }),
        });
        serde_json::from_str(&text).unwrap()
    }

    /// The manifest has to identify the input by its bytes, not its name: the
    /// same filename holding different sequences is the ordinary way a
    /// reproduction quietly fails.
    #[test]
    fn the_input_is_identified_by_its_contents() {
        let m = manifest(Some(Path::new("out.csv")));
        assert_eq!(m["input"]["sha256"], sha256(b">a\nMQIF\n"));
        assert_eq!(m["input"]["bytes"], 8);
        assert_eq!(m["input"]["sequences"], 1);
        assert_eq!(m["input"]["residues"], 4);
        assert_ne!(sha256(b">a\nMQIF\n"), sha256(b">a\nMQIFV\n"), "different bytes, different digest");
    }

    #[test]
    fn the_head_is_identified_by_version_build_and_hash() {
        let m = manifest(Some(Path::new("out.csv")));
        assert_eq!(m["head"]["version"], "1.0.0");
        assert_eq!(m["head"]["build"], "adbdc038cbd171ff");
        assert_eq!(m["head"]["weights_sha256"], "ff".repeat(32));
    }

    /// An unhashed embedding model must say so rather than leaving a field
    /// out, because a reader cannot tell a missing field from a checked one.
    #[test]
    fn an_unhashed_embedding_model_says_why() {
        let meta = metadata();
        let files = HeadFiles::in_dir(Path::new("models"), "disorder");
        let emb = Embedder {
            repo: ModelSize::Large150M.repo(),
            source: "Hugging Face cache",
            weights: None,
            bytes: None,
            digest: None,
            digest_source: "not recorded: hashing 568 MB costs more than most runs",
        };
        let text = render(&Run {
            task: "disorder",
            device: "CPU",
            started: chrono::Utc::now(),
            elapsed: std::time::Duration::from_secs(1),
            batch_tokens: 8192,
            chunk: Some(64),
            input: Input { label: "-", bytes: b">a\nM\n", sequences: 1, residues: 1 },
            output: None,
            rendered: "",
            rows: 0,
            mutations: None,
            head: &meta,
            files: &files,
            embedder: &emb,
            limitations: &[],
            license: None,
        });
        let m: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert!(m["embedding_model"]["sha256"].is_null());
        assert!(
            m["embedding_model"]["sha256_source"].as_str().unwrap().contains("not recorded"),
            "{}",
            m["embedding_model"]["sha256_source"]
        );
        // Standard input has no path, and the manifest says so rather than
        // inventing one.
        assert_eq!(m["input"]["path"], "-");
        assert!(m["output"].is_null());
    }

    /// Paths are stated with forward slashes so a manifest written on
    /// Windows reads the same as one written anywhere else.
    #[test]
    fn paths_do_not_mix_separators() {
        let m = manifest(Some(Path::new("out.csv")));
        for field in [&m["head"]["weights_path"], &m["output"]["path"], &m["input"]["path"]] {
            let text = field.as_str().unwrap_or_default();
            assert!(!text.contains('\\'), "{text} contains a backslash");
        }
        assert_eq!(shown(Path::new("a/b")), "a/b");
    }

    #[test]
    fn the_manifest_goes_beside_the_output() {
        assert_eq!(path_for(Path::new("out.csv")), PathBuf::from("out.manifest.json"));
        assert_eq!(path_for(Path::new("a/b/ddg.tsv")), PathBuf::from("a/b/ddg.manifest.json"));
        // A name with no extension still gets one rather than overwriting the
        // output it describes.
        assert_eq!(path_for(Path::new("results")), PathBuf::from("results.manifest.json"));
    }

    /// A run is not identified unless everything that changed its numbers is
    /// written down, and the batch size changes them.
    /// A share-alike licence that does not reach the manifest is one the
    /// reader of a result file has no way to discover.
    #[test]
    fn a_share_alike_licence_is_recorded_and_flagged() {
        let m = manifest(Some(Path::new("out.csv")));
        assert_eq!(m["head"]["license"]["name"], "cc-by-sa-2.1-jp");
        assert_eq!(m["head"]["license"]["share_alike"], true);
    }

    #[test]
    fn the_batch_size_is_recorded() {
        let m = manifest(Some(Path::new("out.csv")));
        assert_eq!(m["run"]["batch_tokens"], 0);
    }

    #[test]
    fn the_time_is_utc_and_the_duration_is_seconds() {
        let m = manifest(Some(Path::new("out.csv")));
        let started = m["run"]["started_utc"].as_str().unwrap();
        assert!(started.ends_with('Z'), "{started}");
        assert_eq!(m["run"]["seconds"], 2.5);
        assert_eq!(m["manifest_version"], 1);
    }
}
