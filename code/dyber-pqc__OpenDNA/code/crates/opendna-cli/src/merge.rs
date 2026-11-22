//! Putting the parts back together, and checking they belong together.
//!
//! `cat part-*.csv` almost works. It gets the rows right and it cannot notice
//! the thing most worth noticing: that the parts were not produced by the same
//! software.
//!
//! On a cluster that is not a hypothetical. Nodes get upgraded one at a time,
//! so half a run can land on machines with a GPU and half on machines without;
//! a head can be re-downloaded between tasks; somebody can tune
//! `--batch-tokens` on a queue that is busy. Each of those changes the numbers
//! in the last digit, and the resulting file is a mixture nobody can describe.
//! The parts carry run manifests that record exactly these things, so the
//! merge reads them and refuses a mixture.
//!
//! Refuses, rather than warns. A warning on a cluster goes into a log nobody
//! reads, and the file it was warning about goes into the analysis.

use anyhow::{bail, Context, Result};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

/// The parts of a manifest that must agree for two result files to be merged.
///
/// Anything that changes a number belongs here. The input digest deliberately
/// does not: the whole point is that the parts had different inputs.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Agreement {
    /// Which property was predicted.
    task: String,
    /// The head's build id.
    build: String,
    /// The embedding model.
    embedding_model: String,
    /// What it ran on.
    device: String,
    /// Padded tokens per batch.
    batch_tokens: String,
    /// Sequences per checkpoint chunk, or `unrecorded` for a part that was
    /// run without one. Chunking decides which sequences are embedded
    /// together, so two parts chunked differently differ in the last digit of
    /// some values, which is the mixture this refuses.
    chunk: String,
}

impl Agreement {
    /// Read the fields that must match out of a manifest.
    fn read(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let doc: serde_json::Value =
            serde_json::from_str(&text).with_context(|| format!("parsing {}", path.display()))?;
        let field = |pointer: &str| -> String {
            doc.pointer(pointer).map_or_else(
                || "unrecorded".to_string(),
                |v| v.as_str().map_or_else(|| v.to_string(), str::to_owned),
            )
        };
        Ok(Self {
            task: field("/task"),
            build: field("/head/build"),
            embedding_model: field("/embedding_model/repo"),
            device: field("/run/device"),
            batch_tokens: field("/run/batch_tokens"),
            chunk: field("/run/chunk"),
        })
    }

    /// How to show it in a disagreement report.
    fn describe(&self) -> String {
        format!(
            "{} head {} on {} via {} at batch-tokens {}, chunk {}",
            self.task, self.build, self.device, self.embedding_model, self.batch_tokens, self.chunk
        )
    }
}

/// Split a result file into its `#` comments, header and data rows.
fn parts_of(text: &str) -> (Vec<&str>, Option<&str>, Vec<&str>) {
    let mut comments = Vec::new();
    let mut header = None;
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            if header.is_none() {
                comments.push(line);
            }
        } else if line.trim().is_empty() {
        } else if header.is_none() {
            header = Some(line);
        } else {
            rows.push(line);
        }
    }
    (comments, header, rows)
}

/// Merge result files, refusing a mixture.
///
/// Inputs are merged in the order given, so a caller sorting its parts by name
/// gets the original order back.
///
/// # Errors
/// If a file cannot be read, if the headers differ, or if the manifests show
/// the parts were produced by different software or on different devices.
pub fn run(inputs: &[PathBuf], output: Option<&Path>, require_manifests: bool) -> Result<()> {
    if inputs.is_empty() {
        bail!("no input files given");
    }

    // Every part's provenance, grouped so a disagreement can name who differs
    // rather than just that somebody does.
    let mut groups: BTreeMap<Agreement, Vec<String>> = BTreeMap::new();
    let mut missing = Vec::new();
    for input in inputs {
        let manifest = crate::manifest::path_for(input);
        if manifest.is_file() {
            groups.entry(Agreement::read(&manifest)?).or_default().push(input.display().to_string());
        } else {
            missing.push(input.display().to_string());
        }
    }
    if !missing.is_empty() && require_manifests {
        bail!(
            "these parts have no manifest beside them, so there is no way to tell whether they \
             were produced by the same software as the others:\n  {}\n\
             Re-run them without --no-manifest, or pass --allow-unverified to merge anyway.",
            missing.join("\n  ")
        );
    }
    if groups.len() > 1 {
        let mut report = String::from(
            "refusing to merge: these parts were not produced the same way, so the result would \
             be a mixture nobody can describe.\n",
        );
        for (agreement, files) in &groups {
            let _ = writeln!(report, "  {}", agreement.describe());
            for file in files {
                let _ = writeln!(report, "      {file}");
            }
        }
        bail!("{}", report.trim_end());
    }

    let mut out = String::new();
    let mut header: Option<String> = None;
    let mut total_rows = 0usize;
    for input in inputs {
        let text = std::fs::read_to_string(input).with_context(|| format!("reading {}", input.display()))?;
        let (comments, this_header, rows) = parts_of(&text);
        let Some(this_header) = this_header else {
            bail!("{} has no header row", input.display());
        };
        match &header {
            // The first part's comments describe the run, and every part's
            // are the same once the manifests agree, so they are kept once.
            None => {
                for comment in &comments {
                    let _ = writeln!(out, "{comment}");
                }
                let _ = writeln!(out, "{this_header}");
                header = Some(this_header.to_string());
            }
            Some(first) if first != this_header => {
                bail!("{} has columns `{this_header}` but the first part has `{first}`", input.display());
            }
            Some(_) => {}
        }
        for row in rows {
            let _ = writeln!(out, "{row}");
            total_rows += 1;
        }
    }

    if let Some(path) = output {
        std::fs::write(path, &out).with_context(|| format!("writing {}", path.display()))?;
        eprintln!("Merged {} parts, {total_rows} rows, into {}", inputs.len(), path.display());
    } else {
        use std::io::Write as _;
        std::io::stdout().write_all(out.as_bytes()).context("writing to standard output")?;
        eprintln!("Merged {} parts, {total_rows} rows", inputs.len());
    }
    if !missing.is_empty() {
        eprintln!("{} of {} parts had no manifest and could not be checked.", missing.len(), inputs.len());
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("opendna-merge-{name}-{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// Write a part with a manifest describing how it was produced.
    fn part(dir: &Path, name: &str, rows: &[&str], device: &str, build: &str) -> PathBuf {
        let path = dir.join(format!("{name}.csv"));
        let mut body = String::from("# head_version: 1.1.0\nname,position,probability\n");
        for row in rows {
            body.push_str(row);
            body.push('\n');
        }
        std::fs::write(&path, body).unwrap();
        let manifest = serde_json::json!({
            "task": "disorder",
            "head": { "build": build },
            "embedding_model": { "repo": "facebook/esm2_t30_150M_UR50D" },
            "run": { "device": device, "batch_tokens": 0 },
        });
        std::fs::write(crate::manifest::path_for(&path), serde_json::to_string_pretty(&manifest).unwrap())
            .unwrap();
        path
    }

    #[test]
    fn parts_produced_the_same_way_merge_in_order() {
        let dir = scratch("ok");
        let a = part(&dir, "part-0000", &["a,1,0.10", "a,2,0.20"], "CPU", "abc");
        let b = part(&dir, "part-0001", &["b,1,0.30"], "CPU", "abc");
        let out = dir.join("all.csv");
        run(&[a, b], Some(&out), true).unwrap();

        let text = std::fs::read_to_string(&out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "# head_version: 1.1.0", "comments kept once");
        assert_eq!(lines[1], "name,position,probability");
        assert_eq!(&lines[2..], ["a,1,0.10", "a,2,0.20", "b,1,0.30"]);
        assert_eq!(text.matches("name,position").count(), 1, "one header only");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// The case `cat` cannot catch: half the run on a GPU, half not.
    #[test]
    fn parts_from_different_devices_are_refused() {
        let dir = scratch("device");
        let a = part(&dir, "part-0000", &["a,1,0.10"], "CPU", "abc");
        let b = part(&dir, "part-0001", &["b,1,0.30"], "CUDA GPU", "abc");
        let err = run(&[a, b], None, true).unwrap_err().to_string();
        assert!(err.contains("not produced the same way"), "{err}");
        assert!(err.contains("CPU") && err.contains("CUDA GPU"), "{err}");
        assert!(err.contains("part-0000") && err.contains("part-0001"), "the report must name them: {err}");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// And a head re-downloaded or recalibrated between tasks.
    #[test]
    fn parts_from_different_head_builds_are_refused() {
        let dir = scratch("build");
        let a = part(&dir, "part-0000", &["a,1,0.10"], "CPU", "abc");
        let b = part(&dir, "part-0001", &["b,1,0.30"], "CPU", "zzz");
        let err = run(&[a, b], None, true).unwrap_err().to_string();
        assert!(err.contains("not produced the same way"), "{err}");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_part_without_a_manifest_stops_the_merge_unless_allowed() {
        let dir = scratch("nomanifest");
        let a = part(&dir, "part-0000", &["a,1,0.10"], "CPU", "abc");
        let b = dir.join("part-0001.csv");
        std::fs::write(&b, "name,position,probability\nb,1,0.30\n").unwrap();

        let err = run(&[a.clone(), b.clone()], None, true).unwrap_err().to_string();
        assert!(err.contains("no manifest"), "{err}");
        assert!(err.contains("--allow-unverified"), "the error should say how to proceed: {err}");

        let out = dir.join("all.csv");
        run(&[a, b], Some(&out), false).unwrap();
        assert_eq!(std::fs::read_to_string(&out).unwrap().lines().count(), 4);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn parts_with_different_columns_are_refused() {
        let dir = scratch("columns");
        let a = part(&dir, "part-0000", &["a,1,0.10"], "CPU", "abc");
        let b = dir.join("part-0001.csv");
        std::fs::write(&b, "name,position,probability,confidence\nb,1,0.30,0.1\n").unwrap();
        std::fs::copy(crate::manifest::path_for(&a), crate::manifest::path_for(&b)).unwrap();
        let err = run(&[a, b], None, true).unwrap_err().to_string();
        assert!(err.contains("columns"), "{err}");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn merging_nothing_is_an_error_rather_than_an_empty_file() {
        assert!(run(&[], None, true).is_err());
    }
}
