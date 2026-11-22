//! Sequence clustering via MMseqs2.
//!
//! The pipeline does not implement its own clustering. It shells out to
//! MMseqs2 (`easy-cluster`) and fails with an actionable message when the
//! binary cannot be found, rather than falling back to something weaker that
//! would silently make the splits leak homologs.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, Context, Result};
use serde::Serialize;

use crate::fasta;

/// MMseqs2 clustering parameters. Recorded verbatim in manifests and cards.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct ClusterParams {
    /// `--min-seq-id`
    pub min_seq_id: f64,
    /// `-c`
    pub coverage: f64,
    /// `--cov-mode` (0 = bidirectional, 1 = target, 2 = query)
    pub cov_mode: u8,
    /// `--cluster-mode` (0 = set cover, 1 = connected component, 2 = greedy incremental)
    pub cluster_mode: u8,
    /// `-s` sensitivity
    pub sensitivity: f64,
}

impl ClusterParams {
    /// 30% identity over 80% bidirectional coverage.
    pub const ID30: ClusterParams =
        ClusterParams { min_seq_id: 0.3, coverage: 0.8, cov_mode: 0, cluster_mode: 0, sensitivity: 7.5 };

    // There was an ID90 here, for antibody repertoires, on the reasoning that
    // 30% identity collapses a whole species into one cluster while 90%
    // separates clonal lineages. The first half is right and the second is
    // not: measured on 176,543 real heavy chains, clustering at 90% and
    // joining donors through shared clusters put all 37 donors in a single
    // component, because the V(D)J region is mostly germline-constant and 90%
    // over its whole length groups V-gene families. The antibody dataset
    // groups by identical rearrangement instead and does not cluster at all.

    pub fn describe(&self) -> String {
        format!(
            "mmseqs easy-cluster --min-seq-id {} -c {} --cov-mode {} --cluster-mode {} -s {}",
            self.min_seq_id, self.coverage, self.cov_mode, self.cluster_mode, self.sensitivity
        )
    }
}

/// Locates the MMseqs2 executable.
///
/// Order: explicit path, `MMSEQS` environment variable, then `PATH`.
pub fn find_mmseqs(explicit: Option<&Path>) -> Result<PathBuf> {
    if let Some(p) = explicit {
        if p.is_file() {
            return Ok(p.to_owned());
        }
        bail!("--mmseqs {} does not exist", p.display());
    }
    if let Some(p) = std::env::var_os("MMSEQS") {
        let p = PathBuf::from(p);
        if p.is_file() {
            return Ok(p);
        }
        bail!("MMSEQS={} does not exist", p.display());
    }
    let names: &[&str] = if cfg!(windows) { &["mmseqs.exe", "mmseqs.bat", "mmseqs"] } else { &["mmseqs"] };
    if let Some(path) = std::env::var_os("PATH") {
        for dir in std::env::split_paths(&path) {
            for n in names {
                let cand = dir.join(n);
                if cand.is_file() {
                    return Ok(cand);
                }
            }
        }
    }
    bail!(
        "MMseqs2 was not found. This dataset needs real sequence clustering to build \
         leakage-free splits, and the pipeline does not substitute a weaker method.\n\
         Install MMseqs2 (https://github.com/soedinglab/MMseqs2, e.g. `conda install -c bioconda mmseqs2` \
         or a release binary) and put it on PATH, set MMSEQS, or pass --mmseqs <path>."
    )
}

/// Result of clustering: every input id mapped to its cluster representative.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Clustering {
    pub member_to_rep: BTreeMap<String, String>,
}

impl Clustering {
    pub fn cluster_of(&self, id: &str) -> Option<&str> {
        self.member_to_rep.get(id).map(String::as_str)
    }

    pub fn n_clusters(&self) -> usize {
        self.member_to_rep.values().collect::<BTreeSet<_>>().len()
    }
}

/// Parses MMseqs2 `*_cluster.tsv` output (`representative<TAB>member`).
pub fn parse_cluster_tsv(r: impl BufRead) -> Result<Clustering> {
    let mut c = Clustering::default();
    for (i, line) in r.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let mut parts = line.split('\t');
        let (Some(rep), Some(member), None) = (parts.next(), parts.next(), parts.next()) else {
            bail!("cluster tsv line {} is not `rep<TAB>member`: {line:?}", i + 1);
        };
        if let Some(prev) = c.member_to_rep.insert(member.to_owned(), rep.to_owned()) {
            if prev != rep {
                bail!("member {member} assigned to two clusters ({prev}, {rep})");
            }
        }
    }
    Ok(c)
}

/// Clusters `(id, sequence)` pairs. Ids must be unique and must not contain whitespace.
pub fn cluster(
    mmseqs: &Path,
    seqs: &BTreeMap<String, String>,
    params: ClusterParams,
    threads: usize,
) -> Result<Clustering> {
    if seqs.is_empty() {
        return Ok(Clustering::default());
    }
    for id in seqs.keys() {
        if id.is_empty() || id.chars().any(char::is_whitespace) {
            bail!("invalid clustering id {id:?}");
        }
    }
    let work = tempfile::Builder::new().prefix("opendna-mmseqs-").tempdir()?;
    let input = work.path().join("input.fasta");
    {
        let w = BufWriter::new(File::create(&input)?);
        fasta::write(w, seqs.iter().map(|(k, v)| (k.as_str(), v.as_str())))?;
    }
    let out_prefix = work.path().join("clu");
    let tmp = work.path().join("tmp");
    let mut cmd = Command::new(mmseqs);
    cmd.arg("easy-cluster")
        .arg(&input)
        .arg(&out_prefix)
        .arg(&tmp)
        .args(["--min-seq-id", &params.min_seq_id.to_string()])
        .args(["-c", &params.coverage.to_string()])
        .args(["--cov-mode", &params.cov_mode.to_string()])
        .args(["--cluster-mode", &params.cluster_mode.to_string()])
        .args(["-s", &params.sensitivity.to_string()])
        .args(["--threads", &threads.max(1).to_string()])
        .args(["-v", "1"]);
    eprintln!("clustering {} sequences: {}", seqs.len(), params.describe());
    let status = cmd.status().with_context(|| format!("running {}", mmseqs.display()))?;
    if !status.success() {
        bail!("mmseqs easy-cluster failed with {status}");
    }
    let tsv = work.path().join("clu_cluster.tsv");
    let clustering = parse_cluster_tsv(BufReader::new(
        File::open(&tsv).with_context(|| format!("reading {}", tsv.display()))?,
    ))?;
    let missing: Vec<_> =
        seqs.keys().filter(|k| !clustering.member_to_rep.contains_key(*k)).take(5).collect();
    if !missing.is_empty() {
        bail!("mmseqs output is missing input ids, e.g. {missing:?}");
    }
    Ok(clustering)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_cluster_tsv() {
        let tsv = "a\ta\na\tb\nc\tc\n";
        let c = parse_cluster_tsv(tsv.as_bytes()).unwrap();
        assert_eq!(c.cluster_of("b"), Some("a"));
        assert_eq!(c.cluster_of("c"), Some("c"));
        assert_eq!(c.n_clusters(), 2);
    }

    #[test]
    fn rejects_malformed_and_conflicting_lines() {
        assert!(parse_cluster_tsv("a b\n".as_bytes()).is_err());
        assert!(parse_cluster_tsv("a\tb\tc\n".as_bytes()).is_err());
        assert!(parse_cluster_tsv("a\tx\nb\tx\n".as_bytes()).is_err());
    }

    #[test]
    fn explicit_missing_binary_is_an_error() {
        let err = find_mmseqs(Some(Path::new("/definitely/not/here/mmseqs"))).unwrap_err();
        assert!(err.to_string().contains("does not exist"));
    }
}
