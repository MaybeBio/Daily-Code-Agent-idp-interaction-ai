//! Shared build context and the common "write a finished dataset" step.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use serde::Serialize;

use crate::cache::Cache;
use crate::card::{self, CardInput};
use crate::cluster::{self, ClusterParams, Clustering};
use crate::parquet_out::{write_parquet, Table, WrittenFile};
use crate::split::{Fractions, Split, SplitSummary};

pub struct BuildContext {
    pub cache: Cache,
    pub out_root: PathBuf,
    pub mmseqs: Option<PathBuf>,
    pub seed: u64,
    pub threads: usize,
    /// Cap on domains sampled by the burial dataset, for smoke tests.
    pub max_domains: Option<usize>,
}

impl BuildContext {
    pub fn cluster(&self, seqs: &BTreeMap<String, String>, params: ClusterParams) -> Result<Clustering> {
        let bin = cluster::find_mmseqs(self.mmseqs.as_deref())?;
        cluster::cluster(&bin, seqs, params, self.threads)
    }

    pub fn dataset_dir(&self, name: &str) -> PathBuf {
        self.out_root.join(name)
    }
}

#[derive(Serialize)]
struct Manifest<'a> {
    dataset: &'a str,
    pipeline: &'static str,
    pipeline_version: &'static str,
    seed: u64,
    /// Absent for a dataset that is not clustered; see [`FinishedDataset`].
    clustering: Option<String>,
    /// Absent for a dataset with no train/test split.
    fractions: Option<Fractions>,
    sources: BTreeMap<&'a str, &'a str>,
    files: &'a [WrittenFile],
    splits: &'a SplitSummary,
    statistics: &'a [(String, String)],
}

pub struct FinishedDataset<'a> {
    pub name: &'static str,
    pub tables: &'a [(Split, Table)],
    /// Tables written beside the splits, under their own names.
    ///
    /// For rows that are part of the dataset but must not be mistaken for
    /// labelled data: the antibody dataset's out-of-distribution species and
    /// its transgenic-mouse diagnostic. Putting them in `test` with a null
    /// label would let anything that averages over the test split quietly
    /// average over rows that have no answer.
    pub extra: &'a [(&'static str, Table)],
    pub card: CardInput<'a>,
    /// How rows were divided, for a dataset that has a split.
    ///
    /// `None` where there is nothing to divide: the antibody cohorts cannot be
    /// split by group at all, because inbred mice share exact rearrangements
    /// across animals and so form a single indivisible group. Recording a
    /// fraction there would put a number in the manifest that describes
    /// nothing.
    pub fractions: Option<Fractions>,
    /// `None` where no clustering was run.
    pub cluster_params: Option<ClusterParams>,
}

/// Writes `data/<split>.parquet`, `README.md` and `manifest.json` for one dataset.
pub fn write_dataset(ctx: &BuildContext, d: FinishedDataset<'_>) -> Result<Vec<WrittenFile>> {
    let dir = ctx.dataset_dir(d.name);
    let data_dir = dir.join("data");
    if data_dir.exists() {
        // Remove stale split files so an empty split cannot leave an old file behind.
        fs::remove_dir_all(&data_dir)?;
    }
    fs::create_dir_all(&data_dir)?;
    let mut files = Vec::new();
    for (split, table) in d.tables {
        if table.n_rows() == 0 {
            continue;
        }
        let rel = format!("data/{split}.parquet");
        files.push(write_parquet(table, &dir.join(&rel), &rel)?);
    }
    for (name, table) in d.extra {
        if table.n_rows() == 0 {
            continue;
        }
        let rel = format!("data/{name}.parquet");
        files.push(write_parquet(table, &dir.join(&rel), &rel)?);
    }
    fs::write(dir.join("README.md"), card::render(&d.card))?;
    let sources = d
        .card
        .sources
        .iter()
        .map(|s| (s.id, d.card.lock.sources.get(s.id).map_or("", |e| e.sha256.as_str())))
        .collect();
    let manifest = Manifest {
        dataset: d.name,
        pipeline: "opendna-data",
        pipeline_version: env!("CARGO_PKG_VERSION"),
        seed: ctx.seed,
        clustering: d.cluster_params.as_ref().map(ClusterParams::describe),
        fractions: d.fractions,
        sources,
        files: &files,
        splits: d.card.splits,
        statistics: &d.card.statistics,
    };
    fs::write(dir.join("manifest.json"), serde_json::to_string_pretty(&manifest)? + "\n")?;
    if d.card.cohorts.is_empty() {
        print_summary(&dir, d.name, d.card.splits);
    } else {
        println!("{}: wrote {}", d.name, dir.display());
        for (name, rows) in d.card.cohorts {
            println!("  {name:<22} {rows:>9} rows");
        }
    }
    Ok(files)
}

fn print_summary(dir: &Path, name: &str, s: &SplitSummary) {
    println!("{name}: wrote {}", dir.display());
    for sp in Split::ALL {
        let k = s.get(sp);
        println!("  {sp:<10} {:>9} rows  {:>7} groups", k.rows, k.groups);
    }
    println!("  {:<10} {:>9} rows", "total", s.total_rows());
}
