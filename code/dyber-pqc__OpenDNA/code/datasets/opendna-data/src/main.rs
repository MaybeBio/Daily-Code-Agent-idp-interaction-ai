//! `opendna-data`: build the OpenDNA training datasets.
//!
//! Run from the `datasets/` directory; paths default to `cache/`, `out/` and
//! `sources.lock` relative to the working directory.

use std::path::PathBuf;
use std::process::ExitCode;

use anyhow::Result;
use clap::{Parser, Subcommand};

use opendna_data::cache::{hash_file, Cache, LockCheck};
use opendna_data::context::BuildContext;
use opendna_data::datasets;
use opendna_data::sources::{ALL_SOURCES, EXCLUDED_SOURCES};

/// Default split seed. Changing it changes every split; record it if you do.
const DEFAULT_SEED: u64 = 20_261_114;

#[derive(Parser)]
#[command(version, about = "Download, clean, cluster, split and write the OpenDNA datasets")]
struct Cli {
    /// Directory for downloaded upstream files.
    #[arg(long, default_value = "cache", global = true)]
    cache_dir: PathBuf,
    /// Directory for built datasets (one subdirectory per dataset).
    #[arg(long, default_value = "out", global = true)]
    out_dir: PathBuf,
    /// Lockfile with the SHA-256 of every upstream file.
    #[arg(long, default_value = "sources.lock", global = true)]
    lockfile: PathBuf,
    /// Accept upstream files whose checksum differs from the lockfile and rewrite it.
    #[arg(long, global = true)]
    update_lock: bool,
    /// Never touch the network; fail if a file is not cached.
    #[arg(long, global = true)]
    offline: bool,
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// List upstream sources, their licenses, and excluded sources.
    Sources,
    /// Download (or verify cached copies of) the sources of one dataset or `all`.
    Download { dataset: String },
    /// Re-hash every cached file and compare it with the lockfile.
    Verify,
    /// Build one dataset or `all`.
    Build {
        dataset: String,
        /// Path to the MMseqs2 executable (otherwise MMSEQS or PATH).
        #[arg(long)]
        mmseqs: Option<PathBuf>,
        /// Seed for split assignment.
        #[arg(long, default_value_t = DEFAULT_SEED)]
        seed: u64,
        /// Threads passed to MMseqs2.
        #[arg(long, default_value_t = default_threads())]
        threads: usize,
        /// Cap the domains the burial dataset samples, for a quick check
        /// against a handful of real structures before fetching thousands.
        /// Ignored by the other datasets.
        #[arg(long)]
        max_domains: Option<usize>,
    },
    /// Measure a built dataset against what it claims, and write the report.
    ///
    /// Only `humanness` has one. It is not a validation score -- the measure
    /// is a count, so there is nothing to validate -- but whether the table
    /// separates the antibody cohorts it should is a checkable claim, and an
    /// uncheckable one would not be worth making.
    Evaluate { dataset: String },
}

fn default_threads() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get())
}

fn selected(dataset: &str) -> Vec<&str> {
    if dataset == "all" {
        datasets::ALL.to_vec()
    } else {
        vec![dataset]
    }
}

fn run(cli: Cli) -> Result<()> {
    match cli.cmd {
        Cmd::Sources => {
            println!("Included sources:");
            for s in ALL_SOURCES {
                println!("  {:<30} {:<16} {}", s.id, s.license.name(), s.description);
            }
            println!("\nExcluded sources:");
            for e in EXCLUDED_SOURCES {
                println!("  {} ({}): {}", e.name, e.dataset, e.reason);
            }
        }
        Cmd::Download { dataset } => {
            let mut cache = Cache::open(&cli.cache_dir, &cli.lockfile, cli.update_lock, cli.offline)?;
            for name in selected(&dataset) {
                for src in datasets::sources_for(name)? {
                    let p = cache.fetch(src)?;
                    println!("{:<30} {}", src.id, p.display());
                }
            }
        }
        Cmd::Evaluate { dataset } => match dataset.as_str() {
            "humanness" => opendna_data::evaluate::humanness(&cli.out_dir)?,
            other => anyhow::bail!("no evaluation for {other:?}; only `humanness` has one"),
        },
        Cmd::Verify => {
            let cache = Cache::open(&cli.cache_dir, &cli.lockfile, false, true)?;
            let mut bad = 0;
            for src in ALL_SOURCES {
                let p = cache.path_for(src);
                if !p.exists() {
                    println!("{:<30} not cached", src.id);
                    continue;
                }
                let (sha, _) = hash_file(&p)?;
                let state = opendna_data::cache::check_against_lock(cache.lock(), src, &sha);
                let label = match state {
                    LockCheck::Matches => "ok",
                    LockCheck::New => "not in lockfile",
                    LockCheck::Mismatch { .. } => {
                        bad += 1;
                        "MISMATCH"
                    }
                };
                println!("{:<30} {label}", src.id);
            }
            if bad > 0 {
                anyhow::bail!("{bad} cached file(s) do not match the lockfile");
            }
        }
        Cmd::Build { dataset, mmseqs, seed, threads, max_domains } => {
            let cache = Cache::open(&cli.cache_dir, &cli.lockfile, cli.update_lock, cli.offline)?;
            let mut ctx = BuildContext { cache, out_root: cli.out_dir, mmseqs, seed, threads, max_domains };
            for name in selected(&dataset) {
                datasets::build(name, &mut ctx)?;
            }
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}
