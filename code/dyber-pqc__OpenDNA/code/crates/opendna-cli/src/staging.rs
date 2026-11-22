//! Carrying a working installation to a machine with no network.
//!
//! The heads are small and a lab can copy them anywhere. The embedding model
//! is not: it is several hundred megabytes, it lives in the Hugging Face cache
//! rather than beside the heads, and the first prediction on a fresh machine
//! silently downloads it. On an instrument computer, a validated workstation or
//! anything behind an air gap that is not a slow first run — it is a failure,
//! and one that happens the first time somebody actually needs a result.
//!
//! So `opendna stage` writes one directory holding everything: the four heads,
//! the embedding model's own files, a lockfile and a checksum manifest. Copy
//! that directory to the offline machine, point `--models` at it, and pass
//! `--offline` to make any attempt to reach the network an error rather than a
//! stall.
//!
//! Staging deliberately copies what this machine already resolved rather than
//! re-fetching on the target. The bundle is then exactly what was tested here,
//! which is the only version of "the same software" that means anything.

use anyhow::{bail, Context, Result};
use opendna_ml::esm2::{Esm2, ModelSize};
use opendna_ml::heads::{sha256_file, HeadFiles};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

/// Subdirectory of a models directory holding the staged embedding model.
///
/// A subdirectory rather than loose files, so a staged bundle is still a
/// perfectly ordinary models directory: `config.json` beside the heads would
/// read as a head's metadata.
const EMBEDDER: &str = "embedder";

/// The checksum manifest, in `sha256sum` format so the receiving machine can
/// check it with tools it already has.
const SUMS: &str = "sha256sums.txt";

/// The licence notice, named the way a legal review expects to find it.
const NOTICE: &str = "LICENSES.md";

/// The embedding model's files inside a staged bundle, if they are both there.
///
/// Both files or neither: half a bundle that loads is worse than one that
/// does not, because it would fall through to the network on a machine chosen
/// for having none.
pub fn staged_embedder(dir: &Path) -> Option<(PathBuf, PathBuf)> {
    let config = dir.join(EMBEDDER).join("config.json");
    let weights = dir.join(EMBEDDER).join("model.safetensors");
    (config.is_file() && weights.is_file()).then_some((config, weights))
}

/// Load the embedding model, from the bundle if it is staged there.
///
/// # Errors
/// If the model cannot be loaded, or if `offline` is set and nothing is staged
/// — in which case the error says how to stage it rather than describing a
/// network failure the machine was never going to recover from.
pub fn open_embedder(
    dir: &Path,
    size: ModelSize,
    device: &candle_core::Device,
    offline: bool,
) -> Result<Esm2> {
    if let Some((config, weights)) = staged_embedder(dir) {
        return Esm2::from_files(&config, &weights, device)
            .with_context(|| format!("loading the staged embedding model from {}", dir.display()));
    }
    if offline {
        bail!(
            "--offline was given and no embedding model is staged in {}.\n\
             Expected {}/{EMBEDDER}/config.json and {}/{EMBEDDER}/model.safetensors.\n\
             On a machine with network access, run:\n  \
             opendna stage --to {}\n\
             then copy that directory here. OpenDNA will not reach the network with --offline \
             set, and it will not predict without the model.",
            dir.display(),
            dir.display(),
            dir.display(),
            dir.display()
        );
    }
    // The Hugging Face cache is not safe to read from two processes at once.
    // On Windows the second one is refused the weights file outright, with an
    // operating-system error that says nothing about why, so the reason is
    // attached here. Measured: three staged runs side by side all succeed,
    // while two runs through the cache leave one of them failing.
    Esm2::from_pretrained(size, device).map_err(anyhow::Error::from).context(
        "loading the embedding model through the Hugging Face cache. If another OpenDNA run is \
         going at the same time, that is the likely cause: the cache is not safe for two \
         processes to read at once. Run `opendna stage --to <dir>` once and point --models at \
         that directory; several runs can share a staged bundle.",
    )
}

/// Summarise the licences of the heads being staged.
///
/// Read out of the cards rather than written here, for the same reason the
/// limitations are: the card is what is published, so a notice generated from
/// it cannot claim terms the file beside the weights does not.
///
/// A head whose card declares nothing is listed as unstated. That is a true
/// and useful thing to tell somebody about to redistribute a file; defaulting
/// it to the permissive licence the other heads use would be a false one.
fn licences(dir: &Path, installed: &[crate::inventory::Installed]) -> String {
    let mut out = String::from(
        "# Licences in this bundle\n\n\
         Generated by `opendna stage` from each head's model card. The cards\n\
         themselves are in this directory and are the authoritative statement;\n\
         this file exists so the obligations can be seen without reading all\n\
         four of them.\n\n",
    );

    // Grouped by licence so a reader sees how many distinct sets of terms
    // apply, which is the question, rather than a list of four files.
    let mut by_licence: BTreeMap<String, Vec<&str>> = BTreeMap::new();
    let mut links: BTreeMap<String, Option<String>> = BTreeMap::new();
    let mut obliging: Vec<String> = Vec::new();
    for head in installed {
        let files = HeadFiles::in_dir(dir, head.task);
        let (label, link, sa) = match files.license() {
            Some(l) => (l.label().to_owned(), l.link.clone(), l.share_alike()),
            None => ("unstated".to_owned(), None, true),
        };
        by_licence.entry(label.clone()).or_default().push(head.task);
        links.entry(label.clone()).or_insert(link);
        if sa && !obliging.contains(&label) {
            obliging.push(label);
        }
    }

    for (label, tasks) in &by_licence {
        let _ = writeln!(out, "## {label}\n");
        let _ = writeln!(out, "Applies to: {}\n", tasks.join(", "));
        if let Some(Some(link)) = links.get(label) {
            let _ = writeln!(out, "Terms: {link}\n");
        }
        if label == "unstated" {
            let _ = writeln!(
                out,
                "No licence is declared in the card for these. Treat them as \
                 not redistributable until that is resolved.\n"
            );
        }
    }

    if obliging.is_empty() {
        let _ = writeln!(out, "None of these licences requires a redistributor to keep the same terms.");
    } else {
        let _ = writeln!(
            out,
            "**Redistributing this bundle carries obligations.** {} {} share-alike or \
             unstated: passing the bundle on means attributing the source and keeping \
             the same terms for anything derived from it. Each card names the \
             attribution its licence requires.",
            obliging.join(" and "),
            if obliging.len() == 1 { "is" } else { "are" }
        );
    }
    out
}

/// What the manifest should say about the embedding model.
///
/// Nothing here hashes the weights. At 568 MB that costs about as long as a
/// small run, so a staged bundle's recorded digest is used when there is one
/// and the absence is stated plainly when there is not.
pub fn describe_embedder(dir: &Path, size: ModelSize) -> crate::manifest::Embedder {
    let staged = staged_embedder(dir).map(|(_, weights)| weights);
    let bytes = staged.as_ref().and_then(|p| std::fs::metadata(p).ok()).map(|m| m.len());
    let digest = staged.as_ref().and_then(|_| recorded_digest(dir, "embedder/model.safetensors"));
    crate::manifest::Embedder {
        repo: size.repo(),
        source: if staged.is_some() { "staged bundle" } else { "Hugging Face cache" },
        weights: staged,
        bytes,
        digest_source: if digest.is_some() {
            "recorded in the staged bundle's sha256sums.txt"
        } else {
            "not recorded: hashing 568 MB costs about as long as a small run"
        },
        digest,
    }
}

/// The digest a bundle already wrote down for one of its files.
///
/// Read rather than recomputed: the manifest's claim is "the bundle says this
/// file is X", which is what `sha256sum -c` checks, and recomputing it here
/// would cost the time this avoids while proving nothing extra.
fn recorded_digest(dir: &Path, name: &str) -> Option<String> {
    let text = std::fs::read_to_string(dir.join(SUMS)).ok()?;
    text.lines().find_map(|line| {
        let (digest, file) = line.split_once("  ")?;
        (file.trim() == name).then(|| digest.trim().to_owned())
    })
}

/// Build a self-contained bundle in `to`.
///
/// # Errors
/// If a head is missing or corrupt, if the embedding model cannot be obtained,
/// or if anything cannot be copied. A partial bundle is not written: every
/// head is checked before the first byte is copied, because discovering a
/// problem after half a gigabyte has moved is no use to anyone.
pub fn stage(from: &Path, to: &Path, size: ModelSize) -> Result<()> {
    let found = crate::inventory::scan(from, size);
    let unusable: Vec<String> =
        found.iter().filter(|h| !h.ready()).map(|h| format!("  {:11} {}", h.task, h.describe())).collect();
    if !unusable.is_empty() {
        bail!(
            "refusing to stage an installation that does not verify. In {}:\n{}",
            from.display(),
            unusable.join("\n")
        );
    }

    std::fs::create_dir_all(to.join(EMBEDDER))
        .with_context(|| format!("creating {}", to.join(EMBEDDER).display()))?;

    // Recorded as it is copied, so the manifest describes the bytes that
    // landed in the bundle rather than the ones that were meant to.
    let mut sums = String::new();
    let mut copy = |source: &Path, name: &str| -> Result<u64> {
        let target = to.join(name);
        let bytes = std::fs::copy(source, &target)
            .with_context(|| format!("copying {} to {}", source.display(), target.display()))?;
        let digest = sha256_file(&target)?;
        // Two spaces and a relative path: `sha256sum -c sha256sums.txt`.
        let _ = writeln!(sums, "{digest}  {}", name.replace('\\', "/"));
        Ok(bytes)
    };

    let mut total = 0_u64;
    for head in &found {
        let files = HeadFiles::in_dir(from, head.task);
        total += copy(&files.weights, &format!("{}.safetensors", head.task))?;
        total += copy(&files.metadata, &format!("{}.metadata.json", head.task))?;
        // The card travels with the weights: the CLI reads its limitations out
        // of the card, so a bundle without it would state fewer caveats than
        // the model has.
        let card = from.join(format!("{}.README.md", head.task));
        if card.is_file() {
            total += copy(&card, &format!("{}.README.md", head.task))?;
        }
    }

    // On this machine this is a cache hit and copies locally; on a machine
    // that has never run a prediction it downloads first, which is the point
    // of doing it here rather than there.
    eprintln!("Resolving {} ...", size.repo());
    let checkpoint = opendna_ml::weights::fetch(size.repo())?;
    total += copy(&checkpoint.config, &format!("{EMBEDDER}/config.json"))?;
    total += copy(&checkpoint.weights, &format!("{EMBEDDER}/model.safetensors"))?;

    let notice = to.join(NOTICE);
    std::fs::write(&notice, licences(from, &found))
        .with_context(|| format!("writing {}", notice.display()))?;
    std::fs::write(to.join(SUMS), &sums).with_context(|| format!("writing {}", to.join(SUMS).display()))?;
    let lock = to.join("models.lock");
    std::fs::write(&lock, crate::inventory::lockfile(&found, size))
        .with_context(|| format!("writing {}", lock.display()))?;

    #[allow(clippy::cast_precision_loss)]
    let megabytes = total as f64 / 1_048_576.0;
    println!("Staged {} heads and {} into {}", found.len(), size.repo(), to.display());
    println!("{megabytes:.0} MB, {SUMS} and models.lock written");
    println!(
        "\nCopy that directory to the offline machine, then:\n  \
         opendna --models <dir> --offline selftest\n  \
         opendna --models <dir> --offline predict --task disorder -i seqs.fasta -o out.csv"
    );
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("opendna-stage-{name}-{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// Half a bundle must not read as a bundle: on the machine this is for,
    /// falling through to the network is a hang rather than a slow start.
    #[test]
    fn an_embedder_is_staged_only_when_both_files_are_present() {
        let dir = scratch("half");
        assert!(staged_embedder(&dir).is_none(), "an empty directory is not a bundle");

        std::fs::create_dir_all(dir.join(EMBEDDER)).unwrap();
        std::fs::write(dir.join(EMBEDDER).join("config.json"), "{}").unwrap();
        assert!(staged_embedder(&dir).is_none(), "config alone is not a bundle");

        std::fs::write(dir.join(EMBEDDER).join("model.safetensors"), "not really").unwrap();
        assert!(staged_embedder(&dir).is_some(), "both files present");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// The whole promise of `--offline` is that it does not touch the network.
    /// With nothing staged it has to say so, and say how to fix it.
    #[test]
    fn offline_without_a_staged_embedder_explains_how_to_stage_one() {
        let dir = scratch("offline");
        let device = candle_core::Device::Cpu;
        let Err(err) = open_embedder(&dir, ModelSize::Large150M, &device, true) else {
            panic!("--offline loaded a model from an empty directory");
        };
        let err = err.to_string();
        assert!(err.contains("--offline"), "{err}");
        assert!(err.contains("opendna stage --to"), "{err}");
        assert!(err.contains("will not reach the network"), "{err}");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// `describe_embedder` reads the digests back out of the manifest the
    /// staging wrote, so the two formats have to agree.
    #[test]
    fn a_recorded_digest_is_read_back_out_of_the_manifest() {
        let dir = scratch("sums");
        std::fs::write(dir.join(SUMS), "aa11  disorder.safetensors\nbb22  embedder/model.safetensors\n")
            .unwrap();
        assert_eq!(recorded_digest(&dir, "embedder/model.safetensors").as_deref(), Some("bb22"));
        assert_eq!(recorded_digest(&dir, "disorder.safetensors").as_deref(), Some("aa11"));
        assert_eq!(recorded_digest(&dir, "absent.safetensors"), None);
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A bundle whose notice omitted a share-alike licence would be worse
    /// than one with no notice at all: somebody would read it and conclude
    /// there was nothing to comply with.
    #[test]
    fn the_notice_names_every_licence_and_flags_the_obliging_ones() {
        let dir = scratch("licences");
        for (task, front) in [
            ("stability", "---\nlicense: apache-2.0\n---\n"),
            ("solubility", "---\nlicense: other\nlicense_name: cc-by-sa-2.1-jp\n---\n"),
        ] {
            std::fs::write(dir.join(format!("{task}.README.md")), front).unwrap();
        }
        let installed: Vec<crate::inventory::Installed> = ["stability", "solubility"]
            .iter()
            .map(|task| crate::inventory::Installed {
                task,
                metadata: None,
                problem: None,
                integrity: opendna_ml::heads::Integrity::Unrecorded,
            })
            .collect();

        let text = licences(&dir, &installed);
        assert!(text.contains("apache-2.0"), "{text}");
        assert!(text.contains("cc-by-sa-2.1-jp"), "{text}");
        assert!(text.contains("carries obligations"), "{text}");
        assert!(text.contains("stability"), "{text}");
        assert!(text.contains("solubility"), "{text}");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A head with no card must not be quietly dropped from the notice.
    #[test]
    fn an_undeclared_licence_is_reported_as_unstated() {
        let dir = scratch("nolicence");
        let installed = vec![crate::inventory::Installed {
            task: "disorder",
            metadata: None,
            problem: None,
            integrity: opendna_ml::heads::Integrity::Unrecorded,
        }];
        let text = licences(&dir, &installed);
        assert!(text.contains("unstated"), "{text}");
        assert!(text.contains("not redistributable"), "{text}");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// Staging an installation that does not verify would produce a bundle
    /// nobody can trust, carried to a machine where nobody can re-download.
    #[test]
    fn staging_refuses_an_installation_that_does_not_verify() {
        let (from, to) = (scratch("empty-from"), scratch("empty-to"));
        let err = stage(&from, &to, ModelSize::Large150M).unwrap_err().to_string();
        assert!(err.contains("refusing to stage"), "{err}");
        assert!(err.contains("not installed"), "{err}");
        assert!(!to.join(EMBEDDER).exists(), "nothing should have been written");
        std::fs::remove_dir_all(&from).ok();
        std::fs::remove_dir_all(&to).ok();
    }
}
