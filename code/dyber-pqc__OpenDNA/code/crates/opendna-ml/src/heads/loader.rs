//! Loading trained heads from disk.
//!
//! A head is a pair of files written by the training workspace:
//! `<task>.safetensors` and `<task>.metadata.json`. The metadata records
//! which embedding model the head was trained against, and loading checks it,
//! so a head can never be paired with embeddings it was not trained on.

use super::{BurialHead, DisorderHead, HeadMetadata, SolubilityHead, StabilityHead};
use crate::esm2::ModelSize;
use crate::MlError;
use candle_core::{DType, Device};
use candle_nn::VarBuilder;
use std::path::{Path, PathBuf};

/// SHA-256 of a file, streamed.
///
/// # Errors
/// [`MlError::File`] if the file cannot be read.
pub fn sha256_file(path: &Path) -> Result<String, MlError> {
    use sha2::{Digest, Sha256};
    let file = std::fs::File::open(path)
        .map_err(|source| MlError::File { path: path.display().to_string(), source })?;
    let mut reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut hasher = Sha256::new();
    // Heap-allocated: a 64 KiB stack buffer is a lot to ask of a thread with
    // a small stack, and the lint is right.
    let mut buf = vec![0u8; 1 << 16].into_boxed_slice();
    loop {
        let n = std::io::Read::read(&mut reader, &mut buf)
            .map_err(|source| MlError::File { path: path.display().to_string(), source })?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

/// The licence a model card declares.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct License {
    /// The `license` field. `other` when the licence is named separately.
    pub identifier: String,
    /// The `license_name` field, present when `identifier` is `other`.
    pub name: Option<String>,
    /// The `license_link` field.
    pub link: Option<String>,
}

impl License {
    /// How to name this licence to a reader.
    pub fn label(&self) -> &str {
        match self.name.as_deref() {
            Some(name) if self.identifier == "other" => name,
            _ => &self.identifier,
        }
    }

    /// Whether the terms oblige a redistributor to keep them.
    ///
    /// Matched on the identifier rather than parsed, because the only thing
    /// this is used for is deciding whether to say "this one is share-alike"
    /// in a notice -- and a wrong "no" there is the expensive direction.
    /// Anything unrecognised is treated as possibly share-alike.
    pub fn share_alike(&self) -> bool {
        let text = format!("{} {}", self.identifier, self.name.as_deref().unwrap_or_default());
        let text = text.to_ascii_lowercase();
        text.contains("-sa") || text.contains("share") || text.contains("gpl")
    }
}

/// What checking a head's files against its metadata found.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Integrity {
    /// The weights hash the metadata records.
    Verified,
    /// The metadata records no hash, so nothing could be checked. Metadata
    /// written before hashes existed reports this rather than failing.
    Unrecorded,
    /// The weights are not the file the metadata describes.
    Mismatch {
        /// What the metadata says the weights should hash to.
        expected: String,
        /// What they actually hash to.
        actual: String,
    },
}

impl Integrity {
    /// Whether the head may be used. `Unrecorded` is usable: it is an older
    /// head, not a corrupt one.
    pub fn usable(&self) -> bool {
        !matches!(self, Self::Mismatch { .. })
    }
}

/// Where a head's files live.
#[derive(Debug, Clone)]
pub struct HeadFiles {
    /// Weights.
    pub weights: PathBuf,
    /// Description of what the head expects.
    pub metadata: PathBuf,
    /// Model card. Not required to load the head, but it is where the
    /// limitations live.
    pub card: PathBuf,
}

impl HeadFiles {
    /// The conventional file names for `task` inside `dir`.
    pub fn in_dir(dir: &Path, task: &str) -> Self {
        Self {
            weights: dir.join(format!("{task}.safetensors")),
            metadata: dir.join(format!("{task}.metadata.json")),
            card: dir.join(format!("{task}.README.md")),
        }
    }

    /// The head's limitations, as its model card states them.
    ///
    /// Read from the card rather than written into the application, so that
    /// the desktop panel, the CLI's output files and the published card
    /// cannot state three different sets of caveats. Before this existed the
    /// UI carried its own hardcoded copy, which is precisely the drift this
    /// avoids.
    ///
    /// An absent or unparseable card yields none. The prediction is still
    /// valid without them, and a caller that has nothing to show should point
    /// at the card rather than fail.
    pub fn limitations(&self) -> Vec<String> {
        let Ok(card) = std::fs::read_to_string(&self.card) else {
            return Vec::new();
        };
        let mut out: Vec<String> = Vec::new();
        let mut inside = false;
        for line in card.lines() {
            if line.starts_with("## ") {
                inside = line.trim() == "## Limitations";
                continue;
            }
            if !inside {
                continue;
            }
            if let Some(bullet) = line.strip_prefix("- ") {
                out.push(bullet.trim().to_string());
            } else if !line.trim().is_empty() {
                // A continuation of the bullet above. Dropping it would
                // truncate a limitation mid-sentence, which is worse than
                // showing none: the reader gets half a caveat and no sign
                // that the rest was cut.
                if let Some(last) = out.last_mut() {
                    last.push(' ');
                    last.push_str(line.trim());
                }
            }
        }
        out.into_iter().map(|l| l.split_whitespace().collect::<Vec<_>>().join(" ")).collect()
    }

    /// The licence the card declares, if it declares one.
    ///
    /// Hugging Face model cards carry it in YAML front matter as `license`,
    /// optionally refined by `license_name` and `license_link` when the
    /// licence is not one of the recognised identifiers -- which is the case
    /// for the share-alike licence the solubility data carries.
    ///
    /// Returns `None` rather than a guess when the card is absent or silent.
    /// "Unstated" is a true and useful thing to tell somebody about to
    /// redistribute a file; a default of "Apache-2.0" would be a false one.
    pub fn license(&self) -> Option<License> {
        let card = std::fs::read_to_string(&self.card).ok()?;
        let mut front = card.lines();
        // Front matter only: a `license:` line in the prose below is part of a
        // sentence, not a declaration.
        if front.next()?.trim() != "---" {
            return None;
        }
        let mut found = License::default();
        for line in front {
            if line.trim() == "---" {
                break;
            }
            let Some((key, value)) = line.split_once(':') else { continue };
            let value = value.trim().trim_matches('"').to_string();
            match key.trim() {
                "license" => found.identifier = value,
                "license_name" => found.name = Some(value),
                "license_link" => found.link = Some(value),
                _ => {}
            }
        }
        (!found.identifier.is_empty()).then_some(found)
    }

    /// Check the weights against the hash in `metadata`.
    ///
    /// # Errors
    /// [`MlError::File`] if the weights cannot be read.
    pub fn verify(&self, metadata: &HeadMetadata) -> Result<Integrity, MlError> {
        let Some(expected) = metadata.weights_sha256.as_deref() else {
            return Ok(Integrity::Unrecorded);
        };
        let actual = sha256_file(&self.weights)?;
        Ok(if actual == expected {
            Integrity::Verified
        } else {
            Integrity::Mismatch { expected: expected.to_owned(), actual }
        })
    }

    /// Whether both files are present.
    pub fn exist(&self) -> bool {
        self.weights.is_file() && self.metadata.is_file()
    }
}

/// Read and check a head's metadata.
///
/// # Errors
/// [`MlError::File`] if the metadata cannot be read, or [`MlError::Config`]
/// if it does not describe `task` or does not match `embedder`.
pub fn read_metadata(files: &HeadFiles, task: &str, embedder: ModelSize) -> Result<HeadMetadata, MlError> {
    let text = std::fs::read_to_string(&files.metadata)
        .map_err(|source| MlError::File { path: files.metadata.display().to_string(), source })?;
    let metadata: HeadMetadata = serde_json::from_str(&text)
        .map_err(|e| MlError::Config(format!("{}: {e}", files.metadata.display())))?;

    if metadata.task != task {
        return Err(MlError::Config(format!(
            "{} describes the {:?} head, not {task:?}",
            files.metadata.display(),
            metadata.task
        )));
    }
    if metadata.embedding_model != embedder.repo() {
        return Err(MlError::Config(format!(
            "the {task} head was trained on {} but {} is loaded; predictions would be meaningless",
            metadata.embedding_model,
            embedder.repo()
        )));
    }
    if metadata.embedding_dim != embedder.hidden_size() {
        return Err(MlError::Config(format!(
            "the {task} head expects {}-wide embeddings but {} produces {}",
            metadata.embedding_dim,
            embedder.repo(),
            embedder.hidden_size()
        )));
    }
    Ok(metadata)
}

fn var_builder<'a>(path: &Path, device: &Device) -> Result<VarBuilder<'a>, MlError> {
    let tensors = candle_core::safetensors::load(path, device)?;
    Ok(VarBuilder::from_tensors(tensors, DType::F32, device))
}

/// Every head trained against one embedding model, as far as they are present.
///
/// A missing head is `None` rather than an error: the application reports the
/// prediction as unavailable instead of inventing one.
pub struct Heads {
    /// ΔΔG head.
    pub stability: Option<StabilityHead>,
    /// Solubility head.
    pub solubility: Option<SolubilityHead>,
    /// Per-residue disorder head.
    pub disorder: Option<DisorderHead>,
    /// Per-residue solvent accessibility head.
    pub burial: Option<BurialHead>,
}

impl Heads {
    /// Load whichever heads are present in `dir`.
    ///
    /// # Errors
    /// [`MlError`] if a head is present but unreadable or mismatched. A head
    /// that is simply absent is not an error.
    pub fn load(dir: &Path, embedder: ModelSize, device: &Device) -> Result<Self, MlError> {
        let dim = embedder.hidden_size();
        let load_one = |task: &str| -> Result<Option<(VarBuilder<'_>, HeadMetadata)>, MlError> {
            let files = HeadFiles::in_dir(dir, task);
            if !files.exist() {
                return Ok(None);
            }
            let metadata = read_metadata(&files, task, embedder)?;
            Ok(Some((var_builder(&files.weights, device)?, metadata)))
        };

        Ok(Self {
            // The metadata says whether the stability head takes burial, so
            // the right input width is used without guessing from the file.
            stability: load_one("stability")?
                .map(|(vb, m)| StabilityHead::new_with(dim, m.with_burial, vb))
                .transpose()?,
            solubility: load_one("solubility")?.map(|(vb, _)| SolubilityHead::new(dim, vb)).transpose()?,
            // The calibration travels with the weights: a head loaded without
            // its recorded map would quietly emit the uncalibrated numbers.
            disorder: load_one("disorder")?
                .map(|(vb, m)| DisorderHead::new(dim, vb).map(|h| h.with_calibration(m.calibration)))
                .transpose()?,
            burial: load_one("burial")?.map(|(vb, _)| BurialHead::new(dim, vb)).transpose()?,
        })
    }

    /// Names of the heads that loaded, for display.
    pub fn available(&self) -> Vec<&'static str> {
        let mut out = Vec::new();
        if self.stability.is_some() {
            out.push("stability");
        }
        if self.solubility.is_some() {
            out.push("solubility");
        }
        if self.disorder.is_some() {
            out.push("disorder");
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CARD: &str = "\
# Head

## Results

| a | b |

## Limitations

- First thing.
- Second thing,
  wrapped over two lines.

## Intended use

- Not a limitation.
";

    /// A head whose card holds `body`.
    ///
    /// `case` keys the directory, because `cargo test` runs these in parallel
    /// and a directory shared between them means one test overwrites the card
    /// another is reading — which shows up as a failure that disappears when
    /// the test is run on its own.
    fn card_with(case: &str, body: &str) -> HeadFiles {
        let dir = std::env::temp_dir().join(format!("opendna-card-{}-{case}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temp dir");
        std::fs::write(dir.join("disorder.README.md"), body).expect("write card");
        HeadFiles::in_dir(&dir, "disorder")
    }

    #[test]
    fn limitations_come_from_the_limitations_section_only() {
        let got = card_with("section", CARD).limitations();
        assert_eq!(got.len(), 2, "{got:?}");
        assert_eq!(got[0], "First thing.");
        assert!(!got.iter().any(|l| l.contains("Not a limitation")), "{got:?}");
    }

    /// A bullet wrapped across lines used to lose its tail, leaving the
    /// reader half a caveat with no sign the rest was cut.
    #[test]
    fn a_wrapped_bullet_keeps_its_continuation() {
        let got = card_with("wrapped", CARD).limitations();
        assert_eq!(got[1], "Second thing, wrapped over two lines.", "{got:?}");
    }

    #[test]
    fn a_card_without_limitations_yields_none() {
        assert_eq!(card_with("empty", "# Head\n\n## Results\n").limitations(), [] as [String; 0]);
    }

    /// The share-alike licence on the solubility data is the reason this
    /// exists, so it has to be read correctly from a card shaped like that
    /// one: `license: other` with the real name beside it.
    #[test]
    fn a_licence_named_separately_is_read_from_the_front_matter() {
        let files = card_with(
            "sa",
            "---
license: other
license_name: cc-by-sa-2.1-jp
             license_link: https://creativecommons.org/licenses/by-sa/2.1/jp/deed.en
---

# Card
",
        );
        let found = files.license().expect("a licence");
        assert_eq!(found.identifier, "other");
        assert_eq!(found.label(), "cc-by-sa-2.1-jp");
        assert!(found.share_alike(), "a -sa licence must be recognised as share-alike");
        assert!(found.link.unwrap().contains("by-sa/2.1/jp"));
    }

    #[test]
    fn a_plain_identifier_needs_no_name() {
        let files = card_with(
            "apache",
            "---
license: apache-2.0
---

# Card
",
        );
        let found = files.license().expect("a licence");
        assert_eq!(found.label(), "apache-2.0");
        assert!(!found.share_alike());
    }

    /// Telling somebody the terms are unstated is true and useful. Defaulting
    /// to a permissive licence would be neither.
    #[test]
    fn a_card_without_front_matter_declares_nothing() {
        assert_eq!(
            card_with(
                "none",
                "# Card

license: apache-2.0 appears in prose
"
            )
            .license(),
            None
        );
        assert_eq!(
            card_with(
                "empty",
                "---
tags:
  - protein
---
"
            )
            .license(),
            None
        );
    }

    /// Anything unrecognised is treated as possibly share-alike, because a
    /// wrong "no" is the expensive direction when the question is whether a
    /// redistribution carries obligations.
    #[test]
    fn an_unfamiliar_licence_is_assumed_to_carry_obligations() {
        let files = card_with(
            "odd",
            "---
license: other
license_name: some-copyleft-thing
---
",
        );
        assert!(!files.license().expect("a licence").share_alike(), "no marker, so not flagged");
        let gpl = card_with(
            "gpl",
            "---
license: gpl-3.0
---
",
        );
        assert!(gpl.license().expect("a licence").share_alike());
    }

    /// A head can be used without its card; the caller points at the card
    /// instead of failing.
    #[test]
    fn a_missing_card_yields_none() {
        let files = HeadFiles::in_dir(Path::new("definitely/not/here"), "disorder");
        assert_eq!(files.limitations(), [] as [String; 0]);
    }
}
