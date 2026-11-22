//! Fetching model weights from the Hugging Face Hub.
//!
//! Files are cached by the `hf-hub` client (under `HF_HOME`, or the user's
//! cache directory), so a repository is downloaded once per machine.

use crate::MlError;
use hf_hub::HFClientSync;
use std::path::PathBuf;

/// Local paths to a downloaded checkpoint.
#[derive(Debug, Clone)]
pub struct Checkpoint {
    /// `config.json`
    pub config: PathBuf,
    /// `model.safetensors`
    pub weights: PathBuf,
}

/// Download `repo`'s `config.json` and `model.safetensors`, or return the
/// cached copies. `repo` is `owner/name`, e.g. `facebook/esm2_t6_8M_UR50D`.
///
/// # Errors
/// [`MlError::Download`] if `repo` is malformed, the Hub cannot be reached,
/// or a file is absent.
pub fn fetch(repo: &str) -> Result<Checkpoint, MlError> {
    let fail = |message: String| MlError::Download { repo: repo.to_string(), message };
    let (owner, name) = repo
        .split_once('/')
        .ok_or_else(|| fail("expected an \"owner/name\" repository identifier".to_string()))?;
    let client = HFClientSync::new().map_err(|e| fail(e.to_string()))?;
    let model = client.model(owner, name);
    let get =
        |file: &str| model.download_file().filename(file).send().map_err(|e| fail(format!("{file}: {e}")));
    Ok(Checkpoint { config: get("config.json")?, weights: get("model.safetensors")? })
}

/// The files [`fetch`] needs from a checkpoint repository.
const REQUIRED: [&str; 2] = ["config.json", "model.safetensors"];

/// Whether `repo` is already in the local cache, so loading it will not need
/// the network.
///
/// This exists so the application can tell the user that the first prediction
/// will download several hundred megabytes, *before* they trigger it. A
/// failure to scan the cache is reported as "not cached": the only cost of
/// being wrong that way is an unnecessary warning, where the opposite would
/// be a silent, long stall.
///
/// Both files are checked rather than just the repository folder, because an
/// interrupted download leaves the folder behind.
pub fn is_cached(repo: &str) -> bool {
    let Ok(client) = HFClientSync::new() else {
        return false;
    };
    let Ok(cache) = client.scan_cache().send() else {
        return false;
    };
    cache.repos.iter().filter(|r| r.repo_id == repo && r.repo_type == "model").any(|r| {
        r.revisions.iter().any(|rev| {
            REQUIRED
                .iter()
                .all(|needed| rev.files.iter().any(|f| f.file_name == *needed && f.blob_path.is_file()))
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_repository_that_was_never_fetched_is_not_cached() {
        // No network: a miss is decided by scanning the local cache only.
        assert!(!is_cached("opendna-test/definitely-not-a-real-checkpoint"));
    }

    #[test]
    fn a_malformed_identifier_is_not_cached() {
        assert!(!is_cached("no-slash-here"));
    }
}
