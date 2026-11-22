//! Developer tasks for the desktop app.
//!
//! `cargo run -p xtask -- fetch-molstar` downloads the pinned Mol* release
//! from the npm registry, verifies it against the registry's published
//! SHA-512 integrity hash, and extracts the prebuilt viewer bundle into
//! `ui/vendor/molstar/`. The bundle is not committed (≈5 MB, MIT-licensed).

#![allow(clippy::print_stdout)]

use base64::Engine;
use sha2::{Digest, Sha512};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

const MOLSTAR_VERSION: &str = "5.8.0";
/// `dist.integrity` from <https://registry.npmjs.org/molstar/5.8.0>.
const MOLSTAR_INTEGRITY: &str =
    "sha512-1SLIB5SxHEpeSieJmK7gV7MHYZkKgc2gvd3KL6IpIexvoG2u+3RPoqYe3CwJxuE7aDSnBx+frThzGjUA7FcFTA==";
/// Files taken from the tarball, relative to `package/`.
const MOLSTAR_FILES: &[&str] = &["build/viewer/molstar.js", "build/viewer/molstar.css", "LICENSE"];

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn main() -> ExitCode {
    if std::env::args().nth(1).as_deref() != Some("fetch-molstar") {
        eprintln!("usage: cargo run -p xtask -- fetch-molstar");
        return ExitCode::from(2);
    }
    match fetch_molstar() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn app_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().map_or_else(|| PathBuf::from("."), Path::to_path_buf)
}

fn fetch_molstar() -> Result<()> {
    let out = app_root().join("ui/vendor/molstar");
    let stamp = out.join(".version");
    if std::fs::read_to_string(&stamp).is_ok_and(|v| v.trim() == MOLSTAR_VERSION) {
        println!("Mol* {MOLSTAR_VERSION} already present in {}", out.display());
        return Ok(());
    }

    let url = format!("https://registry.npmjs.org/molstar/-/molstar-{MOLSTAR_VERSION}.tgz");
    println!("downloading {url}");
    let tarball = ureq::get(&url).call()?.body_mut().with_config().limit(200 * 1024 * 1024).read_to_vec()?;

    let expected = MOLSTAR_INTEGRITY.strip_prefix("sha512-").ok_or("integrity must be sha512")?;
    let actual = base64::engine::general_purpose::STANDARD.encode(Sha512::digest(&tarball));
    if actual != expected {
        return Err(format!("integrity mismatch: expected {expected}, got {actual}").into());
    }

    std::fs::create_dir_all(&out)?;
    let mut archive = tar::Archive::new(flate2::read::GzDecoder::new(tarball.as_slice()));
    let mut written = 0;
    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()?.into_owned();
        let Some(rel) =
            path.strip_prefix("package").ok().and_then(Path::to_str).map(|s| s.replace('\\', "/"))
        else {
            continue;
        };
        if let Some(wanted) = MOLSTAR_FILES.iter().find(|f| **f == rel) {
            let name = Path::new(wanted).file_name().ok_or("bad file name")?;
            let mut bytes = Vec::new();
            entry.read_to_end(&mut bytes)?;
            std::fs::write(out.join(name), bytes)?;
            written += 1;
        }
    }
    if written != MOLSTAR_FILES.len() {
        return Err(format!("expected {} files in the tarball, found {written}", MOLSTAR_FILES.len()).into());
    }
    std::fs::write(stamp, MOLSTAR_VERSION)?;
    println!("Mol* {MOLSTAR_VERSION} verified and extracted to {}", out.display());
    Ok(())
}
