//! File-format readers and writers.
//!
//! All readers take `&str` and are pure functions: no filesystem access, so
//! they work the same in the desktop app, the CLI and WebAssembly.
//!
//! | Format | Read | Write |
//! |---|---|---|
//! | PDB    | [`pdb::read`]   | [`pdb::write`]   |
//! | mmCIF  | [`mmcif::read`] | —                |
//! | FASTA  | [`fasta::read`] | [`fasta::write`] |
#![cfg_attr(test, allow(clippy::unwrap_used))]

mod builder;
mod error;
pub mod fasta;
pub mod mmcif;
pub mod pdb;

pub use error::{IoError, Location};

use opendna_core::Structure;

/// Detect PDB vs mmCIF from content and parse accordingly.
///
/// # Errors
/// Whatever the selected reader returns.
pub fn read_structure(content: &str) -> Result<Structure, IoError> {
    let first = content.lines().map(str::trim_start).find(|l| !l.is_empty()).unwrap_or("");
    if first.starts_with("data_") {
        mmcif::read(content)
    } else {
        pdb::read(content)
    }
}
