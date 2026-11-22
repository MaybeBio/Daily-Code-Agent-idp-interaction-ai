//! Local persistence.
//!
//! - [`Database`]: a `SQLite` file holding named proteins.
//! - [`History`]: git-like version history of a protein's sequence with
//!   branches and typed point mutations.
//!
//! Encryption at rest and project workspaces arrive with the platform phase.
#![cfg_attr(test, allow(clippy::unwrap_used))]

mod database;
mod history;

pub use database::{Database, StoreError, StoredProtein, StructureFile, StructureFormat};
pub use history::{History, HistoryError, Mutation, Version, VersionId};
