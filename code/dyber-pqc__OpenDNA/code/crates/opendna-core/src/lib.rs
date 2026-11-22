//! Core domain types shared by every `OpenDNA` crate.
//!
//! This crate is deliberately small and dependency-light. It defines *what*
//! a protein is (sequence, structure, residues, atoms) and leaves parsing,
//! analysis and storage to other crates.
#![cfg_attr(test, allow(clippy::unwrap_used))]

mod amino_acid;
mod error;
mod geometry;
mod protein;
mod sequence;
mod structure;

pub use amino_acid::AminoAcid;
pub use error::CoreError;
pub use geometry::Vec3;
pub use protein::Protein;
pub use sequence::Sequence;
pub use structure::{Atom, Chain, Residue, ResidueId, Structure};
