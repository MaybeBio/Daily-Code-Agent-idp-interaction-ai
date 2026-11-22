//! Deterministic sequence and structure analyses.
//!
//! Everything here is a well-defined algorithm with a literature reference
//! and is tested against published values. Nothing in this crate estimates a
//! quantity it cannot actually compute; predictive models live in
//! `opendna-ml`.
//!
//! - [`properties`]: molecular weight, pI, GRAVY, extinction coefficient (ProtParam-compatible)
//! - [`backbone`]: φ/ψ/ω torsions
//! - [`dssp`]: secondary structure (Kabsch & Sander 1983)
//! - [`sasa`]: solvent-accessible surface area (Shrake & Rupley 1973)
//! - [`superpose`]: optimal superposition and RMSD (Kabsch 1976)
//! - [`interactions`]: disulfides, salt bridges, Cα contacts
//! - [`humanness`]: peptide prevalence in human antibody repertoires (`OASis`)
//! - [`shape`]: centroid and radius of gyration
#![cfg_attr(test, allow(clippy::unwrap_used))]

pub mod backbone;
pub mod dssp;
mod error;
pub mod humanness;
pub mod interactions;
pub mod properties;
pub mod sasa;
pub mod shape;
pub mod superpose;

pub use error::AnalysisError;
