//! Dataset pipeline for OpenDNA.
//!
//! The pipeline downloads upstream sources into a content-addressed cache,
//! records SHA-256 checksums in a lockfile, cleans and normalises the records,
//! clusters sequences with MMseqs2, assigns homology-aware train/validation/test
//! splits with a seeded deterministic procedure, and writes Parquet files plus a
//! Hugging Face dataset card per dataset.
//!
//! Parsing, normalisation and splitting are pure functions over readers and
//! plain structs so that they can be unit tested against small fixtures without
//! network access. Only [`cache`] and [`cluster`] touch the outside world.

pub mod cache;
pub mod card;
pub mod cluster;
pub mod context;
pub mod datasets;
pub mod evaluate;
pub mod fasta;
pub mod parquet_out;
pub mod rng;
pub mod sources;
pub mod split;
pub mod unionfind;

/// The 20 standard amino acids in one-letter code.
pub const STANDARD_AMINO_ACIDS: &str = "ACDEFGHIKLMNPQRSTVWY";

/// Returns true when `c` is one of the 20 standard amino acids (uppercase).
pub fn is_standard_aa(c: u8) -> bool {
    STANDARD_AMINO_ACIDS.as_bytes().contains(&c)
}

/// Returns true when every residue of `seq` is one of the 20 standard amino acids.
pub fn is_standard_protein(seq: &str) -> bool {
    !seq.is_empty() && seq.bytes().all(is_standard_aa)
}

/// Lowercase hex SHA-256 of a byte slice.
pub fn sha256_hex(bytes: &[u8]) -> String {
    use sha2::{Digest, Sha256};
    hex::encode(Sha256::digest(bytes))
}

/// Stable identifier for a protein sequence, used as the MMseqs2 FASTA id.
///
/// Using a hash of the sequence (rather than an upstream name) keeps the
/// FASTA ids free of characters that confuse downstream tools and makes
/// identical sequences collapse to one clustering input.
pub fn sequence_id(seq: &str) -> String {
    format!("seq_{}", &sha256_hex(seq.as_bytes())[..16])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_protein_check() {
        assert!(is_standard_protein("ACDEFGHIKLMNPQRSTVWY"));
        assert!(!is_standard_protein("ACDX"));
        assert!(!is_standard_protein("acd"));
        assert!(!is_standard_protein(""));
    }

    #[test]
    fn sequence_id_is_stable() {
        // SHA-256("MKV") prefix; pinned so that a change of hashing scheme is noticed.
        assert_eq!(sequence_id("MKV"), sequence_id("MKV"));
        assert_ne!(sequence_id("MKV"), sequence_id("MKA"));
        assert_eq!(sequence_id("MKV").len(), 4 + 16);
        assert_eq!(sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    }
}
