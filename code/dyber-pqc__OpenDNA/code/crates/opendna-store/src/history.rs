use chrono::{DateTime, Utc};
use opendna_core::{AminoAcid, Sequence};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;
use thiserror::Error;

/// Content-addressed version id (hex SHA-256 prefix of parent, sequence and message).
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct VersionId(String);

impl fmt::Display for VersionId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// A point substitution, written in standard notation (e.g. `K48R`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Mutation {
    /// Zero-based position.
    pub position: usize,
    /// Wild-type residue.
    pub from: AminoAcid,
    /// Substituted residue.
    pub to: AminoAcid,
}

impl fmt::Display for Mutation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.from, self.position + 1, self.to)
    }
}

/// One committed state of the sequence.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Version {
    /// Id of this version.
    pub id: VersionId,
    /// Parent version, `None` for the root.
    pub parent: Option<VersionId>,
    /// Sequence at this version.
    pub sequence: Sequence,
    /// Mutations relative to the parent.
    pub mutations: Vec<Mutation>,
    /// Commit message.
    pub message: String,
    /// Author name.
    pub author: String,
    /// Commit time.
    pub timestamp: DateTime<Utc>,
}

/// Errors from [`History`].
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum HistoryError {
    /// Branch name not found.
    #[error("no branch named {0:?}")]
    UnknownBranch(String),
    /// Branch already exists.
    #[error("branch {0:?} already exists")]
    BranchExists(String),
    /// A mutation does not match the current sequence.
    #[error("mutation {mutation} does not apply: {reason}")]
    InvalidMutation {
        /// The offending mutation.
        mutation: Mutation,
        /// Why.
        reason: &'static str,
    },
}

/// Version history of a single protein sequence.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct History {
    versions: BTreeMap<VersionId, Version>,
    branches: BTreeMap<String, VersionId>,
    current: String,
}

fn version_id(
    parent: Option<&VersionId>,
    sequence: &Sequence,
    message: &str,
    timestamp: DateTime<Utc>,
) -> VersionId {
    use sha2::{Digest, Sha256};
    use std::fmt::Write;
    let mut h = Sha256::new();
    h.update(parent.map_or("", |p| p.0.as_str()));
    h.update([0]);
    h.update(sequence.to_string());
    h.update([0]);
    h.update(message);
    h.update([0]);
    h.update(timestamp.to_rfc3339());
    let digest = h.finalize();
    VersionId(digest.iter().take(8).fold(String::new(), |mut s, b| {
        let _ = write!(s, "{b:02x}");
        s
    }))
}

impl History {
    /// Start a history on branch `main` with `sequence` as the root version.
    pub fn new(sequence: Sequence, author: impl Into<String>) -> Self {
        let timestamp = Utc::now();
        let message = "initial sequence".to_string();
        let id = version_id(None, &sequence, &message, timestamp);
        let root = Version {
            id: id.clone(),
            parent: None,
            sequence,
            mutations: Vec::new(),
            message,
            author: author.into(),
            timestamp,
        };
        Self {
            versions: BTreeMap::from([(id.clone(), root)]),
            branches: BTreeMap::from([("main".to_string(), id)]),
            current: "main".to_string(),
        }
    }

    /// Name of the checked-out branch.
    pub fn current_branch(&self) -> &str {
        &self.current
    }

    /// Head version of the checked-out branch.
    pub fn head(&self) -> &Version {
        // Invariant: `current` always names a branch whose head is in `versions`.
        &self.versions[&self.branches[&self.current]]
    }

    /// Apply point mutations to the head and commit the result. Returns the new id.
    ///
    /// # Errors
    /// [`HistoryError::InvalidMutation`] if a position is out of range or the
    /// wild-type residue doesn't match.
    pub fn commit_mutations(
        &mut self,
        mutations: Vec<Mutation>,
        message: impl Into<String>,
        author: impl Into<String>,
    ) -> Result<VersionId, HistoryError> {
        let mut residues = self.head().sequence.residues().to_vec();
        for &m in &mutations {
            match residues.get_mut(m.position) {
                None => {
                    return Err(HistoryError::InvalidMutation {
                        mutation: m,
                        reason: "position out of range",
                    })
                }
                Some(r) if *r != m.from => {
                    return Err(HistoryError::InvalidMutation {
                        mutation: m,
                        reason: "wild-type residue mismatch",
                    })
                }
                Some(r) => *r = m.to,
            }
        }
        Ok(self.commit(Sequence::from_residues(residues), mutations, message.into(), author.into()))
    }

    fn commit(
        &mut self,
        sequence: Sequence,
        mutations: Vec<Mutation>,
        message: String,
        author: String,
    ) -> VersionId {
        let parent = self.head().id.clone();
        let timestamp = Utc::now();
        let id = version_id(Some(&parent), &sequence, &message, timestamp);
        self.versions.insert(
            id.clone(),
            Version { id: id.clone(), parent: Some(parent), sequence, mutations, message, author, timestamp },
        );
        self.branches.insert(self.current.clone(), id.clone());
        id
    }

    /// Create a branch at the current head.
    ///
    /// # Errors
    /// [`HistoryError::BranchExists`].
    pub fn create_branch(&mut self, name: impl Into<String>) -> Result<(), HistoryError> {
        let name = name.into();
        if self.branches.contains_key(&name) {
            return Err(HistoryError::BranchExists(name));
        }
        let head = self.head().id.clone();
        self.branches.insert(name, head);
        Ok(())
    }

    /// Switch branches.
    ///
    /// # Errors
    /// [`HistoryError::UnknownBranch`].
    pub fn checkout(&mut self, name: &str) -> Result<(), HistoryError> {
        if !self.branches.contains_key(name) {
            return Err(HistoryError::UnknownBranch(name.to_string()));
        }
        self.current = name.to_string();
        Ok(())
    }

    /// Versions from head back to the root.
    pub fn log(&self) -> Vec<&Version> {
        std::iter::successors(Some(self.head()), |v| v.parent.as_ref().and_then(|p| self.versions.get(p)))
            .collect()
    }

    /// Branch names.
    pub fn branches(&self) -> impl Iterator<Item = &str> {
        self.branches.keys().map(String::as_str)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ubq() -> Sequence {
        Sequence::parse("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
            .unwrap()
    }

    #[test]
    fn mutations_apply_and_are_logged() {
        let mut h = History::new(ubq(), "test");
        let k48r = Mutation { position: 47, from: AminoAcid::Lys, to: AminoAcid::Arg };
        assert_eq!(k48r.to_string(), "K48R");
        h.commit_mutations(vec![k48r], "remove K48 chain site", "test").unwrap();
        assert_eq!(h.head().sequence.residues()[47], AminoAcid::Arg);
        assert_eq!(h.log().len(), 2);
    }

    #[test]
    fn wrong_wild_type_is_rejected() {
        let mut h = History::new(ubq(), "test");
        let bad = Mutation { position: 0, from: AminoAcid::Ala, to: AminoAcid::Gly };
        assert!(matches!(h.commit_mutations(vec![bad], "x", "t"), Err(HistoryError::InvalidMutation { .. })));
        assert_eq!(h.log().len(), 1);
    }

    #[test]
    fn branches_diverge() {
        let mut h = History::new(ubq(), "t");
        h.create_branch("exp").unwrap();
        h.checkout("exp").unwrap();
        h.commit_mutations(
            vec![Mutation { position: 0, from: AminoAcid::Met, to: AminoAcid::Ala }],
            "m1a",
            "t",
        )
        .unwrap();
        h.checkout("main").unwrap();
        assert_eq!(h.head().sequence.residues()[0], AminoAcid::Met);
        assert_eq!(h.checkout("nope"), Err(HistoryError::UnknownBranch("nope".into())));
        assert_eq!(h.create_branch("exp"), Err(HistoryError::BranchExists("exp".into())));
    }

    #[test]
    fn identical_sequences_get_distinct_versions() {
        let mut h = History::new(ubq(), "t");
        let a = h.commit_mutations(vec![], "no-op 1", "t").unwrap();
        let b = h.commit_mutations(vec![], "no-op 2", "t").unwrap();
        assert_ne!(a, b);
    }
}
