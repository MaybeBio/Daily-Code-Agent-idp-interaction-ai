use crate::{Sequence, Structure};
use serde::{Deserialize, Serialize};

/// A named protein: a sequence and, optionally, a structure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Protein {
    /// Human-readable name (e.g. a FASTA header or PDB title).
    pub name: String,
    /// Amino-acid sequence.
    pub sequence: Sequence,
    /// 3D structure, if known.
    pub structure: Option<Structure>,
}

impl Protein {
    /// A protein with no structure.
    pub fn new(name: impl Into<String>, sequence: Sequence) -> Self {
        Self { name: name.into(), sequence, structure: None }
    }

    /// Content identifier derived from the sequence. See [`Sequence::digest`].
    pub fn id(&self) -> String {
        self.sequence.digest()
    }
}
