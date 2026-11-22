use thiserror::Error;

/// Errors produced when constructing core types from untrusted input.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum CoreError {
    /// A character that is not a recognised one-letter amino-acid code.
    #[error("invalid amino-acid code {code:?} at position {position}")]
    InvalidResidue {
        /// The offending character.
        code: char,
        /// Zero-based position in the input.
        position: usize,
    },
    /// The input contained no residues.
    #[error("sequence is empty")]
    EmptySequence,
}
