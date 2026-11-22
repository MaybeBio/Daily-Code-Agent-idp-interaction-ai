use thiserror::Error;

/// Reasons an analysis cannot produce a result.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum AnalysisError {
    /// The sequence contains `X`, which has no defined physicochemical values.
    #[error("non-standard residue at position {position}; property is undefined")]
    NonStandardResidue {
        /// Zero-based position.
        position: usize,
    },
    /// The input has nothing to analyse.
    #[error("input is empty")]
    Empty,
    /// Two coordinate sets that must correspond have different lengths.
    #[error("coordinate sets differ in length ({0} vs {1})")]
    LengthMismatch(usize, usize),
    /// A parameter is out of its valid range.
    #[error("invalid parameter: {0}")]
    InvalidParameter(&'static str),
    /// A stored file is not what it claims to be, or contradicts itself.
    ///
    /// Carries a message because the useful part is which field was wrong and
    /// what it said, and that is only known at the point of reading.
    #[error("{0}")]
    Malformed(String),
}
