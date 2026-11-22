use std::fmt;
use thiserror::Error;

/// Where in the input an error happened.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Location {
    /// One-based line number.
    pub line: usize,
}

impl fmt::Display for Location {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}", self.line)
    }
}

/// Errors from the readers in this crate.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum IoError {
    /// A record or value could not be parsed.
    #[error("{at}: {message}")]
    Malformed {
        /// Location of the problem.
        at: Location,
        /// What went wrong.
        message: String,
    },
    /// The file parsed but contained nothing usable.
    #[error("no {0} found in input")]
    Empty(&'static str),
    /// A FASTA record body was not a valid sequence.
    #[error("record {name:?}: {source}")]
    InvalidSequence {
        /// FASTA header of the bad record.
        name: String,
        /// Underlying error.
        source: opendna_core::CoreError,
    },
}

impl IoError {
    pub(crate) fn malformed(line: usize, message: impl Into<String>) -> Self {
        Self::Malformed { at: Location { line }, message: message.into() }
    }
}
