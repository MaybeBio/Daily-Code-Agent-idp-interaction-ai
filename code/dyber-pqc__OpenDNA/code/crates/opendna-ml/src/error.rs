use thiserror::Error;

/// Errors from model loading and inference.
#[derive(Debug, Error)]
pub enum MlError {
    /// A tensor operation failed.
    #[error("tensor operation failed: {0}")]
    Tensor(#[from] candle_core::Error),
    /// The model directory is missing a file, or it could not be read.
    #[error("{path}: {source}")]
    File {
        /// File that could not be read.
        path: String,
        /// Underlying error.
        source: std::io::Error,
    },
    /// `config.json` did not parse, or described a model this code cannot run.
    #[error("unsupported model configuration: {0}")]
    Config(String),
    /// Downloading weights failed.
    #[error("could not fetch {repo} from the Hugging Face Hub: {message}")]
    Download {
        /// Repository identifier.
        repo: String,
        /// Underlying error.
        message: String,
    },
    /// The input could not be tokenised.
    #[error("sequence {index} is empty")]
    EmptySequence {
        /// Index within the batch.
        index: usize,
    },
    /// A sequence exceeded the model's position limit.
    #[error("sequence {index} has {length} residues; the model supports at most {limit}")]
    SequenceTooLong {
        /// Index within the batch.
        index: usize,
        /// Residue count.
        length: usize,
        /// Maximum supported.
        limit: usize,
    },
}
