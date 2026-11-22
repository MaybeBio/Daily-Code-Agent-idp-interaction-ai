//! What the exit status means.
//!
//! A pipeline cannot read an error message. It can only read a number, and
//! what it needs from that number is what to do next: retry, stop and fix the
//! installation, or stop and fix the data. "Non-zero" collapses all three into
//! one, which in practice means either retrying something that will never
//! succeed or giving up on something that would have.
//!
//! So failures carry a [`Fault`], and the status is derived from it. The set is
//! deliberately small. Every code has a different action behind it, and a code
//! nobody can act on differently is a code that should not exist.
//!
//! Anything not explicitly classified is [`Fault::Unexpected`] and exits 1.
//! That is the honest default: a failure we did not anticipate is not one we
//! can tell a caller how to handle, and quietly filing it under a known
//! category would be worse than saying so.

use std::fmt;

/// Why a run stopped, grouped by what the caller should do.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Fault {
    /// Something we did not anticipate. Exit 1.
    Unexpected,
    /// The request itself cannot be served: a flag combination that does not
    /// exist, a required argument missing. Exit 2, the same code `clap` uses
    /// for a usage error, because a caller should treat the two the same way.
    Usage,
    /// The installation is wrong: a head is missing or corrupt, the models do
    /// not match the lockfile, or `--offline` was given with nothing staged.
    /// Exit 3. Retrying will not help; the machine needs fixing.
    Install,
    /// The data is wrong: unparseable FASTA, a mutation naming a residue the
    /// sequence does not have, a sequence longer than the model accepts. Exit
    /// 4. Retrying will not help; the input needs fixing.
    Input,
    /// Checks ran and did not hold. Exit 5. The installation is intact and the
    /// data is fine, but this machine does not reproduce the recorded results,
    /// so nothing it produces should be used yet.
    Checks,
}

impl Fault {
    /// The process exit status.
    ///
    /// 2 is left alone: `clap` already uses it for a usage error, and taking
    /// it would make a mistyped flag indistinguishable from a real failure.
    pub fn code(self) -> u8 {
        match self {
            Self::Unexpected => 1,
            Self::Usage => 2,
            Self::Install => 3,
            Self::Input => 4,
            Self::Checks => 5,
        }
    }

    /// One line telling the reader which kind of problem this is.
    fn advice(self) -> &'static str {
        match self {
            Self::Unexpected => "This was not an anticipated failure; please report it.",
            Self::Usage => "The command asked for something that does not exist.",
            Self::Install => "The installation needs fixing. Retrying will not help.",
            Self::Input => "The input needs fixing. Retrying will not help.",
            Self::Checks => "The install is intact but does not reproduce its recorded results.",
        }
    }
}

impl fmt::Display for Fault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (exit {})", self.advice(), self.code())
    }
}

/// A failure and its classification.
pub struct Failure {
    /// What kind of problem it is.
    pub fault: Fault,
    /// The error itself, with its full context chain.
    pub error: anyhow::Error,
}

/// Anything not classified is unexpected.
///
/// This exists so `?` works unchanged throughout the program: a call site has
/// to do something deliberate to claim a specific fault, and forgetting is
/// safe rather than silently wrong.
impl From<anyhow::Error> for Failure {
    fn from(error: anyhow::Error) -> Self {
        Self { fault: Fault::Unexpected, error }
    }
}

/// Device selection and tensor work surface as candle errors.
impl From<candle_core::Error> for Failure {
    fn from(error: candle_core::Error) -> Self {
        Self { fault: Fault::Unexpected, error: error.into() }
    }
}

impl From<opendna_ml::MlError> for Failure {
    fn from(error: opendna_ml::MlError) -> Self {
        Self { fault: Fault::Unexpected, error: error.into() }
    }
}

/// Claim a fault for one fallible step.
pub trait Classify<T> {
    /// Attach `fault` to whatever this returns.
    fn blame(self, fault: Fault) -> Result<T, Failure>;
}

impl<T, E: Into<anyhow::Error>> Classify<T> for Result<T, E> {
    fn blame(self, fault: Fault) -> Result<T, Failure> {
        self.map_err(|e| Failure { fault, error: e.into() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two codes meaning the same thing, or one code meaning two things, would
    /// both make the contract useless to the pipeline it exists for.
    #[test]
    fn every_fault_has_its_own_code_and_none_is_clap_s() {
        let faults = [Fault::Unexpected, Fault::Usage, Fault::Install, Fault::Input, Fault::Checks];
        let mut codes: Vec<u8> = faults.iter().map(|f| f.code()).collect();
        codes.sort_unstable();
        codes.dedup();
        assert_eq!(codes.len(), faults.len(), "two faults share a code");
        assert!(!codes.contains(&0), "0 means success");
        // 2 is shared with clap on purpose: both mean the request was wrong.
        assert_eq!(Fault::Usage.code(), 2, "a usage error should match clap's own status");
    }

    #[test]
    fn an_unclassified_error_is_unexpected() {
        let failure: Failure = anyhow::anyhow!("something").into();
        assert_eq!(failure.fault, Fault::Unexpected);
        assert_eq!(failure.fault.code(), 1);
    }

    #[test]
    fn blaming_a_step_keeps_the_error_and_its_context() {
        let result: Result<(), anyhow::Error> =
            Err(anyhow::anyhow!("no disorder head").context("loading models"));
        let failure = result.blame(Fault::Install).unwrap_err();
        assert_eq!(failure.fault, Fault::Install);
        let chain = format!("{:#}", failure.error);
        assert!(chain.contains("loading models") && chain.contains("no disorder head"), "{chain}");
    }

    #[test]
    fn the_advice_says_whether_retrying_helps() {
        assert!(Fault::Install.to_string().contains("Retrying will not help"));
        assert!(Fault::Input.to_string().contains("Retrying will not help"));
        assert!(Fault::Install.to_string().contains("exit 3"));
    }
}
