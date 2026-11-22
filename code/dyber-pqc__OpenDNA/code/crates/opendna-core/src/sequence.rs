use crate::{AminoAcid, CoreError};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// An amino-acid sequence.
///
/// Construct with [`Sequence::parse`] (strict: only the 20 standard letters
/// and `X`) or [`Sequence::parse_lenient`] (maps any other letter to `X`).
/// Whitespace and digits are ignored by both, so FASTA bodies and numbered
/// alignment rows can be passed in directly.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
#[serde(into = "String", try_from = "String")]
pub struct Sequence(Vec<AminoAcid>);

impl Sequence {
    /// Strict parse. Fails on the first character that isn't a valid code.
    ///
    /// # Errors
    /// [`CoreError::InvalidResidue`] for an unknown letter, or
    /// [`CoreError::EmptySequence`] if no residues remain after stripping
    /// whitespace and digits.
    pub fn parse(s: &str) -> Result<Self, CoreError> {
        let mut out = Vec::with_capacity(s.len());
        for (position, code) in s.chars().enumerate() {
            if code.is_whitespace() || code.is_ascii_digit() {
                continue;
            }
            let aa = AminoAcid::from_one_letter(code).ok_or(CoreError::InvalidResidue { code, position })?;
            out.push(aa);
        }
        if out.is_empty() {
            return Err(CoreError::EmptySequence);
        }
        Ok(Self(out))
    }

    /// Lenient parse: unrecognised letters become [`AminoAcid::Unknown`].
    /// Non-alphabetic characters are skipped.
    pub fn parse_lenient(s: &str) -> Self {
        Self(
            s.chars()
                .filter(char::is_ascii_alphabetic)
                .map(|c| AminoAcid::from_one_letter(c).unwrap_or(AminoAcid::Unknown))
                .collect(),
        )
    }

    /// Build from residues directly.
    pub fn from_residues(residues: Vec<AminoAcid>) -> Self {
        Self(residues)
    }

    /// Residues as a slice.
    pub fn residues(&self) -> &[AminoAcid] {
        &self.0
    }

    /// Number of residues.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// `true` if there are no residues.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// `true` if every residue is one of the 20 standard amino acids.
    pub fn is_standard(&self) -> bool {
        self.0.iter().all(|aa| aa.is_standard())
    }

    /// Count of each standard amino acid, indexed in [`AminoAcid::STANDARD`] order.
    pub fn composition(&self) -> [usize; 20] {
        let mut counts = [0usize; 20];
        for aa in &self.0 {
            if let Some(i) = AminoAcid::STANDARD.iter().position(|s| s == aa) {
                counts[i] += 1;
            }
        }
        counts
    }

    /// Number of occurrences of `aa`.
    pub fn count(&self, aa: AminoAcid) -> usize {
        self.0.iter().filter(|&&r| r == aa).count()
    }

    /// Stable content identifier: the first 16 hex digits of SHA-256 over the
    /// one-letter string. Two sequences with identical residues always share it.
    pub fn digest(&self) -> String {
        use sha2::{Digest, Sha256};
        use std::fmt::Write;
        let hash = Sha256::digest(self.to_string().as_bytes());
        hash.iter().take(8).fold(String::with_capacity(16), |mut s, b| {
            let _ = write!(s, "{b:02x}");
            s
        })
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for aa in &self.0 {
            write!(f, "{}", aa.one_letter())?;
        }
        Ok(())
    }
}

impl FromStr for Sequence {
    type Err = CoreError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

impl From<Sequence> for String {
    fn from(s: Sequence) -> String {
        s.to_string()
    }
}

impl TryFrom<String> for Sequence {
    type Error = CoreError;
    fn try_from(s: String) -> Result<Self, Self::Error> {
        Self::parse(&s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn strict_parse_reports_position() {
        let err = Sequence::parse("MKT ZV").unwrap_err();
        assert_eq!(err, CoreError::InvalidResidue { code: 'Z', position: 4 });
    }

    #[test]
    fn whitespace_and_numbers_are_ignored() {
        let s = Sequence::parse("1 MKTV\n 11 RQER").unwrap();
        assert_eq!(s.to_string(), "MKTVRQER");
    }

    #[test]
    fn empty_is_an_error() {
        assert_eq!(Sequence::parse("  \n"), Err(CoreError::EmptySequence));
    }

    #[test]
    fn lenient_parse_maps_unknown_letters_to_x() {
        assert_eq!(Sequence::parse_lenient("MUB*").to_string(), "MXX");
    }

    #[test]
    fn serde_uses_plain_string() {
        let s = Sequence::parse("ACDE").unwrap();
        assert_eq!(serde_json::to_string(&s).unwrap(), "\"ACDE\"");
    }

    proptest! {
        #[test]
        fn display_then_parse_round_trips(s in "[ACDEFGHIKLMNPQRSTVWY]{1,200}") {
            let seq = Sequence::parse(&s).unwrap();
            prop_assert_eq!(seq.to_string(), s);
        }

        #[test]
        fn composition_sums_to_length(s in "[ACDEFGHIKLMNPQRSTVWY]{1,200}") {
            let seq = Sequence::parse(&s).unwrap();
            prop_assert_eq!(seq.composition().iter().sum::<usize>(), seq.len());
        }
    }
}
