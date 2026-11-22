//! FASTA reader and writer.

use crate::IoError;
use opendna_core::{Protein, Sequence};
use std::fmt::Write;

/// Parse FASTA text into proteins. Lines starting with `;` are comments.
/// Sequences are parsed strictly; see [`Sequence::parse`].
///
/// # Errors
/// [`IoError::Malformed`] if sequence data appears before the first header,
/// [`IoError::InvalidSequence`] for a bad record body, and
/// [`IoError::Empty`] if there are no records.
pub fn read(content: &str) -> Result<Vec<Protein>, IoError> {
    let mut records: Vec<(String, String)> = Vec::new();
    for (idx, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with(';') {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            records.push((header.trim().to_string(), String::new()));
        } else if let Some((_, body)) = records.last_mut() {
            body.push_str(line);
        } else {
            return Err(IoError::malformed(idx + 1, "sequence data before first '>' header"));
        }
    }
    if records.is_empty() {
        return Err(IoError::Empty("FASTA records"));
    }
    records
        .into_iter()
        .map(|(name, body)| match Sequence::parse(&body) {
            Ok(seq) => Ok(Protein::new(name, seq)),
            Err(source) => Err(IoError::InvalidSequence { name, source }),
        })
        .collect()
}

/// Write proteins as FASTA with 60 residues per line.
pub fn write(proteins: &[Protein]) -> String {
    let mut out = String::new();
    for p in proteins {
        let _ = writeln!(out, ">{}", p.name);
        let s = p.sequence.to_string();
        for chunk in s.as_bytes().chunks(60) {
            out.push_str(std::str::from_utf8(chunk).unwrap_or_default());
            out.push('\n');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn reads_multi_line_records() {
        let ps = read(">sp|P1|one\nMKTV\nRQER\n\n>two desc\nACDE\n").unwrap();
        assert_eq!(ps.len(), 2);
        assert_eq!(ps[0].name, "sp|P1|one");
        assert_eq!(ps[0].sequence.to_string(), "MKTVRQER");
    }

    #[test]
    fn rejects_body_without_header() {
        assert!(matches!(read("MKTV\n"), Err(IoError::Malformed { .. })));
    }

    #[test]
    fn names_the_bad_record() {
        let err = read(">ok\nMK\n>bad\nMK1Z\n").unwrap_err();
        assert!(matches!(err, IoError::InvalidSequence { ref name, .. } if name == "bad"), "{err}");
    }

    proptest! {
        #[test]
        fn write_read_round_trip(seqs in prop::collection::vec("[ACDEFGHIKLMNPQRSTVWY]{1,300}", 1..5)) {
            let proteins: Vec<_> = seqs.iter().enumerate()
                .map(|(i, s)| Protein::new(format!("p{i}"), Sequence::parse(s).unwrap()))
                .collect();
            let back = read(&write(&proteins)).unwrap();
            prop_assert_eq!(back, proteins);
        }
    }
}

/// Parse FASTA into `(header, sequence)` pairs exactly as written, without
/// interpreting the residues.
///
/// [`read`] builds [`Protein`]s, which means every letter has to be one of
/// the twenty standard amino acids or `X`. That is the right strictness when
/// importing something a user will work on, and the wrong strictness for
/// handing sequences to a protein language model: real target files contain
/// selenocysteine (`U`), pyrrolysine (`O`) and the ambiguity codes `B` and
/// `Z`, and ESM-2's alphabet has distinct tokens for all of them. Collapsing
/// them to `X` first would feed the model a different residue than the file
/// specifies, and rejecting them outright means refusing real data — CAID's
/// own disorder targets include a `U`.
///
/// Whitespace and digits are stripped, as in [`read`]. Nothing else is
/// validated; the caller's tokenizer decides what a letter means.
///
/// # Errors
/// [`IoError::Malformed`] if sequence data appears before the first header,
/// and [`IoError::Empty`] if there are no records.
pub fn read_raw(content: &str) -> Result<Vec<(String, String)>, IoError> {
    let mut records: Vec<(String, String)> = Vec::new();
    for (idx, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with(';') {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            records.push((header.trim().to_string(), String::new()));
        } else if let Some((_, body)) = records.last_mut() {
            body.extend(line.chars().filter(|c| !c.is_whitespace() && !c.is_ascii_digit()));
        } else {
            return Err(IoError::malformed(idx + 1, "sequence data before first '>' header"));
        }
    }
    if records.is_empty() {
        return Err(IoError::Empty("FASTA records"));
    }
    Ok(records)
}

#[cfg(test)]
mod raw_tests {
    use super::*;

    /// CAID's own disorder targets contain a selenocysteine. The strict
    /// reader refuses them, which made the CLI unable to read the files it
    /// would be submitted against.
    #[test]
    fn non_standard_residues_are_preserved_rather_than_rejected() {
        let got = read_raw(">DP03745 desc\nMSILLSPPSU\nOBZX\n").unwrap();
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].0, "DP03745 desc");
        assert_eq!(got[0].1, "MSILLSPPSUOBZX");
        assert!(read(">DP03745\nMSILLSPPSU\n").is_err(), "strict reader should still refuse U");
    }

    #[test]
    fn whitespace_and_digits_are_stripped() {
        let got = read_raw(">a\n  MK V 12\nLA\n").unwrap();
        assert_eq!(got[0].1, "MKVLA");
    }

    #[test]
    fn a_body_before_a_header_is_refused() {
        assert!(read_raw("MKV\n>a\nMKV\n").is_err());
        assert!(read_raw("").is_err());
    }
}
