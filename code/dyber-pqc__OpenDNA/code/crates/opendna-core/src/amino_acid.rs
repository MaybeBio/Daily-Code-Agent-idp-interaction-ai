use serde::{Deserialize, Serialize};
use std::fmt;

/// One of the 20 standard proteinogenic amino acids, plus `Unknown` (`X`).
///
/// Selenocysteine (`U`) and pyrrolysine (`O`) are intentionally mapped to
/// [`AminoAcid::Unknown`] by the lenient parser and rejected by the strict
/// one: none of the physicochemical tables used downstream define them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(into = "char", try_from = "char")]
pub enum AminoAcid {
    /// Alanine
    Ala,
    /// Arginine
    Arg,
    /// Asparagine
    Asn,
    /// Aspartate
    Asp,
    /// Cysteine
    Cys,
    /// Glutamine
    Gln,
    /// Glutamate
    Glu,
    /// Glycine
    Gly,
    /// Histidine
    His,
    /// Isoleucine
    Ile,
    /// Leucine
    Leu,
    /// Lysine
    Lys,
    /// Methionine
    Met,
    /// Phenylalanine
    Phe,
    /// Proline
    Pro,
    /// Serine
    Ser,
    /// Threonine
    Thr,
    /// Tryptophan
    Trp,
    /// Tyrosine
    Tyr,
    /// Valine
    Val,
    /// Unknown or non-standard residue.
    Unknown,
}

impl AminoAcid {
    /// The 20 standard amino acids in alphabetical one-letter order (ACDEFGHIKLMNPQRSTVWY).
    pub const STANDARD: [AminoAcid; 20] = [
        AminoAcid::Ala,
        AminoAcid::Cys,
        AminoAcid::Asp,
        AminoAcid::Glu,
        AminoAcid::Phe,
        AminoAcid::Gly,
        AminoAcid::His,
        AminoAcid::Ile,
        AminoAcid::Lys,
        AminoAcid::Leu,
        AminoAcid::Met,
        AminoAcid::Asn,
        AminoAcid::Pro,
        AminoAcid::Gln,
        AminoAcid::Arg,
        AminoAcid::Ser,
        AminoAcid::Thr,
        AminoAcid::Val,
        AminoAcid::Trp,
        AminoAcid::Tyr,
    ];

    /// Parse a one-letter code (case-insensitive). Returns `None` for anything
    /// other than the 20 standard letters and `X`.
    pub fn from_one_letter(code: char) -> Option<Self> {
        Some(match code.to_ascii_uppercase() {
            'A' => Self::Ala,
            'R' => Self::Arg,
            'N' => Self::Asn,
            'D' => Self::Asp,
            'C' => Self::Cys,
            'Q' => Self::Gln,
            'E' => Self::Glu,
            'G' => Self::Gly,
            'H' => Self::His,
            'I' => Self::Ile,
            'L' => Self::Leu,
            'K' => Self::Lys,
            'M' => Self::Met,
            'F' => Self::Phe,
            'P' => Self::Pro,
            'S' => Self::Ser,
            'T' => Self::Thr,
            'W' => Self::Trp,
            'Y' => Self::Tyr,
            'V' => Self::Val,
            'X' => Self::Unknown,
            _ => return None,
        })
    }

    /// Parse a PDB three-letter residue name. Common modified residues
    /// (MSE, SEP, TPO, PTR, HSD/HSE/HSP, HID/HIE/HIP, CYX) map to their parent.
    /// Anything else returns `None`, so callers can tell ligands and water
    /// apart from protein residues.
    pub fn from_three_letter(name: &str) -> Option<Self> {
        Some(match name.trim().to_ascii_uppercase().as_str() {
            "ALA" => Self::Ala,
            "ARG" => Self::Arg,
            "ASN" => Self::Asn,
            "ASP" => Self::Asp,
            "CYS" | "CYX" => Self::Cys,
            "GLN" => Self::Gln,
            "GLU" => Self::Glu,
            "GLY" => Self::Gly,
            "HIS" | "HSD" | "HSE" | "HSP" | "HID" | "HIE" | "HIP" => Self::His,
            "ILE" => Self::Ile,
            "LEU" => Self::Leu,
            "LYS" => Self::Lys,
            "MET" | "MSE" => Self::Met,
            "PHE" => Self::Phe,
            "PRO" => Self::Pro,
            "SER" | "SEP" => Self::Ser,
            "THR" | "TPO" => Self::Thr,
            "TRP" => Self::Trp,
            "TYR" | "PTR" => Self::Tyr,
            "VAL" => Self::Val,
            "UNK" => Self::Unknown,
            _ => return None,
        })
    }

    /// One-letter code.
    pub fn one_letter(self) -> char {
        match self {
            Self::Ala => 'A',
            Self::Arg => 'R',
            Self::Asn => 'N',
            Self::Asp => 'D',
            Self::Cys => 'C',
            Self::Gln => 'Q',
            Self::Glu => 'E',
            Self::Gly => 'G',
            Self::His => 'H',
            Self::Ile => 'I',
            Self::Leu => 'L',
            Self::Lys => 'K',
            Self::Met => 'M',
            Self::Phe => 'F',
            Self::Pro => 'P',
            Self::Ser => 'S',
            Self::Thr => 'T',
            Self::Trp => 'W',
            Self::Tyr => 'Y',
            Self::Val => 'V',
            Self::Unknown => 'X',
        }
    }

    /// Canonical PDB three-letter code.
    pub fn three_letter(self) -> &'static str {
        match self {
            Self::Ala => "ALA",
            Self::Arg => "ARG",
            Self::Asn => "ASN",
            Self::Asp => "ASP",
            Self::Cys => "CYS",
            Self::Gln => "GLN",
            Self::Glu => "GLU",
            Self::Gly => "GLY",
            Self::His => "HIS",
            Self::Ile => "ILE",
            Self::Leu => "LEU",
            Self::Lys => "LYS",
            Self::Met => "MET",
            Self::Phe => "PHE",
            Self::Pro => "PRO",
            Self::Ser => "SER",
            Self::Thr => "THR",
            Self::Trp => "TRP",
            Self::Tyr => "TYR",
            Self::Val => "VAL",
            Self::Unknown => "UNK",
        }
    }

    /// `true` for the 20 standard amino acids.
    pub fn is_standard(self) -> bool {
        self != Self::Unknown
    }
}

impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.one_letter())
    }
}

impl From<AminoAcid> for char {
    fn from(aa: AminoAcid) -> char {
        aa.one_letter()
    }
}

impl TryFrom<char> for AminoAcid {
    type Error = String;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        Self::from_one_letter(c).ok_or_else(|| format!("invalid amino-acid code {c:?}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_and_three_letter_codes_round_trip() {
        for aa in AminoAcid::STANDARD {
            assert_eq!(AminoAcid::from_one_letter(aa.one_letter()), Some(aa));
            assert_eq!(AminoAcid::from_three_letter(aa.three_letter()), Some(aa));
        }
    }

    #[test]
    fn standard_table_is_alphabetical_and_complete() {
        let letters: String = AminoAcid::STANDARD.iter().map(|a| a.one_letter()).collect();
        assert_eq!(letters, "ACDEFGHIKLMNPQRSTVWY");
    }

    #[test]
    fn modified_residues_map_to_parent() {
        assert_eq!(AminoAcid::from_three_letter("MSE"), Some(AminoAcid::Met));
        assert_eq!(AminoAcid::from_three_letter("HOH"), None);
    }
}
