use crate::{AminoAcid, Sequence, Vec3};
use serde::{Deserialize, Serialize};
use std::fmt;

/// A macromolecular structure: one model, made of chains.
///
/// Parsers keep only the first model of multi-model files and the first
/// alternate location of each atom, so every atom here has a single position.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct Structure {
    /// Chains in file order.
    pub chains: Vec<Chain>,
}

/// A polymer chain (or a group of ligands/water sharing a chain id).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Chain {
    /// Chain identifier (`auth_asym_id` in mmCIF, column 22 in PDB).
    pub id: String,
    /// Residues in file order.
    pub residues: Vec<Residue>,
}

/// Position of a residue within its chain as numbered by the depositor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct ResidueId {
    /// Residue sequence number.
    pub number: i32,
    /// Insertion code, if any.
    pub insertion: Option<char>,
}

impl fmt::Display for ResidueId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.insertion {
            Some(c) => write!(f, "{}{c}", self.number),
            None => write!(f, "{}", self.number),
        }
    }
}

/// A residue: an amino acid, nucleotide, ligand or water molecule.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Residue {
    /// Three-letter (or longer, for ligands) component name.
    pub name: String,
    /// Depositor numbering.
    pub id: ResidueId,
    /// `true` for HETATM records.
    pub hetero: bool,
    /// Atoms in file order.
    pub atoms: Vec<Atom>,
}

/// A single atom.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    /// Serial number from the file.
    pub serial: u32,
    /// Atom name (e.g. `CA`, `OG1`).
    pub name: String,
    /// Element symbol, upper-case (e.g. `C`, `FE`).
    pub element: String,
    /// Position in Å.
    pub position: Vec3,
    /// Occupancy, 0–1.
    pub occupancy: f32,
    /// B-factor. For predicted models this column usually holds pLDDT.
    pub b_factor: f32,
}

impl Structure {
    /// Total number of atoms.
    pub fn atom_count(&self) -> usize {
        self.chains.iter().flat_map(|c| &c.residues).map(|r| r.atoms.len()).sum()
    }

    /// Iterate over every residue together with its chain.
    pub fn residues(&self) -> impl Iterator<Item = (&Chain, &Residue)> {
        self.chains.iter().flat_map(|c| c.residues.iter().map(move |r| (c, r)))
    }

    /// Iterate over every atom.
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.chains.iter().flat_map(|c| &c.residues).flat_map(|r| &r.atoms)
    }

    /// Look up a chain by id.
    pub fn chain(&self, id: &str) -> Option<&Chain> {
        self.chains.iter().find(|c| c.id == id)
    }
}

impl Chain {
    /// Residues that are standard or modified amino acids (see
    /// [`AminoAcid::from_three_letter`]) and have a Cα atom.
    pub fn amino_acids(&self) -> impl Iterator<Item = &Residue> {
        self.residues.iter().filter(|r| r.amino_acid().is_some() && r.atom("CA").is_some())
    }

    /// Sequence of the observed amino-acid residues (not SEQRES).
    pub fn sequence(&self) -> Sequence {
        Sequence::from_residues(self.amino_acids().filter_map(Residue::amino_acid).collect())
    }
}

impl Residue {
    /// Amino-acid identity, or `None` for ligands, water and nucleotides.
    pub fn amino_acid(&self) -> Option<AminoAcid> {
        AminoAcid::from_three_letter(&self.name)
    }

    /// Find an atom by name.
    pub fn atom(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().find(|a| a.name == name)
    }

    /// `true` for water molecules.
    pub fn is_water(&self) -> bool {
        matches!(self.name.as_str(), "HOH" | "WAT" | "DOD" | "H2O")
    }
}
