//! IPC message types shared by the desktop backend and the WebAssembly UI.
//!
//! Every Tauri command takes and returns types from this crate, so the two
//! sides cannot drift apart: a field rename is a compile error on both.
//! Command names are constants in [`command`].

use serde::{Deserialize, Serialize};

/// Tauri command names.
pub mod command {
    /// `() -> Vec<ProteinSummary>`
    pub const LIST_PROTEINS: &str = "list_proteins";
    /// `{ path: Option<String> } -> Vec<ProteinSummary>`; `None` opens a file dialog.
    pub const IMPORT_FILE: &str = "import_file";
    /// `{ id } -> ProteinDetail`
    pub const GET_PROTEIN: &str = "get_protein";
    /// `{ id } -> ()`
    pub const DELETE_PROTEIN: &str = "delete_protein";
    /// `{ id } -> SequenceReport`
    pub const ANALYZE_SEQUENCE: &str = "analyze_sequence";
    /// `{ id } -> StructureReport`
    pub const ANALYZE_STRUCTURE: &str = "analyze_structure";
    /// `() -> ModelStatus`. Does not download anything.
    pub const MODEL_STATUS: &str = "model_status";
    /// `{ id, position } -> StabilityScan`
    pub const PREDICT_STABILITY: &str = "predict_stability";
    /// `{ id } -> SolubilityPrediction`
    pub const PREDICT_SOLUBILITY: &str = "predict_solubility";
    /// `{ id } -> DisorderPrediction`
    pub const PREDICT_DISORDER: &str = "predict_disorder";
    /// Per-residue solvent accessibility.
    pub const PREDICT_BURIAL: &str = "predict_burial";
    /// Save text the UI has rendered to a file the user chooses.
    pub const SAVE_TEXT_FILE: &str = "save_text_file";
}

/// One row in the project explorer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProteinSummary {
    /// Content id.
    pub id: String,
    /// Display name.
    pub name: String,
    /// Residue count.
    pub length: usize,
    /// Whether 3D coordinates are available.
    pub has_structure: bool,
}

/// Everything the main view needs to show one protein.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProteinDetail {
    /// Summary fields.
    pub summary: ProteinSummary,
    /// One-letter sequence.
    pub sequence: String,
    /// Coordinate file for the viewer.
    pub structure: Option<StructurePayload>,
}

/// A coordinate file as handed to the 3D viewer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StructurePayload {
    /// `pdb` or `mmcif`.
    pub format: String,
    /// Chain the sequence belongs to.
    pub chain: String,
    /// File contents.
    pub content: String,
}

/// Sequence-derived properties. See `opendna_analysis::properties`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SequenceReport {
    /// Residue count.
    pub length: usize,
    /// Average mass, Da.
    pub molecular_weight: f64,
    /// Isoelectric point.
    pub isoelectric_point: f64,
    /// Net charge at pH 7.
    pub charge_at_ph7: f64,
    /// Kyte–Doolittle grand average.
    pub gravy: f64,
    /// Fraction F+W+Y.
    pub aromaticity: f64,
    /// Ikai aliphatic index.
    pub aliphatic_index: f64,
    /// ε280, cysteines reduced (M⁻¹cm⁻¹).
    pub extinction_reduced: f64,
    /// ε280, all cystines formed.
    pub extinction_cystines: f64,
    /// Nine-residue Kyte–Doolittle window means, one per window start.
    pub hydropathy: Vec<f64>,
}

/// Per-residue structural annotation, aligned with the chain's observed residues.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResidueAnnotation {
    /// Depositor residue number (with insertion code), as shown to users.
    pub label: String,
    /// Residue number without insertion code, for viewer selection.
    pub number: i32,
    /// One-letter amino acid.
    pub amino_acid: char,
    /// DSSP code (`H`, `E`, `G`, `I`, `T`, `S`, `B`, `-`).
    pub dssp: char,
    /// Relative solvent accessibility (Tien 2013 maxima), if defined.
    pub relative_sasa: Option<f64>,
    /// φ in degrees.
    pub phi: Option<f64>,
    /// ψ in degrees.
    pub psi: Option<f64>,
}

/// A residue pair interaction for display.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Contact {
    /// e.g. `A:CYS3`
    pub a: String,
    /// e.g. `A:CYS40`
    pub b: String,
    /// Å
    pub distance: f64,
}

/// Structure-derived analysis of one chain.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StructureReport {
    /// Chain analysed.
    pub chain: String,
    /// One entry per observed amino-acid residue in the chain.
    pub residues: Vec<ResidueAnnotation>,
    /// Cα radius of gyration, Å.
    pub radius_of_gyration: f64,
    /// Total SASA of the whole structure (heavy atoms, no waters), Å².
    pub total_sasa: f64,
    /// Disulfide bonds anywhere in the structure.
    pub disulfides: Vec<Contact>,
    /// Salt bridges anywhere in the structure.
    pub salt_bridges: Vec<Contact>,
}

/// Error returned by every command. `message` is safe to show to users.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CommandError {
    /// Human-readable description.
    pub message: String,
}

/// Whether a trained head can be used, and if not, why not.
///
/// A head that is absent or mismatched is reported as such. The application
/// never substitutes a heuristic for a missing model, so the UI needs the
/// reason in order to say what is wrong rather than show an empty panel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum HeadStatus {
    /// Weights and metadata loaded and matched the embedding model.
    Ready {
        /// Version string from the head's metadata.
        version: String,
        /// The head's limitations, exactly as its model card states them.
        ///
        /// Sent from the backend rather than written into the UI, because the
        /// UI used to carry its own copy and a caveat that disagrees with the
        /// published card is worse than no caveat: a reader who checks the
        /// card finds different text and cannot tell which is current.
        limitations: Vec<String>,
    },
    /// The head's files are not in the models directory.
    Missing,
    /// The files are present but unusable — unreadable, or trained against a
    /// different embedding model.
    Unusable {
        /// What went wrong, safe to show to users.
        reason: String,
    },
}

impl HeadStatus {
    /// Whether a prediction can be requested.
    pub fn is_ready(&self) -> bool {
        matches!(self, Self::Ready { .. })
    }

    /// The head's limitations, or nothing if it is not usable.
    pub fn limitations(&self) -> &[String] {
        match self {
            Self::Ready { limitations, .. } => limitations,
            Self::Missing | Self::Unusable { .. } => &[],
        }
    }

    /// The head's version, or nothing if it is not usable.
    pub fn version(&self) -> Option<&str> {
        match self {
            Self::Ready { version, .. } => Some(version),
            Self::Missing | Self::Unusable { .. } => None,
        }
    }
}

/// What the machine-learning side of the application can currently do.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ModelStatus {
    /// Directory searched for trained heads.
    pub models_dir: String,
    /// Hugging Face repository of the embedding model.
    pub embedding_model: String,
    /// Whether the embedding model is already in the local cache. When false,
    /// the first prediction downloads it, which is slow and needs a network.
    pub embedding_cached: bool,
    /// Device predictions will run on, e.g. `CPU`.
    pub device: String,
    /// ΔΔG head.
    pub stability: HeadStatus,
    /// Solubility head.
    pub solubility: HeadStatus,
    /// Per-residue disorder head.
    pub disorder: HeadStatus,
    /// Per-residue solvent accessibility head.
    pub burial: HeadStatus,
}

/// Per-residue solvent accessibility for one protein.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BurialPrediction {
    /// Relative accessibility per residue, in sequence order: roughly 0 for a
    /// residue in the core and roughly 1 for one fully exposed. Can exceed 1,
    /// because the per-residue maxima it is divided by are theoretical.
    pub accessibility: Vec<f64>,
    /// Accessibility at or below which the interface calls a residue buried.
    /// The convention, not a property of the model.
    pub buried_at: f64,
    /// Whether the sequence was longer than the embedding model's limit and
    /// was covered by overlapping windows.
    pub windowed: bool,
}

/// One substitution's predicted effect on folding stability.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SubstitutionScore {
    /// Mutant residue, one-letter.
    pub mutant: char,
    /// Predicted ΔΔG in kcal/mol. Negative is destabilising.
    pub ddg: f64,
    /// The language model's `log P(mutant) - log P(wild type)` at this
    /// position, shown alongside because it is the zero-shot signal the head
    /// is built on and disagreement between the two is informative.
    pub log_odds: f64,
}

/// Every substitution at one position, ranked.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StabilityScan {
    /// 1-based position scanned.
    pub position: usize,
    /// Wild-type residue at that position.
    pub wild_type: char,
    /// The 19 substitutions, most stabilising first.
    pub scores: Vec<SubstitutionScore>,
}

/// Predicted solubility for a whole protein.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolubilityPrediction {
    /// Percentage on the source assay's scale; not clipped to 0-100.
    pub percent: f64,
}

/// Per-residue disorder probabilities.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DisorderPrediction {
    /// One probability per residue, in sequence order.
    pub probabilities: Vec<f64>,
    /// Probability at or above which a residue is called disordered.
    pub threshold: f64,
    /// Whether the sequence was longer than the embedding model's limit and
    /// was covered by overlapping windows.
    pub windowed: bool,
}
