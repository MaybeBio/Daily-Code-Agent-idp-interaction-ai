from pydantic import BaseModel, Field, field_validator
from typing import List, Optional

# --- PfamDomain ---
class PfamDomainBase(BaseModel):
    pfam_accession: str
    pfam_name: str
    pfam_description: Optional[str] = None
    pfam_clan: Optional[str] = None
    seq_start: int
    seq_end: int
    hmm_start: int
    hmm_end: int
    hmm_length: int
    domain_coverage: float
    e_value: float
    bit_score: float
    domain_status: str  # legacy: "present" | "gained" | "lost" | "truncated"
    functional_class: Optional[str] = None  # see app/analysis/functional_classes.py
    alteration_type: Optional[str] = None   # intact / truncated_n / truncated_c / ...
    is_kinase: bool = False                 # legacy; functional_class is preferred

class PfamDomainResponse(PfamDomainBase):
    id: int
    isoform_id: int
    canonical_isoform_id: Optional[int] = None

    class Config:
        from_attributes = True

# --- PtmAccessibility ---
class PtmAccessibilityBase(BaseModel):
    uniprot_id: str
    residue_position: int
    novel_residue_position: Optional[int] = None
    residue_type: str
    ptm_type: str
    source_db: str
    canonical_sasa: Optional[float] = None
    novel_sasa: Optional[float] = None
    delta_sasa: Optional[float] = None
    # Δ / canonical_sasa — relative shift; more biologically meaningful than
    # the absolute Δ.
    relative_delta_sasa: Optional[float] = None
    # [0, 1] magnitude of the relative shift, capped at 1.0
    change_strength: Optional[float] = None
    # min(canonical_plddt, novel_plddt); < 70 → site treated as undetermined
    min_plddt: Optional[float] = None
    # intact | exposed | buried | lost | undetermined  (legacy: "unchanged")
    accessibility_change: str
    structural_context: Optional[str] = None

class PtmAccessibilityResponse(PtmAccessibilityBase):
    id: int
    isoform_id: int
    canonical_isoform_id: Optional[int] = None

    class Config:
        from_attributes = True

# --- TedDomain ---
class TedDomainBase(BaseModel):
    ted_id: str
    domain_type: str
    seq_start: int
    seq_end: int
    score: Optional[float] = None

class TedDomainResponse(TedDomainBase):
    id: int
    isoform_id: int

    class Config:
        from_attributes = True

# --- NmdPrediction ---
class NmdPredictionBase(BaseModel):
    # nmd_target | nmd_escape | non_translatable | predicted_translatable | unknown
    nmd_status: str
    has_ptc: bool
    truncation_ratio: Optional[float] = None
    orf_type: Optional[str] = None
    ptc_position: Optional[int] = None
    last_ejc_distance: Optional[int] = None
    evidence: str
    # Orthogonal length signal (intact | minor | severe | non_translatable)
    truncation_severity: Optional[str] = None
    # Confidence in [0,1] for nmd_target calls; null otherwise
    nmd_confidence: Optional[float] = None
    is_candidate_smorf: bool = False
    is_candidate_uorf: bool = False

class NmdPredictionResponse(NmdPredictionBase):
    id: int
    isoform_id: int
    canonical_isoform_id: Optional[int] = None

    class Config:
        from_attributes = True

# --- ExonStructuralImpact ---
class ExonStructuralImpactBase(BaseModel):
    residue_number: int
    dssp_assignment: str  # "H" | "E" | "C"
    phi_angle: Optional[float] = None
    psi_angle: Optional[float] = None
    residue_plddt: Optional[float] = None
    exon_number: Optional[int] = None
    canonical_dssp: Optional[str] = None
    structural_change: Optional[str] = None  # "same" | "helix_to_coil" | etc.

class ExonStructuralImpactResponse(ExonStructuralImpactBase):
    id: int
    isoform_id: int
    canonical_isoform_id: Optional[int] = None

    class Config:
        from_attributes = True

# --- CellTypeExpression ---
class CellTypeExpressionBase(BaseModel):
    cell_type: str
    expression_ncpm: float
    source: str
    specificity_category: Optional[str] = None

class CellTypeExpressionResponse(CellTypeExpressionBase):
    id: int
    gene_id: int

    class Config:
        from_attributes = True

# --- IsoformExpression ---
class IsoformExpressionBase(BaseModel):
    sample:       str
    condition:    str                   # "normal" | "tumor"
    fl_count:     int
    cpm:          Optional[float] = None
    isoform_frac: Optional[float] = None  # ψ (psi)
    n_cells:      Optional[int]   = None
    cell_frac:    Optional[float] = None

class IsoformExpressionResponse(IsoformExpressionBase):
    id:         int
    isoform_id: int

    class Config:
        from_attributes = True

# --- NeoantigenPrediction ---
class NeoantigenPredictionBase(BaseModel):
    peptide_sequence: str
    peptide_start: int
    peptide_end: int
    hla_allele: str
    binding_affinity_nm: float
    binding_category: str  # "strong_binder" | "weak_binder" | "non_binder"
    percentile_rank: Optional[float] = None
    junction_type: str  # "novel_insertion" | "novel_substitution" | "frameshift"
    # MHCflurry 2.2 presentation model scores
    presentation_score: Optional[float] = None  # combined presentation probability 0–1
    processing_score:   Optional[float] = None  # proteasomal cleavage score 0–1
    # PTM fields
    ptm_impact_status: Optional[str]   = None
    avg_plddt:         Optional[float] = None
    is_ptm_sensitive:  Optional[bool]  = False
    # Multi-tool ensemble (v5 — NetMHCpan-4.2 + BigMHC)
    netmhcpan_affinity_nm: Optional[float] = None  # NetMHCpan-4.2 BA IC50 nM
    netmhcpan_rank_el:     Optional[float] = None  # NetMHCpan-4.2 EL %rank
    netmhcpan_rank_ba:     Optional[float] = None  # NetMHCpan-4.2 BA %rank (rescue signal)
    bigmhc_el_score:       Optional[float] = None  # BigMHC EL presentation score (0–1)
    bigmhc_im_score:       Optional[float] = None  # BigMHC IM immunogenicity score (0–1)
    consensus_rank:        Optional[float] = None  # EXPERIMENTAL — incommensurable scales; do NOT use as primary ranking
    is_ba_rescue:          Optional[bool]  = False # Cysteine-rich peptide where NetMHCpan BA disagrees with EL
    # Audit-3 additions (Option A)
    tools_agreeing:        Optional[int]  = None   # 0–3 count of tools calling the peptide strong
    agreement_tier:        Optional[str]  = None   # 'high' | 'medium' | 'low'
    ptm_canonical_residue_lost: Optional[bool] = False  # Set only on Direct PTM-Neoantigen rows
    # ranking v2 (app/analysis/neoantigen_tiers.py)
    presentation_tier: Optional[str]   = None   # STRONG | MODERATE | WEAK | NONE (multi-tool presentation)
    priority_tier:     Optional[str]   = None   # P1 | P2 | P3 | P4 (presentation × immunogenicity)
    is_immunogenic:    Optional[bool]  = False  # BigMHC-IM ≥ threshold
    is_expressed:      Optional[bool]  = False  # isoform detected in expression data (decorating)
    min_plddt:         Optional[float] = None   # min pLDDT across epitope window
    # proteome-wide novelty of the peptide vs every GENCODE v47 reference protein
    proteome_uniqueness: Optional[str]  = None  # 'unique' | 'self_gene' | 'other_gene' | None (not computed)
    is_proteome_unique:  Optional[bool] = None

class NeoantigenPredictionResponse(NeoantigenPredictionBase):
    id: int
    isoform_id: int
    canonical_isoform_id: Optional[int] = None

    class Config:
        from_attributes = True

# --- IsoformJunction ---
class IsoformJunctionBase(BaseModel):
    junction_type: str
    chrom: Optional[str] = None
    strand: Optional[str] = None
    subcategory: Optional[str] = None
    exon_count_novel: Optional[int] = None
    exon_count_canonical: Optional[int] = None
    transcript_pos: Optional[int] = None
    is_novel: bool = False

class IsoformJunctionResponse(IsoformJunctionBase):
    id: int
    isoform_id: int

    class Config:
        from_attributes = True

# --- PPI Network (Legacy/General) ---
class PPINetworkBase(BaseModel):
    target_uniprot_id: Optional[str] = None
    interaction_score: Optional[float] = None
    interaction_retained: Optional[bool] = None

class PPINetworkResponse(PPINetworkBase):
    id: int
    isoform_id_1: Optional[int] = None

    class Config:
        from_attributes = True

# --- Homodimer Analysis (Direct Structural) ---
class HomodimerAnalysisBase(BaseModel):
    afdb_accession: Optional[str] = None
    total_interface_residues: Optional[int] = None
    count_absent: Optional[int] = None
    count_sasa_changed: Optional[int] = None
    disruption_score: Optional[float] = None
    interface_residues_json: Optional[str] = None

class HomodimerAnalysisResponse(HomodimerAnalysisBase):
    id: int
    isoform_id: int

    class Config:
        from_attributes = True


# --- Heterodimer Analysis (Direct Structural) ---
class HeterodimerAnalysisBase(BaseModel):
    partner_uniprot_id: Optional[str] = None
    partner_name: Optional[str] = None
    afdb_accession: Optional[str] = None
    total_interface_residues: Optional[int] = None
    count_absent: Optional[int] = None
    count_sasa_changed: Optional[int] = None
    disruption_score: Optional[float] = None
    interface_residues_json: Optional[str] = None

class HeterodimerAnalysisResponse(HeterodimerAnalysisBase):
    id: int
    isoform_id: int

    class Config:
        from_attributes = True

# --- STRING Interaction ---
class StringInteractionBase(BaseModel):
    partner_gene_symbol: Optional[str] = None
    partner_ensp_id: Optional[str] = None    # renamed from partner_uniprot_id in audit-5
    string_score: Optional[int] = None
    evidence_json: Optional[str] = None
    domain_loss_possible: Optional[bool] = False

class StringInteractionResponse(StringInteractionBase):
    id: int
    isoform_id: int

    class Config:
        from_attributes = True

# --- ClinicalFrequency ---
class ClinicalFrequencyBase(BaseModel):
    cancer_type: Optional[str] = None
    frequency: Optional[float] = None

class ClinicalFrequencyResponse(ClinicalFrequencyBase):
    id: int
    isoform_id: Optional[int] = None

    class Config:
        from_attributes = True

# --- DisorderScore ---
class DisorderScoreBase(BaseModel):
    disorder_scores: List[float]
    anchor_scores: List[float]

class DisorderScoreResponse(DisorderScoreBase):
    id: int
    isoform_id: int

    class Config:
        from_attributes = True

# --- StructuralMetrics ---
class StructuralMetricsBase(BaseModel):
    plddt_mean: Optional[float] = None
    plddt_min: Optional[float] = None
    plddt_median: Optional[float] = None
    pae_mean: Optional[float] = None
    pae_inter_domain_mean: Optional[float] = None
    pae_json_path: Optional[str] = None
    afdb_entry_id: Optional[str] = None
    afdb_version: int = 4
    plddt_series: Optional[List[float]] = None
    pae_blocks: Optional[str] = None

class StructuralMetricsResponse(StructuralMetricsBase):
    id: int
    isoform_id: Optional[int] = None

    class Config:
        from_attributes = True

# --- Structure ---
class StructureBase(BaseModel):
    storage_uri: Optional[str] = None
    af3_cif_url: Optional[str] = None
    af3_pae_url: Optional[str] = None
    confidence_score: Optional[float] = None
    avg_plddt: Optional[float] = None

class StructureResponse(StructureBase):
    id: int
    isoform_id: Optional[int] = None

    class Config:
        from_attributes = True

# --- Metadata ---
class TranscriptMetadataBase(BaseModel):
    # Pipeline that produced this row (sqanti3 / gencode / afdb_reference / custom / None).
    # Drives which downstream fields the UI should expose for this isoform.
    metadata_kind: Optional[str] = None

    # Universal (all pipelines)
    chrom: Optional[str] = None
    strand: Optional[str] = None
    length: Optional[int] = None
    biotype: Optional[str] = None
    orf_type: Optional[str] = None
    protein_length_aa: Optional[float] = None
    delta_length_aa: Optional[float] = None
    sample_source: Optional[str] = None
    associated_gene: Optional[str] = None
    canonical_gene_name: Optional[str] = None
    canonical_protein_length: Optional[float] = None

    # ccRCC folded-ORF QC (integrity + transparency)
    folded_orf_mismatch: Optional[bool] = None
    folded_best_ref: Optional[str] = None
    folded_best_ref_identity: Optional[float] = None

    # SQANTI3-specific (metadata_kind == "sqanti3")
    structural_category: Optional[str] = None
    subcategory: Optional[str] = None
    all_canonical: Optional[str] = None

    # GENCODE / Ensembl-specific (metadata_kind == "gencode")
    appris_tag: Optional[str] = None
    tsl: Optional[int] = None
    gencode_tags: Optional[list[str]] = None

    # Cross-database identifiers (universal-optional)
    uniprot_id: Optional[str] = None

class TranscriptMetadataResponse(TranscriptMetadataBase):
    id: int
    transcript_id: Optional[str] = None

    class Config:
        from_attributes = True

# --- Gene ---
class GeneBase(BaseModel):
    symbol: str
    hgnc_symbol: Optional[str] = None
    uniprot_id: Optional[str] = None
    ensembl_gene_id: Optional[str] = None
    description: Optional[str] = None
    protein_name: Optional[str] = None
    function_summary: Optional[str] = None
    canonical_sequence: Optional[str] = None
    canonical_length: Optional[int] = None
    chromosome: Optional[str] = None
    strand: Optional[str] = None
    genomic_start: Optional[int] = None
    genomic_end: Optional[int] = None
    canonical_transcript_id: Optional[int] = None

class GeneResponse(GeneBase):
    id: int

    class Config:
        from_attributes = True

class GeneDetailResponse(GeneResponse):
    isoforms: List["IsoformWithStructuresResponse"] = []

    class Config:
        from_attributes = True

# --- Isoform ---
class IsoformBase(BaseModel):
    transcript_id: str
    # Derived display name (MET-201 / MET-206-NOV-a6d00a). Null when a novel isoform has no
    # close GENCODE match — clients fall back to transcript_id. Not an identifier.
    transcript_name: Optional[str] = None
    sequence: Optional[str] = None
    fasta_sequence: Optional[str] = None
    delta_length: Optional[int] = None
    junction_coords: Optional[str] = None
    is_flagged: Optional[bool] = False
    in_primary_tumor: Optional[bool] = False
    in_normal: Optional[bool] = False
    in_metastasis: bool = False
    is_canonical: bool = False
    exons_json: Optional[str] = None

    # TransDecoder keeps the terminal stop codon as a trailing "*" on ccRCC
    # translations (~26k isoforms). It is not a residue: it inflates the
    # displayed length by one over the predicted structure and shows up as a
    # spurious character in sequence/alignment views. Strip it at serialization
    # so every client sees residues only. No internal "*" exist, so rstrip is safe.
    @field_validator("sequence", "fasta_sequence")
    @classmethod
    def _strip_trailing_stop(cls, v: Optional[str]) -> Optional[str]:
        return v.rstrip("*") if isinstance(v, str) else v

class IsoformResponse(IsoformBase):
    id: int
    gene_id: Optional[int] = None
    is_canonical: Optional[bool] = None
    source_id: Optional[int] = None
    source: Optional["SourceResponse"] = None

    class Config:
        from_attributes = True

class IsoformWithStructuresResponse(IsoformResponse):
    """IsoformResponse extended with structures, used inside GeneDetailResponse."""
    structures: List[StructureResponse] = []
    metadata_record: Optional["TranscriptMetadataResponse"] = None
    pfam_domains: List[PfamDomainResponse] = []
    nmd_prediction: Optional[NmdPredictionResponse] = None
    dataset_presence: List["IsoformSourceMembershipResponse"] = []

    class Config:
        from_attributes = True

class IsoformDetailResponse(IsoformResponse):
    structures: List[StructureResponse] = []
    clinical_frequencies: List[ClinicalFrequencyResponse] = []
    ppi_networks: List[PPINetworkResponse] = []
    pfam_domains: List[PfamDomainResponse] = []
    ptm_accessibility: List[PtmAccessibilityResponse] = []
    nmd_prediction: Optional[NmdPredictionResponse] = None
    structural_impacts: List[ExonStructuralImpactResponse] = []
    neoantigen_predictions: List[NeoantigenPredictionResponse] = []
    expression_data: List[IsoformExpressionResponse] = []
    metadata_record: Optional[TranscriptMetadataResponse] = None
    gene: Optional[GeneResponse] = None
    homodimer_analysis: Optional[HomodimerAnalysisResponse] = None
    heterodimer_analyses: Optional[List[HeterodimerAnalysisResponse]] = None
    string_interactions: List[StringInteractionResponse] = []
    ted_domains: List[TedDomainResponse] = []
    junction_events: List[IsoformJunctionResponse] = []
    disorder_score: Optional[DisorderScoreResponse] = None
    structural_metrics: Optional[StructuralMetricsResponse] = None
    canonical_pfam_domains: List[PfamDomainResponse] = []
    dataset_presence: List["IsoformSourceMembershipResponse"] = []

    class Config:
        from_attributes = True


# --- Source ---
class SourceResponse(BaseModel):
    id: int
    name: str

    class Config:
        from_attributes = True


class IsoformSourceMembershipResponse(BaseModel):
    source_id:        int
    source_name:      str
    in_normal:        bool
    in_primary_tumor: bool
    in_metastasis:    bool

    class Config:
        from_attributes = True

# --- Search Responses ---
class SearchResultItem(BaseModel):
    match_type: str  # "gene" or "transcript"
    id: int          # for transcripts this is the GENE id (row links to /gene/{id})
    gene_symbol: Optional[str] = None
    uniprot_id: Optional[str] = None
    transcript_id: Optional[str] = None
    transcript_name: Optional[str] = None   # display name (e.g. F2-201-NOV-48e626)
    isoform_id: Optional[int] = None        # numeric isoform id, for #isoforms=<id> selection
    source_name: Optional[str] = None
    description: Optional[str] = None

class DrugSearchItem(BaseModel):
    """A target–drug pair. `id` is the GENE id, so the row links to /drugs/{id}."""
    match_type: str = "drug"
    id: int
    gene_symbol: Optional[str] = None
    drug_name: str
    pubchem_cid: Optional[str] = None
    moa: Optional[str] = None
    n_complexes: int = 0

class SearchResponse(BaseModel):
    genes: List[SearchResultItem]
    transcripts: List[SearchResultItem]
    drugs: List[DrugSearchItem] = []
    total_genes: int = 0
    total_transcripts: int = 0
    total_drugs: int = 0

# --- Sequence Search (Phase 1: exact protein match by md5) ---
class SequenceSearchRequest(BaseModel):
    sequence: str  # raw protein sequence or FASTA (header + whitespace tolerated)

class SequenceHit(BaseModel):
    match_type: str = "exact"          # phase 2 adds "similar"
    isoform_id: int
    gene_id: Optional[int] = None      # row links to /gene/{gene_id}
    transcript_id: str
    transcript_name: Optional[str] = None
    gene_symbol: Optional[str] = None
    source_name: Optional[str] = None
    seq_length: int
    identity: Optional[float] = None   # reserved for phase 2 similarity (0..1)

class SequenceSearchResponse(BaseModel):
    query_length: int
    query_md5: str
    exact: List[SequenceHit] = []
    message: Optional[str] = None

# --- Ensembl ID resolution (unknown-id → gene, via cached Ensembl REST) ---
class ResolveResponse(BaseModel):
    query: str
    normalized_id: Optional[str] = None   # version-stripped base accession
    id_type: Optional[str] = None         # ensg / enst / ensp
    # in_db | resolved_gene | gene_not_in_db | not_found | not_ensembl | error
    status: str
    gene_id: Optional[int] = None         # our gene id, when resolvable in-DB
    gene_symbol: Optional[str] = None     # our gene symbol, else Ensembl display name
    resolved_ensg: Optional[str] = None
    isoform_count: Optional[int] = None
    message: Optional[str] = None
    # Most-similar isoform, computed LIVE against the current DB (never cached).
    closest_isoform_id: Optional[int] = None
    closest_transcript_id: Optional[str] = None
    closest_transcript_name: Optional[str] = None
    closest_gene_id: Optional[int] = None
    closest_gene_symbol: Optional[str] = None
    closest_identity: Optional[float] = None   # 0..1 (difflib approximate)
    closest_exact: Optional[bool] = None        # True = md5-identical protein

# --- Global Browse & Analytics ---
class BrowseIsoformSummary(BaseModel):
    id: int
    transcript_id: str
    # Derived display name (MET-201 / MET-206-NOV-a6d00a); null → client falls back to id.
    transcript_name: Optional[str] = None
    gene_id: int
    gene_symbol: Optional[str] = None
    uniprot_id: Optional[str] = None
    plddt_mean: Optional[float] = None
    delta_length: Optional[int] = None
    source_name: Optional[str] = None
    has_neoantigen: bool = False
    neoantigen_count: int = 0
    neoantigen_p1_count: int = 0
    num_domains_lost: int = 0
    num_domains_gained: int = 0
    nmd_status: Optional[str] = None
    has_ptm_exposure_change: bool = False
    has_ptm_buried_change: bool = False
    has_ptm_loss: bool = False
    isoform_count: int = 0
    impact_score: float = 0.0

class GlobalStats(BaseModel):
    plddt_distribution: List[dict]  # e.g. [{"bin": "0-10", "count": 100}, ...]
    impact_categories: List[dict]   # e.g. [{"category": "NMD Target", "count": 50}, ...]
    nmd_distribution: List[dict]    # e.g. [{"status": "nmd_target", "count": 50}, ...]
    ptm_distribution: List[dict]    # e.g. [{"impact": "exposed", "count": 20}, ...]
    source_distribution: List[dict] # e.g. [{"source": "AFDB", "count": 500}, ...]
    # Aggregate quick-stats — optional so old cached responses stay compatible
    avg_delta_length: Optional[float] = None
    avg_domains_lost: Optional[float] = None

class BrowseResponse(BaseModel):
    items: List[BrowseIsoformSummary]
    total_count: int
    stats: Optional[GlobalStats] = None

# Resolve forward references
IsoformResponse.model_rebuild()
IsoformWithStructuresResponse.model_rebuild()
IsoformDetailResponse.model_rebuild()
GeneDetailResponse.model_rebuild()


# ─── Drug–isoform layer ──────────────────────────────────────────────────────
class DrugResponse(BaseModel):
    id: int
    pubchem_cid: Optional[str] = None
    broad_id: Optional[str] = None
    name: str
    slug: str
    # Chemistry from PubChem (backfill_drug_chemistry.py). smiles drives the 2D skeletal
    # depiction on the drug page; the 3D view comes from the co-folded CIF, not from here.
    smiles: Optional[str] = None
    inchikey: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    iupac_name: Optional[str] = None
    chembl_id: Optional[str] = None

    class Config:
        from_attributes = True


class DrugComplexResponse(BaseModel):
    id: int
    job_name: str
    drug: DrugResponse
    gene_id: Optional[int] = None
    target_isoform_id: Optional[int] = None
    target_uniprot_id: Optional[str] = None
    seq_md5: Optional[str] = None
    iptm: Optional[float] = None
    holo_grade: Optional[str] = None
    ligand_plddt: Optional[float] = None
    pocket_plddt: Optional[float] = None
    interface_pae_min: Optional[float] = None
    n_contact_residues: Optional[int] = None
    contact_residues_json: Optional[str] = None
    has_structure: bool = False
    cif_url: Optional[str] = None       # /static path to the co-fold model.cif

    class Config:
        from_attributes = True


class DrugPocketFateResponse(BaseModel):
    id: int
    drug_complex_id: int
    isoform_id: int
    drug_name: Optional[str] = None          # hydrated from the complex's drug
    pubchem_cid: Optional[str] = None
    n_pocket_residues: Optional[int] = None
    n_retained: Optional[int] = None
    frac_retained: Optional[float] = None
    fate: Optional[str] = None
    mechanism_class: Optional[str] = None
    decidable: bool = False
    pocket_plddt: Optional[float] = None
    retained_residues_json: Optional[str] = None

    class Config:
        from_attributes = True


class DrugGeneSummary(BaseModel):
    """One row of the drug lens: a drugged gene with its aggregate fate mix."""
    gene_id: int
    symbol: str
    n_drugs: int
    n_complexes: int
    n_isoforms_scored: int
    n_subtractive: int          # isoforms with ≥1 lost/disrupted decidable pocket
    n_preserved: int
    subtractive_frac: float

    class Config:
        from_attributes = True


class DrugStatsResponse(BaseModel):
    n_drugs: int
    n_complexes: int
    n_drugged_genes: int
    fate_counts: dict            # {fate: count} over decidable pocket-fate rows
    n_subtractive_genes: int
    n_total_scored_genes: int


class GeneDrugsResponse(BaseModel):
    gene_id: int
    symbol: str
    complexes: List[DrugComplexResponse]
    fates: List[DrugPocketFateResponse]
