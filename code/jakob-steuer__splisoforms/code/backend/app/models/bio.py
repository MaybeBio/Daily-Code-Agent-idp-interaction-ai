from sqlalchemy import (Column, Integer, SmallInteger, BigInteger, String, Text, Float, Boolean,
                        ForeignKey, DateTime, Date, Enum, Table, Index, text)
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import relationship
import enum
from datetime import datetime
from app.database import Base

class SourceTypeEnum(enum.Enum):
    pacbio_longread = "pacbio_longread"
    ont_longread = "ont_longread"
    shortread = "shortread"
    afdb = "afdb"
    manual = "manual"
    other = "other"

class Source(Base):
    __tablename__ = "sources"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String(255), nullable=False, unique=True)
    source_type = Column(Enum(SourceTypeEnum), nullable=False)
    publication_doi = Column(String(255), nullable=True)
    cancer_type = Column(String(100), nullable=True)
    tissue_type = Column(String(100), nullable=True)
    sequencing_platform = Column(String(100), nullable=True)
    notes = Column(Text, nullable=True)
    loaded_at = Column(DateTime, default=datetime.utcnow)

    isoforms = relationship("Isoform", back_populates="source_dataset")

class GeneSet(Base):
    __tablename__ = "gene_sets"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String(255), nullable=False, unique=True, index=True)
    description = Column(Text, nullable=True)
    source_db = Column(String(100), nullable=True)  # e.g. KEGG, Hallmark, COSMIC

    genes = relationship("GeneSetMembership", back_populates="gene_set", cascade="all, delete-orphan")

class GeneSetMembership(Base):
    __tablename__ = "gene_set_memberships"

    id = Column(Integer, primary_key=True, index=True)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete="CASCADE"), nullable=False, index=True)
    gene_set_id = Column(Integer, ForeignKey("gene_sets.id", ondelete="CASCADE"), nullable=False, index=True)

    gene = relationship("Gene", back_populates="gene_set_memberships")
    gene_set = relationship("GeneSet", back_populates="genes")

class StructuralMetrics(Base):
    __tablename__ = "structural_metrics"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), unique=True)
    
    # Confidence metrics
    plddt_mean = Column(Float, nullable=True, index=True)
    plddt_min = Column(Float, nullable=True)
    plddt_median = Column(Float, nullable=True)
    
    # Alignment error metrics
    pae_mean = Column(Float, nullable=True)
    # DEPRECATED (audit-6): designed for floppy-domain identification but never
    # populated by current backfill scripts.  See PAEService.get_inter_domain_pae
    # if you need this metric live; the stored column is unused.
    pae_inter_domain_mean = Column(Float, nullable=True)
    pae_json_path = Column(String(511), nullable=True)    # Path to local matrix JSON
    
    # v4 Expanded Confidence series
    plddt_series = Column(ARRAY(Float), nullable=True)    # Per-residue pLDDT series for fast UI/Analytics
    pae_blocks = Column(Text, nullable=True)              # JSON summary of rigid structural blocks
    
    # Metadata
    # DEPRECATED (audit-6): the four fields below were declared for an earlier
    # design but are never populated by the current backfill scripts and are
    # not read by the API.  Left in the schema for forward-compatibility;
    # do not write new code that depends on them.
    afdb_entry_id = Column(String(100), nullable=True)    # legacy: e.g. AF-Q9NQ94-F1 (unused)
    afdb_version = Column(Integer, default=4)             # legacy (unused)
    is_mane_select = Column(Boolean, default=False)       # legacy — canonical-isoform truth lives on Gene now
    created_at = Column(DateTime, default=datetime.utcnow)

    isoform = relationship("Isoform", back_populates="structural_metrics")

class Gene(Base):
    __tablename__ = "genes"

    id = Column(Integer, primary_key=True, index=True)
    symbol = Column(String(100), unique=True, nullable=False)
    hgnc_symbol = Column(String(100), index=True)
    # Optimize lookups during batch ingestion
    uniprot_id = Column(String(100), unique=True, index=True)
    description = Column(Text)
    ensembl_gene_id = Column(String(30), nullable=True, index=True)  # ENSG ID for HPA lookups
    canonical_sequence = Column(Text, nullable=True)
    canonical_length = Column(Integer, nullable=True)
    # UniProt-sourced display fields; populated by backfill_gene_descriptions.py.
    protein_name = Column(Text, nullable=True)        # recommendedName.fullName
    function_summary = Column(Text, nullable=True)    # comments[type=FUNCTION].texts[0].value
    afdb_complex_checked = Column(Boolean, default=False)
    
    # Genomic Location
    chromosome = Column(String(50), nullable=True)
    strand = Column(String(5), nullable=True)
    genomic_start = Column(Integer, nullable=True)
    genomic_end = Column(Integer, nullable=True)
    
    # Canonical & AFDB tracking
    canonical_transcript_id = Column(Integer, ForeignKey("isoforms.id", name="fk_genes_canonical_transcript_id", use_alter=True), nullable=True)
    afdb_model_version = Column(String(50), nullable=True)
    afdb_fetched_at = Column(DateTime, nullable=True)
    afdb_not_found = Column(Boolean, default=False)
    afdb_monomer_cached = Column(Boolean, default=False)
    afdb_homodimer_cached = Column(Boolean, default=False)
    created_at = Column(DateTime, default=datetime.utcnow)

    isoforms = relationship("Isoform", back_populates="gene", foreign_keys="[Isoform.gene_id]")
    canonical_transcript = relationship("Isoform", foreign_keys=[canonical_transcript_id], post_update=True)
    cell_type_expression = relationship("CellTypeExpression", back_populates="gene", foreign_keys="[CellTypeExpression.gene_id]")
    gene_set_memberships = relationship("GeneSetMembership", back_populates="gene", cascade="all, delete-orphan")

class IsoformSourceMembership(Base):
    """
    Records which datasets an isoform was observed in, with per-dataset tissue flags.

    An isoform may appear in multiple datasets (e.g. identical sequence in ccRCC and
    breast cancer).  The primary source_id on Isoform records the dataset that first
    ingested the isoform.  This table captures all subsequent datasets and their
    tissue-level detection flags independently.

    The boolean flags on Isoform (in_normal, in_primary_tumor, in_metastasis,
    in_fibroblast) are the OR-merged view across all memberships and remain the fast
    lookup path.
    """
    __tablename__ = "isoform_source_memberships"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id  = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"),
                         nullable=False, index=True)
    source_id   = Column(Integer, ForeignKey("sources.id", ondelete="CASCADE"),
                         nullable=False, index=True)
    in_normal        = Column(Boolean, default=False, nullable=False)
    in_primary_tumor = Column(Boolean, default=False, nullable=False)
    in_metastasis    = Column(Boolean, default=False, nullable=False)
    in_fibroblast    = Column(Boolean, default=False, nullable=False)

    isoform = relationship("Isoform", back_populates="source_memberships")
    source  = relationship("Source")


class Isoform(Base):
    __tablename__ = "isoforms"

    id = Column(Integer, primary_key=True, index=True)
    gene_id = Column(Integer, ForeignKey("genes.id"), index=True)
    source_id = Column(Integer, ForeignKey("sources.id"), nullable=True, index=True)
    transcript_id = Column(String(255), unique=True, nullable=False)
    # Human-readable display name, derived — NOT an identifier. Never key on this.
    #   GENCODE  → Ensembl transcript_name          (MET-201)
    #   novel    → <closest_match_name>-NOV-<md5:6> (MET-201-NOV-a1b2c3)
    #   no match → NULL; the UI falls back to transcript_id
    # Recompute with backfill_transcript_names.py after any re-ingest or re-gene.
    transcript_name = Column(String(255), nullable=True, index=True)
    sequence = Column(Text, nullable=False)
    fasta_sequence = Column(Text)
    delta_length = Column(Integer)
    junction_coords = Column(Text)
    cds_sequence = Column(Text, nullable=True)  # CDS nucleotide sequence for future NMD 50-nt rule
    is_flagged = Column(Boolean, default=False, index=True) # Isoform of interest flag
    
    # Expression flags (Real data from MAS-seq samples)
    in_primary_tumor = Column(Boolean, default=False)
    in_normal = Column(Boolean, default=False)
    in_metastasis = Column(Boolean, default=False)
    in_fibroblast = Column(Boolean, default=False)  # cultured-fibroblast control (BC); kept distinct from in_normal
    is_canonical = Column(Boolean, default=False, index=True) # Distinguishes Reference (AFDB) from Novel (PacBio)
    exons_json = Column(Text, nullable=True) # Stores list of [start, end] pairs for visualization

    # ── GENCODE v47 reference-build columns (migration k7e4b3a2c9d1) ─────────
    # Filled by the Phase 3 rename pass and the future ingestion of fresh AF3
    # predictions.  See scripts/fetch_gencode_reference.py for the source.
    ensembl_transcript_id  = Column(String(20),  nullable=True, index=True)  # ENST00000...
    ensembl_protein_id     = Column(String(20),  nullable=True)              # ENSP00000...
    gencode_version        = Column(String(10),  nullable=True)              # 'v47'
    # Verbatim record of the pre-rename identifier (PB.*, AFDB_*, ...) so any
    # historical analysis keyed on the legacy name remains queryable.
    legacy_transcript_id   = Column(String(100), nullable=True, index=True)
    # For novel rows, the closest GENCODE-anchored isoform by sequence identity.
    canonical_match_id     = Column(
        Integer,
        ForeignKey("isoforms.id", name="fk_isoforms_canonical_match", use_alter=True),
        nullable=True,
    )
    # md5 of the protein sequence — drives dedup + cross-version matching.
    protein_md5            = Column(String(32),  nullable=True, index=True)
    # True for transcripts marked MANE_Select by GENCODE / NCBI MANE.
    is_mane_select         = Column(Boolean,     nullable=True, index=True)
    # For novel rows, the alignment identity to canonical_match (∈ [0, 1]).
    gencode_match_identity = Column(Float,       nullable=True)
    # Which model produced the stored structure (afdb_legacy, af3_local_v1,
    # af3_local_<date>, ...).  Used by downstream consumers to know whether
    # the structure file is on the same model as canonical-relative analyses.
    structure_provenance   = Column(String(50),  nullable=True, index=True)
    # Populated by Phase 4 (PeptideAtlas / TransCODE) — does this isoform have
    # at least one distinguishing MS-detected peptide?
    ms_evidence            = Column(Boolean,     nullable=False, default=False, server_default="false")

    gene = relationship("Gene", back_populates="isoforms", foreign_keys=[gene_id])
    source_dataset = relationship("Source", back_populates="isoforms")
    source_memberships = relationship("IsoformSourceMembership", back_populates="isoform",
                                      cascade="all, delete-orphan")

    @property
    def source(self):
        """Alias for source_dataset to match schema field name."""
        return self.source_dataset
    structures = relationship("Structure", back_populates="isoform")
    clinical_frequencies = relationship("ClinicalFrequency", back_populates="isoform")
    ppi_networks = relationship("PPINetwork", back_populates="isoform", foreign_keys="[PPINetwork.isoform_id_1]")
    metadata_record = relationship("TranscriptMetadata", back_populates="isoform", uselist=False)
    pfam_domains = relationship("PfamDomain", back_populates="isoform", foreign_keys="[PfamDomain.isoform_id]")
    ptm_accessibility = relationship("PtmAccessibility", back_populates="isoform", foreign_keys="[PtmAccessibility.isoform_id]")
    nmd_prediction = relationship("NmdPrediction", back_populates="isoform", uselist=False, foreign_keys="[NmdPrediction.isoform_id]")
    structural_impacts = relationship("ExonStructuralImpact", back_populates="isoform", foreign_keys="[ExonStructuralImpact.isoform_id]")
    neoantigen_predictions = relationship("NeoantigenPrediction", back_populates="isoform", foreign_keys="[NeoantigenPrediction.isoform_id]")
    homodimer_analysis = relationship("HomodimerAnalysis", back_populates="isoform", uselist=False, foreign_keys="[HomodimerAnalysis.isoform_id]")
    heterodimer_analyses = relationship("HeterodimerAnalysis", back_populates="isoform", foreign_keys="[HeterodimerAnalysis.isoform_id]")
    string_interactions = relationship("StringInteraction", back_populates="isoform", foreign_keys="[StringInteraction.isoform_id]")
    ted_domains = relationship("TedDomain", back_populates="isoform", foreign_keys="[TedDomain.isoform_id]")
    junction_events = relationship("IsoformJunction", back_populates="isoform", foreign_keys="[IsoformJunction.isoform_id]")
    structural_metrics = relationship("StructuralMetrics", back_populates="isoform", uselist=False, foreign_keys="[StructuralMetrics.isoform_id]")
    disorder_score = relationship("DisorderScore", back_populates="isoform", uselist=False, foreign_keys="[DisorderScore.isoform_id]")
    expression_data = relationship("IsoformExpression", back_populates="isoform", foreign_keys="[IsoformExpression.isoform_id]")

class TranscriptMetadata(Base):
    __tablename__ = "transcript_metadata"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id"), unique=True)
    transcript_id = Column(String(255), index=True)

    # ---- Source-of-this-row tag (n7c2b5d9a4e8) ----
    # One of: 'sqanti3' | 'gencode' | 'afdb_reference' | 'custom' | None.
    # Drives which downstream fields are populated and how the frontend
    # Transcript Metadata panel renders this isoform.
    metadata_kind = Column(String(32), nullable=True, index=True)

    # ---- Universal: every annotation pipeline can fill these ----
    chrom = Column(String(100))
    strand = Column(String(20))
    length = Column(Integer)
    biotype = Column(String(64), nullable=True, index=True)        # protein_coding, lncRNA, ...
    orf_type = Column(String(255))
    protein_length_aa = Column(Float)
    delta_length_aa = Column(Float)
    sample_source = Column(String(255))
    associated_gene = Column(String(255))
    canonical_gene_name = Column(String(255))
    canonical_protein_length = Column(Float)

    # ---- ccRCC folded-ORF QC (integrity + transparency) ----
    # The AF3 structure depicts the author-curated folded ORF (af3_targets.fa =
    # forAF). For a few % of ccRCC transcripts that ORF's protein identity diverges
    # from the transcript's genomic-locus gene — e.g. a GAPDH-locus transcript whose
    # folded ORF best-matches UniProt A4UCT1, not GAPDH. We keep the genomic gene as
    # the objective anchor but FLAG the divergence and record the folded protein's
    # best UniProt match, so a structure is never silently labelled as a gene it does
    # not depict (see backfill_ccrcc_orf_flag.py for the provenance rules).
    folded_orf_mismatch = Column(Boolean, default=False, index=True)
    folded_best_ref = Column(String(100), nullable=True)          # e.g. tr|A4UCT1
    folded_best_ref_identity = Column(Float, nullable=True)       # % identity of the fold to that ref

    # ---- SQANTI3-specific (metadata_kind == 'sqanti3') ----
    # Historically these were dumping grounds for non-SQANTI values too;
    # the ingestion scripts have been corrected. Read with that caveat
    # against pre-2026-06 data.
    structural_category = Column(String(255))   # FSM / ISM / NIC / NNC / antisense / ...
    subcategory = Column(String(255))           # mono-exon / alternative_3end / ...
    all_canonical = Column(String(100))         # SQANTI3 'all_canonical' splice junction flag

    # ---- GENCODE/Ensembl-specific (metadata_kind == 'gencode') ----
    appris_tag = Column(String(64), nullable=True)                       # principal:1 / alternative:2 / ...
    tsl = Column(SmallInteger, nullable=True)                            # 1..5 transcript support level
    gencode_tags = Column(ARRAY(String(64)), nullable=True)              # basic, CCDS, MANE_Select, ...

    # ---- Cross-database identifiers (universal but optional) ----
    uniprot_id = Column(String(100))
    ensp_id = Column(String(100), nullable=True)
    enst_id = Column(String(100))
    refseq_id = Column(String(100))
    hgnc_id = Column(String(100))
    is_mane_select = Column(Boolean, default=False)
    pdb_ids = Column(ARRAY(String))

    isoform = relationship("Isoform", back_populates="metadata_record")

class Structure(Base):
    __tablename__ = "structures"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id"), index=True)
    storage_uri = Column(String(255), nullable=True) # Optional back-compat
    af3_cif_url = Column(String(500))
    af3_pae_url = Column(String(500))
    af3_confidence = Column(Float)
    af3_plddt_url = Column(String(500))
    af3_full_cif_path = Column(String(500))
    created_at = Column(DateTime, default=datetime.utcnow)

    isoform = relationship("Isoform", back_populates="structures")

class ClinicalFrequency(Base):
    __tablename__ = "clinical_frequencies"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id"))
    cancer_type = Column(String(100))
    frequency = Column(Float)

    isoform = relationship("Isoform", back_populates="clinical_frequencies")

class PPINetwork(Base):
    __tablename__ = "ppi_networks"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id_1 = Column(Integer, ForeignKey("isoforms.id"))
    target_uniprot_id = Column(String(50))
    interaction_score = Column(Float)
    interaction_retained = Column(Boolean)

    isoform = relationship("Isoform", back_populates="ppi_networks", foreign_keys=[isoform_id_1])

class PfamDomain(Base):
    __tablename__ = "pfam_domains"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, index=True)
    pfam_accession = Column(String(20), nullable=False, index=True)
    pfam_name = Column(String(100), nullable=False)
    pfam_description = Column(Text)
    pfam_clan = Column(String(20))
    seq_start = Column(Integer, nullable=False)
    seq_end = Column(Integer, nullable=False)
    hmm_start = Column(Integer, nullable=False)
    hmm_end = Column(Integer, nullable=False)
    hmm_length = Column(Integer, nullable=False)
    domain_coverage = Column(Float, nullable=False)
    e_value = Column(Float, nullable=False)
    bit_score = Column(Float, nullable=False)
    domain_status = Column(String(20), nullable=False)  # present, gained, lost, truncated
    # Curated functional-activity class derived from Pfam clan + accession
    # (see app/analysis/functional_classes.py).  Values:
    #   catalytic_kinase / catalytic_phosphatase / catalytic_protease /
    #   catalytic_other / ub_ligase / dna_binding / rna_binding /
    #   signaling_adapter / oligomerization / other
    functional_class = Column(String(40), nullable=True, index=True)
    # Where the splicing event disrupts the domain, derived from
    # domain_status + HMM coordinates.  Values:
    #   intact / truncated_n / truncated_c / truncated_both /
    #   internal_disruption / lost / gained / truncated_unknown / unknown
    alteration_type = Column(String(30), nullable=True, index=True)
    # Legacy kinase flag — kept for backward compat; new code should query
    # functional_class == 'catalytic_kinase' instead.
    is_kinase = Column(Boolean, default=False)
    canonical_isoform_id = Column(Integer, ForeignKey("isoforms.id"), nullable=True)

    isoform = relationship("Isoform", back_populates="pfam_domains", foreign_keys=[isoform_id])

class TedDomain(Base):
    __tablename__ = "ted_domains"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, index=True)
    ted_id = Column(String(50), nullable=False, index=True)
    domain_type = Column(String(20), nullable=False, index=True) # e.g. TED, TEDLH
    seq_start = Column(Integer, nullable=False)
    seq_end = Column(Integer, nullable=False)
    score = Column(Float, nullable=True) # Confidence or similarity score
    chopping = Column(String(255), nullable=True)
    plddt = Column(Float, nullable=True)
    cath_label = Column(String(50), nullable=True)
    
    isoform = relationship("Isoform", back_populates="ted_domains", foreign_keys=[isoform_id])

class PtmAccessibility(Base):
    __tablename__ = "ptm_accessibility"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, index=True)
    uniprot_id = Column(String(50), nullable=False)
    residue_position = Column(Integer, nullable=False)       # 1-based, canonical numbering
    novel_residue_position = Column(Integer, nullable=True)  # 1-based, novel numbering
    residue_type = Column(String(5), nullable=False)         # S, T, Y, K, etc.
    ptm_type = Column(String(50), nullable=False)            # phosphorylation, ubiquitination, acetylation
    source_db = Column(String(50), nullable=False)           # PhosphoSitePlus or UniProt
    canonical_sasa = Column(Float, nullable=True)            # Å²
    novel_sasa = Column(Float, nullable=True)                # Å²
    delta_sasa = Column(Float)                                # absolute Δ in Å²
    # Δ / canonical_sasa — relative shift, more biologically meaningful
    # than the absolute Δ because baseline SASA varies 10× across residues.
    relative_delta_sasa = Column(Float, nullable=True)
    # Continuous magnitude on [0, 1] of the relative shift (capped at 1.0).
    change_strength = Column(Float, nullable=True)
    # min(canonical_plddt, novel_plddt) at this site.  pLDDT < 70 → SASA
    # is treated as undetermined (Tunyasuvunakool et al., Nature 2021).
    min_plddt = Column(Float, nullable=True)
    # intact (was "unchanged" / "minor"), exposed, buried, lost, undetermined.
    # The legacy "unchanged" value is retained on pre-redesign rows for
    # backward compat.
    accessibility_change = Column(String(20))
    structural_context = Column(String(50))  # "intra-domain" | "inter-domain" | "disrupted"
    
    canonical_isoform_id = Column(Integer, ForeignKey("isoforms.id"), nullable=True)

    isoform = relationship("Isoform", back_populates="ptm_accessibility", foreign_keys=[isoform_id])

class NmdPrediction(Base):
    __tablename__ = "nmd_predictions"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, unique=True, index=True)
    # NMD verdict.  After the d7f3e2a8b4c5 redesign:
    #   nmd_target              50-nt rule, PTC well upstream of last EJC
    #   nmd_escape              50-nt rule, escape (≤50 nt, 5′-proximal, or in last exon)
    #   non_translatable        no complete ORF detected
    #   predicted_translatable  heuristic only — protein looks intact, no PTC signal
    #   unknown                 insufficient data to call (incl. smORF candidates)
    nmd_status = Column(String(30), nullable=False)
    has_ptc = Column(Boolean, nullable=False)
    truncation_ratio = Column(Float, nullable=True)
    orf_type = Column(String(50), nullable=True)
    ptc_position = Column(Integer, nullable=True)
    last_ejc_distance = Column(Integer, nullable=True)
    evidence = Column(Text, nullable=False)
    canonical_isoform_id = Column(Integer, ForeignKey("isoforms.id"), nullable=True)
    # Orthogonal protein-length signal (intact / minor / severe / non_translatable).
    # Surfaces protein truncation WITHOUT conflating it with NMD evidence.
    truncation_severity = Column(String(30), nullable=True, index=True)
    # Confidence on nmd_target calls; linear ramp from 0 (distance=50) to 1
    # (distance ≥ 200) using last_ejc_distance.  NULL for non-target rows.
    nmd_confidence = Column(Float, nullable=True)
    # TransCODE peptidein candidate: short isoform with a short / absent
    # canonical, where the 50-nt rule does not meaningfully apply.
    is_candidate_smorf = Column(Boolean, nullable=False, default=False, index=True)
    # Candidate upstream ORF: short isoform with a long canonical; uORFs are
    # NMD-resistant via ribosome reinitiation.
    is_candidate_uorf = Column(Boolean, nullable=False, default=False)

    isoform = relationship("Isoform", back_populates="nmd_prediction", foreign_keys=[isoform_id])

class ExonStructuralImpact(Base):
    __tablename__ = "exon_structural_impacts"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, index=True)
    residue_number = Column(Integer, nullable=False)           # 1-based protein residue
    dssp_assignment = Column(String(1), nullable=False)        # H (helix), E (sheet), C (coil)
    phi_angle = Column(Float, nullable=True)
    psi_angle = Column(Float, nullable=True)
    residue_plddt = Column(Float, nullable=True)               # per-residue pLDDT from B-factor
    exon_number = Column(Integer, nullable=True)               # for future exon boundary mapping
    canonical_dssp = Column(String(1), nullable=True)          # canonical residue's DSSP at aligned pos
    structural_change = Column(String(30), nullable=True)      # same, helix_to_coil, sheet_to_coil, etc.
    canonical_isoform_id = Column(Integer, ForeignKey("isoforms.id"), nullable=True)

    isoform = relationship("Isoform", back_populates="structural_impacts", foreign_keys=[isoform_id])

class CellTypeExpression(Base):
    __tablename__ = "cell_type_expression"

    id = Column(Integer, primary_key=True, index=True)
    gene_id = Column(Integer, ForeignKey("genes.id", ondelete="CASCADE"), nullable=False, index=True)
    cell_type = Column(String(100), nullable=False)
    expression_ncpm = Column(Float, nullable=False)            # normalized counts per million from HPA
    source = Column(String(50), nullable=False, default="HPA")  # HPA_v25, etc.
    specificity_category = Column(String(50), nullable=True)    # enriched, group_enriched, enhanced, low

    gene = relationship("Gene", back_populates="cell_type_expression", foreign_keys=[gene_id])

class IsoformExpression(Base):
    """Per-sample FL read counts and derived expression metrics from SQANTI3 CustomFiltered files."""
    __tablename__ = "isoform_expression"

    id           = Column(Integer, primary_key=True, index=True)
    isoform_id   = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, index=True)
    sample       = Column(String(50),  nullable=False)   # "Normal" | "ccRCC2" | ...
    condition    = Column(String(20),  nullable=False)   # "normal" | "tumor"
    fl_count     = Column(Integer,     nullable=False)   # raw Full-Length read count
    cpm          = Column(Float,       nullable=True)    # fl_count / total_sample_fl * 1e6
    isoform_frac = Column(Float,       nullable=True)    # ψ: fl_count / gene_total_fl
    n_cells      = Column(Integer,     nullable=True)    # cells with ≥1 UMI (sparse matrix)
    cell_frac    = Column(Float,       nullable=True)    # n_cells / total_cells_in_sample

    isoform = relationship("Isoform", back_populates="expression_data", foreign_keys=[isoform_id])


class NeoantigenPrediction(Base):
    __tablename__ = "neoantigen_predictions"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, index=True)
    peptide_sequence = Column(String(30), nullable=False)
    peptide_start = Column(Integer, nullable=False)             # 1-based position in novel protein
    peptide_end = Column(Integer, nullable=False)
    hla_allele = Column(String(30), nullable=False)             # e.g. HLA-A*02:01
    binding_affinity_nm = Column(Float, nullable=False)         # IC50 in nM
    binding_category = Column(String(20), nullable=False)       # strong_binder, weak_binder, non_binder
    percentile_rank = Column(Float, nullable=True)
    junction_type = Column(String(30), nullable=False)          # novel_insertion, novel_substitution, frameshift
    
    # MHCflurry 2.0 Presentation Model scores
    presentation_score = Column(Float, nullable=True)           # Combined presentation probability (0–1); primary ranking metric
    processing_score   = Column(Float, nullable=True)           # Antigen-processing / proteasomal cleavage score (0–1)

    # PTM-Sensitive Neoantigen Analysis (v4 Logic)
    ptm_impact_status = Column(String(50), nullable=True)       # Direct PTM-Neoantigen, Processing-Blocked, Splice-Enabled, etc.
    avg_plddt = Column(Float, nullable=True)                    # Average pLDDT across epitope window
    min_plddt = Column(Float, nullable=True)                    # Minimum pLDDT across epitope window
    is_ptm_sensitive = Column(Boolean, default=False)

    # Multi-tool ensemble scores (v5 — NetMHCpan-4.2 + BigMHC)
    netmhcpan_affinity_nm = Column(Float, nullable=True)        # NetMHCpan-4.2 BA IC50 (nM)
    netmhcpan_rank_el     = Column(Float, nullable=True)        # NetMHCpan-4.2 EL %rank (primary; lower = better)
    netmhcpan_rank_ba     = Column(Float, nullable=True)        # NetMHCpan-4.2 BA %rank (rescue signal for atypical binders)
    bigmhc_el_score       = Column(Float, nullable=True)        # BigMHC EL presentation score (0–1)
    bigmhc_im_score       = Column(Float, nullable=True)        # BigMHC IM immunogenicity score (0–1; T-cell reactivity)
    consensus_rank        = Column(Float, nullable=True)        # EXPERIMENTAL median %rank across tools — not the primary ranking metric (incommensurable scales; see audit)
    # Audit-3 additions (Option A)
    # Count of tools agreeing the peptide is a strong binder (0-3):
    #   MHCflurry  %rank < 0.5
    #   NetMHCpan  %rank_el < 0.5
    #   BigMHC     EL score > 0.5
    tools_agreeing = Column(Integer, nullable=True, index=True)
    # Derived bucket: 3 → high, 2 → medium, ≤1 → low.
    agreement_tier = Column(String(10), nullable=True, index=True)
    # When ptm_impact_status starts with "Direct PTM-Neoantigen", True iff the
    # canonical modifiable residue is physically absent in the novel sequence
    # (residue type changed by the splice event).
    ptm_canonical_residue_lost = Column(Boolean, nullable=False, default=False)
    is_ba_rescue          = Column(Boolean, default=False)      # True when NetMHCpan-BA predicts binding but EL tools miss (cysteine/atypical)

    # ── ranking v2 (see app/analysis/neoantigen_tiers.py) ──
    # presentation_tier: STRONG / MODERATE / WEAK / NONE — multi-tool MHC-presentation consensus
    presentation_tier = Column(String(10), nullable=True, index=True)
    # priority_tier: P1–P4 conjunction of presentation × immunogenicity (the headline tier)
    priority_tier     = Column(String(4), nullable=True, index=True)
    is_immunogenic    = Column(Boolean, default=False, index=True)  # BigMHC-IM ≥ threshold
    is_expressed      = Column(Boolean, default=False)              # isoform detected in expression data (decorating flag)
    # is_representative: the single best (peptide_sequence, isoform_id) row among presented
    # predictions — materialises the atlas dedup so /neoantigens/top needs no DISTINCT ON.
    # Set by backfill_priority_tier.py; TRUE implies priority_tier IS NOT NULL.
    is_representative  = Column(Boolean, default=False, index=True)

    # Proteome-wide novelty of the PEPTIDE (not just vs the gene canonical), set by
    # backfill_neoantigen_uniqueness.py against every GENCODE v47 reference protein:
    #   'unique'      — occurs in no reference protein (true neoepitope candidate)
    #   'self_gene'   — occurs in another reference isoform of the same gene
    #   'other_gene'  — occurs in a reference protein of another gene (self-peptide)
    # is_proteome_unique = (proteome_uniqueness == 'unique'); NULL = not yet computed.
    proteome_uniqueness = Column(String(20), nullable=True, index=True)
    is_proteome_unique  = Column(Boolean, nullable=True, index=True)

    canonical_isoform_id = Column(Integer, ForeignKey("isoforms.id"), nullable=True)

    # Partial indexes over the representative (atlas) subset that /neoantigens/top
    # relies on (migrations c3f8a1b0d2e5 + r1b6e3d8c5a2). They are declared HERE as
    # well because `reset_db.py` builds the schema with metadata.create_all() and then
    # `alembic stamp head`: anything that exists only as migration SQL is silently lost
    # at every rebuild. The 2026-10 v2 rebuild lost all three and the atlas went from
    # ~0.85 s to 2–4 s per request.
    __table_args__ = (
        # default atlas sort (priority_tier, consensus_rank NULLS LAST) as an index scan
        Index("ix_neo_representative", "priority_tier", "consensus_rank",
              postgresql_where=text("is_representative")),
        # distinct isoforms → unique_genes; also the cheapest COUNT(*) index-only scan
        Index("ix_neo_repr_isoform", "isoform_id",
              postgresql_where=text("is_representative")),
        # covering index for the atlas stats block: COUNT(DISTINCT peptide_sequence)
        # reads it presorted, and the FILTER counts + every filter on an INCLUDEd
        # column run as index-only scans instead of fetching ~2.6M heap tuples.
        Index("ix_neo_repr_stats", "peptide_sequence",
              postgresql_include=["priority_tier", "presentation_tier", "is_immunogenic",
                                  "is_expressed", "is_ba_rescue", "is_proteome_unique",
                                  "netmhcpan_rank_el", "ptm_impact_status"],
              postgresql_where=text("is_representative")),
    )

    isoform = relationship("Isoform", back_populates="neoantigen_predictions", foreign_keys=[isoform_id])
class HomodimerAnalysis(Base):
    __tablename__ = "homodimer_analysis"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, unique=True)
    afdb_accession = Column(String(100))
    total_interface_residues = Column(Integer)
    count_absent = Column(Integer)
    count_sasa_changed = Column(Integer)
    disruption_score = Column(Float)
    # Audit-9 (migration j9d0c4f6b8a1): TRUE when ≥ 50% of canonical interface
    # residues are absent in the novel — i.e. the interface region itself was
    # spliced out.  This is a different biology from "interface
    # conformationally disrupted" (which is what a high disruption_score on a
    # retained interface measures); UI consumers should de-prioritise or
    # separately render the deletion cohort.
    interface_disrupted_by_deletion = Column(Boolean, default=False, index=True, nullable=False)
    interface_residues_json = Column(Text)  # JSON list of residue info

    isoform = relationship("Isoform", back_populates="homodimer_analysis")

class HeterodimerAnalysis(Base):
    __tablename__ = "heterodimer_analysis"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False)
    partner_uniprot_id = Column(String(100))
    partner_name = Column(String(255))
    afdb_accession = Column(String(100))
    total_interface_residues = Column(Integer)
    count_absent = Column(Integer)
    count_sasa_changed = Column(Integer)
    disruption_score = Column(Float)
    interface_disrupted_by_deletion = Column(Boolean, default=False, index=True, nullable=False)
    interface_residues_json = Column(Text)  # JSON list of residue info

    isoform = relationship("Isoform", back_populates="heterodimer_analyses")

class StringInteraction(Base):
    __tablename__ = "string_interactions"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False)
    partner_gene_symbol = Column(String(100))
    # Renamed from partner_uniprot_id in migration i6f1a3d8e5b2 — the value
    # was always an Ensembl protein ID (ENSP*) extracted from STRING's
    # stringId_B / stringId_A field, never a UniProt accession.
    partner_ensp_id = Column(String(100))
    string_score = Column(Integer)
    evidence_json = Column(Text)  # JSON map of per-channel scores
    # Audit-5 caveat: this is set when the isoform has lost ANY Pfam domain,
    # not specifically the domain that mediates the partner interaction.
    # Treat as a gene-wide flag, not a partner-specific signal.
    domain_loss_possible = Column(Boolean, default=False)

    isoform = relationship("Isoform", back_populates="string_interactions")

class IsoformJunction(Base):
    __tablename__ = "isoform_junctions"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, index=True)
    chrom = Column(String(50))
    strand = Column(String(5))
    junction_type = Column(String(100))  # novel_exons, skipped_exons, alternative_splice_sites, etc.
    subcategory = Column(String(255))
    exon_count_novel = Column(Integer)
    exon_count_canonical = Column(Integer)
    transcript_pos = Column(Integer)
    is_novel = Column(Boolean, default=False)

    isoform = relationship("Isoform", back_populates="junction_events")

class DisorderScore(Base):
    __tablename__ = "disorder_scores"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), unique=True, nullable=False, index=True)
    
    # Store per-residue scores as ARRAY(Float) for efficiency
    disorder_scores = Column(ARRAY(Float), nullable=False)
    anchor_scores = Column(ARRAY(Float), nullable=False)
    
    created_at = Column(DateTime, default=datetime.utcnow)

    isoform = relationship("Isoform", back_populates="disorder_score")

class IsoformDiscoveryIndex(Base):
    __tablename__ = "isoform_discovery_index"

    id = Column(Integer, primary_key=True, index=True)
    isoform_id = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"), unique=True, nullable=False, index=True)
    
    # Discovery Metrics (Pre-calculated for performance)
    plddt_mean = Column(Float, nullable=True)
    delta_length = Column(Integer, nullable=True)
    num_domains_lost = Column(Integer, default=0)
    num_domains_gained = Column(Integer, default=0)
    has_neoantigen = Column(Boolean, default=False)
    # Materialised neoepitope counts (representative rows) so the Impact Explorer
    # can show + sort Total / P1 without aggregating neoantigen_predictions at
    # request time. Populated by backfill_discovery_index.py.
    neoantigen_count = Column(Integer, default=0, index=True)     # total representative neoepitopes
    neoantigen_p1_count = Column(Integer, default=0, index=True)  # P1 (high-confidence) neoepitopes
    has_ptm_exposure_change = Column(Boolean, default=False)
    has_ptm_buried_change = Column(Boolean, default=False)
    has_ptm_loss = Column(Boolean, default=False)
    nmd_status = Column(String(50), nullable=True) # nmd_target, nmd_suspect, etc.
    impact_score = Column(Float, default=0.0)
    
    # onupdate is load-bearing: backfill_discovery_index.py UPDATES existing rows rather
    # than re-inserting, so a bare `default=` freezes this at first-insert time and makes a
    # stale index look freshly built. That is exactly how this table sat stale across the
    # 2026-07 re-attribution + neoantigen load without anyone noticing.
    indexed_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    isoform = relationship("Isoform", backref="discovery_index")


class PredictionJob(Base):
    """
    Per-target AF3 prediction lifecycle (migration k7e4b3a2c9d1).

    One row per protein structure prediction we ask Leonardo to compute as
    part of the GENCODE v47 reference build.  Multiple isoform rows can
    share a single PredictionJob via `protein_md5` when their sequences
    are identical.

    States advance:
        queued      → submitted   (script puts the FASTA on the HPC)
        submitted   → predicting  (Slurm job started)
        predicting  → predicted   (CIF + confidence files exist on the HPC)
        predicted   → ingested    (results pulled back + Structure row written)
        anything    → failed      (with error_message)
        anything    → skipped     (sequence excluded by quality filters)
    """
    __tablename__ = "prediction_jobs"

    id = Column(Integer, primary_key=True, index=True)
    # The display name the eventual Isoform row will adopt once ingested.
    # For GENCODE-known transcripts this is the Ensembl transcript_name
    # (e.g. ANK3-219); for novel PacBio rows it's <closest>-NOV-<6-char-md5>.
    target_id        = Column(String(100), nullable=False, unique=True, index=True)

    protein_md5      = Column(String(32),  nullable=False, index=True)
    protein_sequence = Column(Text,        nullable=False)
    protein_length   = Column(Integer,     nullable=False)

    # 1 = MANE Select, 2 = alternative protein-coding GENCODE, 3 = novel PacBio.
    priority_tier    = Column(SmallInteger, nullable=False, default=2, index=True)
    # Free-text tag for what triggered this prediction
    # (e.g. 'fresh_gencode_v47', 'rename_AFDB_replace', 'novel_pacbio_fresh').
    queued_for       = Column(String(50),  nullable=True)

    # Workflow state — see class docstring for the allowed transitions.
    status           = Column(String(20),  nullable=False, default="queued", index=True)

    # HPC bookkeeping
    hpc_job_id       = Column(String(100), nullable=True)
    submitted_at     = Column(DateTime,    nullable=True)
    completed_at     = Column(DateTime,    nullable=True)
    output_path      = Column(String(511), nullable=True)
    error_message    = Column(Text,        nullable=True)

    gencode_version  = Column(String(10),  nullable=True)
    created_at       = Column(DateTime,    default=datetime.utcnow, nullable=False)
    updated_at       = Column(DateTime,    default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)


class ApiKey(Base):
    """
    Credential for the public REST API (/api/v1).

    Only the SHA-256 hash of the raw key is stored; the plaintext is shown once
    at issuance (see scripts/mint_api_key.py) and never persisted.  Per-key
    minute/day counters back the application-level rate limit; they reset
    lazily when the window bucket rolls over.
    """
    __tablename__ = "api_keys"

    id            = Column(Integer, primary_key=True, index=True)
    key_hash      = Column(String(64), nullable=False, unique=True, index=True)
    key_prefix    = Column(String(20), nullable=False, index=True)
    label         = Column(String(255), nullable=True)
    owner_email   = Column(String(255), nullable=True)
    tier          = Column(String(20), nullable=False, default="standard")

    # A key is the higher-volume tier, so its quota must exceed the ~10 r/s anonymous
    # per-IP edge limit; keyed requests also bypass that per-IP zone at nginx.
    rate_per_min  = Column(Integer, nullable=False, default=600)
    rate_per_day  = Column(Integer, nullable=False, default=100000)

    minute_window = Column(DateTime, nullable=True)
    minute_count  = Column(Integer, nullable=False, default=0)
    day_window    = Column(Date, nullable=True)
    day_count     = Column(Integer, nullable=False, default=0)

    request_count = Column(BigInteger, nullable=False, default=0)
    revoked       = Column(Boolean, nullable=False, default=False, index=True)
    created_at    = Column(DateTime, default=datetime.utcnow, nullable=False)
    last_used_at  = Column(DateTime, nullable=True)


# ─── Drug–isoform layer (co-folded drug–protein complexes) ────────────
# Drugs overlay existing genes/isoforms; they are NOT a new isoform Source. Each
# DrugComplex is one AF3 co-fold, anchored to the DB by protein_md5 (100% join;
# ENST is version-unstable — the drug set ≈ Ensembl 115 vs DB GENCODE v47 = Ensembl 113).
class Drug(Base):
    __tablename__ = "drugs"

    id = Column(Integer, primary_key=True, index=True)
    pubchem_cid = Column(String(20), unique=True, index=True, nullable=True)  # primary external key
    broad_id    = Column(String(50), nullable=True)   # DepMap/Broad BRD id
    name        = Column(String(255), nullable=False)
    slug        = Column(String(255), index=True, nullable=False)
    # Enrichment columns (future ChEMBL mapping — nullable for now)
    chembl_id   = Column(String(20), nullable=True)
    moa         = Column(Text, nullable=True)
    # Chemistry, filled from PubChem by backfill_drug_chemistry.py. The drug-dataset manifest
    # carries only cid/broad_id/name, so none of this arrives at ingest.
    smiles            = Column(Text, nullable=True)
    inchikey          = Column(String(27), nullable=True)
    molecular_formula = Column(String(128), nullable=True)
    molecular_weight  = Column(Float, nullable=True)   # g/mol
    iupac_name        = Column(Text, nullable=True)
    created_at  = Column(DateTime, default=datetime.utcnow)

    complexes = relationship("DrugComplex", back_populates="drug",
                             cascade="all, delete-orphan")


class DrugComplex(Base):
    __tablename__ = "drug_complexes"

    id = Column(Integer, primary_key=True, index=True)
    drug_id           = Column(Integer, ForeignKey("drugs.id", ondelete="CASCADE"),
                               nullable=False, index=True)
    gene_id           = Column(Integer, ForeignKey("genes.id"), nullable=True, index=True)
    target_isoform_id = Column(Integer, ForeignKey("isoforms.id"), nullable=True, index=True)
    target_uniprot_id = Column(String(100), nullable=True, index=True)
    target_enst       = Column(String(30), nullable=True)   # as-provided; may be a renamed ENST
    seq_md5           = Column(String(32), index=True, nullable=True)
    # job_name = the co-fold's structure-dir slug (<drug>--<uniprot>--<hash>); the
    # natural key for idempotent ingest.
    job_name          = Column(String(255), unique=True, nullable=False)
    structure_dir     = Column(String(500), nullable=True)  # source path in the drug dataset until synced
    target_length     = Column(Integer, nullable=True)
    has_structure     = Column(Boolean, default=False, index=True, nullable=False)

    # ── Structure-derived metrics (filled by backfill_drug_contacts.py) ──────
    iptm               = Column(Float, nullable=True)    # top-level interface pTM
    ptm                = Column(Float, nullable=True)
    interface_pae_min  = Column(Float, nullable=True)    # min PAE protein↔ligand
    ligand_plddt       = Column(Float, nullable=True)    # mean pLDDT over ligand atoms
    pocket_plddt       = Column(Float, nullable=True)    # mean pLDDT over contact residues
    holo_grade         = Column(String(1), nullable=True, index=True)  # A/B/C from iptm
    n_contact_residues = Column(Integer, nullable=True)
    contact_residues_json = Column(Text, nullable=True)  # JSON list of contact residue numbers
    created_at         = Column(DateTime, default=datetime.utcnow)

    drug    = relationship("Drug", back_populates="complexes")
    gene    = relationship("Gene")
    isoform = relationship("Isoform", foreign_keys=[target_isoform_id])
    pocket_fates = relationship("DrugPocketFate", back_populates="complex",
                                cascade="all, delete-orphan")


class DrugPocketFate(Base):
    """Our own Axis-1: the drugged transcript's contact pocket carried onto each
    sibling isoform of the same gene, scored by residue retention. Anchored on the
    co-folded transcript (not a gene-level 'canonical'), so it is unaffected by the
    is_canonical vs MANE ambiguity."""
    __tablename__ = "drug_pocket_fate"

    id = Column(Integer, primary_key=True, index=True)
    drug_complex_id   = Column(Integer, ForeignKey("drug_complexes.id", ondelete="CASCADE"),
                               nullable=False, index=True)
    isoform_id        = Column(Integer, ForeignKey("isoforms.id", ondelete="CASCADE"),
                               nullable=False, index=True)
    n_pocket_residues = Column(Integer, nullable=True)
    n_retained        = Column(Integer, nullable=True)
    frac_retained     = Column(Float, nullable=True)
    fate              = Column(String(20), nullable=True, index=True)   # lost/disrupted/preserved/no_call
    mechanism_class   = Column(String(30), nullable=True)              # domain_deletion/pocket_disruption/preserved/neo_pocket_gain
    decidable         = Column(Boolean, default=False, nullable=False)
    pocket_plddt      = Column(Float, nullable=True)
    retained_residues_json = Column(Text, nullable=True)
    created_at        = Column(DateTime, default=datetime.utcnow)

    complex = relationship("DrugComplex", back_populates="pocket_fates")
    isoform = relationship("Isoform", foreign_keys=[isoform_id])


class EnsemblIdCache(Base):
    """Cache for resolving Ensembl IDs (ENST/ENSP/ENSG) that are NOT in our DB.

    When a user searches an Ensembl ID we didn't ingest (a different release, or a
    transcript we don't carry), we resolve it to its gene via the Ensembl REST API
    and cache the answer here — including negative results (`status='not_found'`) so
    we don't re-hit Ensembl for the same miss. `ensembl_id` is stored VERSION-STRIPPED
    (the base accession) since our own IDs are unversioned.
    """
    __tablename__ = "ensembl_id_cache"

    id = Column(Integer, primary_key=True, index=True)
    ensembl_id     = Column(String(40), unique=True, nullable=False, index=True)  # base, unversioned
    id_type        = Column(String(10), nullable=True)     # enst / ensp / ensg
    resolved_ensg  = Column(String(40), nullable=True, index=True)  # the gene stable id
    gene_symbol    = Column(String(100), nullable=True)    # Ensembl display name
    # our gene, iff the resolved ENSG/symbol matches a row in `genes` (else NULL)
    resolved_gene_id = Column(Integer, ForeignKey("genes.id", ondelete="SET NULL"),
                              nullable=True, index=True)
    # found (gene resolved) / not_found (Ensembl has no such id) / error (lookup failed)
    status         = Column(String(20), nullable=False, default="found", index=True)
    checked_at     = Column(DateTime, default=datetime.utcnow)

    # EXTERNAL, stable Ensembl data cached to power the "most similar isoform" line
    # without re-fetching. DB-relative facts (closest isoform + identity, presence)
    # are recomputed live each request, never persisted, so future ingests can't be
    # shadowed by a stale cache.
    query_protein_seq   = Column(Text, nullable=True)
    query_protein_md5   = Column(String(32), nullable=True, index=True)
    protein_seq_checked = Column(Boolean, nullable=False, default=False)

    resolved_gene = relationship("Gene", foreign_keys=[resolved_gene_id])
