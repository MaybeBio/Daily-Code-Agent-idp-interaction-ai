"""
Neoantigen prediction module.

Extracts novel junction peptides from splice-variant isoforms and predicts
MHC class I binding affinity AND antigen presentation likelihood using MHCflurry 2.0.

Pipeline:
  1. Align novel vs canonical protein sequences (BioPython PairwiseAligner)
  2. Identify regions unique to the novel isoform (insertions, substitutions)
  3. Extract overlapping k-mer peptides (8–11 AA) spanning novel junctions,
     including N- and C-terminal flanking sequences for the processing model
  4. Predict via MHCflurry Class1PresentationPredictor:
       • Binding affinity (IC50 nM)
       • Antigen-processing score (proteasomal cleavage probability)
       • Combined presentation score (primary ranking metric)
  5. Classify: strong binder (<50 nM), weak binder (<500 nM), non-binder
     Presentation tier: presented (score > 0.5)

MHCflurry requires model weights downloaded via:
    mhcflurry-downloads fetch

The presentation predictor is MHCflurry 2.0+. Older versions only have
Class1AffinityPredictor; we fall back to affinity-only mode automatically.
"""
import contextlib
import logging
import os
import statistics
import sys
from typing import Any, Dict, List, Optional, Set, Tuple
from dataclasses import dataclass

from sqlalchemy.orm import Session
from sqlalchemy import func

from app.models.bio import Isoform, NeoantigenPrediction, Structure
from app.analysis.pfam_domains import find_canonical_isoform, get_canonical_sequence
from app.analysis.ptm_data import get_ptm_sites, PTMSite
from app.analysis.structural_impact import resolve_structure_path
from app.services.pae_service import PAEService

logger = logging.getLogger(__name__)


@contextlib.contextmanager
def _silence_output():
    """Suppress all stdout and stderr output (Python + C-level).

    MHCflurry emits 'Predicting processing.' / 'Predicting affinities.' via
    print() → stdout, and tqdm progress bars → stderr.  We redirect both
    file descriptors so nothing leaks to the terminal.
    """
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    saved_out_fd = os.dup(1)
    saved_err_fd = os.dup(2)
    saved_py_out = sys.stdout
    saved_py_err = sys.stderr
    os.dup2(devnull_fd, 1)
    os.dup2(devnull_fd, 2)
    sys.stdout = open(os.devnull, "w")
    sys.stderr = open(os.devnull, "w")
    try:
        yield
    finally:
        sys.stdout.close()
        sys.stderr.close()
        os.dup2(saved_out_fd, 1)
        os.dup2(saved_err_fd, 2)
        os.close(saved_out_fd)
        os.close(saved_err_fd)
        os.close(devnull_fd)
        sys.stdout = saved_py_out
        sys.stderr = saved_py_err


# ─── Multi-tool ensemble helpers ─────────────────────────────────────────────

_BA_RESCUE_MIN_CYS     = 2    # min cysteines to consider cysteine-rich
# Audit-3 rewrite (commit f5c2b8d6a3e4 / g8e1d4c7b2a9):
# The original rule required NetMHCpan BA rank < 2.0 AND BOTH EL tools >2.0 —
# never fired (0/1.6M rows) because EL tools agree with BA on the high-rank
# end and the both-EL-miss conjunction is empirically empty.  Replaced with
# a tool-disagreement signal: BA strongly predicts binding (rank < 1.0)
# AND the NetMHCpan EL rank is at least 2× looser than BA's call.
_BA_RESCUE_RANK_THRESH = 1.0  # NetMHCpan-BA %rank below which BA calls binder
_EL_LOOSER_MULTIPLIER  = 2.0  # EL rank must be > multiplier × BA rank for "EL disagrees"

# Agreement thresholds — must match the SQL backfill in the f5c2b8d6a3e4 migration.
_MHCFLURRY_STRONG_RANK = 0.5
_NETMHCPAN_STRONG_RANK = 0.5
_BIGMHC_STRONG_SCORE   = 0.5


def _compute_consensus_rank(
    mhcflurry_percentile: Optional[float],
    netmhcpan_rank_el:    Optional[float],
    bigmhc_el_score:      Optional[float],
) -> Optional[float]:
    """Median %rank across available tools. BigMHC EL score converted via (1-score)*100."""
    ranks: List[float] = []
    if mhcflurry_percentile is not None:
        ranks.append(mhcflurry_percentile)
    if netmhcpan_rank_el is not None:
        ranks.append(netmhcpan_rank_el)
    if bigmhc_el_score is not None:
        ranks.append(max(0.0, min(100.0, (1.0 - bigmhc_el_score) * 100.0)))
    if not ranks:
        return None
    return round(statistics.median(ranks), 3)


def _flag_ba_rescue(
    peptide:           str,
    percentile_rank:   Optional[float],
    netmhcpan_rank_el: Optional[float],
    netmhcpan_rank_ba: Optional[float],
) -> bool:
    """
    True when NetMHCpan-BA confidently predicts binding (rank < 1.0) AND
    NetMHCpan-EL gives a looser rank (> 2× BA), for cysteine-rich peptides.

    Rationale: EL training data (MS-detected presented peptides) is biased
    against cysteine-rich peptides because Cys alkylation and low recovery
    skew the LC-MS sample preparation.  BA-mode training data (in vitro
    binding assays) is unaffected.  When BA strongly predicts binding but
    EL disagrees on a cysteine-rich peptide, the candidate is worth a
    second look.  The `percentile_rank` argument (MHCflurry %rank) is
    accepted for parity with prior versions but is not part of the rule
    after the audit-3 rewrite.
    """
    if netmhcpan_rank_ba is None or netmhcpan_rank_ba >= _BA_RESCUE_RANK_THRESH:
        return False
    if peptide.count("C") < _BA_RESCUE_MIN_CYS:
        return False
    if netmhcpan_rank_el is None:
        return False
    return netmhcpan_rank_el > _EL_LOOSER_MULTIPLIER * netmhcpan_rank_ba


def _compute_agreement(
    mhcflurry_percentile: Optional[float],
    netmhcpan_rank_el:    Optional[float],
    bigmhc_el_score:      Optional[float],
) -> tuple[int, str]:
    """
    Count how many of MHCflurry / NetMHCpan-EL / BigMHC-EL call this
    peptide a strong binder under each tool's native threshold.

    Returns (tools_agreeing, agreement_tier).
    """
    n = 0
    if mhcflurry_percentile is not None and mhcflurry_percentile < _MHCFLURRY_STRONG_RANK:
        n += 1
    if netmhcpan_rank_el is not None and netmhcpan_rank_el < _NETMHCPAN_STRONG_RANK:
        n += 1
    if bigmhc_el_score is not None and bigmhc_el_score > _BIGMHC_STRONG_SCORE:
        n += 1
    tier = "high" if n >= 3 else "medium" if n == 2 else "low"
    return n, tier


# Peptide lengths for MHC-I binding (8-11mers are the standard range)
PEPTIDE_LENGTHS = [8, 9, 10, 11]

# Flanking sequence length for the antigen-processing model
# MHCflurry uses upstream/downstream context to score proteasomal cleavage
FLANK_LEN = 10  # residues on each side of the peptide

# Binding affinity thresholds (IC50 in nM)
STRONG_BINDER_THRESHOLD = 50.0
WEAK_BINDER_THRESHOLD   = 500.0

# Presentation score threshold (MHCflurry 2.0 combined score, 0–1)
PRESENTED_THRESHOLD = 0.5

# ─── HLA Allele Panel ────────────────────────────────────────────────────────
# 35 alleles selected for broad population coverage across European, Asian,
# African, and Latino ancestry groups. Based on the IEDB reference panel
# extended with additional common HLA-C alleles.
#
# Frequency annotations are approximate; see:
#   Sette & Sidney (1999) Immunogenetics 50:201-212
#   Gonzalez-Galarza et al. (2015) Nucleic Acids Res 43:D784-8
DEFAULT_HLA_ALLELES = [
    # ── HLA-A (12 alleles) ──────────────────────────────────────────────────
    "HLA-A*01:01",  # ~16% European
    "HLA-A*02:01",  # ~28% European, ~10% Asian
    "HLA-A*03:01",  # ~14% European
    "HLA-A*11:01",  # ~22% Asian, ~5% European
    "HLA-A*23:01",  # ~9% African
    "HLA-A*24:02",  # ~18% Asian, ~5% European
    "HLA-A*26:01",  # ~5% European
    "HLA-A*29:02",  # ~6% African, ~4% European
    "HLA-A*30:01",  # ~8% African
    "HLA-A*31:01",  # ~7% Asian, ~5% European
    "HLA-A*32:01",  # ~6% European (added for coverage)
    "HLA-A*33:01",  # ~8% Asian
    "HLA-A*68:01",  # ~7% African, ~3% European

    # ── HLA-B (15 alleles) ──────────────────────────────────────────────────
    "HLA-B*07:02",  # ~13% European
    "HLA-B*08:01",  # ~10% European
    "HLA-B*13:02",  # ~7% Asian
    "HLA-B*15:01",  # ~9% European
    "HLA-B*15:02",  # ~8% Asian (SJS-associated)
    "HLA-B*15:03",  # ~10% sub-Saharan African (added for coverage)
    "HLA-B*18:01",  # ~8% European
    "HLA-B*27:05",  # ~8% European (ankylosing spondylitis)
    "HLA-B*35:01",  # ~8% European/Latino
    "HLA-B*38:01",  # ~4% European
    "HLA-B*40:01",  # ~9% Asian
    "HLA-B*44:02",  # ~10% European
    "HLA-B*44:03",  # ~7% African
    "HLA-B*51:01",  # ~7% Asian/Mediterranean
    "HLA-B*53:01",  # ~10% sub-Saharan African (added for coverage)
    "HLA-B*58:01",  # ~8% Asian (allopurinol-HSR allele)

    # ── HLA-C (6 alleles) ───────────────────────────────────────────────────
    # Note: MHCflurry HLA-C coverage is sparser than HLA-A/B.
    # Predictions will fall back to affinity-only mode if the presentation
    # model lacks training data for a given HLA-C allele.
    "HLA-C*03:04",  # ~10% European/Asian
    "HLA-C*04:01",  # ~12% African/European
    "HLA-C*06:02",  # ~9% European
    "HLA-C*07:01",  # ~18% European
    "HLA-C*07:02",  # ~14% European
    "HLA-C*12:03",  # ~5% European
]


@dataclass
class NovelRegion:
    """A region in the novel isoform that differs from canonical."""
    start: int        # 1-based position in novel protein
    end:   int        # 1-based position in novel protein (inclusive)
    region_type: str  # novel_insertion | novel_substitution | frameshift


def find_novel_regions(
    novel_sequence: str,
    canonical_sequence: str,
) -> List[NovelRegion]:
    """
    Identify regions in the novel isoform that are absent or different
    from the canonical isoform.

    Uses the splice-aware aligner (k-mer-anchored chaining, same algorithm
    as the frontend MSA panel).  Previously this routine used BioPython
    global Needleman-Wunsch, which force-aligns biologically-novel exon
    sequences as substitutions and scatters splice events across many
    short mismatches.  That produced wrong novel-region boundaries and,
    therefore, wrong junction-peptide extraction.

    Returns NovelRegion entries with 1-based positions in the novel
    sequence and a region_type of:
      * novel_insertion       — isoform residues with no canonical counterpart
      * novel_substitution    — short bilateral mismatch (≤ SNP_WIDTH in the
                                splice-aware aligner)
      * frameshift            — set heuristically when the novel sequence is
                                substantially shorter and the region lives
                                in the last 10% of the novel
    """
    if not novel_sequence or not canonical_sequence:
        return []

    from app.analysis.splice_aware_align import splice_aware_align

    can_aligned, nov_aligned = splice_aware_align(canonical_sequence, novel_sequence)
    if not can_aligned or not nov_aligned:
        return []

    regions: List[NovelRegion] = []

    try:
        novel_pos = 0         # 0-based in novel sequence
        current_region_start  = None
        current_region_type   = None

        for i in range(min(len(can_aligned), len(nov_aligned))):
            can_char = can_aligned[i]
            nov_char = nov_aligned[i]

            if nov_char == "-":
                if current_region_start is not None:
                    regions.append(NovelRegion(
                        start=current_region_start + 1,
                        end=novel_pos,
                        region_type=current_region_type,
                    ))
                    current_region_start = None
                continue

            novel_pos += 1

            if can_char == "-":
                if current_region_start is None:
                    current_region_start = novel_pos - 1
                    current_region_type  = "novel_insertion"
            elif can_char != nov_char:
                if current_region_start is None:
                    current_region_start = novel_pos - 1
                    current_region_type  = "novel_substitution"
            else:
                if current_region_start is not None:
                    regions.append(NovelRegion(
                        start=current_region_start + 1,
                        end=novel_pos - 1,
                        region_type=current_region_type,
                    ))
                    current_region_start = None

        if current_region_start is not None:
            regions.append(NovelRegion(
                start=current_region_start + 1,
                end=novel_pos,
                region_type=current_region_type,
            ))

    except Exception as e:
        logger.error(f"Alignment parsing failed: {e}")

    # Mark trailing regions as frameshift when novel is substantially shorter
    if len(novel_sequence) < len(canonical_sequence) * 0.7:
        for region in regions:
            if region.end >= len(novel_sequence) * 0.9:
                region.region_type = "frameshift"

    return regions


def extract_junction_peptides(
    sequence: str,
    novel_regions: List[NovelRegion],
    peptide_lengths: List[int] = PEPTIDE_LENGTHS,
    flank_len: int = FLANK_LEN,
) -> List[Tuple[str, int, int, str, str, str]]:
    """
    Extract overlapping k-mer peptides spanning or within novel regions,
    together with their N- and C-terminal flanking sequences.

    Flanking sequences are required for MHCflurry's antigen-processing model,
    which scores whether the proteasome can cleave the peptide out of the source
    protein — critical for splicing-derived neoantigens where the cleavage site
    sits right at the novel junction boundary.

    Args:
        sequence:       Full novel isoform protein sequence.
        novel_regions:  List of novel regions to extract peptides from.
        peptide_lengths: List of peptide lengths to generate (default 8-11).
        flank_len:      Number of flanking residues to extract on each side.

    Returns:
        List of (peptide_seq, start_1based, end_1based, junction_type,
                 n_flank, c_flank) tuples.
        n_flank / c_flank are empty strings when the peptide sits at a terminus.
    """
    seq_len = len(sequence)
    peptides: List[Tuple[str, int, int, str, str, str]] = []
    seen: Set[Tuple[str, str]] = set()   # (peptide_seq, junction_type) dedup

    for region in novel_regions:
        for k in peptide_lengths:
            window_start = max(0, region.start - k)           # 0-based inclusive
            window_end   = min(seq_len - k, region.end - 1)   # 0-based inclusive

            for pos in range(window_start, window_end + 1):
                peptide = sequence[pos : pos + k]
                if len(peptide) < k:
                    continue

                # Skip non-standard residues (X, B, Z, U …)
                if not all(c in "ACDEFGHIKLMNPQRSTVWY" for c in peptide):
                    continue

                key = (peptide, region.region_type)
                if key in seen:
                    continue
                seen.add(key)

                # Extract flanking context for the processing model
                n_flank = sequence[max(0, pos - flank_len) : pos]
                c_flank = sequence[pos + k : min(seq_len, pos + k + flank_len)]

                peptides.append((
                    peptide,
                    pos + 1,          # 1-based start
                    pos + k,          # 1-based end
                    region.region_type,
                    n_flank,
                    c_flank,
                ))

    return peptides


# ─────────────────────────────────────────────────────────────────────────────
# Prediction result type:
#   (affinity_nm, percentile_rank, presentation_score, processing_score)
#   presentation_score / processing_score are None when the presentation
#   predictor is unavailable for this allele (falls back to affinity-only).
# ─────────────────────────────────────────────────────────────────────────────
PredResult = Tuple[float, Optional[float], Optional[float], Optional[float]]


def predict_presentation(
    peptides: List[str],
    alleles: List[str],
    n_flanks: Dict[str, str],
    c_flanks: Dict[str, str],
) -> Dict[Tuple[str, str], PredResult]:
    """
    Predict MHC-I binding affinity AND antigen presentation using MHCflurry 2.0.

    Uses Class1PresentationPredictor (combines affinity + antigen-processing
    model) when available. Falls back to Class1AffinityPredictor per allele
    when the presentation model lacks training data (common for HLA-C).

    Args:
        peptides:  List of unique peptide sequences.
        alleles:   List of HLA alleles to predict against.
        n_flanks:  Dict mapping peptide sequence → N-terminal flank (≤10 aa).
        c_flanks:  Dict mapping peptide sequence → C-terminal flank (≤10 aa).

    Returns:
        Dict mapping (peptide, allele) →
            (affinity_nm, percentile_rank, presentation_score, processing_score)
    """
    if not peptides:
        return {}

    results: Dict[Tuple[str, str], PredResult] = {}

    # Build parallel flank lists aligned to the peptide list
    pep_n_flanks = [n_flanks.get(p, "") for p in peptides]
    pep_c_flanks = [c_flanks.get(p, "") for p in peptides]

    # Try Class1PresentationPredictor (MHCflurry 2.0)
    presentation_predictor = None
    affinity_predictor      = None
    try:
        from mhcflurry import Class1PresentationPredictor
        presentation_predictor = Class1PresentationPredictor.load()
        logger.info("MHCflurry: Using Class1PresentationPredictor (affinity + processing)")
    except (ImportError, Exception) as e:
        logger.warning(f"Class1PresentationPredictor unavailable ({e}), falling back to affinity-only")

    # Also load affinity predictor as fallback for unsupported alleles
    try:
        from mhcflurry import Class1AffinityPredictor
        affinity_predictor = Class1AffinityPredictor.load()
    except (ImportError, Exception) as e:
        logger.error(f"Class1AffinityPredictor also unavailable: {e}")
        if presentation_predictor is None:
            logger.error("MHCflurry not installed. Run: pip install mhcflurry && mhcflurry-downloads fetch")
            return {}

    for allele in alleles:
        # ── Try presentation predictor first ──────────────────────────────
        if presentation_predictor is not None:
            try:
                # MHCflurry 2.0–2.1 used predict_to_dataframe(); 2.2+ uses predict().
                # Try the modern API first, fall back to the legacy name.
                _predict_fn = getattr(
                    presentation_predictor, "predict",
                    getattr(presentation_predictor, "predict_to_dataframe", None),
                )
                if _predict_fn is None:
                    raise AttributeError("No predict method found on Class1PresentationPredictor")

                # Class1PresentationPredictor.predict() takes `alleles` as a
                # *genotype* (list of ≤6 alleles applied to all peptides), NOT
                # a per-peptide list like Class1AffinityPredictor.predict().
                # Pass a single-element genotype so every peptide gets scored
                # against this allele.
                with _silence_output():
                    df = _predict_fn(
                        peptides=peptides,
                        alleles=[allele],
                        n_flanks=pep_n_flanks,
                        c_flanks=pep_c_flanks,
                    )

                for _, row in df.iterrows():
                    pep = row["peptide"]

                    # Affinity (nM) — column name is "affinity" in v2.2+,
                    # "mhcflurry_affinity" in v2.0–2.1, "prediction" in v1.x
                    affinity  = float(row.get("affinity",
                                  row.get("mhcflurry_affinity",
                                  row.get("prediction", 999999))))
                    # Percentile: "presentation_percentile" in v2.2+,
                    # "affinity_percentile" / "prediction_percentile" in older versions
                    percentile = row.get("presentation_percentile",
                                  row.get("affinity_percentile",
                                  row.get("mhcflurry_affinity_percentile",
                                  row.get("prediction_percentile", None))))
                    pres_score = row.get("presentation_score", None)
                    proc_score = row.get("processing_score", None)

                    results[(pep, allele)] = (
                        affinity,
                        float(percentile) if percentile is not None else None,
                        float(pres_score) if pres_score is not None else None,
                        float(proc_score) if proc_score is not None else None,
                    )
                continue  # next allele — done with this one

            except Exception as e:
                logger.debug(
                    f"Presentation predictor failed for {allele} "
                    f"(likely unsupported allele); falling back to affinity. Error: {e}"
                )

        # ── Fallback: affinity-only ────────────────────────────────────────
        if affinity_predictor is not None:
            try:
                with _silence_output():
                    df = affinity_predictor.predict_to_dataframe(
                        peptides=peptides,
                        alleles=[allele] * len(peptides),
                    )

                for _, row in df.iterrows():
                    pep       = row["peptide"]
                    affinity  = float(row.get("mhcflurry_affinity",
                                              row.get("prediction", 999999)))
                    percentile = row.get("mhcflurry_affinity_percentile",
                                         row.get("prediction_percentile", None))

                    results[(pep, allele)] = (
                        affinity,
                        float(percentile) if percentile is not None else None,
                        None,   # presentation_score not available
                        None,   # processing_score not available
                    )

            except Exception as e:
                logger.warning(f"Affinity prediction also failed for {allele}: {e}")

    return results


def classify_binding(affinity_nm: float) -> str:
    """Classify binding affinity into categories."""
    if affinity_nm <= STRONG_BINDER_THRESHOLD:
        return "strong_binder"
    elif affinity_nm <= WEAK_BINDER_THRESHOLD:
        return "weak_binder"
    else:
        return "non_binder"


# Module-level PTM cache: uniprot_id → List[PTMSite].
_ptm_cache: Dict[str, List] = {}


def analyze_isoform(
    db: Session,
    isoform: Isoform,
    force: bool = False,
    hla_alleles: Optional[List[str]] = None,
) -> List[NeoantigenPrediction]:
    """
    Run neoantigen prediction for a single isoform.

    Steps:
        1. Align novel vs canonical to find novel regions
        2. Extract junction-spanning peptides (8-11mers) with flanking sequences
        3. Predict MHC-I binding AND presentation via MHCflurry 2.0
        4. Store results (only binders, or top candidates)

    Args:
        db:          SQLAlchemy session.
        isoform:     The novel isoform to analyze.
        force:       If True, delete existing results and re-analyze.
        hla_alleles: Custom HLA allele panel (default: DEFAULT_HLA_ALLELES).

    Returns:
        List of NeoantigenPrediction ORM objects (binders only).
    """
    if not isoform.sequence:
        logger.warning(f"Isoform {isoform.transcript_id} has no sequence, skipping neoantigen")
        return []

    existing = db.query(NeoantigenPrediction).filter(
        NeoantigenPrediction.isoform_id == isoform.id
    ).first()
    if existing and not force:
        logger.info(f"Neoantigen predictions already exist for {isoform.transcript_id}, skipping")
        return db.query(NeoantigenPrediction).filter(
            NeoantigenPrediction.isoform_id == isoform.id
        ).all()

    if force:
        db.query(NeoantigenPrediction).filter(
            NeoantigenPrediction.isoform_id == isoform.id
        ).delete()
        db.flush()

    canonical_seq      = get_canonical_sequence(db, isoform.gene_id)
    canonical_isoform  = find_canonical_isoform(db, isoform.gene_id)
    canonical_id: Optional[int] = canonical_isoform.id if canonical_isoform else None

    if canonical_isoform and canonical_isoform.id == isoform.id:
        canonical_seq = None

    if not canonical_seq:
        logger.info(f"No canonical to compare for {isoform.transcript_id}, skipping")
        return []

    # 1. Find novel regions
    novel_regions = find_novel_regions(isoform.sequence, canonical_seq)
    if not novel_regions:
        logger.info(f"No novel regions found for {isoform.transcript_id}")
        return []

    # 2. Extract peptides + flanking sequences
    peptides_info = extract_junction_peptides(isoform.sequence, novel_regions)

    # ── Canonical self-peptide filter ──────────────────────────────────────
    pre_filter_count = len(peptides_info)
    peptides_info = [
        entry for entry in peptides_info
        if entry[0] not in canonical_seq
    ]
    filtered_count = pre_filter_count - len(peptides_info)
    if filtered_count:
        logger.info(
            f"Removed {filtered_count} canonical self-peptides from {isoform.transcript_id} "
            f"({pre_filter_count} → {len(peptides_info)} peptides remain)"
        )

    if not peptides_info:
        logger.info(f"No novel junction peptides after canonical filter for {isoform.transcript_id}")
        return []

    # 3. Build per-peptide flank lookup (first occurrence wins for deduped peptides)
    n_flanks: Dict[str, str] = {}
    c_flanks: Dict[str, str] = {}
    for pep, _, _, _, n_flank, c_flank in peptides_info:
        if pep not in n_flanks:
            n_flanks[pep] = n_flank
            c_flanks[pep] = c_flank

    alleles          = hla_alleles or DEFAULT_HLA_ALLELES
    unique_peptides  = list(n_flanks.keys())   # already deduplicated

    logger.info(
        f"Predicting presentation for {len(unique_peptides)} unique peptides "
        f"against {len(alleles)} HLA alleles for {isoform.transcript_id}"
    )

    # 4. MHCflurry 2.2 — binding affinity + antigen-processing + presentation
    binding_results = predict_presentation(unique_peptides, alleles, n_flanks, c_flanks)

    if not binding_results:
        logger.warning(f"No binding predictions returned for {isoform.transcript_id}")
        return []

    # 4b. Multi-tool enrichment — NetMHCpan-4.2 + BigMHC (if available)
    netmhcpan_results: dict = {}
    bigmhc_results:    dict = {}

    try:
        from app.analysis.netmhcpan_runner import is_available as _nmhc_ok, run_netmhcpan_multi_allele
        if _nmhc_ok():
            logger.info(f"  Running NetMHCpan-4.2 for {isoform.transcript_id}…")
            netmhcpan_results = run_netmhcpan_multi_allele(unique_peptides, alleles)
    except Exception as _e:
        logger.debug(f"  NetMHCpan skipped: {_e}")

    try:
        from app.analysis.bigmhc_runner import is_available as _bmhc_ok, run_bigmhc
        if _bmhc_ok():
            logger.info(f"  Running BigMHC EL+IM for {isoform.transcript_id}…")
            bigmhc_results = run_bigmhc(unique_peptides, alleles)
    except Exception as _e:
        logger.debug(f"  BigMHC skipped: {_e}")

    # 5. PTM & structural confidence (unchanged logic, v4)
    from app.analysis.ptm_accessibility import _get_uniprot_id_for_gene
    uniprot_id = _get_uniprot_id_for_gene(db, isoform.gene_id)
    canonical_ptms: List[PTMSite] = []
    if uniprot_id:
        if uniprot_id not in _ptm_cache:
            _ptm_cache[uniprot_id] = get_ptm_sites(uniprot_id)
        canonical_ptms = _ptm_cache[uniprot_id]

    plddt_series: List[float] = []
    pae_blocks:   List[dict]  = []
    if isoform.structural_metrics:
        plddt_series = isoform.structural_metrics.plddt_series or []
        if isoform.structural_metrics.pae_blocks:
            try:
                import json
                pae_blocks = json.loads(isoform.structural_metrics.pae_blocks)
            except Exception:
                pae_blocks = []

    # Build peptide → [(start, end, junction_type)] lookup
    peptide_info_map: Dict[str, List[Tuple[int, int, str]]] = {}
    for pep, start, end, jtype, _, _ in peptides_info:
        peptide_info_map.setdefault(pep, []).append((start, end, jtype))

    # 6. Assemble results — only store binders
    results: List[NeoantigenPrediction] = []
    _seen_keys: Set[Tuple[str, str, int]] = set()   # (peptide, allele, start) dedup guard

    for (peptide, allele), (affinity, percentile, pres_score, proc_score) in binding_results.items():
        category = classify_binding(affinity)

        # Primary filter: must be a binder (affinity) OR presented (presentation_score)
        is_presented = (pres_score is not None and pres_score >= PRESENTED_THRESHOLD)
        if category == "non_binder" and not is_presented:
            continue

        # Multi-tool lookup for this (peptide, allele) pair
        _key  = (peptide, allele)
        _nmhc = netmhcpan_results.get(_key)   # NetMHCpanResult namedtuple or None
        _bmhc = bigmhc_results.get(_key)      # BigMHCResult namedtuple or None

        _consensus = _compute_consensus_rank(
            percentile,
            _nmhc.rank_el if _nmhc else None,
            _bmhc.el_score if _bmhc else None,
        )

        for start, end, junction_type in peptide_info_map.get(peptide, []):
            _key = (peptide, allele, start)
            if _key in _seen_keys:
                continue
            _seen_keys.add(_key)

            # Structural confidence window
            window_plddt = plddt_series[start - 1 : end] if plddt_series and len(plddt_series) >= end else []
            avg_plddt    = sum(window_plddt) / len(window_plddt) if window_plddt else None
            min_plddt    = min(window_plddt) if window_plddt else None

            # PAE rigidity check
            is_rigid = any(
                block["start"] <= start and end <= block["end"]
                for block in pae_blocks
            )

            # Tiered PTM impact classification
            impact_status = "Standard"
            is_sensitive  = False

            ptm_canonical_residue_lost = False
            if canonical_ptms:
                direct_ptms   = [p for p in canonical_ptms if start <= p.residue_position <= end]
                flanking_ptms = [
                    p for p in canonical_ptms
                    if (start - 3 <= p.residue_position < start)
                    or (end < p.residue_position <= end + 3)
                ]

                if direct_ptms:
                    # Distinguish: was the modifiable residue preserved in the
                    # novel (PTM may still occur in a new sequence context) or
                    # is it physically gone (modification removed)?
                    any_lost = False
                    for p in direct_ptms:
                        if 1 <= p.residue_position <= len(isoform.sequence):
                            novel_aa = isoform.sequence[p.residue_position - 1]
                            if novel_aa != p.residue_type:
                                any_lost = True
                        else:
                            any_lost = True
                    impact_status = (
                        "Direct PTM-Neoantigen (residue lost)"
                        if any_lost else
                        "Direct PTM-Neoantigen (residue retained)"
                    )
                    is_sensitive  = True
                    ptm_canonical_residue_lost = any_lost
                elif flanking_ptms:
                    found_flanking = False
                    lost_flanking  = False
                    for p in flanking_ptms:
                        if 1 <= p.residue_position <= len(isoform.sequence):
                            if isoform.sequence[p.residue_position - 1] == p.residue_type:
                                found_flanking = True
                            else:
                                lost_flanking = True
                        else:
                            lost_flanking = True

                    if found_flanking:
                        impact_status = "Processing-Blocked"
                        is_sensitive  = True
                    elif lost_flanking:
                        # ranking v2: keep the mechanistic label only; pLDDT/rigidity are
                        # now the orthogonal high-fidelity-fold flag + priority_tier.
                        impact_status = "Splice-Enabled"
                        is_sensitive  = True

            _ba_rescue = _flag_ba_rescue(
                peptide,
                percentile,
                _nmhc.rank_el if _nmhc else None,
                _nmhc.rank_ba if _nmhc else None,
            )

            _tools_agreeing, _agreement_tier = _compute_agreement(
                percentile,
                _nmhc.rank_el if _nmhc else None,
                _bmhc.el_score if _bmhc else None,
            )

            record = NeoantigenPrediction(
                isoform_id           = isoform.id,
                peptide_sequence     = peptide,
                peptide_start        = start,
                peptide_end          = end,
                hla_allele           = allele,
                binding_affinity_nm  = round(affinity, 2),
                binding_category     = category,
                percentile_rank      = round(percentile, 4) if percentile is not None else None,
                junction_type        = junction_type,
                canonical_isoform_id = canonical_id,
                # Presentation model scores (MHCflurry 2.0)
                presentation_score   = round(pres_score, 4) if pres_score is not None else None,
                processing_score     = round(proc_score, 4) if proc_score is not None else None,
                # PTM / structural fields
                ptm_impact_status    = impact_status,
                avg_plddt            = round(avg_plddt, 2) if avg_plddt is not None else None,
                min_plddt            = round(min_plddt, 2) if min_plddt is not None else None,
                is_ptm_sensitive     = is_sensitive,
                # Multi-tool ensemble (NetMHCpan-4.2 + BigMHC)
                netmhcpan_affinity_nm = round(_nmhc.affinity_nm, 2) if _nmhc else None,
                netmhcpan_rank_el     = round(_nmhc.rank_el, 4)     if _nmhc else None,
                netmhcpan_rank_ba     = round(_nmhc.rank_ba, 4)     if _nmhc else None,
                bigmhc_el_score       = round(_bmhc.el_score, 4)    if (_bmhc and _bmhc.el_score is not None) else None,
                bigmhc_im_score       = round(_bmhc.im_score, 4)    if (_bmhc and _bmhc.im_score is not None) else None,
                consensus_rank        = _consensus,
                is_ba_rescue          = _ba_rescue,
                # Audit-3 additions (Option A)
                tools_agreeing             = _tools_agreeing,
                agreement_tier             = _agreement_tier,
                ptm_canonical_residue_lost = ptm_canonical_residue_lost,
            )
            results.append(record)

    if results:
        db.bulk_save_objects(results)
        db.commit()

    strong    = sum(1 for r in results if r.binding_category == "strong_binder")
    presented = sum(1 for r in results if r.presentation_score is not None and r.presentation_score >= PRESENTED_THRESHOLD)
    splice    = sum(1 for r in results if (r.ptm_impact_status or "").startswith("Splice-Enabled"))
    logger.info(
        f"Neoantigen predictions for {isoform.transcript_id}: "
        f"{len(results)} binders ({strong} strong, {presented} presented, {splice} splice-enabled) "
        f"from {len(novel_regions)} novel regions across {len(alleles)} alleles"
    )

    return results
