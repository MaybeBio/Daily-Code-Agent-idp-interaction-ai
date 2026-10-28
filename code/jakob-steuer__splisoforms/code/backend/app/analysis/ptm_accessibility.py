"""
PTM site accessibility analysis — canonical vs. novel SASA comparison.

Pipeline
--------
For each gene we retrieve known phosphorylation / ubiquitination / acetylation
sites from UniProt + the EBI proteomics PTM API.  For each site we then:

  1. Map the canonical residue position onto the novel isoform via a
     **splice-aware** alignment (the same algorithm the frontend MSA uses).
  2. Compute SASA at that residue in both the canonical and the novel
     structure (FreeSASA via gemmi).
  3. Extract per-residue pLDDT from both structures (CA B-factor).  Sites
     whose minimum pLDDT is < LOW_PLDDT_CUTOFF are flagged `undetermined`
     because SASA at low-confidence residues is unreliable.
  4. Compare SASA using BOTH an absolute ±SASA_CHANGE_THRESHOLD Å² floor
     AND a relative ±RELATIVE_CHANGE_THRESHOLD (Δ/canonical_sasa) check.
     A site must clear both to be called "exposed" or "buried".

Statuses
--------
  intact        baseline retained (Δ within threshold)        — was 'unchanged'
  exposed       Δ above thresholds and positive               — exposed in novel
  buried        Δ above thresholds and negative               — buried in novel
  lost          canonical residue absent in novel (deleted)
  undetermined  could not compute reliably (missing SASA, low pLDDT, or
                pre-redesign rows where we couldn't compute Δ)

`unchanged` is retained as a synonym for `intact` for backward compat with
any caller that still reads the old label, but new writes use `intact`.
"""
import logging
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

from sqlalchemy.orm import Session

from app.models.bio import Isoform, PtmAccessibility, Structure, TranscriptMetadata, TedDomain
from app.analysis.pfam_domains import find_canonical_isoform
from app.analysis.ptm_data import get_ptm_sites, PTMSite
from app.analysis.sasa import (
    calculate_residue_sasa, calculate_residue_plddt, resolve_structure_path,
)
from app.analysis.splice_aware_align import map_canonical_to_novel

logger = logging.getLogger(__name__)

# Absolute SASA change floor — Δ must clear this AND the relative threshold.
SASA_CHANGE_THRESHOLD = 10.0
# Relative SASA change threshold (Δ / canonical_sasa).  ±0.25 = ±25% change.
RELATIVE_CHANGE_THRESHOLD = 0.25
# Per-residue pLDDT below this → treat SASA as undetermined.  Conventional
# AF confidence cutoff (Tunyasuvunakool et al., Nature 2021).
LOW_PLDDT_CUTOFF = 70.0

# Canonical-structure SASA/pLDDT memo. The backfill walks isoforms gene by gene,
# so consecutive isoforms share a canonical; recomputing its FreeSASA each time
# was the dominant cost of the PTM step (review PTM-7).
_CANON_STRUCT_CACHE: "OrderedDict[str, Tuple[Dict[int, float], Dict[int, float]]]" = OrderedDict()
_CANON_STRUCT_CACHE_MAX = 64


def _canonical_sasa_plddt(cif_path) -> Tuple[Dict[int, float], Dict[int, float]]:
    key = str(cif_path)
    hit = _CANON_STRUCT_CACHE.get(key)
    if hit is not None:
        _CANON_STRUCT_CACHE.move_to_end(key)
        return hit
    val = (calculate_residue_sasa(cif_path), calculate_residue_plddt(cif_path))
    _CANON_STRUCT_CACHE[key] = val
    while len(_CANON_STRUCT_CACHE) > _CANON_STRUCT_CACHE_MAX:
        _CANON_STRUCT_CACHE.popitem(last=False)
    return val


def map_canonical_to_novel_positions(
    canonical_seq: str,
    novel_seq: str,
) -> Dict[int, Optional[int]]:
    """
    Map canonical residue positions to novel-isoform positions via the
    splice-aware aligner (same algorithm as the frontend MSA panel).

    Backwards-compatible name; delegates to `splice_aware_align.map_canonical_to_novel`.
    Replaced the previous BioPython global-NW implementation because NW
    force-aligns biologically-novel exon sequences as mismatches, randomly
    mapping canonical PTM positions onto unrelated isoform residues.

    Returns:
        {canonical 1-based pos → novel 1-based pos or None}.  None means the
        residue is deleted in the novel; an absent canonical key (rare) means
        the canonical position lies in a region we couldn't align at all.
    """
    return map_canonical_to_novel(canonical_seq, novel_seq)


def _get_uniprot_id_for_gene(db: Session, gene_id: int) -> Optional[str]:
    """Get the UniProt ID for a gene, checking both Gene and TranscriptMetadata."""
    from app.models.bio import Gene

    gene = db.query(Gene).filter(Gene.id == gene_id).first()
    if gene and gene.uniprot_id:
        return gene.uniprot_id

    # Fallback: check TranscriptMetadata for any isoform of this gene
    meta = (
        db.query(TranscriptMetadata)
        .join(Isoform, Isoform.id == TranscriptMetadata.isoform_id)
        .filter(
            Isoform.gene_id == gene_id,
            TranscriptMetadata.uniprot_id.isnot(None),
        )
        .first()
    )
    return meta.uniprot_id if meta else None


def _get_structure_cif_path(db: Session, isoform: Isoform):
    """Get the CIF file path for an isoform's first structure."""
    structure = db.query(Structure).filter(Structure.isoform_id == isoform.id).first()
    if not structure:
        return None
    return resolve_structure_path(structure.af3_cif_url)

def is_in_ted_domain(residue_position: int, domain: TedDomain) -> bool:
    if not domain.chopping or domain.chopping == "-":
        return domain.seq_start <= residue_position <= domain.seq_end
    
    segments = domain.chopping.split('_')
    for seg in segments:
        try:
            start, end = map(int, seg.split('-'))
            if start <= residue_position <= end:
                return True
        except ValueError:
            pass
    return False

def get_ted_status(position: int, domains: List[TedDomain]) -> bool:
    return any(is_in_ted_domain(position, d) for d in domains)



def analyze_isoform(
    db: Session,
    isoform: Isoform,
    force: bool = False,
) -> List[PtmAccessibility]:
    """
    Run PTM site accessibility analysis for a single isoform.

    Steps:
        1. Check for existing results (skip unless force=True)
        2. Find canonical isoform and UniProt ID for the gene
        3. Retrieve known PTM sites (PhosphoSitePlus → UniProt fallback)
        4. Calculate SASA on both canonical and novel structures
        5. Map canonical PTM positions to novel positions via alignment
        6. Compare SASA, classify changes, store results

    Args:
        db: SQLAlchemy session.
        isoform: The isoform to analyze.
        force: If True, delete existing results and re-analyze.

    Returns:
        List of PtmAccessibility ORM objects created.
    """
    if not isoform.sequence:
        logger.warning(f"Isoform {isoform.transcript_id} has no sequence, skipping PTM analysis")
        return []

    # Check for existing results
    existing = db.query(PtmAccessibility).filter(PtmAccessibility.isoform_id == isoform.id).first()
    if existing and not force:
        logger.info(f"PTM results already exist for {isoform.transcript_id}, skipping")
        return db.query(PtmAccessibility).filter(PtmAccessibility.isoform_id == isoform.id).all()

    if force:
        db.query(PtmAccessibility).filter(PtmAccessibility.isoform_id == isoform.id).delete()

    # Get UniProt ID for PTM lookup
    uniprot_id = _get_uniprot_id_for_gene(db, isoform.gene_id)
    if not uniprot_id:
        logger.info(f"No UniProt ID for gene of {isoform.transcript_id}, skipping PTM analysis")
        return []

    # Get known PTM sites
    ptm_sites = get_ptm_sites(uniprot_id)
    if not ptm_sites:
        logger.info(f"No PTM sites found for {uniprot_id}, skipping")
        return []

    # Find canonical isoform
    canonical = find_canonical_isoform(db, isoform.gene_id)
    canonical_id: Optional[int] = None
    canonical_sasa: Dict[int, float] = {}
    canonical_plddt: Dict[int, float] = {}
    position_map: Dict[int, Optional[int]] = {}
    can_domains = []

    if canonical and canonical.id != isoform.id and canonical.sequence:
        canonical_id = canonical.id
        canonical_cif = _get_structure_cif_path(db, canonical)
        if canonical_cif:
            canonical_sasa, canonical_plddt = _canonical_sasa_plddt(canonical_cif)
        # Splice-aware mapping — same algorithm as the frontend MSA panel.
        position_map = map_canonical_to_novel_positions(canonical.sequence, isoform.sequence)
        can_domains  = db.query(TedDomain).filter(TedDomain.isoform_id == canonical_id).all()

    # The sequence the site positions are supposed to index. UniProt/EBI report
    # positions on the UniProt canonical sequence version, which is not always the
    # protein we hold (older sequence version, a different isoform, a TrEMBL entry).
    # Review PTM-1: ~11% of rows annotated a residue whose letter did not match.
    # Any site whose residue letter (or position) disagrees with our reference is
    # dropped — it cannot be mapped onto this coordinate frame.
    ref_seq = canonical.sequence if canonical_id else (isoform.sequence or "")
    n_frame_mismatch = 0

    # SASA + pLDDT on the novel structure
    novel_cif = _get_structure_cif_path(db, isoform)
    novel_sasa: Dict[int, float] = {}
    novel_plddt: Dict[int, float] = {}
    if novel_cif:
        novel_sasa  = calculate_residue_sasa(novel_cif)
        novel_plddt = calculate_residue_plddt(novel_cif)

    nov_domains = db.query(TedDomain).filter(TedDomain.isoform_id == isoform.id).all()

    results: List[PtmAccessibility] = []

    for site in ptm_sites:
        can_pos = site.residue_position
        if ref_seq:
            if can_pos < 1 or can_pos > len(ref_seq):
                n_frame_mismatch += 1
                continue
            if site.residue_type and site.residue_type != "?" \
                    and ref_seq[can_pos - 1] != site.residue_type:
                n_frame_mismatch += 1
                continue
        can_sasa_val  = canonical_sasa.get(can_pos)  if canonical_sasa  else None
        can_plddt_val = canonical_plddt.get(can_pos) if canonical_plddt else None

        # Canonical → novel position via splice-aware aligner
        if position_map:
            nov_pos = position_map.get(can_pos)
        elif not canonical_id:
            nov_pos = can_pos    # no canonical at all → assume identity (self-comparison)
        else:
            nov_pos = can_pos    # mapping unavailable → assume identity

        nov_sasa_val  = novel_sasa.get(nov_pos)  if (nov_pos and novel_sasa)  else None
        nov_plddt_val = novel_plddt.get(nov_pos) if (nov_pos and novel_plddt) else None

        # min pLDDT across the two sites — drives the undetermined gate below
        min_plddt: Optional[float] = None
        if can_plddt_val is not None and nov_plddt_val is not None:
            min_plddt = min(can_plddt_val, nov_plddt_val)
        elif can_plddt_val is not None:
            min_plddt = can_plddt_val
        elif nov_plddt_val is not None:
            min_plddt = nov_plddt_val

        # ── Classification (5-way) ────────────────────────────────────────────
        delta: Optional[float] = None
        relative_delta: Optional[float] = None
        change_strength: Optional[float] = None

        if nov_pos is None:
            # Site deleted by splicing — unambiguous biology
            accessibility_change = "lost"
        elif can_sasa_val is None or nov_sasa_val is None:
            # Couldn't measure
            accessibility_change = "undetermined"
        elif min_plddt is not None and min_plddt < LOW_PLDDT_CUTOFF:
            # Structures are too low-confidence at this residue
            accessibility_change = "undetermined"
            delta = round(nov_sasa_val - can_sasa_val, 2)
            if can_sasa_val:
                relative_delta = round(delta / can_sasa_val, 4)
        else:
            delta = round(nov_sasa_val - can_sasa_val, 2)
            relative_delta = round(delta / can_sasa_val, 4) if can_sasa_val > 0 else None
            change_strength = (
                min(abs(relative_delta), 1.0) if relative_delta is not None else None
            )

            # Require BOTH the absolute floor AND the relative threshold to flip a call.
            # A fully buried canonical residue (SASA 0) has no defined relative change;
            # treat any absolute gain above the floor as exposure (review PTM-6).
            abs_clears = abs(delta) > SASA_CHANGE_THRESHOLD
            rel_clears = (
                (relative_delta is not None and abs(relative_delta) > RELATIVE_CHANGE_THRESHOLD)
                or (can_sasa_val == 0 and delta > SASA_CHANGE_THRESHOLD)
            )
            if abs_clears and rel_clears:
                accessibility_change = "exposed" if delta > 0 else "buried"
            else:
                accessibility_change = "intact"

        # Structural context
        can_in_domain = get_ted_status(can_pos, can_domains) if canonical_id else False
        nov_in_domain = (get_ted_status(nov_pos, nov_domains)
                         if (nov_pos is not None) else False)
        if nov_in_domain:
            structural_context = "intra-domain"
        elif can_in_domain and not nov_in_domain:
            structural_context = "disrupted"
        else:
            structural_context = "inter-domain"

        ptm_record = PtmAccessibility(
            isoform_id           = isoform.id,
            uniprot_id           = uniprot_id,
            residue_position     = can_pos,
            novel_residue_position = nov_pos,
            residue_type         = site.residue_type,
            ptm_type             = site.ptm_type,
            source_db            = site.source_db,
            canonical_sasa       = round(can_sasa_val, 2) if can_sasa_val is not None else None,
            novel_sasa           = round(nov_sasa_val, 2) if nov_sasa_val is not None else None,
            delta_sasa           = delta,
            relative_delta_sasa  = relative_delta,
            change_strength      = change_strength,
            min_plddt            = round(min_plddt, 1) if min_plddt is not None else None,
            accessibility_change = accessibility_change,
            structural_context   = structural_context,
            canonical_isoform_id = canonical_id,
        )
        db.add(ptm_record)
        results.append(ptm_record)

    db.commit()

    counts = {}
    for r in results:
        counts[r.accessibility_change] = counts.get(r.accessibility_change, 0) + 1
    logger.info(
        f"PTM analysis for {isoform.transcript_id}: "
        f"{len(results)} sites ({counts}); "
        f"{n_frame_mismatch} site(s) dropped (residue/position off the reference frame)"
    )

    return results
