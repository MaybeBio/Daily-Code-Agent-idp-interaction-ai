"""
Pfam domain scanning and canonical-vs-novel diff analysis.

Uses pyhmmer for HMM scanning against the Pfam-A database. Produces PfamDomain
records with domain_status classification:

  - "present":   domain found in both canonical and novel, with >= 80% HMM coverage
  - "gained":    domain found in novel but NOT in canonical
  - "lost":      domain found in canonical but NOT in novel (stored on novel isoform
                 with canonical coordinates for reference)
  - "truncated": domain found in novel but covers < 80% of the HMM model length,
                 suggesting a splice event disrupted part of the domain

Truncation threshold: domain_coverage = (hmm_end - hmm_start + 1) / hmm_length.
Values below 0.80 are classified as truncated.

The canonical isoform for comparison is identified via TranscriptMetadata.is_mane_select,
with a fallback to AFDB_* transcript IDs.
"""
import hashlib
import logging
from collections import OrderedDict
from typing import List, Optional, Tuple
from dataclasses import dataclass, replace

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from sqlalchemy import func
from sqlalchemy.orm import Session

from app.models.bio import Isoform, PfamDomain, TranscriptMetadata
from app.analysis.pfam_db import load_hmm_database, load_clan_map, is_kinase_domain
from app.analysis.functional_classes import classify_domain, alteration_type

logger = logging.getLogger(__name__)

TRUNCATION_THRESHOLD = 0.80          # absolute HMM coverage floor
RELATIVE_TRUNCATION_THRESHOLD = 0.85  # isoform/canonical coverage ratio


@dataclass
class DomainHit:
    """Intermediate representation of a single domain hit from hmmscan."""
    pfam_accession: str
    pfam_name: str
    pfam_description: str
    pfam_clan: Optional[str]
    seq_start: int       # 1-based position on the query sequence
    seq_end: int         # 1-based position on the query sequence
    hmm_start: int       # 1-based position on the HMM model
    hmm_end: int         # 1-based position on the HMM model
    hmm_length: int      # total length of the HMM model
    domain_coverage: float
    e_value: float
    bit_score: float


# Per-process scan memo keyed by sequence: analyze_isoform scans the gene's
# canonical for EVERY isoform of the gene, so with gene-ordered targets each
# worker would otherwise redo the same hmmscan dozens of times (the backfill
# spent ~half its time on canonical rescans). Bounded; hits are returned as copies.
_SCAN_MEMO: "OrderedDict[str, List[DomainHit]]" = OrderedDict()
_SCAN_MEMO_MAX = 512


def scan_sequence(
    sequence: str,
    name: str = "query",
    clan_map: Optional[dict] = None,
) -> List[DomainHit]:
    """Memoised front door for `_scan_sequence_uncached` (same contract)."""
    key = hashlib.md5(sequence.encode()).hexdigest() + ("+clan" if clan_map is not None else "")
    hit = _SCAN_MEMO.get(key)
    if hit is not None:
        _SCAN_MEMO.move_to_end(key)
        return [replace(h) for h in hit]
    hits = _scan_sequence_uncached(sequence, name, clan_map)
    _SCAN_MEMO[key] = [replace(h) for h in hits]
    while len(_SCAN_MEMO) > _SCAN_MEMO_MAX:
        _SCAN_MEMO.popitem(last=False)
    return hits


def _scan_sequence_uncached(
    sequence: str,
    name: str = "query",
    clan_map: Optional[dict] = None,
) -> List[DomainHit]:
    """
    Scan a single protein sequence against the full Pfam-A HMM database.

    Uses pyhmmer.hmmscan() which searches one query sequence against many HMM
    profiles. Filtering uses Pfam gathering thresholds (GA) when available.

    Args:
        sequence:  Amino acid string (single-letter codes, no gaps or stops).
        name:      Identifier for the query sequence (used in pyhmmer internals).
        clan_map:  Optional pre-loaded clan map from load_clan_map().  When
                   provided, each DomainHit.pfam_clan is populated from it.
                   Pass None to skip clan annotation (pfam_clan will be None).

    Returns:
        List of DomainHit for all significant domain hits, sorted by seq_start.
    """
    hmms = load_hmm_database()
    alphabet = Alphabet.amino()

    # Build a digital sequence from the protein string
    text_seq = TextSequence(name=name.encode(), sequence=sequence)
    digital_seq = text_seq.digitize(alphabet)

    hits: List[DomainHit] = []

    # pyhmmer.hmmscan: one sequence vs many HMMs
    for top_hits in pyhmmer.hmmscan(
        queries=[digital_seq],
        profiles=hmms,
        cpus=0, # Use all available CPUs (16 on this machine)
        bit_cutoffs="gathering",
    ):
        for hit in top_hits:
            if not hit.reported:
                continue
            for domain in hit.domains:
                if not domain.reported:
                    continue

                aln = domain.alignment
                hmm_length = aln.hmm_length
                # pyhmmer's Alignment mirrors HMMER's P7_ALIDISPLAY: hmm_from/
                # hmm_to/target_from/target_to are ALREADY 1-based, inclusive —
                # identical to `hmmscan --domtblout` "hmm from/to" and "ali from/to".
                # Do NOT add 1 (review DOM-2: every stored coordinate was +1,
                # 192,478 rows had hmm_end = hmm_length + 1 and none started at 1).
                hmm_from = aln.hmm_from
                hmm_to = aln.hmm_to
                seq_from = aln.target_from
                seq_to = aln.target_to

                coverage = (hmm_to - hmm_from + 1) / hmm_length

                def _decode(val):
                    return val.decode() if isinstance(val, bytes) else val

                hmm_accession = _decode(hit.accession) if hit.accession else _decode(hit.name)
                hmm_name = _decode(hit.name)
                hmm_desc = _decode(hit.description) if hit.description else ""

                # Strip version suffix for clan lookup (e.g. "PF00069.29" → "PF00069")
                acc_no_version = hmm_accession.split(".")[0]
                pfam_clan: Optional[str] = None
                if clan_map is not None:
                    clan_val = clan_map.get(acc_no_version, "")
                    pfam_clan = clan_val if clan_val else None

                hits.append(DomainHit(
                    pfam_accession=acc_no_version,
                    pfam_name=hmm_name,
                    pfam_description=hmm_desc,
                    pfam_clan=pfam_clan,
                    seq_start=seq_from,   # 1-based from pyhmmer, no offset
                    seq_end=seq_to,
                    hmm_start=hmm_from,
                    hmm_end=hmm_to,
                    hmm_length=hmm_length,
                    domain_coverage=round(coverage, 4),
                    e_value=domain.i_evalue,
                    bit_score=domain.score,
                ))

    hits = resolve_clan_overlaps(hits)
    hits.sort(key=lambda h: h.seq_start)
    return hits


def resolve_clan_overlaps(hits: List[DomainHit]) -> List[DomainHit]:
    """
    Pfam clan-overlap resolution (the pfam_scan.pl convention).

    Families in the same clan are related models that frequently match the same
    region (e.g. zf-C2H2 / zf-C2H2_4 / zf-C2H2_6 on one zinc finger). Pfam's own
    rule keeps only the best-scoring hit among same-clan hits that overlap on the
    sequence. Without it a 23-finger protein reports 35 "domains" and every
    fate/instance count is inflated (review DOM-3: 130,123 redundant pairs,
    24,199 isoforms). Hits without a clan are never filtered.
    """
    kept: List[DomainHit] = []
    for h in sorted(hits, key=lambda h: (h.e_value, -h.bit_score, h.seq_start)):
        if h.pfam_clan and any(
            k.pfam_clan == h.pfam_clan
            and k.seq_start <= h.seq_end and h.seq_start <= k.seq_end
            for k in kept
        ):
            continue
        kept.append(h)
    return kept


def diff_domains(
    novel_hits: List[DomainHit],
    canonical_hits: List[DomainHit],
) -> Tuple[List[Tuple[DomainHit, str]], List[DomainHit]]:
    """
    Compare domain hits between a novel isoform and its canonical counterpart.

    Algorithm:
        1. Build set of Pfam accessions found in canonical.
        2. For each novel hit:
           - accession NOT in canonical -> "gained"
           - accession in canonical AND coverage >= 0.80 -> "present"
           - accession in canonical AND coverage < 0.80 -> "truncated"
        3. Canonical accessions NOT in any novel hit -> "lost"

    Args:
        novel_hits: Domain hits from the novel isoform scan.
        canonical_hits: Domain hits from the canonical isoform scan.

    Returns:
        Tuple of (novel_classified, lost_canonical_hits):
        - novel_classified: list of (DomainHit, status_string) for each novel domain
        - lost_canonical_hits: list of DomainHit from canonical that are absent in novel
    """
    canonical_accessions = {h.pfam_accession for h in canonical_hits}
    novel_accessions     = {h.pfam_accession for h in novel_hits}

    # Index canonical hits by accession for relative-coverage comparison
    canonical_by_acc: dict = {}
    for h in canonical_hits:
        canonical_by_acc.setdefault(h.pfam_accession, []).append(h)

    novel_classified: List[Tuple[DomainHit, str]] = []

    for hit in novel_hits:
        if hit.pfam_accession not in canonical_accessions:
            novel_classified.append((hit, "gained"))
        else:
            # Relative truncation: compare isoform HMM coverage against the best
            # (highest-coverage) canonical hit for the same domain family.
            # This handles domains that are naturally partial even in canonical
            # (e.g. a zinc-finger with 70% canonical coverage would need only
            #  70% * 0.85 ≈ 60% isoform coverage to still count as "present").
            best_canon = max(canonical_by_acc[hit.pfam_accession],
                             key=lambda c: c.domain_coverage)
            relative = (hit.domain_coverage / best_canon.domain_coverage
                        if best_canon.domain_coverage > 0 else 0.0)
            is_truncated = (hit.domain_coverage < TRUNCATION_THRESHOLD
                            or relative < RELATIVE_TRUNCATION_THRESHOLD)
            novel_classified.append((hit, "truncated" if is_truncated else "present"))

    lost_canonical = [
        h for h in canonical_hits
        if h.pfam_accession not in novel_accessions
    ]

    return novel_classified, lost_canonical


def get_canonical_sequence(db: Session, gene_id: int) -> Optional[str]:
    """
    Get the best canonical protein sequence for a gene.

    This is the primary function for **sequence-based** canonical comparisons
    (Pfam domain diff, NMD truncation ratio, neoantigen junction detection).

    Strategy (best → worst):
        1. Gene.canonical_sequence — the DESIGNATED canonical isoform's sequence
           (set by finalize_reattribution.py from gene.canonical_transcript_id;
           MANE → APPRIS/GENCODE-canonical → longest). It is no longer a UniProt
           sequence: a UniProt canonical can correspond to no isoform in the DB
           (review ING-3).
        2. Longest AFDB isoform sequence — legacy fallback (no AFDB_ rows today).
        3. Longest isoform sequence for the gene — last resort.

    Returns:
        The canonical protein sequence string, or None.
    """
    from app.models.bio import Gene

    gene = db.query(Gene).filter(Gene.id == gene_id).first()
    if not gene:
        return None

    # 1. Designated canonical isoform's sequence (best source)
    if gene.canonical_sequence:
        return gene.canonical_sequence

    # 2. Best AFDB isoform (prefer longest — more likely to be the reviewed entry)
    afdb_isoform = (
        db.query(Isoform)
        .filter(
            Isoform.gene_id == gene_id,
            Isoform.transcript_id.like("AFDB_%"),
            Isoform.sequence.isnot(None),
        )
        .order_by(func.length(Isoform.sequence).desc())
        .first()
    )
    if afdb_isoform and afdb_isoform.sequence:
        return afdb_isoform.sequence

    # 3. Longest isoform overall
    longest = (
        db.query(Isoform)
        .filter(Isoform.gene_id == gene_id, Isoform.sequence.isnot(None))
        .order_by(func.length(Isoform.sequence).desc())
        .first()
    )
    return longest.sequence if longest else None


def find_canonical_isoform(db: Session, gene_id: int) -> Optional[Isoform]:
    """
    Find the canonical isoform **ORM object** for a given gene.

    Use this when you need the actual Isoform record — primarily for
    **structure-based** comparisons (DSSP) where a CIF file is needed,
    or to obtain the canonical_isoform_id FK.

    For sequence-only comparisons, prefer get_canonical_sequence() which
    returns the UniProt reviewed canonical when available.

    Strategy:
        1. Check Gene.canonical_transcript (the established reference)
        2. Look for TranscriptMetadata.is_mane_select = True
        3. Fallback: longest AFDB isoform (prefer SwissProt-sourced, then longest)
        4. Last resort: Return the longest isoform for the gene
    """
    from app.models.bio import Gene

    # 1. Main canonical link
    gene = db.query(Gene).filter(Gene.id == gene_id).first()
    if gene and gene.canonical_transcript:
        return gene.canonical_transcript

    # 2. MANE Select metadata
    canonical = (
        db.query(Isoform)
        .join(TranscriptMetadata, Isoform.id == TranscriptMetadata.isoform_id)
        .filter(
            Isoform.gene_id == gene_id,
            TranscriptMetadata.is_mane_select == True,
        )
        .first()
    )
    if canonical:
        return canonical

    # 3. Longest AFDB isoform (more likely to be the real canonical than a fragment)
    canonical = (
        db.query(Isoform)
        .filter(
            Isoform.gene_id == gene_id,
            Isoform.transcript_id.like("AFDB_%"),
        )
        .order_by(func.length(Isoform.sequence).desc())
        .first()
    )
    if canonical:
        return canonical

    # 4. Longest isoform as last resort
    canonical = (
        db.query(Isoform)
        .filter(Isoform.gene_id == gene_id)
        .order_by(func.length(Isoform.sequence).desc())
        .first()
    )

    return canonical


def analyze_isoform(
    db: Session,
    isoform: Isoform,
    force: bool = False,
) -> List[PfamDomain]:
    """
    Run full Pfam domain analysis for a single isoform.

    Steps:
        1. Skip if results already exist (unless force=True)
        2. Scan the isoform protein sequence against Pfam-A
        3. If a canonical isoform exists for the same gene, scan it and diff
        4. Store PfamDomain records in the database

    Args:
        db: SQLAlchemy session.
        isoform: The isoform to analyze (must have .sequence populated).
        force: If True, delete existing results and re-analyze.

    Returns:
        List of PfamDomain ORM objects that were created.
    """
    if not isoform.sequence:
        logger.warning(f"Isoform {isoform.transcript_id} has no sequence, skipping Pfam scan")
        return []

    # Check for existing results
    existing = db.query(PfamDomain).filter(PfamDomain.isoform_id == isoform.id).first()
    if existing and not force:
        logger.info(f"Pfam results already exist for {isoform.transcript_id}, skipping")
        return db.query(PfamDomain).filter(PfamDomain.isoform_id == isoform.id).all()

    if force:
        db.query(PfamDomain).filter(PfamDomain.isoform_id == isoform.id).delete()

    # Load clan map once — shared by both scans to avoid double-loading
    clan_map = load_clan_map()

    # Scan the novel/target isoform
    novel_hits = scan_sequence(isoform.sequence, isoform.transcript_id, clan_map=clan_map)

    # Get canonical sequence for domain comparison (the designated canonical isoform).
    canonical_seq = get_canonical_sequence(db, isoform.gene_id)
    canonical_isoform = find_canonical_isoform(db, isoform.gene_id)
    canonical_hits: List[DomainHit] = []
    canonical_id: Optional[int] = canonical_isoform.id if canonical_isoform else None

    # Don't diff against self (when the isoform IS the canonical)
    if canonical_id == isoform.id:
        canonical_seq = None
        canonical_id = None

    if canonical_seq:
        canonical_hits = scan_sequence(
            canonical_seq,
            f"canonical_gene_{isoform.gene_id}",
            clan_map=clan_map,
        )

    results: List[PfamDomain] = []

    if canonical_hits:
        novel_classified, lost_in_novel = diff_domains(novel_hits, canonical_hits)

        # Store novel isoform domains with their classified status
        for hit, status in novel_classified:
            domain = _hit_to_model(hit, isoform.id, status, canonical_id)
            db.add(domain)
            results.append(domain)

        # Store "lost" domains — present in canonical but absent in novel.
        # Recorded on the novel isoform so querying pfam_domains for an isoform
        # returns the complete picture. Uses canonical coordinates for reference.
        for hit in lost_in_novel:
            domain = _hit_to_model(hit, isoform.id, "lost", canonical_id)
            db.add(domain)
            results.append(domain)
    else:
        # No canonical available — apply the absolute HMM coverage threshold so
        # partial hits are still flagged as truncated rather than silently "present".
        for hit in novel_hits:
            status = "present" if hit.domain_coverage >= TRUNCATION_THRESHOLD else "truncated"
            domain = _hit_to_model(hit, isoform.id, status, None)
            db.add(domain)
            results.append(domain)

    db.commit()

    counts = {}
    for r in results:
        counts[r.domain_status] = counts.get(r.domain_status, 0) + 1
    logger.info(
        f"Pfam analysis for {isoform.transcript_id}: "
        f"{len(results)} domains ({counts})"
    )

    return results


def _hit_to_model(
    hit: DomainHit,
    isoform_id: int,
    status: str,
    canonical_isoform_id: Optional[int],
) -> PfamDomain:
    """Convert a DomainHit dataclass to a PfamDomain ORM object."""
    fclass = classify_domain(hit.pfam_accession, hit.pfam_clan)
    return PfamDomain(
        isoform_id=isoform_id,
        pfam_accession=hit.pfam_accession,
        pfam_name=hit.pfam_name,
        pfam_description=hit.pfam_description,
        pfam_clan=hit.pfam_clan,
        seq_start=hit.seq_start,
        seq_end=hit.seq_end,
        hmm_start=hit.hmm_start,
        hmm_end=hit.hmm_end,
        hmm_length=hit.hmm_length,
        domain_coverage=hit.domain_coverage,
        e_value=hit.e_value,
        bit_score=hit.bit_score,
        domain_status=status,
        functional_class=fclass,
        alteration_type=alteration_type(
            status,
            hit.hmm_start, hit.hmm_end, hit.hmm_length, hit.domain_coverage,
        ),
        is_kinase=is_kinase_domain(hit.pfam_accession, hit.pfam_clan),
        canonical_isoform_id=canonical_isoform_id,
    )
