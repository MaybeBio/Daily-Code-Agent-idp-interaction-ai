"""
NMD (Nonsense-Mediated Decay) prediction module.

Classification hierarchy (applied in order):

  1. smORF / uORF candidate guard.  Short-ORF isoforms — TransCODE peptidein
     candidates (Nature 2026) and putative upstream ORFs — are tagged via the
     is_candidate_smorf / is_candidate_uorf flags and routed to
     nmd_status = 'unknown'.  The 50-nt rule is not validated for these
     short, often-translated regulatory ORFs.

  2. Gold-standard 50-nucleotide rule (requires exons_json + TransDecoder
     CDS coords).  Distance is measured between the stop codon and the last
     exon-exon junction in mRNA coordinates.
       - distance > 50 nt AND not 5'-proximal           → nmd_target
       - 5'-proximal exception: stop within 150 nt of the CDS start codon
         → nmd_escape (Karam et al., Hum Mut 2008; ribosome reinitiation
         frees the EJC before NMD machinery engages)
       - distance ≤ 50 nt OR in last exon               → nmd_escape
       - single-exon transcripts have no EJC            → nmd_escape

  3. Heuristic fallback (no CDS data).  Produces an honest verdict:
       - orf_type indicates 5'/3' partial / no_orf      → non_translatable
       - protein looks intact (≥90% of canonical len)   → predicted_translatable
       - any other case                                  → unknown
     The heuristic does NOT set has_ptc=True and does NOT use protein-length
     truncation as evidence of NMD — that signal is exposed separately via
     the orthogonal `truncation_severity` column.

Status values (see also nmd_predictions table):
  nmd_target              50-nt rule, PTC > 50 nt upstream of last EJC,
                          not in 5'-proximal zone.  nmd_confidence is set.
  nmd_escape              50-nt rule, PTC ≤ 50 nt OR in 5'-proximal zone
                          OR in last exon OR single-exon transcript.
  non_translatable        no complete ORF detected.
  predicted_translatable  heuristic only — no PTC evidence found, but we
                          cannot rule out NMD without CDS data.
  unknown                 insufficient data to call (incl. smORF candidates).
"""
import json
import logging
import re
from typing import Dict, List, Optional, Tuple

from sqlalchemy.orm import Session

from app.models.bio import Isoform, NmdPrediction, TranscriptMetadata
from app.analysis.pfam_domains import find_canonical_isoform, get_canonical_sequence

logger = logging.getLogger(__name__)

# truncation_severity thresholds (novel_len / canonical_len)
SEVERE_TRUNCATION = 0.70   # ratio < 0.70 → severe
MINOR_TRUNCATION  = 0.90   # 0.70 <= ratio < 0.90 → minor

# orf_type values that indicate non-translatable transcripts
NON_TRANSLATABLE_PATTERNS = [
    "5prime_partial",
    "3prime_partial",
    "no_orf",
    "ntr",
]

# Short-ORF / peptidein candidate thresholds
SMORF_NOVEL_MAX_AA       = 100   # candidate microprotein cutoff
SMORF_CANONICAL_MAX_AA   = 150   # paired-canonical cutoff for smORF tag
UORF_NOVEL_MAX_AA        = 50    # tiny ORF cutoff for uORF tag
UORF_CANONICAL_MIN_AA    = 150   # only flag uORF when canonical is long

# 5'-proximal escape rule: PTC within this many nt of the CDS start escapes NMD.
# Karam et al., Hum Mut 2008.
PROXIMAL_PTC_ESCAPE_NT   = 150

# nmd_confidence ramp: linear from 0 at distance=50 (threshold) to 1 at
# distance >= NMD_CONFIDENCE_PLATEAU.
NMD_CONFIDENCE_PLATEAU   = 200

# ─────────────────────────────────────────────────────────────────────────────
# TransDecoder pep.fa parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_transdecoder_pep_files(pep_fa_paths: List[str]) -> Dict[str, Tuple[int, int]]:
    """
    Parse one or more TransDecoder *.pep.fa files and return a dict mapping
    each transcript ID to its CDS positions (1-based, in spliced mRNA space).

    Header format (example):
      >PB.100006.4|PB.100006.4:2127764-2145448(+)|molecule/9075797.p2
        ... ORF type:complete (+) ... len:151
        PB.100006.4|PB.100006.4:2127764-2145448(+)|molecule/9075797:344-799(+)

    The last whitespace-separated token encodes
    ``transcript_id|...|molecule/XXXXXXX:CDS_start-CDS_end(strand)``.
    CDS_start and CDS_end are 1-based positions in the transcript.

    When a transcript appears in multiple files (multiple samples), the first
    occurrence wins — coordinates are identical across samples.

    Args:
        pep_fa_paths: List of file paths to TransDecoder protein FASTA files.

    Returns:
        Dict mapping transcript_id → (cds_start_1based, cds_end_1based).
    """
    cds_re = re.compile(r':(\d+)-(\d+)\([+-]\)\s*$')
    result: Dict[str, Tuple[int, int]] = {}
    dups = 0

    for path in pep_fa_paths:
        try:
            with open(path) as fh:
                for line in fh:
                    if not line.startswith(">"):
                        continue
                    line = line.rstrip()
                    # Transcript ID is everything before the first "|"
                    header_id = line[1:].split()[0]          # e.g. "PB.100006.4|..."
                    transcript_id = header_id.split("|")[0]  # e.g. "PB.100006.4"

                    if transcript_id in result:
                        dups += 1
                        continue

                    # CDS position is in the last space-separated token
                    last_token = line.split()[-1]            # e.g. "PB.100006.4|...:344-799(+)"
                    m = cds_re.search(last_token)
                    if m:
                        result[transcript_id] = (int(m.group(1)), int(m.group(2)))
        except OSError as e:
            logger.warning(f"Could not open TransDecoder file {path}: {e}")

    logger.info(
        f"Loaded CDS positions for {len(result):,} transcripts from "
        f"{len(pep_fa_paths)} pep.fa files ({dups:,} cross-sample duplicates skipped)"
    )
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Gold-standard 50-nt rule (exons + TransDecoder CDS)
# ─────────────────────────────────────────────────────────────────────────────

def classify_nmd_50nt_rule(
    exons_json: str,
    strand: str,
    cds_start_1based: int,
    cds_end_1based: int,
) -> Tuple[str, bool, int, Optional[int], Optional[float], str]:
    """
    Apply the 50-nucleotide rule + 5′-proximal escape exception.

    The stop codon occupies mRNA positions [cds_end-2 .. cds_end] (1-based).
    The last exon-exon junction is at mRNA offset =
    total_mRNA_length - last_exon_length.  NMD is predicted when the stop
    is > 50 nt upstream of the last EJC AND not within
    PROXIMAL_PTC_ESCAPE_NT of the CDS start codon.

    Returns:
        (nmd_status, has_ptc, ptc_position_0based, last_ejc_distance,
         nmd_confidence, evidence).
        nmd_confidence is set only for nmd_target rows; NULL otherwise.
    """
    try:
        exons: List[List[int]] = json.loads(exons_json)
    except (json.JSONDecodeError, TypeError):
        return ("unknown", False, cds_start_1based - 1, None, None, "could not parse exons_json")

    if not exons:
        return ("unknown", False, cds_start_1based - 1, None, None, "empty exons_json")

    # Sort exons in transcript (5'→3') order
    if strand == "-":
        exons_ordered = sorted(exons, key=lambda x: x[0], reverse=True)
    else:
        exons_ordered = sorted(exons, key=lambda x: x[0])

    exon_lens = [abs(e[1] - e[0]) + 1 for e in exons_ordered]
    total_mrna_len = sum(exon_lens)

    # Stop codon: last 3 nt of CDS.  cds_end_1based is the last nt of the stop codon.
    stop_start_0based = cds_end_1based - 3
    ptc_position      = stop_start_0based

    # A CDS that runs past the end of the transcript has no annotated stop codon
    # (GENCODE `cds_end_NF` / `mRNA_end_NF`, or a long-read ORF whose exon model
    # belongs to another transcript). The rule cannot be applied — and it must NOT
    # default to "escape" (review NMD-2: 16,360 GENCODE escape calls were of this kind).
    if cds_end_1based > total_mrna_len:
        return (
            "unknown", False, ptc_position, None, None,
            f"CDS end (nt {cds_end_1based}) lies beyond the transcript end "
            f"({total_mrna_len} nt): 3'-incomplete CDS with no annotated stop codon "
            f"→ 50-nt rule not applicable"
        )

    # 5'-proximal escape: PTC within PROXIMAL_PTC_ESCAPE_NT of the CDS start
    # codon escapes NMD (Karam et al., Hum Mut 2008).
    cds_start_0based = cds_start_1based - 1
    distance_from_start = stop_start_0based - cds_start_0based
    if distance_from_start < PROXIMAL_PTC_ESCAPE_NT:
        return (
            "nmd_escape", False, ptc_position, None, None,
            f"5'-proximal escape: stop at mRNA nt {ptc_position} is only "
            f"{distance_from_start} nt past the CDS start "
            f"(<{PROXIMAL_PTC_ESCAPE_NT} → escape via ribosome reinitiation)"
        )

    # Single-exon: no EJC → always NMD escape
    if len(exons_ordered) == 1:
        return (
            "nmd_escape", False, ptc_position, None, None,
            f"single-exon transcript; stop at mRNA nt {ptc_position}, no EJC → NMD escape"
        )

    last_exon_len    = exon_lens[-1]
    last_ejc_0based  = total_mrna_len - last_exon_len

    # Distance: positive = stop is upstream of last EJC (potential NMD)
    distance = last_ejc_0based - stop_start_0based

    if distance > 50:
        # nmd_confidence: linear ramp from 0 at distance=50 to 1 at plateau
        conf = max(0.0, min(1.0, (distance - 50) / (NMD_CONFIDENCE_PLATEAU - 50)))
        return (
            "nmd_target", True, ptc_position, distance, round(conf, 3),
            f"50-nt rule: stop at mRNA nt {ptc_position}, last EJC at nt {last_ejc_0based}, "
            f"distance={distance} nt (>50 → NMD target, confidence={conf:.2f})"
        )
    elif distance > 0:
        return (
            "nmd_escape", False, ptc_position, distance, None,
            f"50-nt rule: stop at mRNA nt {ptc_position}, last EJC at nt {last_ejc_0based}, "
            f"distance={distance} nt (≤50 → NMD escape)"
        )
    else:
        return (
            "nmd_escape", False, ptc_position, None, None,
            f"50-nt rule: stop at mRNA nt {ptc_position} is in or after last exon "
            f"(last EJC at nt {last_ejc_0based}) → NMD escape"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Heuristic fallback (protein truncation + SQANTI3 orf_type)
# ─────────────────────────────────────────────────────────────────────────────

def _parse_orf_type(orf_type: Optional[str]) -> Tuple[str, bool]:
    """
    Parse the SQANTI3 orf_type field.

    Returns:
        Tuple of (cleaned_type, is_complete).
    """
    if not orf_type:
        return ("unknown", False)
    cleaned = orf_type.strip().lower()
    cleaned = re.sub(r"\s*\([+-]\)\s*$", "", cleaned)
    return (cleaned, cleaned == "complete")


def _compute_truncation_ratio(
    isoform: Isoform,
    canonical_sequence: Optional[str],
    metadata: Optional[TranscriptMetadata],
) -> Optional[float]:
    """
    Protein-length ratio novel / canonical, clamped to [0, 1].
    (>1.0 means longer than canonical = not a truncation.)
    """
    novel_len = len(isoform.sequence) if isoform.sequence else None
    if not novel_len:
        return None

    canonical_len = None
    if canonical_sequence:
        canonical_len = len(canonical_sequence)
    if canonical_len is None and metadata and metadata.canonical_protein_length:
        canonical_len = int(metadata.canonical_protein_length)
    if not canonical_len:
        return None

    return round(min(novel_len / canonical_len, 1.0), 4)


def _truncation_severity(truncation_ratio: Optional[float],
                         non_translatable: bool) -> Optional[str]:
    """Map the protein-length ratio to a discrete severity bucket."""
    if non_translatable:
        return "non_translatable"
    if truncation_ratio is None:
        return "intact"
    if truncation_ratio < SEVERE_TRUNCATION:
        return "severe"
    if truncation_ratio < MINOR_TRUNCATION:
        return "minor"
    return "intact"


def classify_nmd_heuristic(
    isoform: Isoform,
    canonical_sequence: Optional[str],
    metadata: Optional[TranscriptMetadata],
    orf_type_override: Optional[str] = None,
) -> Tuple[str, bool, Optional[float], str]:
    """
    Heuristic classifier — used when no CDS coords are available.

    Critically: this no longer uses protein-length truncation as evidence
    of NMD.  Protein-length divergence is a separate orthogonal signal,
    exposed via the `truncation_severity` column.  The heuristic only
    distinguishes:
      - non_translatable          (ORF type indicates incomplete ORF)
      - predicted_translatable    (ORF looks complete; no PTC signal we
                                    can act on)
      - unknown                   (insufficient data to even make that call)

    Returns:
        (nmd_status, has_ptc, truncation_ratio, evidence_string).
        has_ptc is ALWAYS False here — without CDS coords we genuinely
        don't know where the stop codon is.
    """
    evidence_parts = []

    raw_orf_type = orf_type_override or (metadata.orf_type if metadata else None)
    orf_type_clean, is_complete = _parse_orf_type(raw_orf_type)
    evidence_parts.append(f"orf_type={orf_type_clean}")

    for pattern in NON_TRANSLATABLE_PATTERNS:
        if pattern in orf_type_clean:
            evidence_parts.append(f"matches non-translatable pattern '{pattern}'")
            return ("non_translatable", False, None, "; ".join(evidence_parts))

    truncation_ratio = _compute_truncation_ratio(isoform, canonical_sequence, metadata)
    if truncation_ratio is not None:
        evidence_parts.append(f"truncation_ratio={truncation_ratio:.2f}")

    # Canonical / AFDB imports → predicted_translatable (these ARE the canonical)
    if metadata and metadata.is_mane_select:
        evidence_parts.append("is_mane_select=True")
        return ("predicted_translatable", False, truncation_ratio, "; ".join(evidence_parts))
    if isoform.transcript_id and isoform.transcript_id.startswith("AFDB_"):
        evidence_parts.append("AFDB canonical import")
        return ("predicted_translatable", False, truncation_ratio, "; ".join(evidence_parts))

    if is_complete:
        evidence_parts.append("complete ORF, no PTC signal from heuristic")
        return ("predicted_translatable", False, truncation_ratio, "; ".join(evidence_parts))

    if not is_complete and orf_type_clean != "unknown":
        evidence_parts.append("non-complete ORF type")
        return ("non_translatable", False, truncation_ratio, "; ".join(evidence_parts))

    evidence_parts.append("insufficient data for NMD classification")
    return ("unknown", False, truncation_ratio, "; ".join(evidence_parts))


# ─────────────────────────────────────────────────────────────────────────────
# Legacy CDS-sequence path (kept for potential future use)
# ─────────────────────────────────────────────────────────────────────────────

def classify_nmd_from_cds(
    cds_sequence: str,
    exon_junction_positions: List[int],
) -> Tuple[str, bool, Optional[int], Optional[int], str]:
    """
    Apply the 50-nucleotide rule from a stored CDS nucleotide sequence.

    Args:
        cds_sequence:            Nucleotide CDS (ATG … stop codon).
        exon_junction_positions: 0-based CDS positions of exon-exon junctions.

    Returns:
        Tuple of (nmd_status, has_ptc, ptc_position, last_ejc_distance, evidence).
    """
    stop_codons = {"TAA", "TAG", "TGA"}
    cds_upper = cds_sequence.upper()

    ptc_position = None
    for i in range(0, len(cds_upper) - 2, 3):
        if cds_upper[i:i+3] in stop_codons:
            ptc_position = i
            break

    if ptc_position is None:
        return ("non_translatable", False, None, None, "no in-frame stop codon found")

    if not exon_junction_positions:
        return ("nmd_escape", False, ptc_position, None,
                f"stop at nt {ptc_position}, no exon junctions → last exon")

    last_ejc  = max(exon_junction_positions)
    distance  = last_ejc - ptc_position
    if distance > 50:
        return (
            "nmd_target", True, ptc_position, distance,
            f"stop at nt {ptc_position}, last EJC at nt {last_ejc}, "
            f"distance={distance} nt (>50 → NMD target)",
        )
    return (
        "nmd_escape", False, ptc_position, distance,
        f"stop at nt {ptc_position}, last EJC at nt {last_ejc}, "
        f"distance={distance} nt (≤50 → NMD escape)",
    )


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def analyze_isoform(
    db: Session,
    isoform: Isoform,
    force: bool = False,
    cds_map: Optional[Dict[str, Tuple[int, int]]] = None,
    nmd_biotype_map: Optional[Dict[str, str]] = None,
    orf_type_map: Optional[Dict[str, str]] = None,
) -> Optional[NmdPrediction]:
    """
    Run NMD prediction for a single isoform.

    Classification priority:
      0. Authoritative biotype — if ``nmd_biotype_map`` marks the transcript as
         ``nonsense_mediated_decay`` (Ensembl/GENCODE's own curated call), trust
         it directly → nmd_target.  Only annotated (reference) transcripts have
         a biotype; novel isoforms fall through to the computed rule.
      1. 50-nt rule — requires both exons_json AND a CDS entry in ``cds_map``
         (transcript-space CDS coords from TransDecoder pep / projected GTF CDS).
      2. Heuristic — protein truncation ratio + SQANTI3 orf_type.

    Args:
        db:               SQLAlchemy session.
        isoform:          The isoform to analyze.
        force:            Delete existing prediction and re-analyze if True.
        cds_map:          {transcript_id: (cds_start_1based, cds_end_1based)} in
                          spliced-transcript coordinates.  cds_end is the last nt
                          of the stop codon.  None → skip the 50-nt path.
        nmd_biotype_map:  {transcript_id: transcript_biotype}.  A value of
                          'nonsense_mediated_decay' short-circuits to nmd_target.
        orf_type_map:     {transcript_id: orf_type} used to fill orf_type when the
                          stored metadata lacks it (improves the heuristic path).

    Returns:
        NmdPrediction ORM object, or None if isoform has no sequence.
    """
    if not isoform.sequence:
        logger.warning(f"Isoform {isoform.transcript_id} has no sequence, skipping NMD")
        return None

    existing = db.query(NmdPrediction).filter(
        NmdPrediction.isoform_id == isoform.id
    ).first()
    if existing and not force:
        return existing
    if force and existing:
        db.delete(existing)
        db.flush()

    metadata = isoform.metadata_record
    if metadata is None:
        metadata = db.query(TranscriptMetadata).filter(
            TranscriptMetadata.isoform_id == isoform.id
        ).first()

    canonical_seq      = get_canonical_sequence(db, isoform.gene_id)
    canonical_isoform  = find_canonical_isoform(db, isoform.gene_id)
    canonical_id       = (
        canonical_isoform.id
        if (canonical_isoform and canonical_isoform.id != isoform.id)
        else None
    )

    nmd_confidence: Optional[float] = None

    # ── Path 0: smORF / uORF candidate guard ─────────────────────────────────
    # Short ORFs (peptidein candidates per the TransCODE Consortium,
    # Nature 2026) and uORFs are not reliably evaluable by the 50-nt rule.
    novel_len     = len(isoform.sequence) if isoform.sequence else 0
    canonical_len = None
    if canonical_seq:
        canonical_len = len(canonical_seq)
    elif metadata and metadata.canonical_protein_length:
        canonical_len = int(metadata.canonical_protein_length)

    is_candidate_smorf = (
        novel_len > 0 and novel_len < SMORF_NOVEL_MAX_AA
        and (canonical_len is None or canonical_len < SMORF_CANONICAL_MAX_AA)
    )
    is_candidate_uorf = (
        novel_len > 0 and novel_len < UORF_NOVEL_MAX_AA
        and (canonical_len is not None and canonical_len >= UORF_CANONICAL_MIN_AA)
    )

    cds_coords        = (cds_map or {}).get(isoform.transcript_id)
    strand            = metadata.strand if metadata else None
    biotype           = (nmd_biotype_map or {}).get(isoform.transcript_id)
    # orf_type: prefer the injected map, then stored metadata.
    orf_type_injected = (orf_type_map or {}).get(isoform.transcript_id)
    orf_type_stored   = metadata.orf_type if metadata else None
    orf_type_clean    = orf_type_injected or orf_type_stored

    # ── Path 0: authoritative Ensembl/GENCODE NMD biotype ────────────────────
    if biotype == "nonsense_mediated_decay":
        nmd_status, has_ptc, ptc_position, last_ejc_distance, nmd_confidence = \
            "nmd_target", True, None, None, 1.0
        truncation_ratio = _compute_truncation_ratio(isoform, canonical_seq, metadata)
        evidence         = "Ensembl transcript_biotype=nonsense_mediated_decay (curated NMD target)"
        method           = "biotype"

    # ── Path 1: 50-nt rule using exons_json + CDS coords ─────────────────────
    elif cds_coords and isoform.exons_json and strand:
        cds_start, cds_end = cds_coords
        nmd_status, has_ptc, ptc_position, last_ejc_distance, nmd_confidence, evidence = \
            classify_nmd_50nt_rule(isoform.exons_json, strand, cds_start, cds_end)
        truncation_ratio = _compute_truncation_ratio(isoform, canonical_seq, metadata)
        method           = "50nt_rule"

    # ── Path 2: legacy stored CDS sequence (fallback, rarely used) ───────────
    elif isoform.cds_sequence:
        nmd_status, has_ptc, ptc_position, last_ejc_distance, evidence = \
            classify_nmd_from_cds(isoform.cds_sequence, [])
        truncation_ratio = _compute_truncation_ratio(isoform, canonical_seq, metadata)
        method           = "cds_sequence"

    # ── Path 3: heuristic ────────────────────────────────────────────────────
    else:
        nmd_status, has_ptc, truncation_ratio, evidence = \
            classify_nmd_heuristic(isoform, canonical_seq, metadata, orf_type_override=orf_type_injected)
        ptc_position      = None
        last_ejc_distance = None
        method            = "heuristic"

    # ── GENCODE transcript tags (reference set only; stored at ingest) ───────
    # cds_end_NF / mRNA_end_NF: the annotated CDS has no stop codon → the rule is
    # not applicable whatever the arithmetic said. NMD_exception: GENCODE's curated
    # call that this transcript escapes NMD although it has a PTC-like stop (e.g.
    # 3'-UTR introns of MANE transcripts) — honour it over the rule (review NMD-3).
    raw_tags = getattr(metadata, "gencode_tags", None) if metadata else None
    if isinstance(raw_tags, str):            # tolerate a comma-joined string
        raw_tags = raw_tags.split(",")
    gtags = {t for t in (raw_tags or []) if t}
    if method != "biotype" and gtags & {"cds_end_NF", "mRNA_end_NF"}:
        if nmd_status != "unknown":
            evidence = (f"GENCODE tag {'/'.join(sorted(gtags & {'cds_end_NF', 'mRNA_end_NF'}))}: "
                        f"3'-incomplete transcript, no annotated stop codon → "
                        f"50-nt rule not applicable (rule verdict was {nmd_status}: {evidence})")
        nmd_status, has_ptc, nmd_confidence = "unknown", False, None
    elif method != "biotype" and "NMD_exception" in gtags and nmd_status == "nmd_target":
        evidence = (f"GENCODE NMD_exception tag: curated escape overrides the rule "
                    f"(rule verdict was nmd_target: {evidence})")
        nmd_status, nmd_confidence = "nmd_escape", None

    # Apply the candidate guard last so it overrides whichever computed path ran.
    # An authoritative biotype call is never demoted to a smORF guess.
    if (is_candidate_smorf or is_candidate_uorf) and method != "biotype":
        tag = "smORF" if is_candidate_smorf else "uORF"
        evidence = (
            f"{tag} candidate ({novel_len} aa vs canonical "
            f"{canonical_len if canonical_len is not None else '?'} aa); "
            f"50-nt rule not applied — previous verdict was {nmd_status}: {evidence}"
        )
        nmd_status     = "unknown"
        has_ptc        = False
        nmd_confidence = None

    evidence = f"[{method}] {evidence}"

    truncation_severity = _truncation_severity(
        truncation_ratio,
        non_translatable=(nmd_status == "non_translatable"),
    )

    prediction = NmdPrediction(
        isoform_id            = isoform.id,
        nmd_status            = nmd_status,
        has_ptc               = has_ptc,
        truncation_ratio      = truncation_ratio,
        orf_type              = orf_type_clean,
        ptc_position          = ptc_position,
        last_ejc_distance     = last_ejc_distance,
        evidence              = evidence,
        canonical_isoform_id  = canonical_id,
        truncation_severity   = truncation_severity,
        nmd_confidence        = nmd_confidence,
        is_candidate_smorf    = is_candidate_smorf,
        is_candidate_uorf     = is_candidate_uorf,
    )
    db.add(prediction)
    db.commit()

    logger.info(
        f"NMD [{method}] {isoform.transcript_id}: {nmd_status} "
        f"(trunc={truncation_severity}, smORF={is_candidate_smorf}, uORF={is_candidate_uorf})"
    )
    return prediction
