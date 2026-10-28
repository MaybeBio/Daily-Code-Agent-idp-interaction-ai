#!/usr/bin/env python3
"""
backfill_canonical_match.py — link novel isoforms to their closest GENCODE isoform.

Post-ingest step (replaces rename_to_gencode.py Pass 2). For every non-reference
isoform, find the GENCODE (reference) isoform *within the same gene* with the
highest protein-sequence identity and record:

  * Isoform.canonical_match_id       → that GENCODE isoform's id
  * Isoform.gencode_match_identity   → identity ∈ [0, 1]

Exact protein matches never reach here — at ingest they collapse onto the GENCODE
isoform as an IsoformSourceMembership (md5 dedup). This backfill is for the *near*
matches (a splice variant close-but-not-identical to a reference isoform), which
drive the /compare "closest reference" view.

Identity is difflib's SequenceMatcher ratio (same fast, alignment-free measure the
old rename pass used) — captures indels + substitutions without a BLOSUM matrix.

Usage:
    python backend/scripts/backfill_canonical_match.py
    python backend/scripts/backfill_canonical_match.py --min-identity 0.5
"""
from __future__ import annotations

import argparse
import logging
import sys
from difflib import SequenceMatcher
from pathlib import Path
from typing import List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sqlalchemy.orm import Session

from app.database import SessionLocal
from app.models.bio import Gene, Isoform, Source
from _backfill_io import Progress, print_summary, setup_logging
from _sources_registry import ordered_sources

log = logging.getLogger("backfill_canonical_match")

# Below this identity we record no match (leave canonical_match_id NULL).
DEFAULT_MIN_IDENTITY = 0.50


def _identity(a: str, b: str) -> float:
    """Normalised similarity ∈ [0,1] (difflib). Length-prefilter for speed."""
    if not a or not b:
        return 0.0
    min_len, max_len = min(len(a), len(b)), max(len(a), len(b))
    if max_len == 0:
        return 0.0
    if min_len / max_len < 0.2:
        return min_len / max_len
    return SequenceMatcher(a=a, b=b, autojunk=False).ratio()


def _best_match(seq: str, candidates: List[Isoform]) -> Tuple[Optional[int], float]:
    best_id, best_score = None, 0.0
    for c in candidates:
        if not c.sequence:
            continue
        s = _identity(seq, c.sequence)
        if s > best_score:
            best_score, best_id = s, c.id
    return best_id, best_score


def reference_source_ids(db: Session) -> set:
    ref_names = {s.name for s in ordered_sources() if s.is_reference}
    return {s.id for s in db.query(Source).filter(Source.name.in_(ref_names)).all()}


def run(min_identity: float) -> None:
    setup_logging("backfill_canonical_match")
    db = SessionLocal()
    try:
        ref_ids = reference_source_ids(db)
        if not ref_ids:
            log.error("no reference source found — run the GENCODE ingest first")
            return

        # Process gene-by-gene so candidate sets stay small.
        gene_ids = [g.id for g in db.query(Gene.id).all()]
        stats = {"matched": 0, "no_candidate": 0, "below_threshold": 0,
                 "exact": 0, "skipped_ref": 0}

        with Progress(total=len(gene_ids), desc="canonical-match", unit="gene") as p:
            for i, gid in enumerate(gene_ids, 1):
                isos = db.query(Isoform).filter(Isoform.gene_id == gid).all()
                ref = [x for x in isos if x.source_id in ref_ids and x.sequence]
                nonref = [x for x in isos if x.source_id not in ref_ids]
                for iso in nonref:
                    if not ref:
                        stats["no_candidate"] += 1
                        continue
                    mid, score = _best_match(iso.sequence or "", ref)
                    if mid is None or score < min_identity:
                        stats["below_threshold"] += 1
                        iso.canonical_match_id = None
                        iso.gencode_match_identity = round(score, 4) if mid else None
                        continue
                    iso.canonical_match_id = mid
                    iso.gencode_match_identity = round(score, 4)
                    if score >= 0.999:
                        stats["exact"] += 1
                    stats["matched"] += 1
                if i % 500 == 0:
                    db.commit()
                p.update(1, status=f"matched={stats['matched']}")
        db.commit()

        print_summary(
            "backfill_canonical_match",
            [("genes scanned", len(gene_ids)),
             ("isoforms matched", stats["matched"]),
             ("  near-identical (>=0.999)", stats["exact"]),
             ("below threshold (no link)", stats["below_threshold"]),
             ("no gencode candidate in gene", stats["no_candidate"])],
        )
    finally:
        db.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-identity", type=float, default=DEFAULT_MIN_IDENTITY)
    args = ap.parse_args()
    run(args.min_identity)


if __name__ == "__main__":
    main()
