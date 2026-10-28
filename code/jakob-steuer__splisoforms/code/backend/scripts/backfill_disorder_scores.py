#!/usr/bin/env python3
"""
backfill_disorder_scores.py
===========================
Bulk-populate `disorder_scores` for every isoform with a protein sequence using
**AIUPred** (Erdős & Dosztányi, NAR 2024 — the deep-learning successor of IUPred3)
for per-residue disorder and **AIUPred-binding** for binding propensity. The DB
column `anchor_scores` is a legacy name: it holds AIUPred-binding output, not
ANCHOR2 (review DIS-1 — the method must be cited as AIUPred, not IUPred3/ANCHOR2).

Why this exists
---------------
`disorder_scores` currently has 1 row across the whole DB. The lazy
on-demand path in `/isoform/{id}` only fires for sequences < 1500 aa and
silently swallows failures, so most isoforms never get computed.

`StructuralLandscapePanel` falls back to `DATA_PENDING` when ALL of:
  - structural_metrics.plddt_series  (gated on having an AF3 structure)
  - exon_structural_impacts          (same gate)
  - disorder_scores                  (almost always empty)
are absent. The 11,664 isoforms without structures (mostly tiny
fragments) therefore always show DATA_PENDING.

Backfilling disorder makes the bottom two tracks of the panel
(disorder + binding propensity) populated for *every* isoform with a
sequence — independent of structure availability.

Runtime: ~95 ms per isoform at avg sequence length (CPU, single-thread).
For 140 k isoforms ≈ 3.5-4 h. Skip-existing makes resume safe.
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.database import SessionLocal
from app.models.bio import DisorderScore, Isoform

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("backfill_disorder")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--limit", type=int,
                        help="Cap isoforms processed (testing)")
    parser.add_argument("--max-seq-len", type=int, default=4000,
                        help="Skip isoforms with seq length above this "
                             "(AIUPred slows quadratically; default 4000)")
    parser.add_argument("--min-seq-len", type=int, default=25,
                        help="Skip isoforms shorter than this (AIUPred's "
                             "smoothing window needs ~21 aa; default 25)")
    parser.add_argument("--commit-every", type=int, default=200,
                        help="Commit every N isoforms (default 200)")
    parser.add_argument("--force", action="store_true",
                        help="Re-compute even where rows exist")
    args = parser.parse_args()

    db = SessionLocal()

    # Eligible: has a non-empty sequence + (optionally) no existing row
    log.info("collecting candidate isoforms ...")
    q = db.query(Isoform.id, Isoform.transcript_id, Isoform.sequence)\
          .filter(Isoform.sequence.isnot(None))

    if not args.force:
        existing = set(
            row[0] for row in db.query(DisorderScore.isoform_id).all()
        )
        # SQL-side filter via NOT IN would explode with 140k members;
        # fetch all candidates and filter in Python instead.
        candidates = [(i, t, s) for (i, t, s) in q.all() if i not in existing]
    else:
        candidates = [(i, t, s) for (i, t, s) in q.all()]

    # Drop too-long (AIUPred is quadratic for very long; cap to bound runtime)
    # and too-short (AIUPred's SG smoothing window needs ~21 aa).
    too_long  = [c for c in candidates if c[2] and len(c[2]) > args.max_seq_len]
    too_short = [c for c in candidates if c[2] and len(c[2]) < args.min_seq_len]
    candidates = [c for c in candidates
                  if c[2] and args.min_seq_len <= len(c[2]) <= args.max_seq_len]
    log.info("candidates eligible:    %d", len(candidates))
    log.info("skipping (too long):    %d (>%d aa)", len(too_long), args.max_seq_len)
    log.info("skipping (too short):   %d (<%d aa)", len(too_short), args.min_seq_len)
    if too_long:
        log.info("  long examples:  %s",
                 ", ".join(f"{t}={len(s)}" for _, t, s in too_long[:3]))
    if too_short:
        log.info("  short examples: %s",
                 ", ".join(f"{t}={len(s)}" for _, t, s in too_short[:3]))

    if args.limit:
        candidates = candidates[:args.limit]
        log.info("limited to:             %d", len(candidates))

    if not candidates:
        log.info("nothing to do.")
        db.close()
        return 0

    # Load AIUPred once
    log.info("loading AIUPred ...")
    t_load = time.time()
    from aiupred import AIUPred
    predictor = AIUPred()
    log.info("  loaded in %.1fs", time.time() - t_load)

    # Sort by sequence length so the bar shows realistic ETA early (long
    # ones late = back-loaded slowdown is honest about remaining time)
    candidates.sort(key=lambda c: len(c[2]))

    n_ok = 0
    n_err = 0
    err_samples = []
    t0 = time.time()

    for idx, (iso_id, transcript_id, sequence) in enumerate(candidates, 1):
        try:
            d = predictor.predict_disorder(sequence)
            b = predictor.predict_binding(sequence)

            if args.force:
                db.query(DisorderScore).filter(
                    DisorderScore.isoform_id == iso_id
                ).delete(synchronize_session=False)

            db.add(DisorderScore(
                isoform_id=iso_id,
                disorder_scores=d.tolist(),
                anchor_scores=b.tolist(),
            ))
            n_ok += 1
        except Exception as e:
            n_err += 1
            if len(err_samples) < 5:
                err_samples.append(f"{transcript_id}: {type(e).__name__}: {e}")

        if idx % args.commit_every == 0:
            db.commit()
            el = time.time() - t0
            rate = idx / el if el > 0 else 0
            eta = (len(candidates) - idx) / rate if rate > 0 else 0
            log.info("  [%d/%d] ok=%d err=%d (%.1f iso/s · eta %.0fs)",
                     idx, len(candidates), n_ok, n_err, rate, eta)

    db.commit()

    log.info("✅ summary:")
    log.info("  isoforms processed:  %d", len(candidates))
    log.info("  rows added:          %d", n_ok)
    log.info("  errors:              %d", n_err)
    log.info("  total time:          %.0fs", time.time() - t0)
    if err_samples:
        log.info("  error samples:")
        for e in err_samples:
            log.info("    %s", e)
    db.close()
    return 0 if n_err == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
