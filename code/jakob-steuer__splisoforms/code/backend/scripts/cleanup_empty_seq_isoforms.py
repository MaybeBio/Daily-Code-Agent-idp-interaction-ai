#!/usr/bin/env python3
"""
cleanup_empty_seq_isoforms.py
=============================
One-off: remove all isoforms whose sequence is empty / NULL plus their
dependent rows. Also drops orphaned Gene rows (no isoforms left).

Rationale: the upstream long-read pipelines (ccRCC, BC) deliberately
pre-filtered which protein sequences should go to AF3 prediction. The
FASTA (`af3_targets.fa` per cohort) is the authoritative set. Phase-2/3
ingest still creates Isoform rows for metadata-only PBids whose MD5 was
filtered out — those carry no sequence, no structure, no downstream
analyzability and bloat the database without adding signal.

This script deletes those rows in dependency order.

Usage:
    .venv/bin/python backend/scripts/cleanup_empty_seq_isoforms.py
    .venv/bin/python backend/scripts/cleanup_empty_seq_isoforms.py --dry-run
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sqlalchemy import or_

from app.database import SessionLocal
from app.models.bio import (
    Gene,
    Isoform,
    IsoformSourceMembership,
    Structure,
    TranscriptMetadata,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("cleanup_empty_seq_isoforms")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--dry-run", action="store_true",
                        help="Counts only; no DB writes.")
    args = parser.parse_args()

    db = SessionLocal()
    try:
        # 1. Identify isoforms with no sequence.
        empty_iso_ids = [
            iso.id for iso in db.query(Isoform.id)
            .filter(or_(Isoform.sequence == "", Isoform.sequence.is_(None)))
            .all()
        ]
        log.info("isoforms with empty/null sequence: %d", len(empty_iso_ids))

        if not empty_iso_ids:
            log.info("Nothing to do.")
            return 0

        # 2. Count dependents for logging.
        n_meta = (db.query(TranscriptMetadata)
                  .filter(TranscriptMetadata.isoform_id.in_(empty_iso_ids))
                  .count())
        n_struct = (db.query(Structure)
                    .filter(Structure.isoform_id.in_(empty_iso_ids))
                    .count())
        n_member = (db.query(IsoformSourceMembership)
                    .filter(IsoformSourceMembership.isoform_id.in_(empty_iso_ids))
                    .count())
        log.info("  dependent TranscriptMetadata: %d", n_meta)
        log.info("  dependent Structure:          %d", n_struct)
        log.info("  dependent SourceMembership:   %d (auto-cascade)", n_member)

        if args.dry_run:
            # 3. Also identify orphaned genes that WOULD become empty.
            from sqlalchemy import func
            iso_genes = (db.query(Isoform.gene_id, func.count())
                         .filter(Isoform.id.in_(empty_iso_ids))
                         .group_by(Isoform.gene_id).all())
            affected_gene_ids = [g for g, _ in iso_genes if g is not None]
            orphan_count = 0
            for gid in affected_gene_ids:
                total_isos = db.query(Isoform).filter_by(gene_id=gid).count()
                empty_isos = sum(1 for _gid, _n in iso_genes if _gid == gid)
                if total_isos == empty_isos:
                    orphan_count += 1
            log.info("  genes that would become orphan: %d", orphan_count)
            log.info("[dry-run] no changes written.")
            return 0

        # 4. Delete in dependency order. SQLAlchemy bulk-delete with
        #    .delete(synchronize_session=False) is fastest for ~35k rows.
        log.info("deleting dependent TranscriptMetadata ...")
        (db.query(TranscriptMetadata)
           .filter(TranscriptMetadata.isoform_id.in_(empty_iso_ids))
           .delete(synchronize_session=False))

        log.info("deleting dependent Structure ...")
        (db.query(Structure)
           .filter(Structure.isoform_id.in_(empty_iso_ids))
           .delete(synchronize_session=False))

        log.info("deleting dependent IsoformSourceMembership ...")
        (db.query(IsoformSourceMembership)
           .filter(IsoformSourceMembership.isoform_id.in_(empty_iso_ids))
           .delete(synchronize_session=False))

        # 5. Track gene_ids before isoform deletion so we can identify
        #    orphaned genes afterwards.
        affected_gene_ids = set(
            row[0] for row in
            db.query(Isoform.gene_id)
              .filter(Isoform.id.in_(empty_iso_ids))
              .distinct().all()
            if row[0] is not None
        )

        # Null any gene canonical_transcript_id that points to an isoform we are
        # about to delete — that FK is NO ACTION, so the DELETE below would fail
        # otherwise. (These are all-empty genes; they get dropped as orphans in
        # step 6, but the FK must be cleared first.)
        n_canon = (db.query(Gene)
                   .filter(Gene.canonical_transcript_id.in_(empty_iso_ids))
                   .update({"canonical_transcript_id": None}, synchronize_session=False))
        log.info("nulled canonical_transcript_id on %d genes (pointed at empty isoforms)", n_canon)

        log.info("deleting %d empty-sequence isoforms ...", len(empty_iso_ids))
        (db.query(Isoform)
           .filter(Isoform.id.in_(empty_iso_ids))
           .delete(synchronize_session=False))

        db.commit()
        log.info("isoforms deleted.")

        # 6. Drop genes that no longer have any isoforms attached.
        if affected_gene_ids:
            orphans = []
            for gid in affected_gene_ids:
                remaining = db.query(Isoform).filter_by(gene_id=gid).count()
                if remaining == 0:
                    orphans.append(gid)
            log.info("dropping %d orphan genes (no isoforms left) ...", len(orphans))
            if orphans:
                # First null out canonical_transcript_id FK on any such gene
                # (it points to a now-deleted isoform — would FK-constraint anyway).
                (db.query(Gene).filter(Gene.id.in_(orphans))
                 .update({"canonical_transcript_id": None},
                         synchronize_session=False))
                (db.query(Gene).filter(Gene.id.in_(orphans))
                 .delete(synchronize_session=False))
            db.commit()

        # 7. Final counts.
        log.info("─── post-cleanup DB ───")
        log.info(f"  isoforms:     {db.query(Isoform).count():>7}")
        log.info(f"  genes:        {db.query(Gene).count():>7}")
        log.info(f"  structures:   {db.query(Structure).count():>7}")
        log.info(f"  metadata:     {db.query(TranscriptMetadata).count():>7}")
        log.info(f"  memberships:  {db.query(IsoformSourceMembership).count():>7}")
        return 0
    finally:
        db.close()


if __name__ == "__main__":
    sys.exit(main())
