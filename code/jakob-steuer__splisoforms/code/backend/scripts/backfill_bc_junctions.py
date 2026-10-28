#!/usr/bin/env python3
"""
backfill_bc_junctions.py — Populate isoform_junctions for breast-cancer PacBio isoforms.
==========================================================================================

The ccRCC isoforms received their junction events via ingest_junctions.py which reads
SQANTI3 classification files.  The BC dataset uses a different format
(data/lr_breast_cancer_data/isoform_metadata.tsv) and never had a corresponding
ingest step, so the Transcript Splicing Architecture panel is blank for all new BC
isoforms.

This script:
  1. Reads isoform_metadata.tsv (structural_class + tissue per pb_id)
  2. Identifies PB. isoforms that have no isoform_junctions row yet
  3. Derives junction fields from:
       • structural_class   → junction_type  (FSM → "full-splice_match",
                                               novel → "novel_in_catalog")
       • exons_json (DB)    → exon_count_novel
       • canonical exons_json → exon_count_canonical  (via find_canonical_isoform)
       • tissue column      → chrom placeholder (real chrom from exons_json coords)
  4. Inserts one IsoformJunction row per isoform (matching the one-row-per-isoform
     convention already used by the ccRCC data).

Usage
-----
    cd backend
    python scripts/backfill_bc_junctions.py
    python scripts/backfill_bc_junctions.py --data-dir data/lr_breast_cancer_data
    python scripts/backfill_bc_junctions.py --force        # overwrite existing rows
    python scripts/backfill_bc_junctions.py --limit 100    # smoke-test on a subset
"""
import os
import sys
import json
import argparse
import logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
from sqlalchemy import text
from app.database import SessionLocal
from app.models.bio import Isoform, IsoformJunction
from app.analysis.pfam_domains import find_canonical_isoform

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

DEFAULT_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "sources" / "breast_cancer"
METADATA_FILENAME = "isoform_metadata.tsv"

# structural_class → IsoformJunction.junction_type mapping
# Values are chosen to match the ccRCC ingest and the frontend icon map in SequenceMiniMap.
STRUCTURAL_CLASS_MAP = {
    "FSM":   "full-splice_match",
    "novel": "novel_in_catalog",
}


def _exon_count(exons_json_str) -> int | None:
    """Return number of exon blocks from a stored exons_json string, or None."""
    if not exons_json_str:
        return None
    try:
        exons = json.loads(exons_json_str)
        return len(exons) if isinstance(exons, list) else None
    except (json.JSONDecodeError, TypeError):
        return None


def _chrom_from_exons(exons_json_str) -> str | None:
    """Exons are stored as [[start, end], ...] — chrom is not embedded, return None."""
    # The BC assembly GTFs have chrom info but it's not in exons_json.
    # We leave chrom NULL here; it is not rendered by any current UI component.
    return None


def run(data_dir: Path, force: bool, limit: int | None) -> None:
    metadata_path = data_dir / METADATA_FILENAME
    if not metadata_path.exists():
        logger.error(f"Metadata file not found: {metadata_path}")
        return

    logger.info(f"Reading {metadata_path}…")
    df = pd.read_csv(metadata_path, sep="\t", usecols=["pb_id", "structural_class"])

    # De-duplicate: keep one row per pb_id (multiple tissue rows exist per isoform)
    df = df.drop_duplicates(subset="pb_id").set_index("pb_id")
    logger.info(f"{len(df):,} unique pb_ids in metadata")

    db = SessionLocal()
    try:
        # ── Pre-load all PB. isoforms from DB ─────────────────────────────────
        logger.info("Loading PB. isoforms from DB…")
        isoforms = (
            db.query(Isoform)
            .filter(Isoform.transcript_id.like("PB.%"))
            .all()
        )

        # Build lookup: transcript_id → Isoform
        iso_map: dict[str, Isoform] = {iso.transcript_id: iso for iso in isoforms}
        logger.info(f"{len(iso_map):,} PB. isoforms in DB")

        # ── Find those that are missing junction rows ──────────────────────────
        if not force:
            with_junctions_rows = db.execute(
                text("SELECT DISTINCT isoform_id FROM isoform_junctions")
            ).fetchall()
            ids_with_junctions = {r[0] for r in with_junctions_rows}
        else:
            ids_with_junctions = set()

        # ── Build canonical exon-count cache per gene ──────────────────────────
        # Avoid calling find_canonical_isoform() in a tight loop.
        logger.info("Building canonical exon-count cache…")
        gene_ids = {iso.gene_id for iso in isoforms if iso.transcript_id in df.index}
        canon_exon_count: dict[int, int | None] = {}
        for gene_id in gene_ids:
            canon = find_canonical_isoform(db, gene_id)
            if canon and canon.exons_json:
                canon_exon_count[gene_id] = _exon_count(canon.exons_json)
            else:
                canon_exon_count[gene_id] = None
        logger.info(f"Canonical exon counts cached for {len(canon_exon_count):,} genes")

        # ── Insert junction rows ───────────────────────────────────────────────
        inserted  = 0
        skipped   = 0
        not_in_db = 0
        processed = 0

        for pb_id, row in df.iterrows():
            if limit and processed >= limit:
                break

            isoform = iso_map.get(pb_id)
            if isoform is None:
                not_in_db += 1
                continue

            if isoform.id in ids_with_junctions and not force:
                skipped += 1
                continue

            if force:
                db.query(IsoformJunction).filter(
                    IsoformJunction.isoform_id == isoform.id
                ).delete()

            struct_class = str(row["structural_class"]) if pd.notna(row["structural_class"]) else "unknown"
            junction_type = STRUCTURAL_CLASS_MAP.get(struct_class, struct_class)
            is_novel = struct_class != "FSM"
            subcategory = struct_class  # use structural_class as subcategory (ccRCC parity)

            exon_count_novel    = _exon_count(isoform.exons_json)
            exon_count_canonical = canon_exon_count.get(isoform.gene_id)

            junc = IsoformJunction(
                isoform_id          = isoform.id,
                chrom               = None,   # not available in BC metadata
                strand              = None,
                junction_type       = junction_type,
                subcategory         = subcategory,
                exon_count_novel    = exon_count_novel,
                exon_count_canonical = exon_count_canonical,
                is_novel            = is_novel,
            )
            db.add(junc)
            inserted += 1
            processed += 1

            if inserted % 2000 == 0:
                db.flush()
                logger.info(f"  …{inserted:,} inserted so far")

        db.commit()
        logger.info(
            f"Done — inserted: {inserted:,}  |  already had junctions (skipped): {skipped:,}  "
            f"|  pb_id not in DB: {not_in_db:,}"
        )

    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Backfill isoform_junctions for breast-cancer PacBio isoforms",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--data-dir", type=Path, default=DEFAULT_DATA_DIR,
        help=f"Directory containing {METADATA_FILENAME} (default: {DEFAULT_DATA_DIR})",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Delete and re-create junction rows even if they already exist",
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Process at most N isoforms (for smoke-testing)",
    )
    args = parser.parse_args()
    run(data_dir=args.data_dir, force=args.force, limit=args.limit)


if __name__ == "__main__":
    main()
