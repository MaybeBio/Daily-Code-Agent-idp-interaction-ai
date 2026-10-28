#!/usr/bin/env python3
"""
backfill_pfam_clans.py — Patch existing pfam_domains rows with clan IDs.
========================================================================

All pfam_domains rows written before clan support was added have
pfam_clan=NULL.  This script:

  1. Downloads Pfam-A.hmm.dat (if not already present) and builds the
     accession → clan map.
  2. Fetches every distinct pfam_accession in the pfam_domains table.
  3. Bulk-updates pfam_clan and is_kinase for each accession in one
     UPDATE … WHERE pfam_accession = ? statement.

Runtime estimate: a few seconds (only as many UPDATE statements as there
are distinct Pfam families in the DB, not one per row).

Usage
-----
    cd backend
    python scripts/backfill_pfam_clans.py [--dry-run]

Options
-------
    --dry-run   Print what would be updated without writing to the DB.
    --stats     Show per-clan family counts after backfill.
"""

import sys
import argparse
import logging
from pathlib import Path

# ── path setup ────────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent))

from sqlalchemy import text
from app.database import SessionLocal
from app.analysis.pfam_db import load_clan_map, is_kinase_domain

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def run(dry_run: bool = False, show_stats: bool = False) -> None:
    db = SessionLocal()
    try:
        # ── 1. Load clan map ─────────────────────────────────────────────────
        clan_map = load_clan_map()
        if not clan_map:
            logger.error(
                "Clan map is empty — Pfam-A.hmm.dat may not have been downloaded. "
                "Run: python -c \"from app.analysis.pfam_db import download_pfam_dat; "
                "download_pfam_dat()\""
            )
            return

        # ── 2. Collect distinct accessions in the DB ─────────────────────────
        rows = db.execute(
            text("SELECT DISTINCT pfam_accession FROM pfam_domains WHERE pfam_accession IS NOT NULL")
        ).fetchall()
        accessions = [r[0].split(".")[0] for r in rows]  # strip version suffix if any
        logger.info(f"Found {len(accessions):,} distinct Pfam accessions in pfam_domains.")

        # ── 3. Build update plan ──────────────────────────────────────────────
        updates: list[dict] = []
        no_clan = 0
        for acc in accessions:
            clan = clan_map.get(acc, "") or None   # "" → None
            kinase = is_kinase_domain(acc, clan)
            updates.append({"acc": acc, "clan": clan, "kinase": kinase})
            if clan is None:
                no_clan += 1

        with_clan = len(updates) - no_clan
        logger.info(
            f"Update plan: {len(updates):,} accessions — "
            f"{with_clan:,} with clan, {no_clan:,} without."
        )

        if dry_run:
            logger.info("DRY RUN — showing first 20 updates:")
            for u in updates[:20]:
                logger.info(f"  {u['acc']} → clan={u['clan']}, is_kinase={u['kinase']}")
            logger.info("DRY RUN complete, no changes written.")
            return

        # ── 4. Apply updates ──────────────────────────────────────────────────
        updated_rows = 0
        for u in updates:
            result = db.execute(
                text(
                    "UPDATE pfam_domains "
                    "SET pfam_clan = :clan, is_kinase = :kinase "
                    "WHERE pfam_accession = :acc OR pfam_accession LIKE :acc_dot"
                ),
                {"acc": u["acc"], "acc_dot": u["acc"] + ".%", "clan": u["clan"], "kinase": u["kinase"]},
            )
            updated_rows += result.rowcount

        db.commit()
        logger.info(f"Committed — {updated_rows:,} rows updated across {len(updates):,} accessions.")

        # ── 5. Optional stats ─────────────────────────────────────────────────
        if show_stats:
            clan_rows = db.execute(
                text(
                    "SELECT pfam_clan, COUNT(*) AS n "
                    "FROM pfam_domains "
                    "WHERE pfam_clan IS NOT NULL "
                    "GROUP BY pfam_clan "
                    "ORDER BY n DESC "
                    "LIMIT 20"
                )
            ).fetchall()
            logger.info("Top 20 clans by row count:")
            for clan, n in clan_rows:
                logger.info(f"  {clan}: {n:,}")

            kinase_count = db.execute(
                text("SELECT COUNT(*) FROM pfam_domains WHERE is_kinase = TRUE")
            ).scalar()
            logger.info(f"Total is_kinase=TRUE rows: {kinase_count:,}")

    finally:
        db.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Backfill pfam_clan and is_kinase in pfam_domains from Pfam-A.hmm.dat",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print planned updates without modifying the database",
    )
    parser.add_argument(
        "--stats", action="store_true",
        help="Print per-clan family counts after backfill",
    )
    args = parser.parse_args()
    run(dry_run=args.dry_run, show_stats=args.stats)


if __name__ == "__main__":
    main()
