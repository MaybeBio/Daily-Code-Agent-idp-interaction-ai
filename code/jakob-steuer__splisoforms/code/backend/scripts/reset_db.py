#!/usr/bin/env python3
"""
reset_db.py
===========
Drops the entire SPLISOFORMS database schema and recreates it from the
current SQLAlchemy models. Use BEFORE any clean-rebuild ingestion run.

Usage:
    # Interactive (asks for confirmation):
    .venv/bin/python backend/scripts/reset_db.py

    # Non-interactive (CI / scripting):
    .venv/bin/python backend/scripts/reset_db.py --yes

Notes
-----
  * Drops `public` schema with CASCADE — destroys every table, sequence,
    index, view, and trigger. Make a `pg_dump` backup beforehand if you
    care about the current state.
  * Recreates tables via `Base.metadata.create_all()` against the current
    model file. Equivalent to running every migration from scratch, but
    bypasses Alembic (alembic table is NOT created automatically — run
    `alembic stamp head` afterwards if you want migrations to track).
  * Seeds the Source table with the three top-level dataset rows the
    Phase 1/2/3 ingest scripts expect to find or create.

After this script
-----------------
    # 1. stamp Alembic to head so future migrations track from here
    .venv/bin/alembic stamp head

    # 2. Phase 1 — GENCODE reference + AF3 structures
    .venv/bin/python backend/scripts/ingest_gencode_phase1.py \\
        --gencode-tsv data/gencode/v47/gencode_v47.tsv \\
        --structures-dir data/gencode_isoforms/af3_export/structures
"""

import argparse
import os
import sys
from pathlib import Path

# Make `app.*` importable regardless of cwd by adding the backend/ root, and the
# scripts/ dir so sibling modules (`_sources_registry`) import directly — matching
# the `from _backfill_io import ...` convention used by the backfill scripts.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sqlalchemy import text

from app.database import engine, Base, SessionLocal
from app.models import bio as bio_models


def confirm(skip: bool) -> bool:
    if skip:
        return True
    print("⚠️  WARNING: This will permanently DROP ALL DATA in the SPLISOFORMS database.")
    print(f"    Engine URL: {engine.url}")
    print("    Type 'yes' to continue, anything else to abort:")
    return input("> ").strip().lower() == "yes"


def reset_db(skip_confirm: bool = False) -> int:
    if not confirm(skip_confirm):
        print("Aborted.")
        return 1

    print("Dropping public schema (CASCADE)...")
    with engine.connect() as conn:
        conn.execute(text("DROP SCHEMA public CASCADE;"))
        conn.execute(text("CREATE SCHEMA public;"))
        # Best-effort: re-grant to splisoforms_user if that role exists.
        # On dev machines the role may differ — skip silently rather than
        # fail the whole reset.
        for stmt in (
            "GRANT ALL ON SCHEMA public TO splisoforms_user;",
            "GRANT ALL ON SCHEMA public TO public;",
        ):
            try:
                conn.execute(text(stmt))
            except Exception as e:  # noqa: BLE001
                print(f"  (note) grant skipped: {e}")
        conn.commit()

    print("Recreating tables from app.models.bio ...")
    Base.metadata.create_all(bind=engine)

    print("Seeding sources from the _sources_registry (single source of truth) ...")
    db = SessionLocal()
    try:
        from _sources_registry import ordered_sources

        sources = [
            bio_models.Source(
                name=s.name,
                source_type=bio_models.SourceTypeEnum(s.source_type),
                cancer_type=s.cancer_type or None,
                tissue_type=s.tissue_type or None,
                sequencing_platform=s.sequencing_platform or None,
                notes=s.notes or None,
            )
            for s in ordered_sources()
        ]
        db.add_all(sources)
        db.commit()
        print(f"  Seeded {len(sources)} sources: " + ", ".join(s.name for s in sources))
    except Exception as e:  # noqa: BLE001
        db.rollback()
        print(f"❌ Error seeding sources: {e}")
        return 2
    finally:
        db.close()

    print("✅ Database reset complete.")
    print("Next: `.venv/bin/alembic stamp head` then Phase 1 ingest.")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--yes", action="store_true",
                        help="Skip the interactive confirmation prompt.")
    args = parser.parse_args()
    sys.exit(reset_db(skip_confirm=args.yes))
