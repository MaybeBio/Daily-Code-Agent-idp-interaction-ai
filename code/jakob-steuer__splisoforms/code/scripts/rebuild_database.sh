#!/usr/bin/env bash
# rebuild_database.sh — clean, full rebuild of the SPLISOFORMS database from the
# verified per-source metadata.tsv files. Replaces the old phase-script pipeline.
#
#   data/sources/<src>/metadata.tsv  (built by scripts/build_source_metadata.py)
#        │  ingest_from_metadata.py   (md5 dedup, GENCODE-first, inclusion policy)
#        ▼
#   genes · isoforms · memberships · structures · transcript_metadata
#        │  analysis backfills (re-run from scratch — no AF3 re-prediction)
#        ▼
#   fully populated DB
#
# ⚠️  DESTRUCTIVE + ID CHURN: this DROPs the target schema and re-ingests. Isoform
#     numeric IDs are renumbered, so any /isoform/<id> bookmarks change. The public
#     v1 API (keyed by transcript_id / gene symbol) stays stable.
#
# Stepwise is the recommended cadence for the one-time foundational rebuild
# (GENCODE is the dedup anchor — verify it before layering). Pass --all to do all
# three sources in one ingest call.
#
# Usage:
#   DATABASE_URL=postgresql://…/splisoforms_rebuild  bash scripts/rebuild_database.sh
#   bash scripts/rebuild_database.sh --all            # non-stepwise ingest
#   bash scripts/rebuild_database.sh --skip-backfills # core build only
set -euo pipefail

cd "$(dirname "$0")/.."                 # repo root (splisoforms/)
PY=".venv/bin/python"
DB_URL="${DATABASE_URL:-$($PY -c 'import sys;sys.path.insert(0,"backend");from app.database import SQLALCHEMY_DATABASE_URL as u;print(u)')}"

ALL=0; SKIP_BACKFILLS=0
for a in "$@"; do
  case "$a" in
    --all) ALL=1 ;;
    --skip-backfills) SKIP_BACKFILLS=1 ;;
  esac
done

echo "════════════════════════════════════════════════════════════════"
echo "  SPLISOFORMS full rebuild"
echo "  TARGET DB: ${DB_URL%%\?*}"
echo "  ⚠️  This DROPS all data in the target DB and renumbers isoform IDs."
echo "════════════════════════════════════════════════════════════════"
read -r -p "Type 'rebuild' to proceed: " ans
[ "$ans" = "rebuild" ] || { echo "Aborted."; exit 1; }

echo "── [1/5] backup (best-effort pg_dump) ──"
ts=$(date +%Y%m%d_%H%M%S)
mkdir -p backups
if command -v pg_dump >/dev/null 2>&1; then
  pg_dump -Fc "$DB_URL" -f "backups/prerebuild_${ts}.dump" \
    && echo "  → backups/prerebuild_${ts}.dump" \
    || echo "  (pg_dump failed — continuing)"
else
  echo "  (pg_dump not on PATH — skipping backup; dump manually if this is the live DB)"
fi

echo "── [2/5] reset schema from models + seed sources ──"
$PY backend/scripts/reset_db.py --yes
( cd backend && ../.venv/bin/alembic stamp head ) || echo "  (alembic stamp skipped)"

echo "── [3/6] ingest sources ──"
if [ "$ALL" = "1" ]; then
  $PY backend/scripts/ingest_from_metadata.py --all
else
  for src in gencode_v47 ccrcc breast_cancer; do
    echo "    → $src"
    $PY backend/scripts/ingest_from_metadata.py --source "$src" --no-finalize
  done
  echo "── [4/6] finalize (gene canonical + delta_length) ──"
  $PY backend/scripts/ingest_from_metadata.py --finalize-only
fi

# Ground-truth gene reconciliation is part of ingest ("canonical assignment"), NOT
# an analysis backfill: the metadata-driven ingest drops each source's own gene call,
# so novel isoforms must be re-homed to their true gene before anything downstream
# keys off gene identity (canonical baselines, drug pockets, neoantigens). Runs even
# with --skip-backfills so the core DB is always correctly grouped.
echo "── [5/6] reconcile genes (ground-truth re-attribution) ──"
bash scripts/reconcile_genes.sh

if [ "$SKIP_BACKFILLS" = "1" ]; then
  echo "Core build complete (reconciled; analysis backfills skipped)."; exit 0
fi

echo "── [6/6] analysis backfills + drug layer ──"
# Single source of truth for the analysis order lives in run_local_backfills.sh
# (continue-on-error, per-step logs). Neoantigen is excluded there by design — it
# runs on the HPC and is loaded/scored afterwards. See docs/dataset_ingestion_tutorial.md.
bash scripts/run_local_backfills.sh

echo "✅ Rebuild complete."
