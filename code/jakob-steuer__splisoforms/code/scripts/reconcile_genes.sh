#!/usr/bin/env bash
# reconcile_genes.sh — ground-truth gene reconciliation ("canonical assignment"
# stage of the construction pipeline). Runs AFTER core ingest+finalize and BEFORE
# the analysis backfills / drug layer.
#
# The metadata-driven ingest deliberately drops each source's own gene call (BC
# gene_symbol empty, ccRCC mostly empty), so novel isoforms land on UniProt best-hit
# stubs at ingest. This stage re-homes them to their TRUE gene using each dataset's
# ground truth (BC = assembly-GTF gene_id ENSG; ccRCC = SQANTI associated_gene) and
# the full GENCODE v47 annotation as the registry, then re-normalises canonicals and
# drops the emptied stub genes. Without it the gene count roughly doubles.
#
# Order-critical → fail-fast (set -e). Every step is idempotent and re-runnable.
#
# Usage:  bash scripts/reconcile_genes.sh          # apply
#         DRY=1 bash scripts/reconcile_genes.sh    # dry-run every step
set -euo pipefail
cd "$(dirname "$0")/.."
PY=".venv/bin/python"
APPLY="--apply"; [ "${DRY:-0}" = "1" ] && APPLY=""

step() { echo "── reconcile: $1 ──"; $PY "backend/scripts/$1" $APPLY; }

# 1. registry + target map (creates missing ENSG/novel genes, sets coords)
step backfill_gene_registry.py
# 2. re-home novel isoforms per the ground-truth target map
step backfill_regene_from_ground_truth.py
# 3. one canonical per gene (MANE → GENCODE-canonical ENST → longest; never a novel
#    PB.* when an ENST exists)
step backfill_canonical_designation.py
# 4/5. recompute canonical_length + delta_length, drop orphaned stub genes
step finalize_reattribution.py

echo "✅ gene reconciliation complete."
