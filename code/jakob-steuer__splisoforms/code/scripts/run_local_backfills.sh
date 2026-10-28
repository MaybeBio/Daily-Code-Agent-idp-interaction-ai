#!/usr/bin/env bash
# run_local_backfills.sh — the ANALYSIS stage of the construction pipeline. Runs on
# the local DB AFTER core ingest + gene reconciliation (rebuild_database.sh, which
# calls scripts/reconcile_genes.sh). Neoantigen is EXCLUDED (runs on the HPC and is
# loaded/scored afterwards — see docs/dataset_ingestion_tutorial.md).
#
# Assumes genes are already reconciled: every isoform sits on its true gene, so the
# gene-keyed steps (uniprot/descriptions/complex/string/drugs) and canonical-relative
# steps (canonical_match/ptm/structural_impact) all key off correct identities.
#
# Continue-on-error: each step logs to backend/logs/backfill_<ts>/NN_<name>.log and
# the chain proceeds even if one step fails. A PASS/FAIL summary is written at the end.
#
# Uses the default DATABASE_URL (the active `splisoforms` DB). Override by exporting
# DATABASE_URL before running.
#
# Usage:  bash scripts/run_local_backfills.sh
#         bash scripts/run_local_backfills.sh --from backfill_ptm.py   # resume
set -uo pipefail
cd "$(dirname "$0")/.."
PY=".venv/bin/python"

# --from <script.py-or-prefix>: skip steps until the named one (resume support).
FROM=""
while [ $# -gt 0 ]; do
  case "$1" in
    --from) FROM="$2"; shift 2 ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
done

TS=$(date +%Y%m%d_%H%M%S)
LOGDIR="backend/logs/backfill_${TS}"
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/summary.txt"
echo "Backfill run ${TS}" | tee "$SUMMARY"

# Ordered chain — REAL script names verified against backend/scripts/. Neoantigen
# omitted by design (remote HPC). Gene reconciliation is NOT here — it runs in the
# ingest stage (reconcile_genes.sh) before this chain.
STEPS=(
  # ── gene-annotation (network: UniProt) — genes are already reconciled ──
  "backfill_gene_uniprot.py"          # ENSG→SwissProt (must precede descriptions)
  # backfill_canonical_uniprot.py is RETIRED from the chain (2026-10): gene.canonical_
  # sequence is now set by finalize_reattribution.py from the designated canonical
  # isoform; overwriting it with a UniProt sequence that may correspond to no
  # isoform broke the Pfam/NMD/neoantigen baselines (review ING-3). It remains as a
  # QC tool (`--report`) only.
  "backfill_gene_descriptions.py"     # protein_name + function_summary
  "backfill_canonical_match.py"       # novel→closest GENCODE isoform WITHIN its gene
  "backfill_transcript_names.py"      # readable display names (GENCODE name / <match>-NOV-<md5>)
  "backfill_ccrcc_orf_flag.py"        # ccRCC folded-ORF QC flag (transparency)
  # ── structural (secondary structure via PyDSSP — pure-Python) ──
  "backfill_structural_metrics_v4.py" # pLDDT/PAE + DSSP  (absolute)
  "backfill_pae_blocks.py"            # PAE rigid blocks   (absolute)
  "backfill_structural_impact.py"     # per-residue exon impact (canonical-relative)
  # ── domains (HMMER + Pfam DB; chainsaw results pre-computed) ──
  "backfill_pfam_parallel.py"
  "backfill_pfam_clans.py"
  "backfill_ted.py"
  "generate_novel_ted.py"
  # ── sequence-level (nmd: 50-nt rule + biotype; disorder: aiupred; ptm: UniProt) ──
  "backfill_nmd.py"
  "backfill_disorder_scores.py"
  "backfill_ptm.py"                   # PTM site accessibility (canonical-relative)
  # ── interaction & expression ──
  "backfill_complex.py"               # AFDB homodimer_analysis + heterodimer_analysis
  "backfill_string_only.py"           # STRING interactions (gene uniprot + STRING API)
  "backfill_cell_type.py"             # HPA single-cell expression (ensembl_gene_id)
  "backfill_expression.py"            # per-sample FL/CPM (data/sources/*/expression/)
  # ── drug–isoform layer (data/drugs/*.tsv + drug_complexes structures present) ──
  # build_drug_metadata re-resolves seq_md5 → CURRENT isoform_id/gene_id from the DB;
  # MUST precede ingest_drugs on a rebuild because the tsv bakes in DB IDs that churn.
  "build_drug_metadata.py"            # (re)build data/drugs/*.tsv from the drug dataset + current DB
  "ingest_drugs.py"                   # drugs (CID-keyed) + drug_complexes (md5-anchored)
  "backfill_drug_chemistry.py"        # SMILES/InChI/formula/weight from PubChem
  "backfill_drug_contacts.py"         # 4.5Å ligand-contact residues from the CIFs
  "backfill_drug_pocket_fate.py"      # pocket fate onto sibling isoforms (needs reconciled genes)
  # ── junctions ──
  "ingest_junctions.py"               # ccRCC SQANTI junctions
  "backfill_bc_junctions.py"          # BC GTF junctions
  # ── scoring (LAST). priority_tier + discovery depend on neoantigen → re-run these
  #    two AFTER the HPC neoantigen load; here they run on an empty table (placeholder).
  "backfill_priority_tier.py"          # neoepitope ranking v2 (P1–P4); supersedes gold_medal
  "backfill_neoantigen_uniqueness.py"  # peptide vs every GENCODE reference protein (unique/self_gene/other_gene)
  "backfill_discovery_index.py"        # Structural Impact Score (SIS v2, 0–100; tier-aware)
  "backfill_genomic_coords.py"        # gene chrom/strand/span (registry already set most)
)

n=0
started=0
[ -z "$FROM" ] && started=1
for step in "${STEPS[@]}"; do
  n=$((n+1))
  name="${step%% *}"; name="${name%.py}"
  # resume gate: skip until we reach the --from step
  if [ "$started" -eq 0 ]; then
    case "$step" in
      "$FROM"*|"${FROM}.py"*) started=1 ;;
      *) echo "── [$n/${#STEPS[@]}] $step  (skipped: before --from $FROM)"; continue ;;
    esac
  fi
  log="$LOGDIR/$(printf '%02d' $n)_${name}.log"
  echo "── [$n/${#STEPS[@]}] $step → $log"
  start=$(date +%s)
  if $PY "backend/scripts/$step" >"$log" 2>&1; then
    dur=$(( $(date +%s) - start ))
    echo "  PASS  $step  (${dur}s)" | tee -a "$SUMMARY"
  else
    dur=$(( $(date +%s) - start ))
    echo "  FAIL  $step  (${dur}s) — see $log" | tee -a "$SUMMARY"
  fi
done

echo "" | tee -a "$SUMMARY"
echo "Done. Summary: $SUMMARY"
echo "Neoantigen NOT run here — stage_hpc_input.py → HPC → load_neoantigen_results.py"
echo "  (use --skip-ptm for a fast load, then backfill_neoantigen_ptm.py for the PTM chips),"
echo "then RE-RUN backfill_priority_tier.py + backfill_neoantigen_uniqueness.py + backfill_discovery_index.py."
