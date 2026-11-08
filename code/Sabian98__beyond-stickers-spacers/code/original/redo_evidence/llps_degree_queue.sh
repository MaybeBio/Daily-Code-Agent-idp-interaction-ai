#!/bin/bash
# llps_degree_queue.sh — llps_degree_plan.md queue: waits until the degree_models chain has finished (its task output
# ends with "chain exit"), records the hash of the frozen consensus degree model in gaivi/FROZEN.sha256, then runs the three
# LLPS degree copies once, the Holm family summary, and hashes the outputs.
#   bash llps_degree_queue.sh <chain task output file>
set -u
CHAIN_OUT="$1"
EV="/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence"
PY=~/.venvs/bender/bin/python
until grep -q "chain exit" "$CHAIN_OUT" 2>/dev/null; do sleep 120; done
if ! grep -q "chain exit 0" "$CHAIN_OUT"; then echo "degree chain failed; LLPS queue not started"; exit 1; fi
cd "$EV/gaivi" || exit 1
M="../hub_predictor/degree_consensus_2026-10-13.pkl"
[ -f "$M" ] || { echo "consensus degree model missing"; exit 1; }
echo "$(shasum -a 256 "$M" | cut -d' ' -f1)  $M   # $(date +%Y-%m-%dT%H:%M:%S); frozen consensus degree model (freeze_degree_consensus.py 9b67be94), recorded by llps_degree_queue.sh before any LLPS degree prediction" >> FROZEN.sha256
cd "$EV" || exit 1
export OMP_NUM_THREADS=1
for s in llps_csat_degree llps_unevenness_degree llps_mutations_degree; do
    echo "== $s $(date)"; $PY "$s.py" || { echo "$s failed"; exit 1; }
done
$PY llps_degree_family.py
cd gaivi && TS=$(date +%Y-%m-%dT%H:%M:%S)
for f in ../llps_csat/llps_csat_summary_degree.txt ../llps_csat/llps_features_degree.csv ../llps_csat/llps_unevenness_summary_degree.txt \
         ../llps_csat/llps_unevenness_degree.csv ../llps_csat/llps_mutations_degree.txt ../llps_csat/llps_mutations_degree.csv \
         ../llps_csat/llps_degree_family.txt; do
    echo "$(shasum -a 256 "$f" | cut -d' ' -f1)  $f   # $TS; LLPS degree RESULTS (llps_degree_plan.md), single run" >> FROZEN.sha256
done
echo "llps degree queue done"
