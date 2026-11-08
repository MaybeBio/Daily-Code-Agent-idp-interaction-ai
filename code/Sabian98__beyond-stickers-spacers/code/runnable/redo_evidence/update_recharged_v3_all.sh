#!/bin/bash
# update_recharged_v3_all.sh — everything downstream of the Mpipi-Recharged labels, after label_recharged_v3.sh has
# added new runs (2026-10-04 evening: 5,417 -> about 6,060 labelled). Supersedes update_recharged_v3.sh, which predates
# fig3_arrangement.py and the analyses added on Oct 4. Current outputs are kept with the suffix _<TAG> (default _5417)
# so the old and new numbers can be compared. Truncated or short trajectories never enter: hub_labels_gaivi.py rejects
# them (header vs file size, and frames below the protocol target) and lists them in hub_frames_v3_recharged/failed.csv.
# Frozen predictions and the glycine sub-batch design are never touched. Run from redo_evidence/:
#   bash update_recharged_v3_all.sh [TAG]      (about 1-1.5 h; each step can be rerun on its own)
#   START=3 THREADS=8 bash update_recharged_v3_all.sh [TAG]   resume at step 3; THREADS = OpenMP threads for the
#   single-process steps (the residue-neutral null runs its own 8 workers and always uses 1 thread each)
set -e
cd "$(dirname "$0")"
TAG=${1:-5417}
START=${START:-1}
export OMP_NUM_THREADS=${THREADS:-${OMP_NUM_THREADS:-1}}
step() { (( START <= $1 )); }
PY=~/.venvs/bender/bin/python
HP=hub_predictor
FIG=../methods_paper/figures
DECK=../methods_paper/deck
PAPER=$BSS_DATA/Downloads/RECOMB_paper_template/figs
keep() { for f in "$@"; do [ -f "$f" ] && [ ! -f "${f%.*}_$TAG.${f##*.}" ] && cp -p "$f" "${f%.*}_$TAG.${f##*.}"; done; return 0; }
if step 1; then
echo "[$(date +%H:%M)] 1/9 keep the current outputs (_$TAG)"
keep $HP/cv_full_per_protein.csv $HP/sticker_spacer_by_ff_summary.txt $HP/sticker_definitions_summary.txt \
     $HP/context_model_by_ff_summary.txt $HP/context_proline_full.txt $HP/context_size_control_summary.txt \
     $HP/context_neighbourhood_typed.csv $HP/context_neighbourhood_null_typed.csv sim_null_v3_summary.txt \
     sim_null_v3_by_type.csv $FIG/fig3_arrangement_values.txt arrangement_test_summary.txt hub_overlap_by_ff_summary.txt
fi
if step 2; then
echo "[$(date +%H:%M)] 2/9 dataset with all labelled Mpipi-Recharged runs"
HUB_LABELS=v3 python3 $HP/hub_dataset.py $HP/data_v3
cp $HP/experimental_occupancy.csv $HP/data_v3/
fi
if step 3; then
echo "[$(date +%H:%M)] 3/9 cross-validation and transfer"
(cd $HP && $PY cv_full.py --data data_v3 > cv_full_v3_run.log 2>&1)
fi
if step 4; then
echo "[$(date +%H:%M)] 4/9 stickers and spacers by label set, sticker definitions"
(cd $HP && $PY sticker_spacer_by_ff.py data_v3 > /dev/null)
fi
if step 5; then
echo "[$(date +%H:%M)] 5/9 residue-neutral nulls on the sequences simulated with all three force fields"
OMP_NUM_THREADS=1 $PY sim_null_v3.py run 8 > sim_null_v3_run_$TAG.log 2>&1
OMP_NUM_THREADS=1 $PY sim_null_gg_pairs.py run 8 > sim_null_gg_pairs_run.log 2>&1
$PY sim_null_v3.py summarise > /dev/null
(cd $HP && $PY sticker_definitions.py > /dev/null)
fi
if step 6; then
echo "[$(date +%H:%M)] 6/9 arrangement test tables (the pre-registered result stays in arrangement_test_summary_$TAG.txt)"
$PY arrangement_test.py > arrangement_test_rerun_$TAG.log 2>&1 || echo "  arrangement_test.py failed, see log"
fi
if step 7; then
echo "[$(date +%H:%M)] 7/9 context models, proline, size control, neighbourhoods"
(cd $HP && $PY context_model_by_ff.py > context_model_by_ff_run.log 2>&1 \
        && $PY context_proline_full.py > /dev/null 2>&1 \
        && $PY context_size_control.py > context_size_control_run.log 2>&1 \
        && $PY context_neighbourhood_typed.py both > /dev/null 2>&1)
$PY hub_overlap_by_ff.py > /dev/null   # shared hubs between force fields (Results 3.3 sentence)
fi
if step 8; then
echo "[$(date +%H:%M)] 8/9 figures"
(cd $FIG && $PY fig2_flicker_null.py > fig2_run.log 2>&1 && $PY fig3_arrangement.py > fig3_run.log 2>&1)
cp $FIG/fig2_flicker_null.pdf $FIG/fig3_arrangement.pdf $PAPER/
fi
if step 9; then
echo "[$(date +%H:%M)] 9/9 deck charts (panel crops of Fig. 3 are re-cut separately)"
(cd $DECK && $PY make_enrichment.py > /dev/null 2>&1 && $PY make_sticker_definitions.py > /dev/null 2>&1 \
          && $PY make_context_charts.py > /dev/null 2>&1)
fi
echo "[$(date +%H:%M)] done"
