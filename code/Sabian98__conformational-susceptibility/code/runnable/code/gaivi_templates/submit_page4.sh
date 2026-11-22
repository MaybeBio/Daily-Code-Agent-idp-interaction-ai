#!/bin/bash
# =============================================================================
# submit_page4.sh — PAGE4 phosphoform runs on GAIVI (CALVADOS-2 and Mpipi-Recharged)
#
# Usage (activate the calvados env first; jobs inherit it, as for submit_calvados.sh / submit_openmpipi.sh):
#   bash submit_page4.sh calvados                 # the 36 CALVADOS-2 runs in calvados/runs.tsv
#   bash submit_page4.sh mpipi                    # the 46 Mpipi-Recharged runs in mpipi/runs.tsv
#   bash submit_page4.sh both
#   options: --sets heat            only the 278-320 K series (or smFRET, SAXS; comma-separated)
#            --outroot DIR          default /data/taseefr/page4 (CALVADOS-2), /general/taseefr/page4 (Mpipi);
#                                   with 'both', DIR/calvados and DIR/mpipi so the two models never share a folder
#            --time-factor N        run length in units of the BENDER run (default 3 = 300 ns for PAGE4)
#            --exclude NODES        keep jobs off these nodes (default GPU4,GPU6,GPU8, which misbehave; GPU6 runs jobs without
#                                   writing to /data, GPU8 fails every job at once; comma-separated; --exclude '' for none)
#
# Output : <outroot>/<T>K/<name>/        (<name>.dcd + top.pdb + config/components for CALVADOS-2)
#          <outroot>/<T>K/data/<name>/   (conf_prop.csv)
# Logs   : <outroot>/logs/<model>_a<attempt>_<task>.out/.err, with pending_<model>_a<attempt>_*.tsv giving the run
#          of each task
# Only runs without a finished result are submitted (CALVADOS-2: conf_prop.csv; Mpipi: <name>.dcd). After each
# array finishes, a follow-up job resubmits the leftovers (up to MAX_ATTEMPTS rounds), as in submit_calvados.sh.
# Mpipi runs resume from their last saved segment.
# =============================================================================
set -euo pipefail

MAX_ATTEMPTS=3
SCRIPT_PATH="$(realpath "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
MODEL="${1:-}"; shift || true
SETS="heat,smFRET,SAXS"
OUT_CAL="/data/taseefr/page4"
OUT_MP="/general/taseefr/page4"
TF=3
ATTEMPT=1
EXCLUDE="GPU4,GPU6,GPU8"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sets)        SETS="$2"; shift 2 ;;
        --outroot)     OUT_CAL="${2%/}"; OUT_MP="${2%/}"; OUTROOT_SET=1; shift 2 ;;
        --time-factor) TF="$2"; shift 2 ;;
        --attempt)     ATTEMPT="$2"; shift 2 ;;
        --exclude)     EXCLUDE="$2"; shift 2 ;;
        *) echo "unknown option $1"; exit 1 ;;
    esac
done
[[ "$MODEL" =~ ^(calvados|mpipi|both)$ ]] || { echo "usage: bash submit_page4.sh calvados|mpipi|both [options]"; exit 1; }
if [[ "$MODEL" == both && -n "${OUTROOT_SET:-}" ]]; then
    OUT_CAL="$OUT_CAL/calvados"; OUT_MP="$OUT_MP/mpipi"
fi

submit() {   # $1 model dir, $2 outroot
    local M="$1" OUT="$2" RUNS="$SCRIPT_DIR/$1/runs.tsv"
    mkdir -p "$OUT/logs"
    # Runs that still have a queued or running job, from any earlier submission of this model, are skipped, so a run is
    # never simulated twice (a second job on a running simulation corrupts its trajectory, seen once in GAIVI_PAGE4).
    # Each array stores its run list in --comment, which maps its task numbers back to runs; tasks that are only
    # completing after scancel (CG) do not count. Jobs submitted by older versions carry no comment and are not seen.
    local JN; [[ "$M" == calvados ]] && JN=PAGE4_cal || JN=PAGE4_mp
    local ACTIVE; ACTIVE="$(mktemp)"
    if command -v squeue >/dev/null 2>&1; then
        { squeue -u "$USER" -h -t PENDING,RUNNING,CONFIGURING -n "$JN" -o "%k|%K" 2>/dev/null || true; } |
        while IFS='|' read -r L IDS; do
            [[ -f "$L" ]] || continue
            IDS="${IDS%%\%*}"; IDS="${IDS//[\[\]]/}"
            for P in ${IDS//,/ }; do
                A="${P%-*}"; B="${P#*-}"
                [[ "$A$B" =~ ^[0-9]+$ ]] || continue
                for (( i = A; i <= B; i++ )); do sed -n "${i}p" "$L" | cut -f1,4; done
            done
        done > "$ACTIVE"
    fi
    : > "$ACTIVE.skip"
    local LIST="$OUT/logs/pending_${M}_a${ATTEMPT}_$(date +%Y%m%d_%H%M%S).tsv"
    : > "$LIST"
    tail -n +2 "$RUNS" | while IFS=$'\t' read -r NAME FORM SET TEMP IONIC REST; do
        [[ ",$SETS," == *",$SET,"* ]] || continue
        if grep -qxF "$NAME"$'\t'"$TEMP" "$ACTIVE"; then echo "$NAME" >> "$ACTIVE.skip"; continue; fi
        if [[ "$M" == calvados ]]; then DONE="$OUT/${TEMP}K/data/$NAME/conf_prop.csv"; else DONE="$OUT/${TEMP}K/$NAME/$NAME.dcd"; fi
        [[ -f "$DONE" ]] || printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$NAME" "$FORM" "$SET" "$TEMP" "$IONIC" "$REST" >> "$LIST"
    done
    local NSKIP; NSKIP=$(wc -l < "$ACTIVE.skip" | tr -d ' ')
    rm -f "$ACTIVE" "$ACTIVE.skip"
    local N; N=$(wc -l < "$LIST" | tr -d ' ')
    echo "$M: $N runs pending (sets $SETS, attempt $ATTEMPT/$MAX_ATTEMPTS) -> $OUT"
    if (( NSKIP > 0 )); then echo "  $NSKIP runs skipped because a job for them is still queued or running"; fi
    (( N == 0 )) && { rm -f "$LIST"; return; }
    local JOB_ID
    if [[ "$M" == calvados ]]; then
        JOB_ID=$(sbatch --parsable --job-name="PAGE4_cal" --comment="$LIST" --array="1-${N}" ${EXCLUDE:+--exclude="$EXCLUDE"} \
            --output="$OUT/logs/calvados_a${ATTEMPT}_%a.out" --error="$OUT/logs/calvados_a${ATTEMPT}_%a.err" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
LINE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$LIST")
NAME=\$(echo "\$LINE" | cut -f1); FORM=\$(echo "\$LINE" | cut -f2); SET=\$(echo "\$LINE" | cut -f3)
TEMP=\$(echo "\$LINE" | cut -f4); IONIC=\$(echo "\$LINE" | cut -f5); PH=\$(echo "\$LINE" | cut -f6); RES=\$(echo "\$LINE" | cut -f7)
RUNDIR="$OUT/\${TEMP}K/\$NAME"; DATADIR="$OUT/\${TEMP}K/data/\$NAME"
if [[ -f "\$DATADIR/conf_prop.csv" ]]; then echo "\$NAME @ \${TEMP}K done"; exit 0; fi
echo "\$NAME (\$FORM, \$SET) @ \${TEMP}K, I \$IONIC, pH \$PH on \$(hostname)"
python "$SCRIPT_DIR/calvados/prepare_page4.py" --name "\$NAME" --form "\$FORM" --residues "\$RES" --temp "\$TEMP" --ionic "\$IONIC" --ph "\$PH" --outdir "$OUT" --time-factor $TF
python "\$RUNDIR/run.py" --path "\$RUNDIR"
if [ -d "\$RUNDIR" ]; then
    find "\$RUNDIR" -type f ! \( -name "*.dcd" -o -name "top.pdb" -o -name "*.yaml" -o -name "*.fasta" \) -delete
fi
if [ -d "\$DATADIR" ]; then
    find "\$DATADIR" -type f ! -name "conf_prop.csv" -delete
fi
EOF
        )
    else
        JOB_ID=$(sbatch --parsable --job-name="PAGE4_mp" --comment="$LIST" --array="1-${N}" ${EXCLUDE:+--exclude="$EXCLUDE"} \
            --output="$OUT/logs/mpipi_a${ATTEMPT}_%a.out" --error="$OUT/logs/mpipi_a${ATTEMPT}_%a.err" <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
LINE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$LIST")
NAME=\$(echo "\$LINE" | cut -f1); FORM=\$(echo "\$LINE" | cut -f2); SET=\$(echo "\$LINE" | cut -f3)
TEMP=\$(echo "\$LINE" | cut -f4); IONIC=\$(echo "\$LINE" | cut -f5); PHOS=\$(echo "\$LINE" | cut -f6)
RUNDIR="$OUT/\${TEMP}K/\$NAME"; DATADIR="$OUT/\${TEMP}K/data/\$NAME"
if [[ -f "\$RUNDIR/\$NAME.dcd" ]]; then echo "\$NAME @ \${TEMP}K done"; exit 0; fi
echo "\$NAME (\$FORM, \$SET) @ \${TEMP}K, I \$IONIC, phospho \$PHOS on \$(hostname)"
python "$SCRIPT_DIR/mpipi/run_openmpipi_page4.py" --name "\$NAME" --form "\$FORM" --phospho "\$PHOS" --temp "\$TEMP" --ionic "\$IONIC" --outdir "$OUT" --time-factor $TF --seed "\$SLURM_ARRAY_TASK_ID"
# Clean up only finished runs, so unfinished ones keep state.xml / segments to resume from (as submit_openmpipi.sh)
if [[ -f "\$DATADIR/conf_prop.csv" ]]; then
    find "\$DATADIR" -type f ! -name "conf_prop.csv" -delete
    find "\$RUNDIR" -type f ! \( -name "*.dcd" -o -name "top.pdb" \) -delete
fi
EOF
        )
    fi
    echo "  array $JOB_ID ($N tasks)"

    # Follow-up: once the array has ended, resubmit whatever is still missing (as submit_calvados.sh)
    if (( ATTEMPT < MAX_ATTEMPTS )); then
        local RID
        RID=$(sbatch --parsable --job-name="PAGE4_resub_${M}" \
            --output="$OUT/logs/resubmit_${M}_a${ATTEMPT}_%j.out" --error="$OUT/logs/resubmit_${M}_a${ATTEMPT}_%j.err" \
            --dependency="afterany:${JOB_ID}" --ntasks=1 --mem=1G ${EXCLUDE:+--exclude="$EXCLUDE"} \
            --wrap="bash '${SCRIPT_PATH}' ${M} --sets '${SETS}' --outroot '${OUT}' --time-factor ${TF} --attempt $((ATTEMPT + 1)) --exclude '${EXCLUDE}'")
        echo "  resubmitter $RID (retries leftovers after $JOB_ID)"
    else
        echo "  last attempt, no further resubmission; check $OUT/logs for failures"
    fi
}

[[ "$MODEL" == calvados || "$MODEL" == both ]] && submit calvados "$OUT_CAL"
[[ "$MODEL" == mpipi || "$MODEL" == both ]] && submit mpipi "$OUT_MP"
echo "Done."
