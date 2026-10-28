#!/bin/bash

# NUM_JOBS jobs splitting every system of every dataset (main, sabari, exp,
# known; see DATASETS in run_ddm_het.py) into contiguous chunks. Each job
# loops over its indices and runs run_ddm_het.py once per system.
#
# A system whose results/<dataset>/<base_name>.npz already exists is skipped
# unless its trajectory has grown, so resubmitting only picks up what is
# missing.

PYTHON_SCRIPT_NAME="run_ddm_het.py"
NUM_JOBS=400

if [ ! -f "$PYTHON_SCRIPT_NAME" ]; then
    echo "Error: Python script '$PYTHON_SCRIPT_NAME' not found."
    exit 1
fi

mkdir -p logs

# ask the python script how many systems there are
TOTAL=$(python $PYTHON_SCRIPT_NAME --count)
if [ -z "$TOTAL" ] || [ "$TOTAL" -lt 1 ]; then
    echo "Error: no systems found. Check DATASETS in $PYTHON_SCRIPT_NAME."
    exit 1
fi

echo "Splitting $TOTAL systems across $NUM_JOBS jobs..."

BASE=$((TOTAL / NUM_JOBS))
REM=$((TOTAL % NUM_JOBS))

start=0
for ((i=0; i<NUM_JOBS; i++)); do
    # first REM jobs get one extra system so the chunks cover all TOTAL
    size=$BASE
    if [ $i -lt $REM ]; then
        size=$((size + 1))
    fi
    if [ $size -eq 0 ]; then
        continue
    fi
    end=$((start + size - 1))

    echo "Submitting job $i for systems $start-$end..."

    qsub -N "het_${i}" -o logs/ -e logs/ -v JOB_START=$start,JOB_END=$end << 'EOF'

# --- PBS Directives ---
#PBS -lselect=1:ncpus=1:ngpus=0:mem=16gb
#PBS -lwalltime=24:00:0
#PBS -qcpu

# --- Job Setup ---
echo "Starting job task for systems: $JOB_START-$JOB_END"

# Load required modules
module load anaconda3

# Change to the directory where qsub was run
cd $PBS_O_WORKDIR

# Activate the Conda environment
source activate /storage/brno14-ceitec/shared/softmatter/shakespearem/.conda/envs/omm2

# --- Run Analysis: one system at a time within this job's slice ---
for idx in $(seq $JOB_START $JOB_END); do
    echo "  -> system $idx"
    python run_ddm_het.py $idx
done

echo "Finished job task for systems: $JOB_START-$JOB_END"

EOF

    start=$((end + 1))
done

echo "All jobs have been submitted to the queue."
echo "You can monitor the jobs with the 'qstat' command."
echo "When they have finished, run 'python run_ddm_het.py --collect' to write"
echo "heterotypic_<dataset>.csv. Logs are in logs/"
