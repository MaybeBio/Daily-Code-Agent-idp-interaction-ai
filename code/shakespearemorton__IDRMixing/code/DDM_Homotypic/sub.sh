#!/bin/bash

# One job per sequence directory. Each job runs DDM.py at every voxel size
# (run_ddm.py) and the per-frame two-peak count for Figure 1E
# (frac_two_peak.py). Directories are enumerated by run_ddm.py from
# ../../calv_comp/simulations_{training,validation_IDRome10,validation_IDRome90},
# sorted, so the index is stable between submissions.
#
# A sequence whose results_<S>/<name>.npz or peaks_<S>/<name>.npz already
# exists is skipped, so resubmitting only picks up what is missing.

for f in run_ddm.py frac_two_peak.py; do
    if [ ! -f "$f" ]; then
        echo "Error: Python script '$f' not found."
        exit 1
    fi
done

mkdir -p logs

# ask the python script how many sequence directories it found
NUM_JOBS=$(python run_ddm.py --count)
if [ -z "$NUM_JOBS" ] || [ "$NUM_JOBS" -lt 1 ]; then
    echo "Error: no sequence directories found. Check SIM_ROOT in run_ddm.py."
    exit 1
fi

echo "Submitting $NUM_JOBS individual jobs (indices 0-$((NUM_JOBS-1)))..."

for i in $(seq 0 $((NUM_JOBS-1))); do

    echo "Submitting job for index $i..."

    qsub -N "slab_${i}" -o logs/ -e logs/ -v JOB_INDEX=$i << 'EOF'

# --- PBS Directives ---
#PBS -lselect=1:ncpus=1:ngpus=0:mem=16gb
#PBS -lwalltime=24:00:0
#PBS -qcpu

# --- Job Setup ---
echo "Starting job task for index: $JOB_INDEX"

# Load required modules
module load anaconda3

# Change to the directory where qsub was run
cd $PBS_O_WORKDIR

# Activate the Conda environment
source activate /storage/brno14-ceitec/shared/softmatter/shakespearem/.conda/envs/omm2

# --- Run Analysis ---
python run_ddm.py $JOB_INDEX
python frac_two_peak.py $JOB_INDEX

echo "Finished job task for index: $JOB_INDEX"

EOF

done

echo "All jobs have been submitted to the queue."
echo "You can monitor the jobs with the 'qstat' command."
echo "When they have finished, run 'python run_ddm.py --collect' to write"
echo "slab_results_<S>.csv. Logs are in logs/"
