#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# The SLURM job that ./noodle launch submits. You never run this yourself.
#
# It is submitted as a job array: every array element runs one job of the run (a few
# tasks one after another, or a single long one),
#     job number = <offset> + array index      (offset is 0 unless the run has more
#                                                jobs than the cluster's array limit)
# ./noodle launch adds the job name, time limit, array range and SLURM's own output file.
# Arguments: <run folder> <offset>
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# 10 minutes before the time limit, SLURM sends SIGUSR1: the running task then saves a
# checkpoint and the job stops cleanly, so a job that runs out of time loses nothing.
#SBATCH --signal=B:USR1@600

# One core per job: HOOMD's CPU build runs a simulation on one core, and a run gets its
# speed from many jobs side by side (a brush too). Print progress straight into the
# output file.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1

[ -g . ] && umask 002  # shared group folder: files group-writable (see ./noodle)

JOB=$(( $2 + SLURM_ARRAY_TASK_ID ))
# Everything this job prints goes to jobs/job_0007.out (for job 7)
exec >> "$1/jobs/job_$(printf %04d "$JOB").out" 2>&1
echo "=== job $JOB   SLURM job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}   $(hostname)   $(date) ==="

# The tool's own Python environment (made by setup.sh): no conda activate needed.
# exec: Python replaces this shell, so it is the process that receives SLURM's signals.
exec ./env/bin/python src/simulate.py "$1" "$JOB"
