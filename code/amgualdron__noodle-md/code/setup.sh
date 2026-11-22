#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# One-time setup: builds the tool's own Python environment in ./env from
# src/environment.yml, then checks the physics (src/check_physics.py).
#
#   on the cluster:          sbatch setup.sh     (builds on a compute node; ~15-30 min;
#                                                  progress and result in setup.log)
#   on your own computer:    bash setup.sh
#
# After this you never need to load or activate conda again: ./noodle calls
# ./env/bin/python directly. At the end it prints two lines for your ~/.bashrc, so
# that noodle works from any folder and Tab completes commands and run names.
# To rebuild from scratch, delete the env folder and run this again.
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH --job-name=noodle-setup
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=setup.log

CONDA_MODULE=anaconda   # name of your cluster's conda module (module avail conda)

# sbatch runs a copy of this script, so go to the folder it was submitted from
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}" || exit 1
[ -g . ] && umask 002  # shared group folder: env/ usable by the whole group (see ./noodle)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

if ! command -v conda > /dev/null; then
    module load "$CONDA_MODULE" || { echo "Could not load the module '$CONDA_MODULE'. Edit CONDA_MODULE in setup.sh."; exit 1; }
fi
if [ -d env ]; then
    echo "./env already exists. Delete it first to rebuild: rm -rf env"
    exit 1
fi

echo "Building ./env from src/environment.yml ($(date))"
if command -v mamba > /dev/null; then
    mamba env create --prefix ./env --file src/environment.yml || exit 1
else
    conda env create --prefix ./env --file src/environment.yml || exit 1
fi

echo
./env/bin/python -c "import hoomd; print('HOOMD-blue', hoomd.version.version, '| GPU support:', hoomd.version.gpu_enabled)" || exit 1
echo
./env/bin/python src/check_physics.py || exit 1
echo
echo "Setup finished ($(date))."
echo
echo "To type noodle from any folder, with Tab completing commands and run names, add these"
echo "two lines to the end of your ~/.bashrc (then open a new terminal):"
echo
echo "    alias noodle='$PWD/noodle'"
echo "    source '$PWD/noodle-completion.bash'"
echo
echo "Then start a run with:  noodle new <run name>      (noodle help lists every command)"
