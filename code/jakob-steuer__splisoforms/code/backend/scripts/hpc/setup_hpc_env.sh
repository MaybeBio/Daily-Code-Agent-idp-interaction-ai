#!/usr/bin/env bash
#
# setup_hpc_env.sh
# ----------------
# Self-contained installer for the neoantigen prediction toolchain on a
# no-sudo HPC node (A40/L40 GPUs, SLURM-managed). Safe to run interactively
# on the headnode. Builds:
#
#   $HPC_ROOT/miniforge3/     # conda (bootstrapped if missing)
#   $HPC_ROOT/envs/neo/       # conda env with python, pytorch+CUDA,
#                             #   mhcflurry, pyarrow, ...
#   $HPC_ROOT/tools/bigmhc/   # cloned/copied BigMHC tree (with weights)
#   $HPC_ROOT/tools/netmhcpan-4.2/   # YOU drop the DTU x86_64 tarball here
#
# What this script does:
#   1. Bootstraps Miniforge into $HPC_ROOT/miniforge3 if conda is missing.
#   2. Detects CUDA via nvidia-smi and picks the matching pytorch-cuda version.
#   3. Creates the `neo` env with python 3.12, pytorch (CUDA build), tcsh,
#      numpy/pandas/scipy, pyarrow/biopython, and pip-installs mhcflurry 2.2.1.
#   4. Runs `mhcflurry-downloads fetch` (~3 GB, login-node only — needs egress).
#   5. Smoke-tests: import torch + cuda check, import mhcflurry, import bigmhc.
#
# What it does NOT do:
#   - Install NetMHCpan-4.2 (license-locked DTU tarball; manual scp + extract).
#   - Copy BigMHC source (you scp the local tools/bigmhc/ tree first).
#
# Usage on the HPC headnode:
#
#   1. From laptop, scp this script + bigmhc tree:
#        HPC=user@hpc-login
#        DST=/path/to/neoantigen_predictions
#        scp backend/scripts/hpc/setup_hpc_env.sh "$HPC:$DST/"
#        rsync -av tools/bigmhc/ "$HPC:$DST/tools/bigmhc/"
#
#   2. On the HPC:
#        cd /path/to/neoantigen_predictions
#        bash setup_hpc_env.sh
#
#   3. NetMHCpan (manual): from DTU, pick
#      "Version 4.2c static → Linux" (statically-linked x86_64 binary).
#      Then on the HPC:
#        cd $HPC_ROOT/tools/
#        tar xf ~/netMHCpan-4.2cstatic.Linux.tar
#        mv netMHCpan-4.2 netmhcpan-4.2
#      Edit netmhcpan-4.2/netMHCpan and set the absolute NMHOME.

set -euo pipefail

# ──────────────────────────────────────────────────────────────────────────
#  CONFIG
# ──────────────────────────────────────────────────────────────────────────
HPC_ROOT="${HPC_ROOT:-/path/to/neoantigen_predictions}"
ENV_NAME="${ENV_NAME:-neo}"
PYTHON_VER="${PYTHON_VER:-3.12}"
MHCFLURRY_VER="${MHCFLURRY_VER:-2.2.1}"

MINIFORGE_DIR="${HPC_ROOT}/miniforge3"
ENVS_DIR="${HPC_ROOT}/envs"
TOOLS_DIR="${HPC_ROOT}/tools"
BIGMHC_DIR="${TOOLS_DIR}/bigmhc"
NETMHC_DIR="${TOOLS_DIR}/netmhcpan-4.2"

mkdir -p "$HPC_ROOT" "$ENVS_DIR" "$TOOLS_DIR"

# ──────────────────────────────────────────────────────────────────────────
#  1. Bootstrap miniforge if conda is missing
# ──────────────────────────────────────────────────────────────────────────
if [[ -x "${MINIFORGE_DIR}/bin/conda" ]]; then
    echo "[1/5] conda already at ${MINIFORGE_DIR}"
elif command -v conda &>/dev/null; then
    echo "[1/5] system conda found at $(command -v conda) — using it"
    MINIFORGE_DIR="$(conda info --base)"
else
    echo "[1/5] bootstrapping Miniforge into ${MINIFORGE_DIR} ..."
    INSTALLER="/tmp/miniforge_$$.sh"
    curl -fsSL \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
        -o "$INSTALLER"
    bash "$INSTALLER" -b -p "$MINIFORGE_DIR"
    rm -f "$INSTALLER"
fi

# shellcheck disable=SC1091
source "${MINIFORGE_DIR}/etc/profile.d/conda.sh"
conda config --set channel_priority strict
conda config --add channels conda-forge

# ──────────────────────────────────────────────────────────────────────────
#  2. Detect CUDA via nvidia-smi → pick pytorch-cuda meta version
# ──────────────────────────────────────────────────────────────────────────
echo "[2/5] detecting CUDA ..."
if command -v nvidia-smi &>/dev/null; then
    CUDA_FULL=$(nvidia-smi 2>/dev/null | grep -oE 'CUDA Version:[[:space:]]+[0-9]+\.[0-9]+' | head -1 | awk '{print $3}')
    CUDA_MAJOR="${CUDA_FULL%%.*}"
    CUDA_MINOR="${CUDA_FULL#*.}"
    echo "    detected CUDA ${CUDA_FULL}"
else
    echo "    nvidia-smi not on PATH — install login-node may not see GPUs."
    echo "    Defaulting to CUDA 12.4. Override with CUDA_FULL=X.Y env var."
    CUDA_FULL="${CUDA_FULL:-12.4}"
    CUDA_MAJOR=12
    CUDA_MINOR=4
fi

case "${CUDA_MAJOR}.${CUDA_MINOR}" in
    12.4|12.5|12.6) PT_CUDA="12.4" ;;
    12.1|12.2|12.3) PT_CUDA="12.1" ;;
    11.8|11.9)      PT_CUDA="11.8" ;;
    *)
        echo "    unrecognized CUDA ${CUDA_FULL}; falling back to pytorch-cuda=12.4"
        PT_CUDA="12.4"
        ;;
esac
echo "    will install pytorch-cuda=${PT_CUDA}"

# ──────────────────────────────────────────────────────────────────────────
#  3. Create / update the conda env
# ──────────────────────────────────────────────────────────────────────────
ENV_PREFIX="${ENVS_DIR}/${ENV_NAME}"
echo "[3/5] building env at ${ENV_PREFIX} ..."

if [[ -d "$ENV_PREFIX" ]]; then
    echo "    env already exists — pip-only sync"
else
    conda create -y -p "$ENV_PREFIX" \
        -c pytorch -c nvidia -c conda-forge \
        "python=${PYTHON_VER}" \
        "pytorch>=2.2" "pytorch-cuda=${PT_CUDA}" \
        tcsh \
        numpy pandas scipy scikit-learn \
        biopython requests tqdm psutil \
        pyarrow \
        pip
fi

# pip installs — anything not reliably packaged on conda-forge
# conda activate scripts (e.g. MKL) reference unset vars; disable nounset briefly
set +u
conda activate "$ENV_PREFIX"
set -u
pip install --upgrade pip
pip install "mhcflurry==${MHCFLURRY_VER}"

# ──────────────────────────────────────────────────────────────────────────
#  4. Fetch MHCflurry models (~3 GB; needs internet egress)
# ──────────────────────────────────────────────────────────────────────────
echo "[4/5] fetching MHCflurry models ..."
if mhcflurry-downloads info 2>/dev/null | grep -q "models_class1_presentation.*[*]"; then
    echo "    models already fetched."
else
    mhcflurry-downloads fetch
fi

# ──────────────────────────────────────────────────────────────────────────
#  5. Smoke tests
# ──────────────────────────────────────────────────────────────────────────
echo "[5/5] smoke tests ..."
python - <<'PY'
import sys, torch
print(f"  python:       {sys.version.split()[0]}")
print(f"  torch:        {torch.__version__}")
print(f"  cuda built:   {torch.version.cuda}")
print(f"  cuda visible: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"  gpu count:    {torch.cuda.device_count()}")
    for i in range(torch.cuda.device_count()):
        print(f"    [{i}] {torch.cuda.get_device_name(i)}")
PY

python -c "
from mhcflurry import Class1PresentationPredictor
p = Class1PresentationPredictor.load()
print(f'  mhcflurry:    OK ({len(p.supported_alleles)} alleles)')
"

# BigMHC: verify the tree is present
if [[ -d "${BIGMHC_DIR}/src" && -d "${BIGMHC_DIR}/models" ]]; then
    echo "  bigmhc tree:  OK at ${BIGMHC_DIR}"
    python -c "
import sys, importlib.util
sys.path.insert(0, '${BIGMHC_DIR}/src')
spec = importlib.util.find_spec('bigmhc')
print(f'  bigmhc import: {\"OK\" if spec else \"MISSING\"}')
"
else
    echo "  ⚠ BigMHC NOT at ${BIGMHC_DIR}"
    echo "     rsync it from your laptop:"
    echo "       rsync -av tools/bigmhc/ <hpc-host>:${BIGMHC_DIR}/"
fi

# NetMHCpan: verify
if [[ -x "${NETMHC_DIR}/netMHCpan" ]]; then
    echo "  netMHCpan:    found at ${NETMHC_DIR}/netMHCpan"
else
    echo "  ⚠ NetMHCpan NOT at ${NETMHC_DIR}/netMHCpan"
    echo "     1. Download netMHCpan-4.2c.Linux_x86_64.tar from DTU."
    echo "     2. tar xf netMHCpan-4.2c.Linux_x86_64.tar -C ${TOOLS_DIR}/"
    echo "        and rename the extracted dir to netmhcpan-4.2/"
    echo "     3. Edit ${NETMHC_DIR}/netMHCpan and set:"
    echo "          setenv NMHOME ${NETMHC_DIR}"
fi

# ──────────────────────────────────────────────────────────────────────────
#  Summary
# ──────────────────────────────────────────────────────────────────────────
cat <<EOF

════════════════════════════════════════════════════════════════════════
  Environment ready.
════════════════════════════════════════════════════════════════════════
  Activate:
    source ${MINIFORGE_DIR}/etc/profile.d/conda.sh
    conda activate ${ENV_PREFIX}

  Required env vars at runtime:
    export NETMHCPAN_BIN=${NETMHC_DIR}/netMHCpan
    export BIGMHC_DIR=${BIGMHC_DIR}

  Next: scp hpc_input.tsv + prediction scripts, then sbatch neoantigen.sbatch
EOF
