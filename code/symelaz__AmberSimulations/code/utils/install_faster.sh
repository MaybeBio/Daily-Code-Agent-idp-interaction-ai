#!/bin/bash
# install_faster.sh - build the conda environment for the --faster engine
# (NVIDIA cuPotential OpenMM plugin, https://github.com/NVIDIA/cuPotential).
#
# Usage:  bash utils/install_faster.sh [--env openmm_faster] [--src ~/software/cuPotential]
#                                      [--driver-cuda 13.3] [--arch "80;90"]
#
#   --env          conda environment to create (default: openmm_faster)
#   --src          where to clone cuPotential (default: ~/software/cuPotential)
#   --driver-cuda  CUDA version the GPU driver supports (default: read from nvidia-smi;
#                  REQUIRED on a machine without GPUs, e.g. a cluster login node:
#                  run nvidia-smi on a GPU node and pass its "CUDA Version")
#   --arch         GPU compute capabilities to compile for (default: this machine's
#                  GPUs, else "80;86;89;90"; A100 = 80, RTX 4090 = 89, H100 = 90)
#
# What it does, and why each step is there:
#  1. Clones cuPotential and downloads its prebuilt runtime (libcupo.so is not in
#     git). The runtime fixes the CUDA toolkit version the plugin must be built
#     with (runtime_build.json, "cuda_min"); NVIDIA's blog post says 13.1, but
#     release 0.3.5 needs 13.4. Read, not assumed.
#  2. Creates the conda env: OpenMM 8.4 + that CUDA toolkit + compilers (conda-forge).
#  3. Builds the plugin against the env's OpenMM and installs it as a wheel.
#  4. Driver check. If the driver supports an OLDER CUDA minor than the toolkit
#     (e.g. driver 13.3, toolkit 13.4), stock OpenMM fails at the first Context
#     with CUDA_ERROR_UNSUPPORTED_PTX_VERSION: its kernels are compiled at run
#     time by NVRTC, and the driver rejects PTX newer than itself. cuPotential's
#     own kernels are precompiled and are fine. Fix: an NVRTC matching the driver,
#     preloaded automatically on `conda activate` (activate.d hook). Remove it
#     once the driver is updated: rm $CONDA_PREFIX/etc/conda/activate.d/nvrtc_driver.sh
#  5. Next: python utils/check_faster.py on a GPU node.
set -euo pipefail

ENV=openmm_faster; SRC="$HOME/software/cuPotential"; DRIVER_CUDA=""; ARCH=""
while [ $# -gt 0 ]; do
    case "$1" in
        --env) ENV="$2"; shift 2 ;;
        --src) SRC="$2"; shift 2 ;;
        --driver-cuda) DRIVER_CUDA="$2"; shift 2 ;;
        --arch) ARCH="$2"; shift 2 ;;
        *) echo "unknown option $1"; exit 1 ;;
    esac
done
command -v conda >/dev/null || { echo "ERROR: conda not found (load/initialise Anaconda first)"; exit 1; }
eval "$(conda shell.bash hook)"

# ---- driver ----
trap 'echo "ERROR: install_faster.sh failed at line $LINENO: $BASH_COMMAND" >&2' ERR
if [ -z "$DRIVER_CUDA" ]; then
    # Ask the driver itself (newer nvidia-smi no longer prints "CUDA Version:" in its header)
    DRIVER_CUDA=$(python3 -c "
import ctypes
v = ctypes.c_int()
ctypes.CDLL('libcuda.so.1').cuDriverGetVersion(ctypes.byref(v))
print(f'{v.value // 1000}.{(v.value % 1000) // 10}')" 2>/dev/null || true)
fi
if [ -z "$DRIVER_CUDA" ] && command -v nvidia-smi >/dev/null; then
    DRIVER_CUDA=$(nvidia-smi -q 2>/dev/null | grep -m1 -oP 'CUDA (UMD )?Version\s*:\s*\K[0-9]+\.[0-9]+' || true)
fi
if [ -z "$DRIVER_CUDA" ]; then
    echo "ERROR: no GPU driver visible here. Run nvidia-smi on a GPU node and pass"
    echo "       --driver-cuda <its CUDA Version>, e.g. --driver-cuda 13.3"; exit 1
fi
if [ -z "$ARCH" ]; then
    ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | tr -d '.' | sort -u | paste -sd';' || true)
    ARCH=${ARCH:-"80;86;89;90"}
fi
echo "== driver supports CUDA $DRIVER_CUDA; building for compute capabilities $ARCH"

# ---- 1. source + prebuilt runtime ----
if [ ! -d "$SRC/.git" ]; then
    mkdir -p "$(dirname "$SRC")"
    git clone --depth 1 https://github.com/NVIDIA/cuPotential.git "$SRC"
fi
cd "$SRC"
echo "== cuPotential $(cat VERSION) ($(git log -1 --format=%h))"
[ -f build/runtime_build.json ] || python3 fetch_cupo
CUDA=$(python3 -c "import json; print(json.load(open('build/runtime_build.json'))['cuda_min'])")
echo "== cuPotential runtime needs CUDA toolkit $CUDA"
if [ "${DRIVER_CUDA%%.*}" != "${CUDA%%.*}" ]; then
    echo "ERROR: driver CUDA $DRIVER_CUDA and runtime CUDA $CUDA differ in the major version;"
    echo "       the GPU driver has to be updated (ask the cluster admins)."; exit 1
fi

# ---- 2. environment ----
if conda env list | awk '{print $1}' | grep -qx "$ENV"; then
    echo "== conda env $ENV exists, updating CUDA toolkit to $CUDA"
    conda install -y -n "$ENV" -c conda-forge --override-channels "cuda-version=$CUDA" "cuda-toolkit=$CUDA"
else
    conda create -y -n "$ENV" -c conda-forge --override-channels \
        python=3.12 openmm=8.4.0 "cuda-version=$CUDA" "cuda-toolkit=$CUDA" \
        cmake make gcc_linux-64=13 gxx_linux-64=13 swig pip setuptools wheel numpy
fi
set +u; conda activate "$ENV"; set -u   # conda/CUDA activate.d scripts are not set -u safe
python -m pip install "pybind11>=2.12" nvidia-mathdx

# ---- 3. plugin ----
sed -i 's|^OPENMM=/opt/openmm$|OPENMM="${CONDA_PREFIX:?Activate the conda environment first}"|' plugins/openmm/run_build
export CUDACXX="$CONDA_PREFIX/bin/nvcc" CUDAHOSTCXX="$CXX" CUDA_ARCHITECTURES="$ARCH"
rm -rf plugins/openmm/build plugins/openmm/dist
./plugins/openmm/run_build > plugins/openmm/build_install_faster.log 2>&1 \
    || { echo "ERROR: plugin build failed, see $SRC/plugins/openmm/build_install_faster.log"; exit 1; }
python -m pip wheel --no-deps --no-build-isolation --wheel-dir plugins/openmm/dist ./plugins/openmm
# not --force-reinstall: that would also reinstall OpenMM from PyPI over conda's
python -m pip uninstall -y nv-openmm-plugin >/dev/null 2>&1 || true
python -m pip install plugins/openmm/dist/*.whl
(cd /tmp && python -c "import nv_openmm_plugin") && echo "== plugin installed"

# ---- 4. driver older than the toolkit: NVRTC matching the driver ----
HOOK="$CONDA_PREFIX/etc/conda/activate.d/nvrtc_driver.sh"
UNHOOK="$CONDA_PREFIX/etc/conda/deactivate.d/nvrtc_driver.sh"
rm -f "$HOOK" "$UNHOOK"
if python -c "import sys; d, c = (tuple(map(int, v.split('.'))) for v in sys.argv[1:]); sys.exit(0 if d < c else 1)" "$DRIVER_CUDA" "$CUDA"; then
    echo "== driver CUDA $DRIVER_CUDA < toolkit $CUDA: installing NVRTC $DRIVER_CUDA for stock OpenMM kernels"
    D="$CONDA_PREFIX/opt/nvrtc_driver"; rm -rf "$D"; mkdir -p "$D"
    python -m pip download -q --no-deps --only-binary=:all: -d "$D/whl" "nvidia-cuda-nvrtc==${DRIVER_CUDA}.*"
    python -m zipfile -e "$D"/whl/*.whl "$D/x"
    LIB=$(dirname "$(find "$D/x" -name 'libnvrtc.so.*' | head -1)")
    mkdir -p "$(dirname "$HOOK")" "$(dirname "$UNHOOK")"
    cat > "$HOOK" <<EOF
# install_faster.sh: GPU driver supports CUDA $DRIVER_CUDA, toolkit is $CUDA. OpenMM compiles its
# kernels at run time with NVRTC; use the NVRTC that matches the driver. Delete this file
# (and deactivate.d/nvrtc_driver.sh) once the driver supports CUDA $CUDA.
export _NVRTC_OLD_LD_PRELOAD="\${LD_PRELOAD:-}" _NVRTC_OLD_LD_LIBRARY_PATH="\${LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$LIB:\${LD_LIBRARY_PATH:-}"
export LD_PRELOAD="$LIB/$(basename "$(ls "$LIB"/libnvrtc.so.* | grep -v builtins | head -1)")"
EOF
    cat > "$UNHOOK" <<'EOF'
export LD_PRELOAD="${_NVRTC_OLD_LD_PRELOAD:-}" LD_LIBRARY_PATH="${_NVRTC_OLD_LD_LIBRARY_PATH:-}"
unset _NVRTC_OLD_LD_PRELOAD _NVRTC_OLD_LD_LIBRARY_PATH
EOF
else
    echo "== driver CUDA $DRIVER_CUDA >= toolkit $CUDA: no NVRTC workaround needed"
fi

echo
echo "Done. Environment: $ENV. Next, on a GPU node:"
echo "  conda activate $ENV"
echo "  python utils/check_faster.py                                   # generic water box"
echo "  python utils/check_faster.py step3_input.parm7 step3_input.rst7  # your built system"
