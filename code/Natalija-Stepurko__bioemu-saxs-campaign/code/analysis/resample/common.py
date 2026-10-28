"""Paths and fixed settings of the v1.2 re-sample analysis (docs/RESAMPLE_DESIGN.md).

Standard library only: imported both by the scripts that run in the main environment and by the conversion
script, which runs in a separate environment with the BioEmu package.

RESAMPLE_SAMPLES  the notebook's output folder (one sub-folder per protein, progress.json)
RESAMPLE_WORK     recovered conformers, back-mapped structures and computed curves (not tracked; ~1 GB)
BSC_ROOT          where `bsc fetch` put data/ (measured profiles, the archive's ensembles and curves)
SIM_WORK          the simulation analysis's work folder; PULCHRA 3.06 and Pepsi-SAXS 3.0 are taken from its
                  tools/, where analysis/simulations/02_fetch_sources.py installs them (PULCHRA and
                  PEPSI override them)
RESAMPLE_CPUS     the cores the scripts may use, e.g. "5,6,7" (default: the last three of this machine)
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WORK = Path(os.environ.get("RESAMPLE_WORK", REPO / "analysis" / "resample" / "work"))
SAMPLES = Path(os.environ.get("RESAMPLE_SAMPLES", WORK / "bioemu_resample"))
BSC_ROOT = Path(os.environ.get("BSC_ROOT", REPO))
DATA = BSC_ROOT / "data"
ENSEMBLES = DATA / "Predictions" / "PeptoneDB-SAXS-ensembles" / "bioemu"
ARCHIVE_CURVES = DATA / "Predictions" / "PeptoneDB-SAXS-expt" / "bioemu"
SAXS_DIR = DATA / "PeptoneDB-SAXS" / "sasbdb-clean_data"

V12 = WORK / "v12"            # per protein: every sample recovered from the npz batches, the reproduction
SAXS = WORK / "saxs"          # per side and protein: per-conformer Pepsi-SAXS output, curves.csv, logs
OUT = REPO / "results" / "resample"
RES = REPO / "results"
FIG = RES / "figures"
LIST = REPO / "notebooks" / "resample_list.csv"
VERSIONS = OUT / "versions.json"

# the binaries and flags of analysis/simulations/common.py (not imported: it creates its work folders)
SIM_TOOLS = Path(os.environ.get("SIM_WORK", REPO / "analysis" / "simulations" / "work")) / "tools"
PULCHRA = Path(os.environ.get("PULCHRA", SIM_TOOLS / "pulchra_306" / "pulchra"))
PEPSI = Path(os.environ.get("PEPSI", SIM_TOOLS / "pepsi" / "Pepsi-SAXS"))
PEPSI_FLAGS = ["--maximum_scattering_vector", "10"]   # as in the PeptoneBench wrapper
PEPSI_FALLBACK = ["--nCoeff", "50"]                   # the wrapper's retry after a non-zero exit

MODEL = "bioemu-v1.2"
N_SAMPLES = 100
# the atoms BioEmu's output carries; the physical filter is evaluated on these on both sides
FILTER_ATOMS = ("N", "CA", "C", "O", "CB")
# the largest RMSD two xtc copies of one conformer can differ by: each file rounds every coordinate to
# 0.001 nm, so an atom moves by at most sqrt(3) x 0.005 Å per file
XTC_RMSD_TOL_A = 3**0.5 * 0.01
MIN_FILTERED = 10             # fewer passing conformers on either side: out of the paired filtered comparison
SEED = 20261016
N_BOOT = 2000
MAX_WORKERS = 3
SIDES = ("v12", "archive")


def cpus() -> list[int]:
    env = os.environ.get("RESAMPLE_CPUS")
    if env:
        return [int(c) for c in env.split(",")]
    n = os.cpu_count() or 1
    return list(range(max(0, n - MAX_WORKERS), n))


def pin() -> list[int]:
    """Restrict this process and its children to at most MAX_WORKERS cores (Pepsi-SAXS uses every core it
    can see, and has no option to limit its threads)."""
    c = cpus()[:MAX_WORKERS]
    os.sched_setaffinity(0, c)
    os.environ["OMP_NUM_THREADS"] = "1"
    return c


def commit() -> str:
    try:
        h = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True,
                           check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "analysis/resample"],
                               capture_output=True, text=True).stdout.strip()
        return h + ("+uncommitted" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def record(stage: str, info: dict) -> None:
    """Write this stage's parameter record (software, inputs, wall time, code commit) into versions.json."""
    OUT.mkdir(parents=True, exist_ok=True)
    v = json.loads(VERSIONS.read_text()) if VERSIONS.exists() else {}
    v[stage] = {"run_at": datetime.datetime.now().isoformat(timespec="seconds"), "commit": commit(), **info}
    VERSIONS.write_text(json.dumps(v, indent=1, sort_keys=True) + "\n")


def labels() -> list[str]:
    import csv
    with open(LIST) as f:
        return [r["label"] for r in csv.DictReader(f)]
