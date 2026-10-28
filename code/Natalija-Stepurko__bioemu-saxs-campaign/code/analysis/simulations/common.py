"""Paths and pinned sources of the simulation analysis (docs/SIMULATION_DESIGN.md).

SIM_WORK holds downloads, tools, trajectories, back-mapped frames and computed curves (not tracked;
several GB); the small outputs go to results/simulations/ (tracked)."""

from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WORK = Path(os.environ.get("SIM_WORK", REPO / "analysis" / "simulations" / "work"))
DL = WORK / "downloads"
TOOLS = WORK / "tools"
RUNS = WORK / "runs"
SAXS = WORK / "saxs"
LOGS = WORK / "logs"
OUT = REPO / "results" / "simulations"
RES = REPO / "results"

CALVADOS_REPO = "https://github.com/KULL-Centre/CALVADOS"
CALVADOS_COMMIT = "69815423632cb9b9b08f4b0aea34d8ba43a057c4"  # package version 0.9.0
CALVADOS3_RESIDUES = (
    f"https://raw.githubusercontent.com/KULL-Centre/CALVADOS/{CALVADOS_COMMIT}/"
    "examples/single_MDP/input/residues_CALVADOS3.csv"
)
PAPERS_COMMIT = "85979b1b4123b6b5391b617d16551969eda9f56e"  # github.com/KULL-Centre/papers
CAO_COMMIT = "51564e0edd126c44d6a7299fbe4414c579394ef4"  # github.com/KULL-Centre/_2024_Cao_CALVADOSCOM
PULCHRA_URL = f"https://raw.githubusercontent.com/KULL-Centre/papers/{PAPERS_COMMIT}/2021/CG-IDPs-Tesei-et-al/pulchra_306.tgz"
PEPSI_URL = "https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/Linux/3.0/Pepsi-SAXS-Linux.zip"
PULCHRA = TOOLS / "pulchra_306" / "pulchra"
PEPSI = TOOLS / "pepsi" / "Pepsi-SAXS"
PEPSI_FLAGS = ["--maximum_scattering_vector", "10"]  # as in the PeptoneBench Pepsi-SAXS wrapper
N_SAXS_FRAMES = 100

for _d in (DL, TOOLS, RUNS, SAXS, LOGS, OUT):
    _d.mkdir(parents=True, exist_ok=True)
