"""Computed SAXS curves of the simulated ensembles (design §3, secondary measure).

    python 06_saxs.py [--labels ...] [--workers 6]

For every protein whose run is complete: N_SAXS_FRAMES frames evenly spaced over the production frames (the
first n_discard dropped), each back-mapped to all atoms with PULCHRA 3.06 (default options, as in the CALVADOS
authors' scripts) and passed to Pepsi-SAXS 3.0 with the measured profile and the flags of the PeptoneBench
wrapper; the fitted curve of each frame (fourth column) is kept. Writes SIM_WORK/saxs/<label>/curves.csv
(frames x q, the format of the archive's computed curves) and pepsi_log.csv (d_rho, r0, chi2 per frame).
Per-frame outputs are cached, so a restart continues where it stopped.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd
from common import N_SAXS_FRAMES, OUT, PEPSI, PEPSI_FLAGS, PULCHRA, RUNS, SAXS

BSC_ROOT = Path(os.environ.get("BSC_ROOT", Path(__file__).resolve().parents[2]))
SAXS_DIR = BSC_ROOT / "data" / "PeptoneDB-SAXS" / "sasbdb-clean_data"


def load_traj(label: str) -> md.Trajectory:
    run = RUNS / label
    chunks = sorted(run.glob("chunk_*.dcd"))
    return md.load([str(c) for c in chunks], top=str(run / "top.pdb"))


def frame_indices(n_total: int, n_discard: int, n_pick: int = N_SAXS_FRAMES) -> np.ndarray:
    """n_pick frames evenly spaced over the production frames n_discard .. n_total - 1."""
    return np.unique(np.linspace(n_discard, n_total - 1, n_pick).round().astype(int))


def frames_planned(label: str, default: int) -> int:
    """Frames the compute plan (04_plan.py) gives this protein."""
    f = OUT / "plan.json"
    if not f.exists():
        return default
    return next(
        (int(p["frames_planned"]) for p in json.loads(f.read_text())["proteins"] if p["label"] == label),
        default,
    )


def discard_planned(label: str, default: int) -> int:
    f = OUT / "plan.json"
    if not f.exists():
        return default
    return next(
        (int(p["discard_planned"]) for p in json.loads(f.read_text())["proteins"] if p["label"] == label),
        default,
    )


def one(label: str) -> str:
    r = pd.read_csv(OUT / "inputs.csv").set_index("label").loc[label]
    out = SAXS / label
    if (out / "curves.csv").exists():
        return f"{label}: done"
    pf = RUNS / label / "progress.json"
    prog = json.loads(pf.read_text()) if pf.exists() else {"frames": 0}
    if prog["frames"] < frames_planned(label, int(r["n_frames"])):
        return f"{label}: run incomplete ({prog['frames']} frames)"
    out.mkdir(parents=True, exist_ok=True)
    exp = pd.read_csv(
        SAXS_DIR / f"{label}-bift.dat", sep=r"\s+", comment="#", header=None, names=["q", "I(q)", "sigma"]
    )
    exp_f = out / "exp.dat"
    exp.to_csv(exp_f, sep="\t", index=False)  # as the PeptoneBench wrapper writes it
    t = load_traj(label)
    idx = frame_indices(t.n_frames, discard_planned(label, int(r["n_discard"])))
    curves, logs = [], []
    for i in idx:
        dat = out / f"frame_{i:04d}.dat"
        if not dat.exists():
            pdb = out / f"frame_{i:04d}.pdb"
            t[i].save_pdb(str(pdb))
            subprocess.run([str(PULCHRA), pdb.name], cwd=out, check=True, capture_output=True)
            rebuilt = out / f"frame_{i:04d}.rebuilt.pdb"
            log = out / f"frame_{i:04d}.log"
            cmd = [str(PEPSI), rebuilt.name, exp_f.name, "-o", dat.name, *PEPSI_FLAGS]
            p = subprocess.run(cmd, cwd=out, capture_output=True, text=True)
            if p.returncode != 0:  # the wrapper's fallback: a fixed expansion order of 50
                p = subprocess.run(cmd + ["--nCoeff", "50"], cwd=out, capture_output=True, text=True)
            log.write_text(p.stdout)
            pdb.unlink()
            rebuilt.unlink()
            if p.returncode != 0:
                dat.unlink(missing_ok=True)
                continue
        d = pd.read_csv(dat, sep=r"\s+", header=None, names=["q", "I", "sigma", "I_fit"], comment="#")
        if not np.allclose(d["q"], exp["q"], atol=1e-6):
            raise ValueError(f"{label} frame {i}: q grid differs from the measured profile")
        curves.append(d["I_fit"].to_numpy())
        rec = {"frame": int(i)}
        for line in (out / f"frame_{i:04d}.log").read_text().splitlines():
            if line.startswith("Best d_rho found"):
                rec["d_rho"] = float(line.split()[-2])
            elif line.startswith("Best r0 found"):
                rec["r0"] = float(line.split()[-2])
            elif line.startswith("Chi^2"):
                rec["chi2_frame"] = float(line.split()[-1])
        logs.append(rec)
    frames = [rec["frame"] for rec in logs]
    pd.DataFrame(curves, columns=exp["q"].to_numpy(), index=frames).to_csv(out / "curves.tmp")
    pd.DataFrame(logs).to_csv(out / "pepsi_log.csv", index=False)
    (out / "curves.tmp").replace(out / "curves.csv")
    return f"{label}: {len(curves)} curves"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", nargs="*")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    labels = a.labels or pd.read_csv(OUT / "inputs.csv")["label"].tolist()
    with ProcessPoolExecutor(a.workers) as ex:
        for msg in ex.map(one, labels):
            print(msg, flush=True)


if __name__ == "__main__":
    sys.exit(main())
