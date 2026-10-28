"""Computed SAXS curves of every conformer, v1.2 (path A) and archive (path B), through the same steps
(design §4).

    python 04_saxs.py [--labels ...] [--sides v12 archive] [--workers 3]

Per conformer: the Cα atoms are written to a PDB file, back-mapped to all atoms with PULCHRA 3.06 (default
options) and passed to Pepsi-SAXS 3.0 with the measured profile and PEPSI_FLAGS; after a non-zero exit, once
more with PEPSI_FALLBACK. The fitted curve (fourth output column) is kept. Per-conformer outputs (.dat, .log)
are cached under RESAMPLE_WORK/saxs/<side>/<label>/ (the log is written last, so a frame with a log is
complete), and a restart continues where it stopped; the PDB files are deleted after use. Writes curves.csv
(conformers x q, NaN where Pepsi-SAXS gave no curve), pepsi_log.csv and timing.json per side and protein.
At most three worker processes, pinned to three cores.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import mdtraj as md
import numpy as np
import pandas as pd
import rslib
from common import (
    ENSEMBLES,
    LIST,
    MAX_WORKERS,
    PEPSI,
    PEPSI_FALLBACK,
    PEPSI_FLAGS,
    PULCHRA,
    SAXS,
    SAXS_DIR,
    SIDES,
    V12,
    labels,
    pin,
    record,
    sha256,
)


def traj(side: str, lab: str) -> md.Trajectory:
    if side == "v12":
        return md.load(str(V12 / lab / "samples_all.xtc"), top=str(V12 / lab / "topology.pdb"))
    return md.load(str(ENSEMBLES / f"{lab}.xtc"), top=str(ENSEMBLES / f"{lab}.pdb"))


def read_dat(p) -> pd.DataFrame:
    return pd.read_csv(p, sep=r"\s+", header=None, names=["q", "I", "sigma", "I_fit"], comment="#")


def one(task: tuple[str, str]) -> str:
    side, lab = task
    out = SAXS / side / lab
    if (out / "curves.csv").exists():
        return f"{side} {lab}: done"
    t0 = time.time()
    out.mkdir(parents=True, exist_ok=True)
    exp = pd.read_csv(SAXS_DIR / f"{lab}-bift.dat", sep=r"\s+", comment="#", header=None,
                      names=["q", "I(q)", "sigma"])
    exp_f = out / "exp.dat"
    exp.to_csv(exp_f, sep="\t", index=False)  # as the PeptoneBench wrapper and analysis/simulations write it
    t = traj(side, lab)
    ca = t.atom_slice(t.topology.select("name CA"))
    curves, logs, computed = [], [], 0
    for i in range(ca.n_frames):
        dat, log = out / f"frame_{i:04d}.dat", out / f"frame_{i:04d}.log"
        rec = {"conformer": i}
        if not log.exists():
            computed += 1
            pdb = out / f"frame_{i:04d}.pdb"
            ca[i].save_pdb(str(pdb))
            subprocess.run([str(PULCHRA), pdb.name], cwd=out, check=True, capture_output=True)
            rebuilt = out / f"frame_{i:04d}.rebuilt.pdb"
            cmd = [str(PEPSI), rebuilt.name, exp_f.name, "-o", dat.name, *PEPSI_FLAGS]
            p = subprocess.run(cmd, cwd=out, capture_output=True, text=True)
            fallback = p.returncode != 0
            if fallback:
                p = subprocess.run(cmd + PEPSI_FALLBACK, cwd=out, capture_output=True, text=True)
            log.write_text(f"# returncode {p.returncode} fallback {int(fallback)}\n" + p.stdout)
            pdb.unlink()
            rebuilt.unlink()
            if p.returncode != 0:
                dat.unlink(missing_ok=True)
        head = log.read_text().splitlines()[0].split()
        rec.update({"returncode": int(head[2]), "fallback": bool(int(head[4]))})
        rec.update(rslib.pepsi_log_values(log.read_text()))
        if dat.exists():
            d = read_dat(dat)
            if not np.allclose(d["q"], exp["q"], atol=1e-6):
                raise ValueError(f"{side} {lab} conformer {i}: q grid differs from the measured profile")
            curves.append(d["I_fit"].to_numpy())
        else:
            curves.append(np.full(len(exp), np.nan))
        rec["finite"] = bool(np.isfinite(curves[-1]).all())
        logs.append(rec)
    pd.DataFrame(logs).to_csv(out / "pepsi_log.csv", index=False)
    (out / "timing.json").write_text(json.dumps({"seconds": round(time.time() - t0, 1),
                                                 "conformers_computed": computed}))
    pd.DataFrame(np.array(curves), columns=exp["q"].to_numpy()).to_csv(out / "curves.tmp")
    (out / "curves.tmp").replace(out / "curves.csv")
    n_ok = sum(r["finite"] for r in logs)
    return f"{side} {lab}: {n_ok}/{len(logs)} finite curves, {time.time() - t0:.0f} s"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", nargs="*")
    ap.add_argument("--sides", nargs="*", default=list(SIDES))
    ap.add_argument("--workers", type=int, default=MAX_WORKERS)
    a = ap.parse_args()
    if a.workers > MAX_WORKERS:
        sys.exit(f"at most {MAX_WORKERS} workers")
    t0 = time.time()
    used = pin()
    length = pd.read_csv(LIST).set_index("label")["length"]
    # longest chains first, so the three workers finish together
    tasks = sorted(((s, lab) for lab in (a.labels or labels()) for s in a.sides), key=lambda t: -length[t[1]])
    with ProcessPoolExecutor(a.workers) as ex:
        for msg in ex.map(one, tasks):
            print(msg, flush=True)
    timing = {f"{s}/{lab}": json.loads((SAXS / s / lab / "timing.json").read_text()) for s, lab in tasks}
    ver = subprocess.run([str(PEPSI), "--version"], capture_output=True, text=True).stdout.strip()
    ver = ver.splitlines()
    record("04_saxs", {
        "pulchra": str(PULCHRA), "pulchra_sha256": sha256(PULCHRA), "pepsi": str(PEPSI),
        "pepsi_sha256": sha256(PEPSI), "pepsi_version": ver[-1] if ver else None,
        "pepsi_flags": PEPSI_FLAGS, "pepsi_fallback": PEPSI_FALLBACK, "workers": a.workers, "cpus": used,
        "wall_seconds_this_run": round(time.time() - t0, 1),
        "seconds_per_task_sum": round(sum(v["seconds"] for v in timing.values()), 1), "tasks": len(tasks)})


if __name__ == "__main__":
    main()
