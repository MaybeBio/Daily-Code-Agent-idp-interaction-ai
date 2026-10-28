"""Run the CALVADOS 3 single-chain simulations (design §4) with resume.

    python 05_run.py [--labels L1 L2 ...] [--workers 6] [--frames N]

The driver starts one process per protein (`--one LABEL`), at most --workers at a time, each on one CPU
thread, with its output in SIM_WORK/logs/<label>.log. A protein's system is built by the CALVADOS package
from the input table (results/simulations/inputs.csv); the OpenMM simulation then runs in chunks of CHUNK
frames, each chunk written to its own DCD file and followed by a checkpoint and a progress record, so a
restart continues from the last complete chunk (a partial chunk is deleted). --frames N stops a protein
after N frames (the pilot uses it on the longest protein); a later run without it continues to its full
length. Once results/simulations/plan.json exists (04_plan.py), its frame counts and inclusions apply.
Finished proteins are skipped. Per protein, SIM_WORK/runs/<label>/progress.json records frames written,
steps, wall time and speed."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from common import CALVADOS_COMMIT, DL, LOGS, OUT, RUNS

CHUNK = 50  # frames per DCD file and checkpoint


def plan() -> dict | None:
    """results/simulations/plan.json (04_plan.py) by label, once the pilot has been planned."""
    f = OUT / "plan.json"
    return {p["label"]: p for p in json.loads(f.read_text())["proteins"]} if f.exists() else None


def build(label: str, r: pd.Series, path: Path):
    """CALVADOS configuration of one protein; returns the built calvados.sim.Sim."""
    from calvados import sim
    from calvados.cfg import Components, Config

    path.mkdir(parents=True, exist_ok=True)
    (path / "input.fasta").write_text(f">{label}\n{r['sequence']}\n")
    L = float(r["box_nm"])
    cfg = Config(
        sysname=label,
        box=[L, L, L],
        temp=float(r["temp_k"]),
        ionic=float(r["ionic_m"]),
        pH=float(r["ph"]),
        topol="center",
        wfreq=int(r["n_save"]),
        steps=int(r["steps"]),
        platform="CPU",
        threads=1,
        restart=None,
        random_number_seed=int(r["seed"]),
        verbose=False,
    )
    comp = Components(
        molecule_type="protein",
        nmol=1,
        restraint=False,
        charge_termini="both",
        fresidues=str(DL / "residues_CALVADOS3.csv"),
        ffasta=str(path / "input.fasta"),
    )
    comp.add(name=label)
    cfg.write(str(path), name="config.yaml")
    comp.write(str(path), name="components.yaml")
    import yaml

    s = sim.Sim(
        str(path), yaml.safe_load(open(path / "config.yaml")), yaml.safe_load(open(path / "components.yaml"))
    )
    s.build_system()
    return s


def run_one(label: str, max_frames: int | None) -> None:
    import openmm
    from openmm import app, unit

    r = pd.read_csv(OUT / "inputs.csv").set_index("label").loc[label]
    path = RUNS / label
    prog_f, chk = path / "progress.json", path / "restart.chk"
    prog = (
        json.loads(prog_f.read_text())
        if prog_f.exists()
        else {"frames": 0, "wall_seconds": 0.0, "chunks": []}
    )
    pl = plan()
    if pl is not None and label in pl and not pl[label]["included"]:
        print(f"{label}: left out by the compute plan", flush=True)
        return
    total = int(pl[label]["frames_planned"]) if pl is not None and label in pl else int(r["n_frames"])
    target = min(total, max_frames) if max_frames else total
    if prog["frames"] >= target:
        print(f"{label}: {prog['frames']} frames already written", flush=True)
        return
    s = build(label, r, path)
    integ = openmm.LangevinMiddleIntegrator(
        float(r["temp_k"]) * unit.kelvin, 0.01 / unit.picosecond, 0.01 * unit.picosecond
    )
    integ.setRandomNumberSeed(int(r["seed"]))
    simulation = app.Simulation(
        app.PDBxFile(s.cif_cg).topology,
        s.system,
        integ,
        openmm.Platform.getPlatformByName("CPU"),
        {"Threads": "1"},
    )
    k = prog["frames"] // CHUNK
    for stale in path.glob("chunk_*.dcd"):  # a chunk not covered by the checkpoint is incomplete
        if int(stale.stem.split("_")[1]) >= k:
            stale.unlink()
    if prog["frames"] > 0:
        simulation.loadCheckpoint(str(chk))
        print(f"{label}: resumed at frame {prog['frames']} (step {simulation.currentStep})", flush=True)
    else:
        simulation.context.setPositions(app.PDBxFile(s.cif_cg).positions)
        simulation.minimizeEnergy()
        simulation.currentStep = 0
    n_save = int(r["n_save"])
    while prog["frames"] < target:
        n = min(CHUNK, target - prog["frames"])
        rep = app.DCDReporter(str(path / f"chunk_{k:03d}.dcd"), n_save, enforcePeriodicBox=False)
        simulation.reporters = [rep]
        t0 = time.time()
        simulation.step(n * n_save)
        dt = time.time() - t0
        rep._out.close()
        simulation.saveCheckpoint(str(chk))
        prog["frames"] += n
        prog["wall_seconds"] += dt
        prog["chunks"].append(
            {"chunk": k, "frames": n, "seconds": round(dt, 2), "steps_per_second": round(n * n_save / dt, 1)}
        )
        prog.update(
            {
                "label": label,
                "length": int(r["length"]),
                "n_save": n_save,
                "target_frames": total,
                "steps": prog["frames"] * n_save,
                "steps_per_second": round(prog["frames"] * n_save / prog["wall_seconds"], 1),
                "calvados_commit": CALVADOS_COMMIT,
                "openmm": openmm.__version__,
                "threads": 1,
                "host_cpus": os.cpu_count(),
            }
        )
        tmp = prog_f.with_suffix(".tmp")
        tmp.write_text(json.dumps(prog, indent=1))
        tmp.replace(prog_f)
        print(f"{label}: {prog['frames']}/{total} frames, {prog['steps_per_second']:.0f} steps/s", flush=True)
        k += 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", nargs="*")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--frames", type=int, default=None)
    ap.add_argument("--one")
    a = ap.parse_args()
    if a.one:
        run_one(a.one, a.frames)
        return
    inp = pd.read_csv(OUT / "inputs.csv")
    pl = plan()
    labels = a.labels or [
        x
        for x in inp.sort_values("length", ascending=False)["label"]  # longest first
        if pl is None or pl[x]["included"]
    ]

    def launch(lab):
        cmd = [sys.executable, __file__, "--one", lab] + (["--frames", str(a.frames)] if a.frames else [])
        with open(LOGS / f"{lab}.log", "a") as log:
            log.write(f"=== {time.strftime('%Y-%m-%dT%H:%M:%S')} {' '.join(cmd[2:])}\n")
            log.flush()
            rc = subprocess.run(
                cmd, stdout=log, stderr=subprocess.STDOUT, env={**os.environ, "OPENMM_CPU_THREADS": "1"}
            ).returncode
        print(f"{lab}: exit {rc}", flush=True)
        return lab, rc

    with ThreadPoolExecutor(a.workers) as ex:
        res = list(ex.map(launch, labels))
    bad = [lab for lab, rc in res if rc]
    print("failed:", bad if bad else "none", flush=True)


if __name__ == "__main__":
    main()
