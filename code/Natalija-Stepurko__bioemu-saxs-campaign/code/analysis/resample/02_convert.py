"""Every v1.2 sample recovered from the npz batches with BioEmu's own conversion, and the physical filter
evaluated per conformer on both sides (design §2 and §5). Runs in the environment with the BioEmu package
(analysis/resample/requirements.txt), CPU only:

    /path/to/venv-bioemu/bin/python 02_convert.py [--labels ...]

Per protein:
  1. The npz batches are read in BioEmu's order (sorted file names) and concatenated, as bioemu.sample.main
     does before converting.
  2. Reproduction: bioemu.convert_chemgraph.save_pdb_and_xtc with filter_samples=True, the notebook's setting,
     writes into RESAMPLE_WORK/v12/<label>/reproduced/. The filter is wrapped so that the per-conformer
     criteria are recorded on the exact in-memory trajectory the filter sees. The reproduced files are
     compared with the notebook's: same file names, same frame count, largest coordinate difference.
  3. Recovery: save_pdb_and_xtc with filter_samples=False writes every sample to
     RESAMPLE_WORK/v12/<label>/samples_all.xtc (+ topology.pdb). Each frame of the notebook's xtc file is
     compared with the recovered sample it should be (the passing samples for samples.xtc, every sample for
     samples_unphysical.xtc) by RMSD after superposition.
  4. Archive: the same criteria, with BioEmu's own function, on the atoms BioEmu's output carries
     (N, CA, C, O, CB) taken from the archive's all-atom conformers.

Checks (design §2; any failure stops the analysis after the tables are written): 100 samples recovered;
the recovered atoms are BioEmu's backbone set; the notebook wrote the same files as the reproduction; the
frames of samples.xtc equal the recomputed passes, or every sample is in samples_unphysical.xtc when none
passes; every notebook frame matches its recovered sample within XTC_RMSD_TOL_A.

Outputs: results/resample/conversion.csv (one row per protein), results/resample/filter.csv (one row per
conformer and side: the three criteria and their conjunction).
"""

from __future__ import annotations

import argparse
import inspect
import sys
import time
from importlib.metadata import version
from pathlib import Path

import bioemu.convert_chemgraph as cc
import mdtraj as md
import numpy as np
import pandas as pd
import rslib
import torch
from common import (
    ENSEMBLES,
    FILTER_ATOMS,
    N_SAMPLES,
    OUT,
    SAMPLES,
    V12,
    XTC_RMSD_TOL_A,
    labels,
    pin,
    record,
    sha256,
)

CRITERIA = ("pass_ca_seq", "pass_cn_seq", "pass_clash")
XTC_NAMES = ("samples.xtc", "samples_unphysical.xtc")


def load_batches(folder: Path) -> tuple[np.ndarray, np.ndarray, str]:
    """As bioemu.sample.main reads them: sorted file names, one sequence, positions and orientations."""
    files = sorted(folder.glob("batch_*.npz"))
    seqs = {np.load(f)["sequence"].item() for f in files}
    if len(seqs) != 1:
        raise ValueError(f"{folder.name}: sequences differ between batches")
    pos = np.concatenate([np.load(f)["pos"] for f in files])
    ori = np.concatenate([np.load(f)["node_orientations"] for f in files])
    return pos, ori, seqs.pop()


def masks(traj: md.Trajectory) -> dict[str, np.ndarray]:
    """BioEmu's three criteria per conformer, with its default thresholds."""
    ca, cn, clash = cc._filter_unphysical_traj_masks(traj)
    return {"pass_ca_seq": np.asarray(ca), "pass_cn_seq": np.asarray(cn), "pass_clash": np.asarray(clash)}


def convert(pos, ori, seq, top: Path, xtc: Path, filter_samples: bool) -> None:
    cc.save_pdb_and_xtc(pos_nm=torch.tensor(pos), node_orientations=torch.tensor(ori), sequence=seq,
                        topology_path=top, xtc_path=xtc, filter_samples=filter_samples)


def reproduce(pos, ori, seq, out: Path) -> dict[str, np.ndarray]:
    """save_pdb_and_xtc as the notebook ran it, recording the criteria on the trajectory the filter sees."""
    seen = {}
    real = cc.filter_unphysical_traj

    def recording(traj, *args, **kw):
        seen.update(masks(traj))
        seen["n_atoms"] = traj.n_atoms
        return real(traj, *args, **kw)

    out.mkdir(parents=True, exist_ok=True)
    for f in out.glob("*.xtc"):
        f.unlink()
    cc.filter_unphysical_traj = recording
    try:
        convert(pos, ori, seq, out / "topology.pdb", out / "samples.xtc", filter_samples=True)
    finally:
        cc.filter_unphysical_traj = real
    if "pass_clash" not in seen:
        raise RuntimeError("save_pdb_and_xtc did not call filter_unphysical_traj")
    return seen


def max_rmsd_angstrom(a: md.Trajectory, b: md.Trajectory) -> float:
    """Largest RMSD (Å) between frame i of a and frame i of b, each pair superposed (float64)."""
    if a.n_frames == 0:
        return float("nan")
    return float(10.0 * rslib.kabsch_rmsd(a.xyz, b.xyz).max())


def one(lab: str) -> tuple[dict, list[dict]]:
    src = SAMPLES / lab
    pos, ori, seq = load_batches(src)
    dest = V12 / lab
    rep = reproduce(pos, ori, seq, dest / "reproduced")
    convert(pos, ori, seq, dest / "topology.pdb", dest / "samples_all.xtc", filter_samples=False)
    allt = md.load(str(dest / "samples_all.xtc"), top=str(dest / "topology.pdb"))
    keep = rep["pass_ca_seq"] & rep["pass_cn_seq"] & rep["pass_clash"]
    names = {a.name for a in allt.topology.atoms}
    r = {"label": lab, "n_recovered": allt.n_frames, "n_pass_recomputed": int(keep.sum()),
         "recovered_atoms": allt.n_atoms, "recovered_atom_names": ";".join(sorted(names)),
         "atoms_are_filter_atoms": names <= set(FILTER_ATOMS)}
    for k in CRITERIA:
        r[f"n_{k}"] = int(rep[k].sum())
    # the notebook's files against the reproduction, and each notebook frame against its recovered sample
    for name in XTC_NAMES:
        col, rp = src / name, dest / "reproduced" / name
        r[f"colab_{name}"] = col.exists()
        r[f"reproduced_{name}"] = rp.exists()
        if not col.exists():
            continue
        tc = md.load(str(col), top=str(src / "topology.pdb"))
        r[f"frames_colab_{name}"] = tc.n_frames
        if rp.exists():
            tr = md.load(str(rp), top=str(dest / "reproduced" / "topology.pdb"))
            same = tc.n_frames == tr.n_frames and tc.n_atoms == tr.n_atoms
            r[f"max_abs_diff_A_{name}"] = float(10.0 * np.abs(tc.xyz - tr.xyz).max()) if same else np.nan
        idx = np.where(keep)[0] if name == "samples.xtc" else np.arange(allt.n_frames)
        r[f"frames_expected_{name}"] = len(idx)
        if len(idx) == tc.n_frames and tc.n_atoms == allt.n_atoms:
            r[f"max_rmsd_A_colab_vs_recovered_{name}"] = max_rmsd_angstrom(tc, allt[idx])
            # the same comparison by mdtraj.rmsd (float32), kept to show why it is not used
            r[f"max_mdtraj_rmsd_A_colab_vs_recovered_{name}"] = float(max(
                10.0 * md.rmsd(tc[i], allt[j])[0] for i, j in enumerate(idx))) if len(idx) else np.nan
    r["files_match"] = all(r[f"colab_{n}"] == r[f"reproduced_{n}"] for n in XTC_NAMES)
    fails = []
    if r["n_recovered"] != N_SAMPLES:
        fails.append("samples recovered")
    if not r["atoms_are_filter_atoms"]:
        fails.append("recovered atoms")
    if not r["files_match"]:
        fails.append("file names")
    for name in XTC_NAMES:
        if r[f"colab_{name}"]:
            if r[f"frames_colab_{name}"] != r[f"frames_expected_{name}"]:
                fails.append(f"frames of {name}")
            elif not r.get(f"max_rmsd_A_colab_vs_recovered_{name}", np.inf) <= XTC_RMSD_TOL_A:
                fails.append(f"RMSD of {name}")
    r["check_failures"] = ";".join(fails)
    frames = [{"label": lab, "side": "v12", "conformer": i, **{k: bool(rep[k][i]) for k in CRITERIA},
               "pass_all": bool(keep[i])} for i in range(allt.n_frames)]
    # archive conformers, judged on the same atoms with the same function
    t = md.load(str(ENSEMBLES / f"{lab}.xtc"), top=str(ENSEMBLES / f"{lab}.pdb"))
    ta = t.atom_slice([a.index for a in t.topology.atoms if a.name in FILTER_ATOMS])
    ma = masks(ta)
    ka = ma["pass_ca_seq"] & ma["pass_cn_seq"] & ma["pass_clash"]
    r.update({"archive_conformers": t.n_frames, "archive_residues": t.topology.n_residues,
              "archive_filter_atoms": ta.n_atoms, "archive_n_pass": int(ka.sum())})
    for k in CRITERIA:
        r[f"archive_n_{k}"] = int(ma[k].sum())
    frames += [{"label": lab, "side": "archive", "conformer": i, **{k: bool(ma[k][i]) for k in CRITERIA},
                "pass_all": bool(ka[i])} for i in range(t.n_frames)]
    return r, frames


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", nargs="*")
    a = ap.parse_args()
    t0 = time.time()
    used = pin()
    torch.set_num_threads(1)
    rows, frames = [], []
    for lab in a.labels or labels():
        t1 = time.time()
        r, f = one(lab)
        r["seconds"] = round(time.time() - t1, 1)
        rows.append(r)
        frames += f
        print(f"{lab}: recovered {r['n_recovered']}, pass {r['n_pass_recomputed']}, "
              f"archive pass {r['archive_n_pass']}/{r['archive_conformers']}, "
              f"failures [{r['check_failures']}], {r['seconds']} s", flush=True)
    conv = pd.DataFrame(rows)
    if a.labels:   # a trial on a few proteins: report, leave the tracked tables alone
        print(conv.T.to_string())
        return
    OUT.mkdir(parents=True, exist_ok=True)
    conv.to_csv(OUT / "conversion.csv", index=False)
    pd.DataFrame(frames).to_csv(OUT / "filter.csv", index=False)
    params = inspect.signature(cc._filter_unphysical_traj_masks).parameters
    record("02_convert", {
        "bioemu": version("bioemu"), "torch": version("torch"), "mdtraj": version("mdtraj"),
        "bioemu_convert_chemgraph_sha256": sha256(Path(cc.__file__)), "filter_atoms": list(FILTER_ATOMS),
        "filter_defaults": {k: v.default for k, v in params.items() if k != "traj"},
        "xtc_rmsd_tolerance_A": XTC_RMSD_TOL_A, "cpus": used, "wall_seconds": round(time.time() - t0, 1),
        "n_proteins_failing_checks": int((conv["check_failures"] != "").sum())})
    bad = conv[conv["check_failures"] != ""]
    if len(bad):
        print(bad[["label", "check_failures"]].to_string(), file=sys.stderr)
        sys.exit(f"{len(bad)} proteins fail the conversion checks")


if __name__ == "__main__":
    main()
