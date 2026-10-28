"""Stage `coordrg` (DESIGN §7, B2): the radius of gyration of the BioEmu-1 conformers from their
coordinates, for every clean disordered entry and a random sample of folded ones, to tell whether the
excess ensemble size of disordered proteins sits in the conformers or in the forward model.

Needs the conformer ensembles (xtc + pdb) under data/Predictions/PeptoneDB-SAXS-ensembles/bioemu/.
When the PeptoneBench Predictions archive is available the missing ensembles are extracted from it in
one pass (`--archive PATH`). Output: results/coord_rg.csv (one row per entry).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from bsc import analyse
from bsc import config as C

ENSEMBLES = C.DATA / "Predictions" / "PeptoneDB-SAXS-ensembles"
N_FOLDED_SAMPLE = 40


def selection(d: pd.DataFrame) -> dict[str, str]:
    """Labels to analyse with their reason: every clean disordered entry of the model under study,
    a seeded random sample of folded ones, and the three worked examples."""
    m = d[(d["model"] == C.MODEL) & d["clean"]]
    out = {lab: "disordered" for lab in sorted(m.loc[m["disorder_class"] == "disordered", "label"])}
    folded = sorted(m.loc[m["disorder_class"] == "folded", "label"])
    rng = np.random.default_rng(C.SEED)
    for lab in sorted(rng.choice(folded, min(N_FOLDED_SAMPLE, len(folded)), replace=False)):
        out.setdefault(lab, "folded sample")
    for lab in analyse.choose_examples(d).values():
        out.setdefault(lab, "worked example")
    return out


def archive_members(labels: list[str], examples: list[str]) -> list[str]:
    mem = []
    for lab in labels:
        mem += [f"Predictions/PeptoneDB-SAXS-ensembles/{C.MODEL}/{lab}.{ext}" for ext in ("xtc", "pdb")]
    mem += [f"Predictions/PeptoneDB-SAXS-ensembles/alphafold/{lab}.pdb" for lab in examples]
    return mem


def extract(archive: Path, members: list[str]) -> None:
    """One pass over the archive for every missing member."""
    missing = [m for m in members if not (C.DATA / m).exists()]
    if not missing:
        return
    lst = C.DATA / "members_to_extract.txt"
    lst.write_text("\n".join(missing) + "\n")
    print(f"  extracting {len(missing)} files from {archive.name} (one pass)", flush=True)
    subprocess.run(["tar", "-xzf", str(archive), "-C", str(C.DATA), "-T", str(lst)], check=True)
    lst.unlink()


def ensemble_rg(label: str) -> dict | None:
    """Per-conformer Rg from the C-alpha and from all heavy atoms (Å), summarised over conformers."""
    import mdtraj as md
    xtc, pdb = ENSEMBLES / C.MODEL / f"{label}.xtc", ENSEMBLES / C.MODEL / f"{label}.pdb"
    if not (xtc.exists() and pdb.exists()):
        return None
    tr = md.load(str(xtc), top=str(pdb))
    ca = tr.atom_slice(tr.topology.select("name CA"))
    rg_ca = 10.0 * md.compute_rg(ca)
    rg_heavy = 10.0 * md.compute_rg(tr)
    return {"label": label, "n_conformers": tr.n_frames, "n_residues": tr.topology.n_residues,
            "rg_ca_mean": float(rg_ca.mean()), "rg_ca_median": float(np.median(rg_ca)),
            "rg_ca_sd": float(rg_ca.std()),
            # the ensemble-average curve's Rg is the z-average over conformers: sqrt(mean(Rg^2))
            "rg_ca_rms": float(np.sqrt(np.mean(rg_ca**2))),
            "rg_heavy_mean": float(rg_heavy.mean()), "rg_heavy_rms": float(np.sqrt(np.mean(rg_heavy**2)))}


def main(archive: Path | None = None) -> None:
    _, d = analyse.load()
    sel = selection(d)
    examples = list(analyse.choose_examples(d).values())
    if archive is not None and archive.exists():
        extract(archive, archive_members(sorted(sel), examples))
    rows = []
    for lab, why in sel.items():
        r = ensemble_rg(lab)
        if r is None:
            continue
        r["reason"] = why
        rows.append(r)
    if not rows:
        print("  no ensembles under data/Predictions/PeptoneDB-SAXS-ensembles; stage skipped", flush=True)
        return
    out = pd.DataFrame(rows)
    m = d[(d["model"] == C.MODEL)].set_index("label")
    out["disorder_class"] = m.loc[out["label"], "disorder_class"].to_numpy()
    out["length"] = m.loc[out["label"], "length"].to_numpy()
    out["rg_exp"] = m.loc[out["label"], "rg_exp"].to_numpy()
    out["rg_curve"] = m.loc[out["label"], "rg_model"].to_numpy()
    out["ratio_coord_exp"] = out["rg_ca_rms"] / out["rg_exp"]
    out["ratio_curve_exp"] = out["rg_curve"] / out["rg_exp"]
    out["ratio_curve_coord"] = out["rg_curve"] / out["rg_ca_rms"]
    C.RESULTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(C.RESULTS / "coord_rg.csv", index=False)
    for cls, g in out.groupby("disorder_class"):
        print(f"  {cls}: n = {len(g)}; median coordinate Rg / exp {g['ratio_coord_exp'].median():.3f}; "
              f"curve Rg / exp {g['ratio_curve_exp'].median():.3f}; curve / coordinate "
              f"{g['ratio_curve_coord'].median():.3f}", flush=True)


if __name__ == "__main__":
    main()
