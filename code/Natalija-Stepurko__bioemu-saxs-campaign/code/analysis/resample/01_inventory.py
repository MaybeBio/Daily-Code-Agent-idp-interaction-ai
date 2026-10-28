"""Inventory of the v1.2 samples and of the archive's ensembles for the listed proteins (design §2).

    python 01_inventory.py [--archive Predictions.tar.gz]

Per protein: the npz batches (files, samples, contiguous index ranges, stored sequence), the files the
notebook wrote, the frames in samples.xtc and samples_unphysical.xtc, the entry in progress.json (status,
model, requested samples, seconds), the class and Guinier Rg from the tracked results, and the archive's
conformers and computed curves. With --archive, archive ensembles missing under BSC_ROOT are extracted from
Predictions.tar.gz in one pass (bsc.coordrg.extract). Output: results/resample/inventory.csv.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import mdtraj as md
import pandas as pd
import rslib
from common import (
    ARCHIVE_CURVES,
    ENSEMBLES,
    LIST,
    MODEL,
    N_SAMPLES,
    OUT,
    RES,
    SAMPLES,
    labels,
    pin,
    record,
    sha256,
)

from bsc import coordrg
from bsc.analyse import clean_mask


def n_frames(xtc: Path, top: Path) -> int | None:
    return md.load(str(xtc), top=str(top)).n_frames if xtc.exists() else None


def fasta_sequence(p: Path) -> str:
    return "".join(line.strip() for line in p.read_text().splitlines() if not line.startswith(">"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path)
    a = ap.parse_args()
    t0 = time.time()
    used = pin()
    lst = pd.read_csv(LIST).set_index("label")
    labs = labels()
    if a.archive is not None:
        members = [f"Predictions/PeptoneDB-SAXS-ensembles/bioemu/{lab}.{ext}" for lab in labs
                   for ext in ("xtc", "pdb")]
        coordrg.extract(a.archive, members)
    prog = json.loads((SAMPLES / "progress.json").read_text())
    res = pd.read_csv(RES / "resolvability.csv")
    res = res[res["model"] == "bioemu"].set_index("label")
    ent = pd.read_csv(RES / "entries.csv").set_index("label")
    clean = clean_mask(ent)   # the study's one definition of an entry passing the quality checks
    rows = []
    for lab in labs:
        d = SAMPLES / lab
        p = prog.get(lab, {})
        r = {"label": lab, "length": int(lst.loc[lab, "length"]),
             "disorder_class": res.loc[lab, "disorder_class"],
             "class_in_list": lst.loc[lab, "disorder_class"],
             "clean": bool(clean.loc[lab]),
             "rg_exp": float(ent.loc[lab, "rg_exp"]), "guinier_valid": bool(ent.loc[lab, "guinier_valid"]),
             "progress_ok": bool(p.get("ok", False)), "progress_model": p.get("model"),
             "progress_num_samples": p.get("num_samples"), "colab_seconds": p.get("seconds")}
        c = rslib.count_npz_samples(d)
        r.update({k: v for k, v in c.items() if k != "npz_sequences"})
        r["npz_sequence_matches"] = c["npz_sequences"] == [lst.loc[lab, "sequence"]]
        r["fasta_matches"] = fasta_sequence(d / "sequence.fasta") == lst.loc[lab, "sequence"]
        top = d / "topology.pdb"
        r["has_samples_xtc"] = (d / "samples.xtc").exists()
        r["has_unphysical_xtc"] = (d / "samples_unphysical.xtc").exists()
        r["frames_samples_xtc"] = n_frames(d / "samples.xtc", top)
        r["frames_unphysical_xtc"] = n_frames(d / "samples_unphysical.xtc", top)
        r["topology_residues"] = md.load_topology(str(top)).n_residues
        r["other_files"] = ";".join(sorted(f.name for f in d.iterdir()
                                           if not rslib.BATCH.search(f.name)
                                           and f.name not in ("samples.xtc", "samples_unphysical.xtc",
                                                              "topology.pdb", "sequence.fasta")))
        xtc, pdb = ENSEMBLES / f"{lab}.xtc", ENSEMBLES / f"{lab}.pdb"
        r["archive_conformers"] = n_frames(xtc, pdb) if pdb.exists() else None
        cur = ARCHIVE_CURVES / f"Pepsi-{lab}.csv"
        r["archive_curve_rows"] = len(pd.read_csv(cur, index_col=0)) if cur.exists() else None
        rows.append(r)
        print(f"{lab}: {r['n_samples_npz']} samples, samples.xtc {r['frames_samples_xtc']}, "
              f"unphysical {r['frames_unphysical_xtc']}, archive {r['archive_conformers']}", flush=True)
    inv = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    inv.to_csv(OUT / "inventory.csv", index=False)
    bad = inv[(inv["n_samples_npz"] != N_SAMPLES) | ~inv["npz_contiguous"] | ~inv["npz_sequence_matches"]
              | ~inv["fasta_matches"] | (inv["progress_model"] != MODEL) | ~inv["progress_ok"]
              | (inv["disorder_class"] != inv["class_in_list"]) | ~inv["clean"]
              | (inv["topology_residues"] != inv["length"])]
    record("01_inventory", {
        "samples_dir": str(SAMPLES), "progress_json_sha256": sha256(SAMPLES / "progress.json"),
        "resample_list_sha256": sha256(LIST), "cpus": used, "n_proteins": len(inv),
        "n_proteins_failing_checks": len(bad), "wall_seconds": round(time.time() - t0, 1)})
    if len(bad):
        print(bad.to_string(), file=sys.stderr)
        sys.exit(f"{len(bad)} proteins fail the inventory checks")


if __name__ == "__main__":
    main()
