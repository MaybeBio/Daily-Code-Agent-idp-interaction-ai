"""Cα radius of gyration of every conformer on both sides (design §3).

    python 03_rg.py

v1.2: every recovered sample (RESAMPLE_WORK/v12/<label>/samples_all.xtc). Archive: the archive's BioEmu
conformers under BSC_ROOT. The Cα atoms only, equal weights. Output: results/resample/rg_conformers.csv
(label, side, conformer, rg_ca in Å), in the conformer order of results/resample/filter.csv.
"""

from __future__ import annotations

import time

import mdtraj as md
import pandas as pd
import rslib
from common import ENSEMBLES, OUT, V12, labels, pin, record


def ca_rg(xtc, top) -> tuple[list[float], int]:
    t = md.load(str(xtc), top=str(top))
    ca = t.atom_slice(t.topology.select("name CA"))
    return (10.0 * rslib.rg_per_conformer(ca.xyz)).tolist(), ca.n_atoms


def main() -> None:
    t0 = time.time()
    used = pin()
    rows = []
    for lab in labels():
        for side, xtc, top in (("v12", V12 / lab / "samples_all.xtc", V12 / lab / "topology.pdb"),
                               ("archive", ENSEMBLES / f"{lab}.xtc", ENSEMBLES / f"{lab}.pdb")):
            rg, n_ca = ca_rg(xtc, top)
            rows += [{"label": lab, "side": side, "conformer": i, "rg_ca": v, "n_ca": n_ca}
                     for i, v in enumerate(rg)]
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "rg_conformers.csv", index=False)
    record("03_rg", {"mdtraj": md.__version__, "cpus": used, "wall_seconds": round(time.time() - t0, 1)})


if __name__ == "__main__":
    main()
