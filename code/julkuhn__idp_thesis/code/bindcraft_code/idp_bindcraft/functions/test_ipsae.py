"""Regression test: our ipSAE must reproduce the reference implementation.

Expected values were produced by the reference script (DunbrackLab/IPSAE v4):

    python ipsae.py <pae>.json <design>.pdb 10 15

on two BindCraft designs that have stored PAE matrices. Run with:

    python test_ipsae.py
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipsae import compute_ipsae, compute_lis, compute_pdockq  # noqa: E402

BASE = "/capstor/store/cscs/swissai/ab011/jkuhn"

CASES = {
    "CEACAM19": {
        "pae": f"{BASE}/idp_bindcraft/known_binders_v2/CEACAM19_pae_matrix.npy",
        "pdb": f"{BASE}/Downloads_for_review/ecto_CEACAM19_top_design.pdb",
        "expect": {
            "ipsae_asym_AB": 0.783776, "ipsae_asym_BA": 0.782044, "ipsae": 0.783776,
            "ipsae_d0chn_asym_AB": 0.876252, "ipsae_d0chn_asym_BA": 0.884072,
            "ipsae_d0dom_asym_AB": 0.870033, "ipsae_d0dom_asym_BA": 0.879588,
            "iptm_d0chn_asym_AB": 0.876252, "iptm_d0chn_asym_BA": 0.809025,
            "n0res": 125, "n0chn": 250, "d0res": 4.14, "d0chn": 5.85,
            "pdockq": 0.4849, "pdockq2": 0.8695, "lis": 0.6867,
        },
    },
    "Q8N535": {
        "pae": f"{BASE}/idp_bindcraft/known_binders_v2/Q8N535_pae_matrix.npy",
        "pdb": f"{BASE}/Downloads_for_review/daple_Q8N535_top_design.pdb",
        "expect": {
            "ipsae_asym_AB": 0.838177, "ipsae_asym_BA": 0.249453, "ipsae": 0.838177,
            "ipsae_d0chn_asym_AB": 0.904579, "ipsae_d0chn_asym_BA": 0.713725,
            "ipsae_d0dom_asym_AB": 0.862179, "ipsae_d0dom_asym_BA": 0.667742,
            "iptm_d0chn_asym_AB": 0.904579, "iptm_d0chn_asym_BA": 0.288701,
            "n0res": 90, "n0chn": 157, "d0res": 3.43, "d0chn": 4.67,
            "pdockq": 0.3375, "pdockq2": 0.7104, "lis": 0.6587,
        },
    },
}


def parse_pdb(path):
    """Chain labels, CB coordinates (CA for glycine) and per-residue pLDDT, in PDB order."""
    order, rows = [], {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        key = (line[21], int(line[22:26]))
        if key not in rows:
            rows[key] = {"chain": line[21], "resname": line[17:20].strip()}
            order.append(key)
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        if name == "CA":
            rows[key]["ca"] = xyz
            rows[key]["plddt"] = float(line[60:66])
        elif name == "CB":
            rows[key]["cb"] = xyz

    chains, coords, plddt = [], [], []
    for key in order:
        r = rows[key]
        chains.append(r["chain"])
        coords.append(r["cb"] if r["resname"] != "GLY" and "cb" in r else r["ca"])
        plddt.append(r["plddt"])
    return chains, np.array(coords), np.array(plddt)


def main():
    failures = 0
    for name, case in CASES.items():
        if not (os.path.exists(case["pae"]) and os.path.exists(case["pdb"])):
            print(f"SKIP {name}: input files missing")
            continue

        pae = np.load(case["pae"])
        chains, coords, plddt = parse_pdb(case["pdb"])

        got = dict(compute_ipsae(pae, chains, "A", "B", pae_cutoff=10.0))
        got.update(compute_pdockq(pae, chains, coords, plddt, "A", "B"))
        got.update(compute_lis(pae, chains, "A", "B"))

        print(f"=== {name}")
        for key, want in case["expect"].items():
            # the reference prints to 6 dp for scores, 2 dp for d0, 4 dp for pDockQ/LIS
            tol = 0.5 if isinstance(want, int) else (0.005 if want > 1 else 5e-5)
            if key in ("pdockq", "pdockq2", "lis"):
                tol = 5e-5
            if key.startswith("d0"):
                tol = 5e-3
            ok = abs(float(got[key]) - float(want)) <= tol
            failures += not ok
            print(f"  {'ok  ' if ok else 'FAIL'} {key:24s} got={float(got[key]):.6f} want={float(want):.6f}")

    print("\nPASS" if not failures else f"\n{failures} MISMATCHES")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
