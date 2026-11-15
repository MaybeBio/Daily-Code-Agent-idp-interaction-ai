"""End-to-end tests for ipsae_cli.py.

Covers what test_ipsae.py (which tests the library against hard-coded reference
values) cannot: the file readers, all-pairs looping over >2 chains, and the three
output files.

The expected ipSAE values for the real designs are the reference implementation's
(DunbrackLab/IPSAE v4), so this keeps CLI parity pinned without needing the
reference script on disk. Run with:

    python test_ipsae_cli.py
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(HERE, "ipsae_cli.py")
BASE = "/capstor/store/cscs/swissai/ab011/jkuhn"

# reference values, same source as test_ipsae.py
REAL_CASES = {
    "CEACAM19": {
        "pae": f"{BASE}/idp_bindcraft/known_binders_v2/CEACAM19_pae_matrix.npy",
        "pdb": f"{BASE}/Downloads_for_review/ecto_CEACAM19_top_design.pdb",
        "ipsae_ab": 0.783776, "ipsae_ba": 0.782044,
    },
    "Q8N535": {
        "pae": f"{BASE}/idp_bindcraft/known_binders_v2/Q8N535_pae_matrix.npy",
        "pdb": f"{BASE}/Downloads_for_review/daple_Q8N535_top_design.pdb",
        "ipsae_ab": 0.838177, "ipsae_ba": 0.249453,
    },
}

failures = []


def check(name, cond, detail=""):
    if not cond:
        failures.append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  -- ' + detail) if detail else ''}")


def run_cli(pae, struct, outdir, cutoff=10, dist=15):
    # stdout/stderr spelled out rather than capture_output= so this runs on py3.6
    res = subprocess.run(
        [sys.executable, CLI, pae, struct, str(cutoff), str(dist), "--outdir", outdir],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if res.returncode != 0:
        raise RuntimeError(f"CLI failed:\n{res.stdout}\n{res.stderr}")
    return res.stdout


def summary_rows(path):
    """Parse the summary file into {(chn1, chn2, type): {col: value}}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 10 or p[0] == "Chn1":
                continue
            out[(p[0], p[1], p[4])] = {
                "ipsae": float(p[5]), "ipsae_d0chn": float(p[6]),
                "ipsae_d0dom": float(p[7]), "iptm_d0chn": float(p[9]),
            }
    return out


def build_synthetic(d):
    """A 3-chain complex in PDB, mmCIF, AF2-json, AF3-json and Boltz-npz form."""
    rng = np.random.default_rng(7)
    lens = {"A": 30, "B": 25, "C": 20}
    chains = sum(([c] * n for c, n in lens.items()), [])
    n = len(chains)
    ch = np.array(chains)

    pae = rng.uniform(18, 31, (n, n))
    for c in lens:
        m = ch == c
        pae[np.ix_(m, m)] = rng.uniform(0.5, 3.0, (m.sum(), m.sum()))
    A, B, C = ch == "A", ch == "B", ch == "C"
    pae[np.ix_(A, B)] = rng.uniform(1.0, 8.0, (A.sum(), B.sum()))
    pae[np.ix_(B, A)] = rng.uniform(1.0, 8.0, (B.sum(), A.sum()))
    pae[np.ix_(A, C)] = rng.uniform(8.0, 14.0, (A.sum(), C.sum()))
    pae[np.ix_(C, A)] = rng.uniform(8.0, 14.0, (C.sum(), A.sum()))
    np.fill_diagonal(pae, 0.0)
    plddt = rng.uniform(70, 95, n)

    lines, idx, serial, atoms = [], 0, 1, []
    for c, ln in lens.items():
        off = {"A": 0.0, "B": 12.0, "C": 26.0}[c]
        for i in range(ln):
            x, y, z = off + 1.5 * (i % 5), 1.5 * i, 1.5 * (i % 3)
            for atom, el in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")):
                lines.append(
                    f"ATOM  {serial:5d}  {atom:<3s} ALA {c}{i+1:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{plddt[idx]:6.2f}           {el}\n")
                atoms.append((atom, "ALA", c, i + 1, x, y, z, plddt[idx]))
                serial += 1
            idx += 1
    lines.append("END\n")
    open(f"{d}/tri.pdb", "w").writelines(lines)

    with open(f"{d}/tri.cif", "w") as fh:
        fh.write("data_tri\n#\nloop_\n")
        for f in ("group_PDB", "id", "label_atom_id", "label_comp_id", "label_asym_id",
                  "auth_asym_id", "label_seq_id", "auth_seq_id", "Cartn_x", "Cartn_y",
                  "Cartn_z", "occupancy", "B_iso_or_equiv", "type_symbol"):
            fh.write(f"_atom_site.{f}\n")
        for i, (atom, res, c, rn, x, y, z, b) in enumerate(atoms, 1):
            fh.write(f"ATOM {i} {atom} {res} {c} {c} {rn} {rn} "
                     f"{x:.3f} {y:.3f} {z:.3f} 1.00 {b:.2f} {atom[0]}\n")
        fh.write("#\n")

    json.dump({"predicted_aligned_error": pae.tolist(), "plddt": plddt.tolist()},
              open(f"{d}/tri.json", "w"))
    json.dump({"pae": pae.tolist(), "atom_plddts": np.repeat(plddt, 5).tolist()},
              open(f"{d}/tri_af3.json", "w"))
    np.savez_compressed(f"{d}/pae_tri.npz", pae=pae)
    np.savez_compressed(f"{d}/plddt_tri.npz", plddt=plddt / 100.0)
    shutil.copy(f"{d}/tri.pdb", f"{d}/boltz_tri.pdb")
    return lens


def main():
    with tempfile.TemporaryDirectory() as tmp:
        # ---------------------------------------------------- real designs, AF2 path
        print("\n[AF2 json + pdb] reference-value parity on real designs")
        for name, case in REAL_CASES.items():
            if not (os.path.exists(case["pae"]) and os.path.exists(case["pdb"])):
                print(f"  SKIP {name}: inputs missing")
                continue
            d = os.path.join(tmp, name)
            os.makedirs(d)
            shutil.copy(case["pdb"], f"{d}/{name}.pdb")
            pae = np.load(case["pae"])
            plddt = [float(l[60:66]) for l in open(case["pdb"])
                     if l.startswith("ATOM") and l[12:16].strip() == "CA"]
            json.dump({"predicted_aligned_error": pae.tolist(), "plddt": plddt},
                      open(f"{d}/{name}.json", "w"))

            run_cli(f"{d}/{name}.json", f"{d}/{name}.pdb", d)
            rows = summary_rows(f"{d}/{name}_10_15.txt")
            ab = rows[("A", "B", "asym")]["ipsae"]
            ba = rows[("B", "A", "asym")]["ipsae"]
            mx = rows[("A", "B", "max")]["ipsae"]
            check(f"{name} A->B", abs(ab - case["ipsae_ab"]) < 5e-6,
                  f"{ab:.6f} vs {case['ipsae_ab']:.6f}")
            check(f"{name} B->A", abs(ba - case["ipsae_ba"]) < 5e-6,
                  f"{ba:.6f} vs {case['ipsae_ba']:.6f}")
            check(f"{name} max row = max of directions",
                  abs(mx - max(ab, ba)) < 1e-9)
            for suffix in (".txt", "_byres.txt", ".pml"):
                p = f"{d}/{name}_10_15{suffix}"
                check(f"{name} wrote {suffix}",
                      os.path.exists(p) and os.path.getsize(p) > 0)

        # ------------------------------------------------- all pairs over 3 chains
        print("\n[all-pairs] three chains -> three unordered pairs")
        d = os.path.join(tmp, "multi")
        os.makedirs(d)
        lens = build_synthetic(d)
        run_cli(f"{d}/tri.json", f"{d}/tri.pdb", d)
        rows = summary_rows(f"{d}/tri_10_15.txt")
        pairs = {(a, b) for (a, b, t) in rows if t == "max"}
        check("all three chain pairs scored", pairs == {("A", "B"), ("A", "C"), ("B", "C")},
              str(sorted(pairs)))
        check("both directions present for every pair",
              all((a, b, "asym") in rows and (b, a, "asym") in rows for a, b in pairs))
        ab = rows[("A", "B", "max")]["ipsae"]
        ac = rows[("A", "C", "max")]["ipsae"]
        bc = rows[("B", "C", "max")]["ipsae"]
        check("score ordering tracks the built-in interaction strength (A-B > A-C > B-C)",
              ab > ac > bc, f"AB={ab:.4f} AC={ac:.4f} BC={bc:.4f}")
        check("non-interacting pair scores exactly 0", bc == 0.0)

        nres = sum(lens.values())
        byres = [l for l in open(f"{d}/tri_10_15_byres.txt").readlines()[1:] if l.strip()]
        # one row per aligned residue per ordered chain pair
        expected = sum(lens[c1] for c1 in lens for c2 in lens if c1 != c2)
        check("byres has one row per residue per ordered pair",
              len(byres) == expected, f"{len(byres)} rows, expected {expected}")
        aliases = [l for l in open(f"{d}/tri_10_15.pml") if l.startswith("alias")]
        check("one PyMOL alias per ordered chain pair", len(aliases) == 6,
              f"{len(aliases)} aliases")
        check("no invalid 'resi None' in the PyMOL script",
              not any("resi None" in a for a in aliases))

        # ------------------------------------------- format equivalence
        print("\n[format readers] every input form must give the same numbers")
        ref_txt = open(f"{d}/tri_10_15.txt").read()
        for label, pae_file, struct in (
            ("mmCIF", f"{d}/tri.json", f"{d}/tri.cif"),
            ("AF3 per-atom pLDDT", f"{d}/tri_af3.json", f"{d}/tri.cif"),
            ("Boltz npz + 0-1 pLDDT", f"{d}/pae_tri.npz", f"{d}/boltz_tri.pdb"),
        ):
            sub = os.path.join(tmp, label.replace(" ", "_").replace("-", "_"))
            os.makedirs(sub, exist_ok=True)
            for f in (pae_file, struct):
                shutil.copy(f, sub)
            if pae_file.endswith(".npz"):
                shutil.copy(f"{d}/plddt_tri.npz", sub)
            run_cli(os.path.join(sub, os.path.basename(pae_file)),
                    os.path.join(sub, os.path.basename(struct)), sub)
            got = summary_rows(os.path.join(
                sub, os.path.basename(struct).rsplit(".", 1)[0] + "_10_15.txt"))
            want = summary_rows(f"{d}/tri_10_15.txt")
            same = all(abs(got[k]["ipsae"] - want[k]["ipsae"]) < 5e-6 for k in want)
            check(f"{label} matches the AF2 json + pdb path", same)

        # ------------------------------------------------------- cutoff plumbing
        print("\n[cutoffs]")
        run_cli(f"{d}/tri.json", f"{d}/tri.pdb", d, cutoff=15, dist=15)
        r15 = summary_rows(f"{d}/tri_15_15.txt")
        check("a looser PAE cutoff raises the marginal A-C pair",
              r15[("A", "C", "max")]["ipsae"] > ac,
              f"10A: {ac:.4f} -> 15A: {r15[('A','C','max')]['ipsae']:.4f}")
        check("output filename encodes both cutoffs",
              os.path.exists(f"{d}/tri_15_15.txt"))

    print(f"\n{'PASS' if not failures else 'FAILED: ' + ', '.join(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
