#!/usr/bin/env python3
"""
Paired region bootstrap at every matched layer.

Usage:
  python scripts/19_paired_by_layer.py
"""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REG = ROOT / "outputs" / "probe" / "regions"
OUT = ROOT / "outputs" / "probe" / "paired_by_layer.json"

AXES = ["disprot-molten_globule", "disprot-pre_molten_globule",
        "prediction-disorder-alphafold", "prediction-proline_rich-mobidb_lite_sub"]
LAYERS = [1, 3, 6, 9, 10, 12]


def main():
    out, tmp = {}, ROOT / "outputs" / "probe" / ".tmp_paired.json"
    for axis in AXES:
        out[axis] = {}
        for L in LAYERS:
            rp = REG / f"esm2_35m_layer{L}_{axis}.npz"
            cp = REG / f"seqdance_layer{L}_{axis}.npz"
            if not (rp.exists() and cp.exists()):
                print(f"  skip {axis} L{L}: missing dump")
                continue
            r = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "17_probe_paired.py"),
                 "--ref", str(rp), "--cmp", str(cp),
                 "--label", f"{axis} L{L}", "--out", str(tmp)],
                capture_output=True, text=True)
            if r.returncode != 0:
                print(f"  FAILED {axis} L{L}: {r.stderr.strip()[:200]}")
                continue
            d = json.load(open(tmp))
            out[axis][f"layer{L}"] = {
                "gap": d["gap"], "ci95": d["gap_ci95"],
                "p": d.get("boot_p_directional"), "sig": d["gap_excludes_0"],
                "n": d["n_paired_regions"],
            }
            flag = "SIG" if d["gap_excludes_0"] else " ns"
            print(f"{flag} {axis:45} L{L:<3} gap={d['gap']:+.3f} "
                  f"CI[{d['gap_ci95'][0]:+.3f},{d['gap_ci95'][1]:+.3f}] "
                  f"p={d.get('boot_p_directional'):.4f}")
    tmp.unlink(missing_ok=True)
    json.dump(out, open(OUT, "w"), indent=2)
    print(f"\nwrote → {OUT}")


if __name__ == "__main__":
    main()
