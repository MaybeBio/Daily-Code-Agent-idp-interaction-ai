#!/usr/bin/env python3
"""
Regenerate the main paired probe tests from the per-region .npz dumps.

Usage:
  python scripts/18_rebuild_paired.py
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REG = ROOT / "outputs" / "probe" / "regions"
OUT = ROOT / "outputs" / "probe"

MOLTEN = "disprot-molten_globule"
PREMOLTEN = "disprot-pre_molten_globule"
ALPHAFOLD = "prediction-disorder-alphafold"
PROLINE = "prediction-proline_rich-mobidb_lite_sub"

# (label, ref stem, cmp stem, axis, output filename)
COMPARISONS = [
    # conformational: SeqDance expected ahead
    ("MOLTEN matched: ESM-2-35M L10 vs SeqDance L10",
     "esm2_35m_layer10", "seqdance_layer10", MOLTEN, "paired_matched_L10.json"),
    ("MOLTEN paper layers: ESM-2-8M L4 vs SeqDance L10",
     "esm2_8m_layer4", "seqdance_layer10", MOLTEN, "paired_paper_layers.json"),
    ("MOLTEN best layer of each model",
     "esm2_8m_layer1", "seqdance_layer10", MOLTEN, "paired_best_layers.json"),
    ("PRE-MOLTEN matched: ESM-2-35M L10 vs SeqDance L10",
     "esm2_35m_layer10", "seqdance_layer10", PREMOLTEN, "paired_premolten_matched.json"),
    ("PRE-MOLTEN paper layers: ESM-2-8M L4 vs SeqDance L10",
     "esm2_8m_layer4", "seqdance_layer10", PREMOLTEN, "paired_premolten_paper.json"),
    # compositional: ESM-2 expected ahead (the sign flip)
    ("COMPOSITIONAL AlphaFold-disorder: 35M L10 vs SeqDance L10",
     "esm2_35m_layer10", "seqdance_layer10", ALPHAFOLD, "paired_comp_alphafold.json"),
    ("COMPOSITIONAL proline-rich: 35M L10 vs SeqDance L10",
     "esm2_35m_layer10", "seqdance_layer10", PROLINE, "paired_comp_proline.json"),
    # sequence-composition baseline
    ("MOLTEN vs composition baseline: SeqDance L10",
     "esm2_8m_composition", "seqdance_layer10", MOLTEN, "paired_vs_composition_molten.json"),
    ("PRE-MOLTEN vs composition baseline: SeqDance L10",
     "esm2_8m_composition", "seqdance_layer10", PREMOLTEN, "paired_vs_composition_premolten.json"),
    # does ESM-2 beat the composition baseline at all? (it barely does / does not)
    ("MOLTEN ESM-2-35M L10 vs composition",
     "esm2_8m_composition", "esm2_35m_layer10", MOLTEN, "paired_esmVScomp_esm35m_molten.json"),
    ("MOLTEN ESM-2-8M L4 vs composition",
     "esm2_8m_composition", "esm2_8m_layer4", MOLTEN, "paired_esmVScomp_esm8m_molten.json"),
    ("PRE-MOLTEN ESM-2-35M L10 vs composition",
     "esm2_8m_composition", "esm2_35m_layer10", PREMOLTEN, "paired_esmVScomp_esm35m_pre.json"),
    ("MOLTEN ESM-2-8M L1 (best ESM-2 layer) vs composition",
     "esm2_8m_composition", "esm2_8m_layer1", MOLTEN, "paired_esmVScomp_esm8m_L1_molten.json"),
]


def main():
    made, skipped = 0, []
    for label, ref, cmp_, axis, outname in COMPARISONS:
        rp = REG / f"{ref}_{axis}.npz"
        cp = REG / f"{cmp_}_{axis}.npz"
        if not rp.exists() or not cp.exists():
            missing = [p.name for p in (rp, cp) if not p.exists()]
            skipped.append((label, missing))
            continue
        cmd = [sys.executable, str(ROOT / "scripts" / "17_probe_paired.py"),
               "--ref", str(rp), "--cmp", str(cp), "--label", label,
               "--out", str(OUT / outname)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  FAILED {label}\n{r.stderr.strip()}")
            continue
        made += 1

    print(f"\nrebuilt {made} paired result(s) → {OUT}")
    for label, missing in skipped:
        print(f"  skipped (no dumps yet): {label}  [missing: {', '.join(missing)}]")

    import glob, json
    print("\n=== paired summary ===")
    for f in sorted(glob.glob(str(OUT / "paired_*.json"))):
        d = json.load(open(f))
        flag = "SIG" if d["gap_excludes_0"] else " ns"
        print(f"{flag} {d['label'][:52]:52} gap={d['gap']:+.3f} "
              f"CI[{d['gap_ci95'][0]:+.3f},{d['gap_ci95'][1]:+.3f}] "
              f"p={d.get('boot_p_directional', float('nan')):.4f} {d.get('direction','')}")


if __name__ == "__main__":
    main()
