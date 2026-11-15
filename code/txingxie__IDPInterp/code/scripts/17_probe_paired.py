#!/usr/bin/env python3
"""
Paired region-level bootstrap of the probe gap between two models.

Usage:
  python scripts/17_probe_paired.py \
     --ref outputs/probe/regions/esm2_8m_layer4_disprot-molten_globule.npz \
     --cmp outputs/probe/regions/seqdance_layer10_disprot-molten_globule.npz

  # each model's best layer, which is a post-hoc pairing and reported as such
  python scripts/17_probe_paired.py \
     --ref outputs/probe/regions/esm2_8m_layer1_disprot-molten_globule.npz \
     --cmp outputs/probe/regions/seqdance_layer10_disprot-molten_globule.npz \
     --label "best layer of each model"
"""
import argparse
import json
from pathlib import Path

import numpy as np


def _f1_from(hit, tp, fp):
    """Identical to scripts/12b and 16."""
    R = hit.mean() if hit.size else 0.0
    TP = tp.sum()
    P = TP / (TP + fp + 1e-10)
    return 2 * P * R / (P + R + 1e-10)


def load(path):
    z = np.load(path, allow_pickle=False)
    keys = [(u, int(s), int(e)) for u, s, e in zip(z["uid"], z["start"], z["end"])]
    return keys, z["hit"], z["tp"], int(z["fp"][0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="baseline model npz (e.g. ESM-2)")
    ap.add_argument("--cmp", required=True, help="comparison model npz (e.g. SeqDance)")
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--label", default="")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    kr, hr, tr, fr = load(a.ref)
    kc, hc, tc, fc = load(a.cmp)

    # Pair on identical biological regions; both dumps come from the same corpus/axis so
    # this should be a full intersection. Any shortfall is reported rather than hidden.
    idx_r = {k: i for i, k in enumerate(kr)}
    idx_c = {k: i for i, k in enumerate(kc)}
    common = sorted(set(idx_r) & set(idx_c))
    if len(common) < 2:
        raise SystemExit(f"Only {len(common)} shared regions — cannot pair.")
    if len(common) != len(kr) or len(common) != len(kc):
        print(f"  note: pairing {len(common)} regions "
              f"(ref had {len(kr)}, cmp had {len(kc)})")

    ir = np.array([idx_r[k] for k in common])
    ic = np.array([idx_c[k] for k in common])
    hr, tr = hr[ir], tr[ir]
    hc, tc = hc[ic], tc[ic]
    n = len(common)

    f1_r, f1_c = _f1_from(hr, tr, fr), _f1_from(hc, tc, fc)

    rng = np.random.default_rng(a.seed)
    gap = np.empty(a.n_boot)
    for i in range(a.n_boot):
        j = rng.integers(0, n, n)          # ONE resample, applied to both models
        gap[i] = _f1_from(hc[j], tc[j], fc) - _f1_from(hr[j], tr[j], fr)
    lo, hi = np.percentile(gap, [2.5, 97.5])
    p_cmp_le = float((gap <= 0).mean())
    # Two-sided: the compositional axes are EXPECTED to come out negative (ESM-2 ahead),
    # so a one-sided "is cmp > ref" test would report the confirming result as a failure.
    excludes = bool(lo > 0 or hi < 0)
    direction = "cmp>ref" if (f1_c - f1_r) > 0 else "ref>cmp"
    p_dir = min(p_cmp_le, 1.0 - p_cmp_le)      # one-sided in the OBSERVED direction

    res = {
        "label": a.label or f"{Path(a.cmp).stem} vs {Path(a.ref).stem}",
        "ref": Path(a.ref).stem, "cmp": Path(a.cmp).stem,
        "n_paired_regions": n,
        "ref_f1": round(float(f1_r), 4), "cmp_f1": round(float(f1_c), 4),
        "gap": round(float(f1_c - f1_r), 4),
        "gap_ci95": [round(float(lo), 4), round(float(hi), 4)],
        "gap_excludes_0": excludes,
        "direction": direction,
        "boot_p_cmp_le_ref": round(p_cmp_le, 4),
        "boot_p_directional": round(p_dir, 4),
        "method": "paired region bootstrap over out-of-fold probe outcomes "
                  "(same resample applied to both models; bg FP fixed per model)",
    }
    print(json.dumps(res, indent=2))
    print(f"\n  {res['ref']}  F1={f1_r:.3f}")
    print(f"  {res['cmp']}  F1={f1_c:.3f}")
    print(f"  gap = {f1_c - f1_r:+.3f}  95%CI[{lo:+.3f},{hi:+.3f}]  "
          f"excludes 0: {excludes} ({direction})   "
          f"(paired bootstrap p={p_dir:.4f} one-sided in observed direction, n={n})")

    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        json.dump(res, open(a.out, "w"), indent=2)
        print(f"  saved → {a.out}")


if __name__ == "__main__":
    main()
