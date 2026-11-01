"""Score the fits of fit.py, paired per protein against the published form and against s+/N (reference/).

    python experiments/u_terms/compare.py [--against TAG] [--dir runs/u_terms]

--against TAG: also pair each model against that one, e.g. dS when judging dS_no_scd_pos. Standard errors treat one
sequence cluster as one unit, as the c1 rule did. Cross-validation is on training proteins only; validation is shown
only if a standard:-1 fit exists.
"""
import argparse
import glob
import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from idprg import data, distribution as dist  # noqa: E402

P = data.proteins()[["id", "split", "cluster", "scd_pos", "scd_neg"]]
MD_U = data.md_fits().set_index("id").u
SUBSETS = {"all": lambda d: d.N > 0, "s+>.3,N<=50": lambda d: (d.scd_pos > .3) & (d.N <= 50), "s+>.3,N<=100": lambda d: (d.scd_pos > .3) & (d.N <= 100),
           "s+>.3,N>200": lambda d: (d.scd_pos > .3) & (d.N > 200), "s-<-.05": lambda d: d.scd_neg < -.05, "s~0": lambda d: (d.scd_pos + d.scd_neg).abs() < .01}


def load(paths):
    d = pd.concat([pd.read_csv(p) for p in paths]).merge(P, on="id")
    if (d.split == "test").any():
        sys.exit("test-split rows found: the test split is sealed")
    d["nu_md"] = dist.nu_of(MD_U.reindex(d.id).to_numpy(float))
    return d.sort_values("id").reset_index(drop=True)


def paired(a, b):
    """Mean NLL change a - b, its cluster standard error, proteins improved, on the proteins both have."""
    m = a[["id", "nll", "cluster"]].merge(b[["id", "nll"]], on="id", suffixes=("", "_b"))
    x = m.nll - m.nll_b
    g = x.groupby(m.cluster).sum()
    se = np.sqrt((g ** 2).sum() * len(g) / max(len(g) - 1, 1)) / max(len(x), 1)
    return x.mean(), se, int((x < 0).sum()), len(x)


def row(d):
    e = d.nu - d.nu_md
    r = np.corrcoef(d.nu, d.nu_md)[0, 1] if len(d) > 2 else np.nan
    return (f"NLL {d.nll.mean():.5f}  Rg {d.mean_err_pct.abs().mean():.3f}  sigma {d.std_err_pct.abs().mean():.3f}  W1 {d.w1_pct.median():.3f}  "
            f"nu r {r:.3f} |err| {e.abs().mean():.4f} sd {d.nu.std():.3f} (MD {d.nu_md.std():.3f}) max {d.nu.max():.3f}")


def report(tag, d, refs):
    print(f"\n=== {tag}  n={len(d)}")
    for lab, sel in SUBSETS.items():
        m = sel(d)
        if m.sum() < 3:
            continue
        print(f"  {lab:13s} n={m.sum():4d}  {row(d[m])}")
        for ref, b in refs.items():
            ids = d.id[m]
            c, se, better, n = paired(d[m], b[b.id.isin(ids)])
            print(f"  {'':13s} vs {ref:12s} dNLL {c:+.5f} se {se:.5f} ({c / se if se else 0:+.1f} SE), {better}/{n} better")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--against")
    ap.add_argument("--dir", default=str(ROOT / "runs" / "u_terms"))
    args = ap.parse_args()
    R = HERE / "reference"
    refs = {"published": load([R / "published" / "closed_form_cv.csv"]), "s+/N": load(sorted(glob.glob(str(R / "s_over_N.standard.[012].csv"))))}
    print("baselines over all three folds:")
    for k, b in refs.items():
        print(f"  {k:10s} n={len(b)}  {row(b)}")
    files = sorted(glob.glob(f"{args.dir}/*.standard.[012].csv"))
    tags = sorted({re.sub(r"\.standard\.\d\.csv$", "", pathlib.Path(f).name) for f in files})
    if args.against:
        refs[args.against] = load(sorted(glob.glob(f"{args.dir}/{args.against}.standard.[012].csv")))
    for tag in tags:
        fs = sorted(glob.glob(f"{args.dir}/{tag}.standard.[012].csv"))
        js = [json.loads(pathlib.Path(f[:-3] + "json").read_text()) for f in fs]
        report(f"{tag} ({len(fs)} of 3 folds, {js[0]['n_parameters']} parameters)", load(fs), refs)
        print("  objective per fold", [round(j["objective_end"], 6) for j in js], " iterations", [j["iterations"] for j in js])
        if len(js[0].get("shape_params", [])) > 1:
            print("  shape parameters per fold", [[round(v, 3) for v in j["shape_params"]] for j in js])
        for k in js[0]["terms"]:
            print(f"  {k}: delta per fold", [round(j["terms"][k]["delta"], 3) for j in js], " coef", [round(j["terms"][k]["coef"], 3) for j in js])
    print("\nextrapolation sets: the c1 rule allows at most +0.001 NLL against the published form")
    for s in ("length", "sticky"):
        base = load([HERE / "reference" / "published" / f"closed_form_{s}.csv"])
        for f in sorted(glob.glob(f"{args.dir}/*.{s}.-1.csv")):
            d = load([f])
            c, se, better, n = paired(d, base[base.id.isin(d.id)])
            top = (d.merge(base[["id", "nll"]], on="id", suffixes=("", "_b")).eval("nll - nll_b")).max()
            print(f"  {pathlib.Path(f).stem:30s} {s}: n={n} dNLL {c:+.5f} se {se:.5f}, {better}/{n} better, worst protein {top:+.3f}, Rg {d.mean_err_pct.abs().mean():.3f}, sigma {d.std_err_pct.abs().mean():.3f}")
    for f in sorted(glob.glob(f"{args.dir}/*.standard.-1.csv")):
        report(f"{pathlib.Path(f).stem} on VALIDATION", load([f]), {"published": load([R / "published" / "closed_form_val.csv"]), "s+/N": load([R / "s_over_N.standard.-1.csv"])})


if __name__ == "__main__":
    main()
