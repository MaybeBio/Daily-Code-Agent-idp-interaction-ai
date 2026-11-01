"""How often two human disordered chains with the same mean R_g differ in skewness by a set amount.

    python proteome/skew_pairs.py [--tol 0.01] [--dskew 0.2] [--keep-extrapolated]

Pairs are all pairs of distinct chains (one per sequence, proteome/results/predictions.parquet) whose predicted <R_g>
agree within --tol (|ln <R_g>_i - ln <R_g>_j| <= ln(1 + tol)). For each predictor, the fraction of those pairs with
|skew_i - skew_j| >= --dskew, overall and by <R_g> bin, and the fraction of chains with at least one such partner.
Writes proteome/results/skew_pairs.csv.
"""
import argparse
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from rg_groups import chains  # noqa: E402

RES = HERE / "results"
BINS = [1, 1.5, 2, 2.5, 3, 4, 5, 6, 10]


def pairs(rg, skew, tol, dskew):
    """Per chain i: number of partners j > i (in <R_g> order) within tol, and how many of them differ by >= dskew."""
    o = np.argsort(rg)
    lr, sk = np.log(rg[o]), skew[o]
    hi = np.searchsorted(lr, lr + np.log1p(tol), side="right")
    n, big = np.zeros(len(lr), int), np.zeros(len(lr), int)
    any_big = np.zeros(len(lr), bool)
    for i in range(len(lr)):
        d = np.abs(sk[i + 1:hi[i]] - sk[i]) >= dskew
        n[i], big[i] = d.size, d.sum()
        if big[i]:
            any_big[i] = True
            any_big[i + 1:hi[i]] |= d
    inv = np.empty_like(o)
    inv[o] = np.arange(len(o))
    return n[inv], big[inv], any_big[inv]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--dskew", type=float, default=0.2)
    ap.add_argument("--keep-extrapolated", action="store_true")
    args = ap.parse_args(argv)

    D = chains(args.keep_extrapolated, "none")
    rows = []
    for m in ("cf", "net"):
        rg, sk = D[f"{m}_mean_rg"].to_numpy(), D[f"{m}_skew"].to_numpy()
        n, big, any_big = pairs(rg, sk, args.tol, args.dskew)
        b = pd.cut(rg, BINS)
        for lab, g in pd.DataFrame(dict(b=b, n=n, big=big, any_big=any_big)).groupby("b", observed=True):
            rows.append(dict(model=m, rg_bin=str(lab), chains=len(g), pairs=g.n.sum(), pairs_dskew=g.big.sum(),
                             frac_pairs=g.big.sum() / max(g.n.sum(), 1), frac_chains_with_partner=g.any_big.mean()))
        rows.append(dict(model=m, rg_bin="all", chains=len(D), pairs=n.sum(), pairs_dskew=big.sum(),
                         frac_pairs=big.sum() / n.sum(), frac_chains_with_partner=any_big.mean()))
        print(f"{m}: skew 1-99% {np.quantile(sk, .01):.2f}-{np.quantile(sk, .99):.2f}", file=sys.stderr)
    S = pd.DataFrame(rows)
    S.to_csv(RES / "skew_pairs.csv", index=False, float_format="%.4g")
    print(f"pairs with <R_g> within {100 * args.tol:g}%, |d skew| >= {args.dskew}  ({len(D)} chains)")
    print(S.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
