"""P(R_g) of the disordered chains of the human proteome, from the network and the closed form.

    python proteome/disorder.py first (in the metapredict environment), then
    python proteome/characterize.py [--max-n 600] [--min-idr-frac 0.9] [--sample K] [--jobs 8]

Two kinds of chain, from proteome/results/disorder_*.csv:
    protein   whole proteins with at least --min-idr-frac of their residues inside predicted IDRs
    idr       every predicted IDR, cut out of its protein
Both are kept between --min-n and --max-n residues; chains with residues other than the 20 standard ones (U, X) are
dropped. Writes proteome/results/predictions.csv (and the same table as .parquet): one row per chain with its
descriptors (h, the closed form's thresholded contact energy, included), both predictors' moments and distribution
parameters (net_*, cf_*; cf_A = exp(cf_lnA)), and where it sits relative to the training proteins:

    whole_protein  the chain is its entire protein; an IDR spanning a whole protein also appears as a `protein` row

    out_of_range   descriptors outside the training min-max (`;`-separated), N beyond 490 included
    outside_1_99   descriptors outside the training 1st-99th percentiles
    sticky         contact energy >= the training 90th percentile (the weak regime, MODEL.md)
    compact        the network's <Rg> below 0.8 x the Flory baseline (MD-compact chains are under-flagged: the network
                   over-predicts them)
    tier           confident: none of the above; extrapolated: anything out of range; caution: the rest
    md_split       the protein table's split when the identical sequence was simulated (train / val / test)
"""
import argparse
import pathlib
import sys
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import idprg  # noqa: E402,F401  (float64 in JAX)
from idprg import closed_form as cf, data, descriptors, distribution as dist, network  # noqa: E402

RES = HERE / "results"
STANDARD = set(descriptors.RESIDUES)
RANGED = ("N",) + descriptors.COLUMNS          # the descriptors of MODEL.md's domain table
MOMENTS = ("mean_rg", "std_rg", "skew", "exkurt", "nu", "S_m", "lnA")


def chains(min_n, max_n, min_idr_frac):
    P = pd.read_csv(RES / "disorder_proteins.csv")
    I = pd.read_csv(RES / "disorder_idrs.csv")
    whole = P[P.frac_in_idr >= min_idr_frac]
    seq = {r.id.split("|")[1]: str(r.seq) for r in __import__("Bio.SeqIO", fromlist=["x"]).parse(
        HERE.parent / "data" / "UP000005640_9606.fasta", "fasta")}
    a = whole.assign(kind="protein", chain=whole.accession, idr=0, start=0, end=whole.length,
                     sequence=whole.accession.map(seq))
    b = I.merge(P[["accession", "entry", "gene", "description", "reviewed", "length"]], on="accession",
                suffixes=("", "_protein"))
    b = b.assign(kind="idr", chain=b.accession + "_idr" + b.idr.astype(str))
    a["whole_protein"], b["whole_protein"] = True, (b.start == 0) & (b.end == b.length_protein)
    cols = ["chain", "kind", "accession", "entry", "gene", "description", "reviewed", "idr", "start", "end",
            "whole_protein", "mean_disorder", "sequence"]
    c = pd.concat([a[cols], b[cols]], ignore_index=True)
    n = c.sequence.str.len()
    std = c.sequence.map(lambda s: set(s) <= STANDARD)
    print(f"candidates: {(c.kind == 'protein').sum()} whole proteins (>= {min_idr_frac:.0%} in IDRs), "
          f"{(c.kind == 'idr').sum()} IDRs; {(n < min_n).sum()} shorter than {min_n}, {(n > max_n).sum()} longer than "
          f"{max_n}, {(~std).sum()} with non-standard residues", file=sys.stderr)
    return c[(n >= min_n) & (n <= max_n) & std].reset_index(drop=True)


def predictions(X):
    """Moments and parameters of both predictors, in blocks (the quadrature grid is 4,001 points per chain)."""
    net, spec = network.load()
    coef = cf.load()
    out = {}
    for tag, params in (("net", lambda t: network.parameters(net, spec, t)), ("cf", lambda t: cf.parameters(coef, t))):
        cols = []
        for b in range(0, len(X), 2000):
            t = X.iloc[b:b + 2000]
            lnA, lnSm, u = params(t)
            nu = dist.nu_of(u)
            m, s, k, e = dist.moments(lnA, lnSm, nu)
            cols.append(np.column_stack([m, s, k, e, nu, np.exp(lnSm), lnA]))
        out.update({f"{tag}_{c}": v for c, v in zip(MOMENTS, np.concatenate(cols).T)})
    return pd.DataFrame(out, index=X.index)


def domain(X, pred):
    """Flags against the training proteins (MODEL.md, domain of validity)."""
    P = data.proteins()
    tr = P[P.split == "train"]
    lo, hi = tr[list(RANGED)].min(), tr[list(RANGED)].max()
    p01, p99 = tr[list(RANGED)].quantile(0.01), tr[list(RANGED)].quantile(0.99)
    q90 = tr.contact_energy.quantile(0.90)
    names = lambda m: m.apply(lambda r: ";".join(r.index[r]), axis=1)
    V = X[list(RANGED)]
    f = pd.DataFrame(index=X.index)
    f["out_of_range"] = names((V < lo) | (V > hi))
    f["outside_1_99"] = names((V < p01) | (V > p99))
    f["sticky"] = X.contact_energy >= q90
    f["rg_flory_ratio"] = pred.net_mean_rg / np.exp(data.LN_A_F + data.NU_F * np.log(X.N))
    f["compact"] = f.rg_flory_ratio < 0.8
    f["tier"] = np.where(f.out_of_range != "", "extrapolated",
                         np.where((f.outside_1_99 != "") | f.sticky | f.compact, "caution", "confident"))
    f["md_split"] = X.sequence.map(P.set_index("sequence").split.to_dict()) if "sequence" in X else ""
    return f


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-n", type=int, default=20)
    ap.add_argument("--max-n", type=int, default=600, help="training proteins reach N = 490")
    ap.add_argument("--min-idr-frac", type=float, default=0.9, help="whole proteins: fraction of residues in IDRs")
    ap.add_argument("--sample", type=int, help="characterize only this many random chains (seed 0)")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", default=str(RES / "predictions.csv"))
    args = ap.parse_args(argv)

    C = chains(args.min_n, args.max_n, args.min_idr_frac)
    if args.sample:
        C = C.sample(args.sample, random_state=0).reset_index(drop=True)
    with ProcessPoolExecutor(args.jobs, mp_context=multiprocessing.get_context("spawn")) as pool:  # JAX is threaded
        X = pd.DataFrame(list(pool.map(descriptors.descriptors, C.sequence, chunksize=16)))
    X["sequence"] = C.sequence
    X.insert(X.columns.get_loc("mean_hydro") + 1, "h", cf.regressors(X, **cf.load()["ce_threshold"]).h)
    pred = predictions(X)
    pred.insert(pred.columns.get_loc("cf_lnA") + 1, "cf_A", np.exp(pred.cf_lnA))
    flags = domain(X, pred)
    out = pd.concat([C.drop(columns="sequence"), X.drop(columns="sequence"), pred, flags, C[["sequence"]]], axis=1)
    out["cf_vs_net_pct"] = 100 * (out.cf_mean_rg / out.net_mean_rg - 1)
    out.to_csv(args.out, index=False, float_format="%.6g")
    out.to_parquet(pathlib.Path(args.out).with_suffix(".parquet"), index=False)

    for kind, g in out.groupby("kind", sort=False):
        print(f"\n{kind}: {len(g)} chains, N {g.N.min()}-{g.N.max()} (median {g.N.median():.0f}); "
              f"{g.md_split.notna().sum()} already simulated", file=sys.stderr)
        print(g.tier.value_counts().to_string(), file=sys.stderr)
        print(f"closed form vs network <Rg>: median |diff| {g.cf_vs_net_pct.abs().median():.2f}%, "
              f"p95 {g.cf_vs_net_pct.abs().quantile(0.95):.2f}%", file=sys.stderr)
    print(f"-> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
