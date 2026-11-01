"""How far can residue order alone move P(R_g)? Composition-preserving Monte Carlo on long human IDRs rich in both
charge and hydropathy, scored by the network.

    python proteome/shuffle_search.py [--chains Q6ZQX7_idr2 ...] [--n-candidates 6] [--walkers 32] [--steps 400]
                                      [--random 1000] [--jobs 16]

Candidates (default): chains with 200 <= N <= 450, not extrapolated and not domain outliers (proteome/domain_map.py),
ranked by min(z(fcr), z(<lambda>)), i.e. rich in both; the two anticorrelate (charged residues have low lambda in HPS).

Every move keeps the amino-acid composition (so N, fcr, <lambda>, frac_sticky are fixed): swap two residues (1/2), move
a segment of up to N/4 residues elsewhere (1/4), reverse a segment (1/4). Proposals that leave the model's domain are
rejected: every descriptor within the training min-max, and the distance to the nearest training protein (standardized
descriptor space, as in domain_map.py) within the validation 99th percentile. Per candidate and objective (max / min
<R_g>, max / min skewness), --walkers independent walkers anneal from random in-domain shuffles over --steps steps,
Metropolis temperature 0.05 -> 0.001 (geometric) on ln <R_g> or on skewness. --random plain random shuffles give the
unoptimized spread. Native, random and optimized sequences are re-scored by both predictors (characterize.py), with
the usual domain flags. Writes to proteome/results/:

    shuffle_search.csv        native and best sequence per candidate and objective: moments, descriptors, flags
    shuffle_random.parquet    the random shuffles, with both predictors
    shuffle_trace.csv         best objective per step (convergence)
    shuffle_search.png        random shuffles in the (<R_g>, skewness) plane, native and optimized marked
    shuffle_search_pdfs.png   P(R_g) of native and optimized sequences
"""
import os

os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
import argparse  # noqa: E402
import multiprocessing  # noqa: E402
import pathlib  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import idprg  # noqa: E402,F401
from idprg import data, descriptors, distribution as dist  # noqa: E402

RES = HERE / "results"
OBJECTIVES = ("rg_max", "rg_min", "skew_max", "skew_min")
DIMS = ("log_seq_length",) + descriptors.COLUMNS
T0, T1 = 0.05, 0.001


class Domain:
    """Training min-max on every descriptor and the nearest-training-protein distance within the validation q."""

    def __init__(self, q=0.99):
        from scipy.spatial import cKDTree
        P = data.proteins()
        tr, va = P[P.split == "train"], P[P.split == "val"]
        self.lo, self.hi = tr[list(DIMS)].min().to_numpy(), tr[list(DIMS)].max().to_numpy()
        self.mu, self.sd = tr[list(DIMS)].mean().to_numpy(), tr[list(DIMS)].std().to_numpy()
        self.tree = cKDTree((tr[list(DIMS)].to_numpy() - self.mu) / self.sd)
        self.cut = np.quantile(self.tree.query((va[list(DIMS)].to_numpy() - self.mu) / self.sd)[0], q)

    def inside(self, X):
        V = X[list(DIMS)].to_numpy()
        d = self.tree.query((V - self.mu) / self.sd)[0]
        return np.all((V >= self.lo) & (V <= self.hi), axis=1) & (d <= self.cut), d


def propose(seq, rng):
    s, n = list(seq), len(seq)
    r = rng.random()
    if r < 0.5:
        i, j = rng.choice(n, 2, replace=False)
        s[i], s[j] = s[j], s[i]
    else:
        L = int(rng.integers(2, max(3, n // 4) + 1))
        i = int(rng.integers(0, n - L + 1))
        if r < 0.75:
            seg, rest = s[i:i + L], s[:i] + s[i + L:]
            k = int(rng.integers(0, len(rest) + 1))
            s = rest[:k] + seg + rest[k:]
        else:
            s[i:i + L] = s[i:i + L][::-1]
    return "".join(s)


def shuffle(seq, rng):
    a = np.array(list(seq))
    rng.shuffle(a)
    return "".join(a)


class Scorer:
    def __init__(self):
        from idprg import network
        self.net, self.spec = network.load()
        self.network = network
        self.domain = Domain()

    def __call__(self, seqs):
        """Descriptors, in-domain mask and network (mean_rg, skew) for a list of sequences."""
        X = pd.DataFrame([descriptors.descriptors(s) for s in seqs])
        ok, _ = self.domain.inside(X)
        m, k = np.full(len(seqs), np.nan), np.full(len(seqs), np.nan)
        if ok.any():
            lnA, lnSm, u = self.network.parameters(self.net, self.spec, X[ok])
            mm, _, kk, _ = dist.moments(lnA, lnSm, dist.nu_of(u))
            m[ok], k[ok] = mm, kk
        return ok, m, k


def value(objective, m, k):
    return {"rg_max": np.log(m), "rg_min": -np.log(m), "skew_max": k, "skew_min": -k}[objective]


def anneal(task):
    chain, seq, objective, walkers, steps, seed = task
    rng = np.random.default_rng(seed)
    score = Scorer()
    t0 = time.time()
    # starts: random in-domain shuffles (the native sequence if none is found)
    cur = []
    for _ in range(20):
        c = [shuffle(seq, rng) for _ in range(walkers)]
        ok, m, k = score(c)
        cur += [(s, mm, kk) for s, o, mm, kk in zip(c, ok, m, k) if o]
        if len(cur) >= walkers:
            break
    if not cur:
        _, m, k = score([seq])
        cur = [(seq, m[0], k[0])]
    cur = [cur[i % len(cur)] for i in range(walkers)]
    S = [c[0] for c in cur]
    f = np.array([value(objective, c[1], c[2]) for c in cur])
    best = (f.max(), S[int(f.argmax())])
    trace, acc = [], 0
    for t in range(steps):
        T = T0 * (T1 / T0) ** (t / max(steps - 1, 1))
        prop = [propose(s, rng) for s in S]
        ok, m, k = score(prop)
        g = np.where(ok, value(objective, m, k), -np.inf)
        take = ok & (np.log(rng.random(walkers)) < (g - f) / T)
        for w in np.flatnonzero(take):
            S[w], f[w] = prop[w], g[w]
        acc += take.sum()
        if f.max() > best[0]:
            best = (f.max(), S[int(f.argmax())])
        trace.append(best[0])
    print(f"{chain} {objective}: best {best[0]:.4f}, acceptance {acc / (walkers * steps):.2f}, "
          f"{time.time() - t0:.0f} s", file=sys.stderr, flush=True)
    return chain, objective, best[1], trace


def randoms(task):
    chain, seq, n, seed = task
    rng = np.random.default_rng(seed)
    dom = Domain()
    seqs = [shuffle(seq, rng) for _ in range(n)]
    X = pd.DataFrame([descriptors.descriptors(s) for s in seqs])
    X["in_domain"], X["nn_dist"] = dom.inside(X)
    X["sequence"], X["chain"] = seqs, chain
    return X


def candidates(k):
    from rg_groups import chains
    D = chains(False, "none")
    D = D[~D.chain.isin(pd.read_csv(RES / "domain_outliers.csv").chain)]
    L = D[(D.N >= 200) & (D.N <= 450)]
    z = lambda v: (v - v.mean()) / v.std()
    return L.assign(score=np.minimum(z(L.fcr), z(L.mean_hydro))).sort_values("score", ascending=False).head(k)


def rescore(seqs, extra):
    """Both predictors and the domain flags of characterize.py for a list of sequences."""
    from characterize import domain, predictions
    from idprg import closed_form as cf
    X = pd.DataFrame([descriptors.descriptors(s) for s in seqs])
    X["h"] = cf.regressors(X, **cf.load()["ce_threshold"]).h
    X["sequence"] = seqs
    pred = predictions(X)
    flags = domain(X, pred)
    _, X["nn_dist"] = Domain().inside(X)
    return pd.concat([extra.reset_index(drop=True), X, pred, flags.drop(columns="md_split")], axis=1)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chains", nargs="+")
    ap.add_argument("--n-candidates", type=int, default=6)
    ap.add_argument("--walkers", type=int, default=32)
    ap.add_argument("--steps", type=int, default=400)
    ap.add_argument("--random", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    if args.chains:
        D = pd.read_parquet(RES / "predictions.parquet").drop_duplicates("chain").set_index("chain")
        C = D.loc[args.chains].reset_index()
    else:
        C = candidates(args.n_candidates)
    print(C[["chain", "gene", "N", "fcr", "mean_hydro", "net_mean_rg", "net_skew"]].round(3).to_string(index=False),
          file=sys.stderr)

    ctx = multiprocessing.get_context("spawn")
    tasks = [(c.chain, c.sequence, o, args.walkers, args.steps, args.seed + 100 * i + j)
             for i, c in enumerate(C.itertuples()) for j, o in enumerate(OBJECTIVES)]
    rtasks = [(c.chain, c.sequence, args.random // 4, args.seed + 7919 * i + b)
              for i, c in enumerate(C.itertuples()) for b in range(4)]
    with ProcessPoolExecutor(args.jobs, mp_context=ctx) as pool:
        rfut = pool.map(randoms, rtasks)
        out = list(pool.map(anneal, tasks))
        R = pd.concat(list(rfut), ignore_index=True)

    # rescore everything with both predictors
    info = C.set_index("chain")[["gene", "N"]]
    rows = [dict(chain=c.chain, gene=c.gene, kind="native") for c in C.itertuples()]
    seqs = list(C.sequence)
    for chain, obj, s, _ in out:
        rows.append(dict(chain=chain, gene=info.loc[chain, "gene"], kind=obj))
        seqs.append(s)
    S = rescore(seqs, pd.DataFrame(rows))
    S.to_csv(RES / "shuffle_search.csv", index=False, float_format="%.5g")
    Rs = rescore(list(R.sequence), R[["chain", "in_domain"]])
    Rs.to_parquet(RES / "shuffle_random.parquet", index=False)
    pd.DataFrame([dict(chain=c, objective=o, step=t, best=v) for c, o, _, tr in out for t, v in enumerate(tr)]).to_csv(
        RES / "shuffle_trace.csv", index=False, float_format="%.5g")

    cols = ["gene", "kind", "net_mean_rg", "net_skew", "net_std_rg", "net_nu", "cf_mean_rg", "cf_skew", "scd_pos",
            "scd_neg", "h", "nn_dist", "tier"]
    print(S[cols].round(3).to_string(index=False))
    print("\nrandom shuffles (in domain only), per chain: network <R_g> and skew 1-99%")
    for ch, g in Rs[Rs.in_domain].groupby("chain", sort=False):
        print(f"{ch}: in domain {len(g)}/{(Rs.chain == ch).sum()}; <R_g> {g.net_mean_rg.quantile(.01):.2f}-"
              f"{g.net_mean_rg.quantile(.99):.2f}; skew {g.net_skew.quantile(.01):.2f}-{g.net_skew.quantile(.99):.2f}")


if __name__ == "__main__":
    main()
