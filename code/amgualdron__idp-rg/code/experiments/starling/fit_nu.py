#!/usr/bin/env python
"""The two-exponent distribution fitted to each STARLING ensemble, the fit data/md_fits.csv and paper fig24 use for MD
(idprg.mdfit.fit_free: ln A, ln S_m and u free). Each ensemble is first divided by its own mean, so only shape enters.

    PYTHONPATH=. python experiments/starling/fit_nu.py runs/starling/benchmark --jobs 12

Writes nu_fits.csv to the run directory: nu and ln A of the fit, and nu_se, the sd of nu over resamples of the ensemble.
"""
import argparse
import multiprocessing
import os
import pathlib

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from idprg import mdfit
from idprg.distribution import nu_of

N_BOOT = 20


def fit_one(task):
    name, x = task
    x = x / x.mean()
    lnA, _, u, _ = mdfit.fit_free(x)
    rng = np.random.default_rng(abs(hash_name(name)))
    boot = []
    for _ in range(N_BOOT):
        b = rng.choice(x, len(x))
        boot.append(float(nu_of(mdfit.fit_free(b / b.mean())[2])))
    return dict(name=name, n=len(x), nu=float(nu_of(u)), lnA=float(lnA), nu_se=float(np.std(boot, ddof=1)))


def hash_name(name):
    return int.from_bytes(name.encode(), "little") % (2 ** 32)


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("run", type=pathlib.Path)
    a.add_argument("--jobs", type=int, default=4)
    args = a.parse_args()
    s = np.load(args.run / "rg_samples.npz")
    with ProcessPoolExecutor(args.jobs, mp_context=multiprocessing.get_context("spawn")) as pool:
        rows = list(pool.map(fit_one, [(n, np.asarray(s[n], float)) for n in s.files]))
    pd.DataFrame(rows).to_csv(args.run / "nu_fits.csv", index=False)


if __name__ == "__main__":
    main()
