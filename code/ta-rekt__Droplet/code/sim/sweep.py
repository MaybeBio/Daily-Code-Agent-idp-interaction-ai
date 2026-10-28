"""
Python sweep driver, cognate to sim/sweep.jl.

Loads a spec file that defines:
    alphas, betas, gammas    -- 1D parameter axes (alphas is the parallelized axis)
    nparams                   -- len(alphas)
    fixed                     -- dict of fixed Parametebetas kwargs + solver settings
    blurb                     -- str
    target(res, keya, alpha, keyb, R, keyc, phi_inf, fixed, blurb) -> dict

Slices alphas across CPUs and writes data/<...>/output_cpu<cpukey>.json.
"""

import sys
import os
import json
import importlib.util
from threadpoolctl import threadpool_limits
import numpy as np

_blas_threads = int(os.environ.get("BLAS_NUM_THREADS", "1"))
threadpool_limits(_blas_threads)


def _to_jsonable(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(f"not JSON-serializable: {type(o)}")


def _load_spec(path):
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    spec_path = path if os.path.isabs(path) else os.path.join(repo_root, path)
    spec = importlib.util.spec_from_file_location("sweep_spec", spec_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    specfile = sys.argv[1]
    cpukey = int(sys.argv[2])
    ncpus = int(sys.argv[3])
    outdir = sys.argv[4]

    os.makedirs(outdir, exist_ok=True)

    spec = _load_spec(specfile)

    alphas = spec.alphas
    betas = spec.betas
    gammas = spec.gammas
    nparams = spec.nparams
    fixed = spec.fixed
    blurb = spec.blurb
    target = spec.target

    nparams_per_cpu = nparams // ncpus

    results = {}

    for k in range(1, nparams_per_cpu + 1):
        keya = (cpukey - 1) * nparams_per_cpu + k
        alpha = float(alphas[keya - 1])
        results[keya] = {}

        for keyb, beta in enumerate(betas, start=1):
            results[keya][keyb] = {}
            for keyc, gamma in enumerate(gammas, start=1):
                res = {}
                results[keya][keyb][keyc] = target(
                    res, keya, alpha, keyb, float(beta), keyc, float(gamma), fixed, blurb
                )

    outfile = os.path.join(outdir, f"output_cpu{cpukey}.json")
    with open(outfile, "w") as f:
        json.dump(results, f, indent=4, default=_to_jsonable)


if __name__ == "__main__":
    main()
