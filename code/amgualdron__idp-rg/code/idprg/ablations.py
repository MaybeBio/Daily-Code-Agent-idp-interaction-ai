"""The ablation grid (results/ablations.csv). Each row retrains the network or refits the closed form with one change
from the released configuration and scores its evaluation set with metrics.summary.

network, seed 42 unless stated, scored on the validation, length and sticky splits:
    reference        the released configuration; seeds 42-46 give the seed spread
    -<input>         one of the eight descriptors left out (ln N is always kept)
    ln N only        no sequence descriptor at all
    <width>x<depth>  128x2, 1024x2, 512x1, 512x3 (validation only)
closed form, validation only:
    reference, -scd_pos, -scd_neg, -mean_hydro, -h (each descriptor removed everywhere it enters), no ln N interactions,
    constant CE threshold (c1 = 0), free ln A (linear in the descriptors instead of the law)
"""
import json
import multiprocessing
import os
import pathlib
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from . import data
from .descriptors import INPUTS

SETS = ("standard", "length", "sticky")
SET_NAME = {"standard": "val", "length": "length", "sticky": "sticky"}
KEY = ["model", "ablation", "seed", "set"]


def tasks():
    t = []
    for split in SETS:
        t += [dict(model="network", ablation="reference", seed=s, split=split) for s in range(42, 47)]
        t += [dict(model="network", ablation=f"-{x}", seed=42, split=split, inputs=[i for i in INPUTS if i != x])
              for x in INPUTS[1:]]
        t.append(dict(model="network", ablation="ln N only", seed=42, split=split, inputs=["log_seq_length"]))
    t += [dict(model="network", ablation=f"{w}x{d}", seed=42, split="standard", width=w, depth=d)
          for w, d in ((128, 2), (1024, 2), (512, 1), (512, 3))]
    cf = [("reference", {}), ("-scd_pos", dict(drop=["scd_pos"])), ("-scd_neg", dict(drop=["scd_neg"])),
          ("-mean_hydro", dict(drop=["mean_hydro"])), ("-h", dict(drop=["h"])),
          ("no ln N interactions", dict(interactions=False)), ("constant CE threshold", dict(c1=0.0)),
          ("free ln A", dict(law=False))]
    t += [dict(model="closed_form", ablation=a, seed=0, split="standard", fit=kw) for a, kw in cf]
    return t


def _run(task):
    from . import closed_form as cf, metrics, network
    import equinox as eqx
    import jax
    P, M, F = data.proteins(), data.md_fits(), data.frames()
    fit, ev = data.split(P, task["split"])
    t0 = time.time()
    if task["model"] == "network":
        w, d = task.get("width", 512), task.get("depth", 2)
        net, spec = network.train(P, M, fit, F, task["split"], task.get("inputs", INPUTS), w, d, seed=task["seed"],
                                  log=lambda s: None)
        params = network.parameters(net, spec, P.iloc[ev])
        n_par = sum(x.size for x in jax.tree_util.tree_leaves(eqx.filter(net, eqx.is_array)))
        extra = dict(width=w, depth=d)
    else:
        coef, info = cf.fit(P, M, fit, F, **task["fit"])
        params = cf.parameters(coef, P.iloc[ev])
        n_par, extra = info["n_parameters"], dict(objective=info["objective_end"])
    s = metrics.summary(metrics.score(*params, P.iloc[ev], F))
    return dict(model=task["model"], ablation=task["ablation"], seed=task["seed"], set=SET_NAME[task["split"]],
                n_parameters=int(n_par), **{k: s[k] for k in ("n", "nll", "mean_err_pct", "bias_pct", "std_err_pct",
                                                               "skew_mae", "w1_pct_median", "compact_mean_err_pct")},
                **extra, seconds=round(time.time() - t0))


def run(args):
    if args.frames:
        os.environ[data.FRAMES_ENV] = args.frames
    data.frames()
    out = pathlib.Path(args.out)
    done = pd.read_csv(out) if out.exists() else pd.DataFrame(columns=KEY)
    have = set(map(tuple, done[KEY].astype(str).to_numpy()))
    todo = [t for t in tasks() if (args.only is None or t["model"] == args.only.replace("-", "_"))
            and (t["model"], t["ablation"], str(t["seed"]), SET_NAME[t["split"]]) not in have]
    print(f"{len(todo)} runs to do, {len(have)} done", flush=True)
    # spawn, not fork: the parent has already started JAX, and a forked worker deadlocks on its locks
    with ProcessPoolExecutor(args.jobs, mp_context=multiprocessing.get_context("spawn")) as pool:
        for f in as_completed([pool.submit(_run, t) for t in todo]):
            row = f.result()
            done = pd.read_csv(out) if out.exists() else pd.DataFrame()     # re-read: another runner may have added rows
            pd.concat([done, pd.DataFrame([row])], ignore_index=True).to_csv(out, index=False)
            print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
