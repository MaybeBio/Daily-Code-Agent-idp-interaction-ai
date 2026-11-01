"""Closed-form variants: u with terms x (N/100)^-Delta, Delta fixed or fitted jointly by maximum likelihood; ln S_m's
x ln N products removable one at a time. See README.md here. model/ and results/ are not touched.

    python experiments/u_terms/fit.py JOB [JOB ...]        (IDPRG_FRAMES = the frame store)
    JOB = model:split:fold[:product,product]               e.g. dS:standard:0   dH:standard:1:scd_pos,h   dS:length:-1

split: standard (fold 0, 1, 2) | length | sticky (fold -1). A job with standard:-1 scores the validation split and needs
--validation. There is no way to reach the test split: it is sealed. Finished jobs are skipped.
--smoke: a few proteins and iterations, to check the plumbing; writes to runs/u_terms_smoke, not runs/u_terms.
"""
import argparse
import json
import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from idprg import closed_form as cf, data, distribution as dist, metrics  # noqa: E402

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

#         terms (descriptor, Delta at the start, fitted?)      plain terms dropped from u
MODELS = {"check": ([("scd_pos", 1.0, False)], []),                                   # s+/N, must reproduce the reference
          "dS": ([("scd_pos", 1.0, True)], []),
          "dH": ([("scd_pos", 1.0, False), ("h", 0.0, True)], ["h"]),
          "hN": ([("scd_pos", 1.0, False), ("h", 1.0, False)], []),
          "dSH": ([("scd_pos", 1.0, True), ("h", 0.0, True)], ["h"])}
# Saturating s+ terms in place of s+/N: (descriptor, shape, parameters at the start). Parameters are logs, of N0/100 and of s0.
# N0 -> 0 and s0 -> infinity both recover s+/N, so each model nests its parent.
SAT = {"satN": ("scd_pos", lambda x, r, d: x / (1 + jnp.exp(r - d[0])), [np.log(0.2)]),                                   # s+ / (1 + N/N0)
       "satS": ("scd_pos", lambda x, r, d: jnp.exp(d[0]) * jnp.tanh(x * jnp.exp(-d[0])) * jnp.exp(-r), [0.0]),            # s0 tanh(s+/s0) / N
       "satSN": ("scd_pos", lambda x, r, d: jnp.exp(d[0]) * jnp.tanh(x * jnp.exp(-d[0])) / (1 + jnp.exp(r - d[1])), [0.0, np.log(0.2)])}
for _m, (_k, _, _d) in SAT.items():
    MODELS[_m] = ([(_k, _d[0], True)], [])
# Variants of s+/N: u terms or ln A law terms removed, another c1, or a perturbed start. "-" in a c1 name would break tags.
VARIANT = {"noUlnN": dict(drop_u=["log_seq_length"]), "noUs": dict(drop_u=["scd_pos"]), "noAs": dict(drop_law=["scd_pos"]),
           **{f"c1_{t}": dict(c1=v) for t, v in (("m50", -0.5), ("0", 0.0), ("50", 0.5), ("75", 0.75), ("100", 1.0))},
           **{f"start{i}": dict(jitter=i) for i in (1, 2, 3)},
           "a1": dict(drop_law=["scd_pos"], drop_u=["log_seq_length"]),                       # stepwise, each on the one before
           "a2": dict(drop_law=["scd_pos"], drop_u=["log_seq_length", "scd_pos"]),
           "a3": dict(drop_law=["scd_pos"], drop_u=["log_seq_length", "scd_pos"], drop_s=["log_seq_length^2"])}
for _m, _v in VARIANT.items():
    MODELS[_m] = ([("scd_pos", 1.0, False)], _v.get("drop_u", []))
PRODUCTS = ("scd_pos", "scd_neg", "mean_hydro", "h")                                  # removable x ln N products of ln S_m
# L-BFGS-B's first trial step can send the objective to about 1e298, after which the line search reports false
# convergence at iteration 1. Same objective, parameters rescaled so that step is 100 times smaller; restarted until
# it stops improving.
SCALE = 0.01
# The form these tests started from, the published 24-parameter one; idprg.closed_form has since adopted model a1.
cf.LAW = ("1", "nu", "1/N", "mean_hydro", "scd_pos")


def _columns(desc, interactions):
    S = list(desc) + ["log_seq_length"] + ([f"{d}*log_seq_length" for d in desc] if interactions else []) + ["log_seq_length^2"]
    return {"S": S, "u": list(desc) + ["log_seq_length"], "A": list(desc) + ["log_seq_length"]}


def run(name, split, fold, drop_S, smoke, allow_validation):
    terms, drop_u = MODELS[name]
    gone = [f"{d}*log_seq_length" for d in drop_S] + VARIANT.get(name, {}).get("drop_s", [])
    cf._columns = lambda desc, inter: {**_columns(desc, inter), "S": [c for c in _columns(desc, inter)["S"] if c not in gone],
                                       "u": [c for c in _columns(desc, inter)["u"] if c not in drop_u]}
    P, M, F = data.proteins(), data.md_fits(), data.frames()
    if split not in ("standard", "length", "sticky"):
        raise ValueError(f"split {split!r}: the test split is sealed")
    if split == "standard" and fold < 0 and not allow_validation:
        raise ValueError("standard:-1 scores the validation split: pass --validation for the one final check")
    fit, ev = data.split(P, split, fold)
    if smoke:
        fit, ev = fit[:40], ev[:12]
    assert (P.split.iloc[fit] == "train").all() and not (P.split.iloc[ev] == "test").any()
    t = P.iloc[fit]
    variant = VARIANT.get(name, {})
    c1 = variant.get("c1", 0.25)
    c0 = float(np.median(t.contact_energy - c1 * t.log_seq_length))
    a = cf.regressors(t, c0, c1).reset_index(drop=True)
    form = cf._Form(a, list(cf.DESC), True, True)
    form.law = [k for k in form.law if k not in variant.get("drop_law", [])]
    free = M.set_index("id").loc[a.id].reset_index()
    th_lin = form.start(a, free)
    n_lin = len(th_lin)
    keys = [x for x, _, _ in terms]
    d0 = np.array(SAT[name][2] if name in SAT else [d for _, d, _ in terms])
    fitted = list(range(len(d0))) if name in SAT else [i for i, (_, _, f) in enumerate(terms) if f]
    term = (lambda x, r, d: SAT[name][1](x[:, 0], r, d)[:, None]) if name in SAT else (lambda x, r, d: x * jnp.exp(-r[:, None] * d))
    ln_n = lambda df: np.log(df.N.to_numpy(float) / 100)
    X = lambda df: np.column_stack([df[k].to_numpy(float) for k in keys])
    T0 = np.asarray(term(X(a), ln_n(a), d0))
    sd = T0.std(0)
    p = form.unpack(th_lin)
    Z0 = form.Z(a)
    res = free.u.to_numpy(float) - np.asarray(p["u"][0] + Z0["u"] @ p["u"][1])
    g0 = np.linalg.lstsq(np.column_stack([np.ones(len(a)), T0 / sd]), res, rcond=None)[0][1:]
    th0 = np.concatenate([th_lin, g0, d0[fitted]])
    if "jitter" in variant:
        th0 = th0 + 0.1 * np.random.default_rng(variant["jitter"]).standard_normal(len(th0))
    fit_idx = jnp.asarray(fitted, int)

    def params(th, Z, law_raw, Xj, rj):
        q = form.unpack(th[:n_lin])
        g = th[n_lin:n_lin + len(keys)]
        d = jnp.asarray(d0).at[fit_idx].set(th[n_lin + len(keys):])
        lnSm = q["S"][0] + Z["S"] @ q["S"][1]
        u = q["u"][0] + Z["u"] @ q["u"][1] + jnp.sum(g * term(Xj, rj, d) / sd, axis=1)
        nu = 1 / 3 + (2 / 3) * jax.nn.sigmoid(u)
        L = (jnp.column_stack([nu] + [law_raw[:, j] for j in range(law_raw.shape[1])]) - form.law_mu) / form.law_sd
        return q["law"][0] + L @ q["law"][1:], lnSm, u

    law_of = lambda df: jnp.asarray(np.column_stack(form.law_X(np.zeros(len(df)), df)[1:]))
    mask = jnp.asarray(np.concatenate([form.penalty_mask(), np.ones(len(keys)), np.zeros(len(fitted))]))
    lr, c, seg, n = data.packed(F, a.frame_row)
    args = tuple(map(jnp.asarray, (lr, c, seg))) + (mask, Z0, law_of(a), jnp.asarray(X(a)), jnp.asarray(ln_n(a)))

    def objective(th, lr, c, seg, mask, Z, law_raw, Xj, rj):
        lp = dist.jax_logpdf_flat(lr, seg, *params(th, Z, law_raw, Xj, rj))
        return jnp.mean(-jax.ops.segment_sum(c * lp, seg, n) / jax.ops.segment_sum(c, seg, n)) + 1.0 / n * jnp.sum(mask * th ** 2)

    vg = jax.jit(jax.value_and_grad(objective))
    f = lambda th: tuple(np.asarray(v, float) for v in vg(jnp.asarray(th), *args))
    th, nit, last = th0.copy(), 0, np.inf
    options = dict(maxiter=25 if smoke else 3000, maxfun=60 if smoke else 6000, ftol=1e-12, gtol=1e-7 * SCALE)
    for _ in range(1 if smoke else 6):
        origin = th.copy()

        def scaled(z):
            v, g = f(origin + SCALE * z)
            return v, SCALE * g

        r = minimize(scaled, np.zeros_like(th), jac=True, method="L-BFGS-B", options=options)
        th, nit, stop = origin + SCALE * r.x, nit + r.nit, last - r.fun < 1e-8
        last = r.fun
        if stop:
            break
    if nit <= 20:
        raise RuntimeError(f"optimizer stopped after {nit} iterations: {r.message}")
    ev_df = cf.regressors(P.iloc[ev], c0, c1)
    pr = params(jnp.asarray(th), form.Z(ev_df), law_of(ev_df), jnp.asarray(X(ev_df)), jnp.asarray(ln_n(ev_df)))
    d = metrics.score(*[np.asarray(v, float) for v in pr], P.iloc[ev], F)
    deltas = d0.copy()
    deltas[fitted] = th[n_lin + len(keys):]
    info = dict(model=name, shape_params=[float(v) for v in deltas], drop_S=drop_S, split=split, fold=fold, n_fit=len(a), n_parameters=len(th) + 2, objective_end=float(last),
                iterations=int(nit), message=str(r.message), ce_threshold=dict(c0=c0, c1=c1),
                terms={k: dict(delta=float(dl), coef=float(gi / s)) for k, dl, gi, s in zip(keys, deltas, th[n_lin:], sd)},
                **form.to_raw(th[:n_lin]))
    return d, info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobs", nargs="+")
    ap.add_argument("--validation", action="store_true", help="allow standard:-1, the one final check on validation")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    assert jnp.asarray(1.0).dtype == jnp.float64, "idprg should have switched JAX to float64"
    out = ROOT / "runs" / ("u_terms_smoke" if args.smoke else "u_terms")
    out.mkdir(parents=True, exist_ok=True)
    for job in args.jobs:
        parts = job.split(":")
        name, split, fold = parts[0], parts[1], int(parts[2])
        drop_S = parts[3].split(",") if len(parts) > 3 else []
        if set(drop_S) - set(PRODUCTS):
            raise ValueError(f"{job}: products are named by descriptor, from {PRODUCTS}")
        tag = f"{name}{'_no_' + '+'.join(drop_S) if drop_S else ''}.{split}.{fold}"
        if (out / f"{tag}.csv").exists():
            continue
        t = time.time()
        d, info = run(name, split, fold, drop_S, args.smoke, args.validation)
        d.to_csv(out / f"{tag}.csv", index=False)
        (out / f"{tag}.json").write_text(json.dumps(info, indent=1))
        print(f"{tag} parameters={info['n_parameters']} objective={info['objective_end']:.6f} iterations={info['iterations']} "
              f"nll={d.nll.mean():.5f} " + " ".join(f"{k}: delta={v['delta']:.3f} coef={v['coef']:.3f}" for k, v in info["terms"].items())
              + f" {time.time() - t:.0f}s", flush=True)


if __name__ == "__main__":
    main()
