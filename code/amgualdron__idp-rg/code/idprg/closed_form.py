"""The closed form: the network's two-exponent distribution with the network replaced by linear regressions on four
descriptors and ln N, fitted by maximum likelihood to the MD frames (paper/closed_form.tex).

    ln S_m = b_S + beta  . [s+, s-, <lambda>, h, ln N, s- ln N, <lambda> ln N, (ln N)^2]
    u      = b_u + gamma . [s+, s+/N, s-, <lambda>, h],          nu = 1/3 + (2/3) expit(u)
    ln A   = a . [1, nu, 1/N, <lambda>]
    h      = max(CE - c(N), 0),   c(N) = c0 + c1 ln N  (the CE threshold; CE = contact_energy)

Fit (`fit`): the mean over the fit proteins of the NLL per frame on all their stored frames, plus (lam / n) times the
squared slopes of ln S_m and u on standardized regressors (lam = 1; the ln A coefficients are not penalised).
L-BFGS-B in float64, restarted until it stops improving, from a two-step start: ridge regressions of the per-protein maximum-likelihood fits
(data/md_fits.csv) on the regressors, and least squares for the ln A law. c0 is the fit proteins' median of
CE - c1 ln N; c1 = 0.25 was chosen by 3-fold cross-validation on the training split (results/ce_threshold_cv.csv).
The options of `fit` (drop a descriptor, no ln N interactions, a free linear ln A) exist for the ablations.
"""
import json

import jax
import jax.numpy as jnp
import numpy as np
from scipy.optimize import minimize

from . import distribution as dist
from .data import ROOT, packed

DESC = ("scd_pos", "scd_neg", "mean_hydro", "h")
LAW = ("1", "nu", "1/N", "mean_hydro")
PRODUCTS = ("scd_neg", "mean_hydro")                   # the descriptors with an x ln N product in ln S_m
# L-BFGS-B's first trial step can send the objective to about 1e298, after which the line search reports false
# convergence at iteration 1. Same objective, parameters rescaled so that step is 100 times smaller.
SCALE = 0.01
NAME = {"S": "ln_Sm", "u": "u", "A": "ln_A"}           # regression -> key in closed_form.json
FREE = {"S": "lnSm", "u": "u", "A": "lnA"}            # regression -> per-protein fit in md_fits.csv


def load(path=ROOT / "model" / "closed_form.json"):
    return json.loads(open(path).read())


def _columns(desc, interactions):
    S = list(desc) + ["log_seq_length"] + ([f"{d}*log_seq_length" for d in desc if d in PRODUCTS] if interactions else []) + ["log_seq_length^2"]
    return {"S": S, "u": list(desc) + ["scd_pos/N"] * ("scd_pos" in desc), "A": list(desc) + ["log_seq_length"]}


def regressors(table, c0, c1):
    """The table with h, s+/N, the descriptor x ln N products and (ln N)^2 added."""
    d = table.copy()
    d["h"] = np.maximum(d.contact_energy - c1 * d.log_seq_length - c0, 0.0)
    for x in DESC:
        d[f"{x}*log_seq_length"] = d[x] * d.log_seq_length
    d["log_seq_length^2"] = d.log_seq_length ** 2
    d["scd_pos/N"] = d.scd_pos / d.N
    return d


def parameters(coef, table):
    """(ln A, ln S_m, u) for the rows of a descriptor table."""
    c = coef["ce_threshold"]
    d = regressors(table, c["c0"], c["c1"])
    lin = lambda r: coef[r]["intercept"] + sum(v * d[k].to_numpy(float) for k, v in coef[r]["coef"].items())
    lnSm, u = lin("ln_Sm"), lin("u")
    a = coef["ln_A"]
    if "intercept" in a:
        lnA = lin("ln_A")
    else:
        cov = {"1": 1.0, "nu": dist.nu_of(u), "1/N": 1 / d.N.to_numpy(float)}
        lnA = sum(v * (cov[k] if k in cov else d[k].to_numpy(float)) for k, v in a.items())
    return np.asarray(lnA, float), np.asarray(lnSm, float), np.asarray(u, float)


class _Form:
    """Parameter vector <-> regressions, on regressors standardized over the fit proteins."""

    def __init__(self, fit_df, desc, interactions, law):
        self.law = [k for k in LAW if k in ("1", "nu", "1/N") or k in desc] if law else None
        self.cols = _columns(desc, interactions)
        self.reg = ["S", "u"] + ([] if law else ["A"])
        self.mu = {r: fit_df[self.cols[r]].to_numpy(float).mean(0) for r in self.reg}
        self.sd = {r: fit_df[self.cols[r]].to_numpy(float).std(0) for r in self.reg}
        self.k = {r: len(self.cols[r]) for r in self.reg}

    def Z(self, df):
        return {r: jnp.asarray((df[self.cols[r]].to_numpy(float) - self.mu[r]) / self.sd[r]) for r in self.reg}

    def law_X(self, nu, df):
        return [nu, 1 / df.N.to_numpy(float)] + [df[k].to_numpy(float) for k in self.law[3:]]

    def unpack(self, th):
        out, i = {}, 0
        for r in self.reg:
            out[r] = (th[i], th[i + 1:i + 1 + self.k[r]])
            i += 1 + self.k[r]
        if self.law:
            out["law"] = th[i:i + len(self.law)]
        return out

    def penalty_mask(self):
        m = np.concatenate([[0.0] + [1.0] * self.k[r] for r in self.reg])
        return np.concatenate([m, np.zeros(len(self.law))]) if self.law else m

    def params(self, th, Z, law_raw):
        p = self.unpack(th)
        lnSm = p["S"][0] + Z["S"] @ p["S"][1]
        u = p["u"][0] + Z["u"] @ p["u"][1]
        if not self.law:
            return p["A"][0] + Z["A"] @ p["A"][1], lnSm, u
        nu = 1 / 3 + (2 / 3) * jax.nn.sigmoid(u)
        X = (jnp.column_stack([nu] + [law_raw[:, j] for j in range(law_raw.shape[1])]) - self.law_mu) / self.law_sd
        return p["law"][0] + X @ p["law"][1:], lnSm, u

    def start(self, fit_df, free):
        """Ridge (alpha 1) of the per-protein fits on the standardized regressors; ln A law by least squares."""
        th = []
        Z = self.Z(fit_df)
        for r in self.reg:
            t = free[FREE[r]].to_numpy(float)
            Zr = np.asarray(Z[r])
            th += [t.mean()] + list(np.linalg.solve(Zr.T @ Zr + np.eye(Zr.shape[1]), Zr.T @ (t - t.mean())))
        if self.law:
            X = np.column_stack(self.law_X(dist.nu_of(free.u.to_numpy(float)), fit_df))
            self.law_mu, self.law_sd = X.mean(0), X.std(0)
            c, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(X)), (X - self.law_mu) / self.law_sd]),
                                    free.lnA.to_numpy(float), rcond=None)
            th += list(c)
        return np.array(th, float)

    def to_raw(self, th):
        p, out = self.unpack(np.asarray(th)), {}
        for r in self.reg:
            b, s = p[r]
            c = s / self.sd[r]
            out[NAME[r]] = dict(intercept=float(b - np.sum(c * self.mu[r])), coef=dict(zip(self.cols[r], map(float, c))))
        if self.law:
            ck = p["law"][1:] / self.law_sd
            out["ln_A"] = dict(zip(self.law, map(float, [p["law"][0] - np.sum(ck * self.law_mu)] + list(ck))))
        return out


def fit(proteins, md_fits, fit, F, c1=0.25, drop=(), interactions=True, law=True, lam=1.0):
    """Fit on the rows `fit` of the protein table. Returns (coefficients in the closed_form.json layout, fit info)."""
    desc = [d for d in DESC if d not in drop]
    t = proteins.iloc[fit]
    c0 = float(np.median(t.contact_energy - c1 * t.log_seq_length))
    a = regressors(t, c0, c1).reset_index(drop=True)
    form = _Form(a, desc, interactions, law)
    free = md_fits.set_index("id").loc[a.id].reset_index()
    th0 = form.start(a, free)
    Z = form.Z(a)
    law_raw = jnp.asarray(np.column_stack(form.law_X(np.zeros(len(a)), a)[1:])) if law else None
    lr, c, seg, n = packed(F, a.frame_row)
    args = tuple(map(jnp.asarray, (lr, c, seg, form.penalty_mask()))) + (Z, law_raw)

    def objective(th, lr, c, seg, mask, Z, law_raw):         # data as arguments, not constants: XLA would fold them
        lp = dist.jax_logpdf_flat(lr, seg, *form.params(th, Z, law_raw))
        nll = -jax.ops.segment_sum(c * lp, seg, n) / jax.ops.segment_sum(c, seg, n)
        return jnp.mean(nll) + lam / n * jnp.sum(mask * th ** 2)

    vg = jax.jit(jax.value_and_grad(objective))
    f = lambda th: tuple(np.asarray(v, float) for v in vg(jnp.asarray(th), *args))
    th, nit, last = th0.copy(), 0, np.inf
    for _ in range(6):
        origin = th.copy()

        def scaled(z):
            v, g = f(origin + SCALE * z)
            return v, SCALE * g

        r = minimize(scaled, np.zeros_like(th), jac=True, method="L-BFGS-B",
                     options=dict(maxiter=3000, maxfun=6000, ftol=1e-12, gtol=1e-7 * SCALE))
        th, nit, stop = origin + SCALE * r.x, nit + r.nit, last - r.fun < 1e-8
        last = r.fun
        if stop:
            break
    if nit <= 20:
        raise RuntimeError(f"optimizer stopped after {nit} iterations: {r.message}")
    coef = dict(ce_threshold=dict(c0=c0, c1=c1), **form.to_raw(th))
    n_par = len(th) + ("h" in desc) * (1 + (c1 != 0))
    info = dict(n_fit=len(a), n_coefficients=len(th), n_parameters=int(n_par), objective_start=float(f(th0)[0]),
                objective_end=float(last), iterations=int(nit), message=str(r.message))
    return coef, info
