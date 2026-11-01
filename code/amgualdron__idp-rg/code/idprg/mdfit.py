"""Per-protein maximum-likelihood fits of the two-exponent distribution to a protein's stored MD frames
(data/md_fits.csv):

    free             ln A, ln S_m, u all free. Starts: S_m at the mode of a log-normal matched to the moments of ln Rg,
                     A from its curvature there, nu = 0.588 or nu = 0.5 (with S_m at the mean of ln Rg)
    law-constrained  ln A fixed to the law a . [1, nu, nu^2, 1/N, <lambda>, s+] (a = network.law_start on a fit set's
                     free fits); ln S_m and u fitted, started from the free fit

Each start: L-BFGS-B, then Nelder-Mead from its result; the best of all is kept. Box: ln A in [-5, 12], u in [-8, 8].
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import logit

from .data import compress
from .distribution import T, exponents, nu_of

FREE_BOUNDS = [(-5.0, 12.0), (None, None), (-8.0, 8.0)]
LAW_BOUNDS = [(None, None), (-8.0, 8.0)]


def logpdf(r, lnA, lnSm, nu):
    """ln p(Rg) at the frame values r, one protein (normaliser by the trapezoid rule on T)."""
    g, dl = exponents(nu)
    A = np.exp(lnA)
    f = lambda ls: np.exp(np.clip(dl * ls, None, 700)) / dl + np.exp(np.clip(-g * ls, None, 700)) / g
    lg = -A * f(T) + T
    m = lg.max()
    return -A * f(np.log(r) - lnSm) - (m + np.log(np.trapezoid(np.exp(lg - m), T)) + lnSm)


def _clip(x, bounds):
    return np.array([min(max(v, -np.inf if lo is None else lo), np.inf if hi is None else hi) for v, (lo, hi) in zip(x, bounds)])


def _minimize(f, starts, bounds):
    best = None
    for x0 in starts:
        x0 = _clip(np.asarray(x0, float), bounds)
        with np.errstate(all="ignore"):
            s = minimize(f, x0, method="L-BFGS-B", bounds=bounds)
            s2 = minimize(f, s.x if np.isfinite(s.fun) else x0, method="Nelder-Mead", bounds=bounds,
                          options=dict(xatol=1e-8, fatol=1e-12, maxiter=6000, maxfev=6000))
        for c in (s, s2):
            if np.isfinite(c.fun) and (best is None or c.fun < best.fun):
                best = c
    return best.x


def fit_free(frames):
    """(ln A, ln S_m, u, NLL per frame) of one protein's frames."""
    u, c = compress(frames)
    y, n = np.log(u), c.sum()
    m = np.sum(c * y) / n
    s = np.sqrt(np.sum(c * (y - m) ** 2) / n)
    g, dl = exponents(0.588)
    lnA0 = -np.log((g + dl) * s ** 2)
    starts = [[lnA0, m - s ** 2, logit((0.588 - 1 / 3) * 1.5)], [lnA0, m, logit((0.5 - 1 / 3) * 1.5)]]
    nll = lambda t: -np.sum(c * logpdf(u, t[0], t[1], nu_of(t[2]))) / n
    t = _minimize(nll, starts, FREE_BOUNDS)
    with np.errstate(all="ignore"):
        return (*t, nll(t))


def fit_law(frames, N, mean_hydro, scd_pos, law, start):
    """(ln S_m, u) of one protein's frames with ln A on the law; start = (ln S_m, u) of its free fit."""
    u, c = compress(frames)
    law = np.asarray(law, float)

    def nll(t):
        nu = nu_of(t[1])
        lnA = np.stack([np.ones_like(nu), nu, nu ** 2, 1 / N, mean_hydro, scd_pos], -1) @ law
        with np.errstate(all="ignore"):
            return -np.sum(c * logpdf(u, lnA, t[0], nu)) / c.sum()
    t = _minimize(nll, [start], LAW_BOUNDS)
    return t[0], t[1]
