"""Scattering-curve statistics: scale-and-fit, chi-squared, Guinier analysis and maximum-entropy
reweighting of a conformer ensemble against an experimental profile.

Conventions follow PeptoneBench (Invernizzi et al. 2025) so that numbers are comparable with that
benchmark: intensities are standardised by the experimental error, and a computed curve is scaled to
the data with the single factor of Svergun et al. (1995). A constant background is fitted as well,
because Pepsi-SAXS curves and SASBDB data differ in how the buffer was subtracted.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax


def scale_and_background(model: np.ndarray, exp: np.ndarray, sigma: np.ndarray) -> tuple[float, float]:
    """Weighted least-squares scale c and background b minimising sum(((c*model + b) - exp)/sigma)^2."""
    w = 1.0 / sigma**2
    A = np.vstack([model, np.ones_like(model)]).T * np.sqrt(w)[:, None]
    y = exp * np.sqrt(w)
    (c, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(c), float(b)


def chi2(model: np.ndarray, exp: np.ndarray, sigma: np.ndarray, fit_background: bool = True) -> float:
    """Reduced chi-squared of a computed curve against the data after scaling (and background)."""
    if fit_background:
        c, b = scale_and_background(model, exp, sigma)
    else:
        c, b = float(np.sum(exp * model / sigma**2) / np.sum(model**2 / sigma**2)), 0.0
    r = (c * model + b - exp) / sigma
    dof = len(exp) - (2 if fit_background else 1)
    return float(np.sum(r**2) / dof)


def nrmsd_log(model: np.ndarray, exp: np.ndarray, sigma: np.ndarray, snr_min: float = 3.0) -> float:
    """A second fit score that does not weight by the reported errors.

    The computed curve is scaled to the data with the same scale factor and constant background as
    the chi-squared fit. Over the usable q range, the points at which the experimental intensity is
    positive and at least `snr_min` standard errors above zero and the scaled model is positive, the
    root-mean-square deviation of ln I between model and data is divided by the range of ln I of the
    data over those points. 0 is a perfect fit; 0.02 is a deviation of 2% of the curve's dynamic
    range in log intensity. Errors enter only through the point selection and the scale fit."""
    c, b = scale_and_background(model, exp, sigma)
    fit = c * model + b
    use = (exp > 0) & (exp >= snr_min * sigma) & (fit > 0)
    if use.sum() < 10:
        return float("nan")
    d = np.log(fit[use]) - np.log(exp[use])
    rng = np.ptp(np.log(exp[use]))
    return float(np.sqrt(np.mean(d**2)) / rng) if rng > 0 else float("nan")


def ensemble_curve(curves: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    """Weighted mean of per-conformer curves, shape (n_conformers, n_q) -> (n_q,)."""
    if weights is None:
        return curves.mean(axis=0)
    return weights @ curves


def guinier(q: np.ndarray, I: np.ndarray, sigma: np.ndarray, qrg_max: float = 1.3) -> dict:
    """Guinier fit ln I = ln I0 - q^2 Rg^2 / 3 over the lowest q with q*Rg <= qrg_max (iterated).

    The fit uses the upper three quarters of the window, so that an aggregation upturn at the
    lowest angles shows up as a positive excess (`upturn`) of those points above the line and does
    not inflate Rg. Returns Rg, I0, Rg's standard error, the number of points and the validity."""
    pos = I > 0
    q, I, sigma = q[pos], I[pos], sigma[pos]
    n, rg = max(12, len(q) // 10), np.nan
    for _ in range(10):
        lo = n // 4
        qq, ll, ww = q[lo:n] ** 2, np.log(I[lo:n]), (I[lo:n] / sigma[lo:n]) ** 2
        A = np.vstack([np.ones(n - lo), -qq / 3.0]).T * np.sqrt(ww)[:, None]
        (lnI0, rg2), *_ = np.linalg.lstsq(A, ll * np.sqrt(ww), rcond=None)
        if rg2 <= 0:
            return {"rg": np.nan, "rg_err": np.nan, "i0": float(np.exp(lnI0)), "n_points": n,
                    "upturn": np.nan, "valid": False}
        rg_new = float(np.sqrt(rg2))
        n_new = max(12, min(int(np.searchsorted(q * rg_new, qrg_max)), len(q)))
        converged = n_new == n
        rg, n = rg_new, n_new
        if converged:
            break
    cov = np.linalg.inv(A.T @ A)
    rg_err = float(0.5 * np.sqrt(cov[1, 1]) / rg)
    line = np.exp(lnI0 - q[:n] ** 2 * rg**2 / 3.0)
    upturn = float(np.mean(I[:lo] / line[:lo]) - 1.0)
    return {"rg": rg, "rg_err": rg_err, "i0": float(np.exp(lnI0)), "n_points": n, "upturn": upturn,
            "valid": bool(q[n - 1] * rg <= qrg_max + 0.05)}


def conformer_rg(coords: np.ndarray) -> np.ndarray:
    """Radius of gyration of each conformer from heavy-atom coordinates, shape (n, atoms, 3)."""
    c = coords - coords.mean(axis=1, keepdims=True)
    return np.sqrt((c**2).sum(axis=2).mean(axis=1))


# ---- maximum-entropy reweighting ------------------------------------------------------------------
def _standardise(curves: np.ndarray, exp: np.ndarray, sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return curves / sigma, exp / sigma


def reweight_objective(curves: np.ndarray, exp: np.ndarray, sigma: np.ndarray, theta: float):
    """The objective over log-weights, chi2/2 + theta * KL(w || uniform), with the computed curve scaled
    to the data (Svergun factor) and a constant background fitted by least squares inside it; returns a
    function of the log-weights giving (value, gradient)."""
    X, y = _standardise(curves, exp, sigma)
    ones = 1.0 / sigma           # the constant background, standardised like the data

    def f(logw):
        lw = logw - logsumexp(logw)
        w = np.exp(lw)
        avg = w @ X
        # scale and background by least squares on the standardised curve
        A = np.vstack([avg, ones]).T
        (c, b), *_ = np.linalg.lstsq(A, y, rcond=None)
        r = c * avg + b * ones - y
        loss = 0.5 * r @ r + theta * np.sum(w * lw)
        g_w = X @ (c * r) + theta * (1.0 + lw)
        jac = w * g_w - w * (w @ g_w)
        return loss, jac

    return f


def reweight(curves: np.ndarray, exp: np.ndarray, sigma: np.ndarray, theta: float,
             x0: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    """Weights minimising chi2/2 + theta * KL(w || uniform) (see reweight_objective), optimised over
    log-weights; returns (weights, reduced chi2 of the reweighted curve)."""
    n = curves.shape[0]
    f = reweight_objective(curves, exp, sigma, theta)
    if x0 is None:
        x0 = np.full(n, -np.log(n))
    res = minimize(f, x0, jac=True, method="L-BFGS-B", options={"maxiter": 5000, "gtol": 1e-6})
    w = softmax(res.x)
    return w, chi2(ensemble_curve(curves, w), exp, sigma)


def kish_fraction(w: np.ndarray) -> float:
    """Effective sample size as a fraction of the number of conformers."""
    return float(np.sum(w) ** 2 / np.sum(w**2) / len(w))


def reweighting_curve(curves: np.ndarray, exp: np.ndarray, sigma: np.ndarray, thetas: np.ndarray,
                      return_weights: bool = False) -> list[dict] | tuple[list[dict], list]:
    """Chi2 and Kish fraction along a decreasing sequence of theta, from almost no reweighting
    towards a full fit. Each solve starts from uniform weights. With `return_weights` the weight
    vector of every point on the path is returned as well."""
    out, ws = [], []
    for th in sorted(thetas, reverse=True):
        w, c2 = reweight(curves, exp, sigma, th)
        out.append({"theta": float(th), "chi2": c2, "phi": kish_fraction(w)})
        ws.append(w)
    return (out, ws) if return_weights else out


def phi_at_chi2(path: list[dict], target: float) -> float:
    """Largest Kish fraction along the path at which chi2 <= target; nan if never reached."""
    ok = [p["phi"] for p in path if p["chi2"] <= target]
    return max(ok) if ok else float("nan")


def chi2_at_phi(path: list[dict], phi_min: float) -> float:
    """Best chi2 reachable while keeping at least phi_min of the effective sample."""
    ok = [p["chi2"] for p in path if p["phi"] >= phi_min]
    return min(ok) if ok else float("nan")
