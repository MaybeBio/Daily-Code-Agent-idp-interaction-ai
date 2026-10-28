"""
fit_model.py

Least-squares model fitting with parameter bounds, random restarts,
and uncertainty estimation.

Usage:
    from viz_util.fit_model import fit

    def model(x, params):
        return params["a"] * x[0] + params["b"] * x[1]

    result = fit(
        model,
        x,                          # (n, d) array of inputs
        y,                          # (n,)   array of outputs
        params={"a": (0.0, 0.0, 1.0),   # name: (initial, lower, upper)
                "b": (0.5, 0.0, 1.0)},  #   initial=None → midpoint of bounds
        n_restarts=10,
    )

    print(result["params"])       # dict of fitted values
    print(result["stderr"])       # dict of standard errors
    print(result["loss"])         # scalar MSE at optimum
    print(result["residuals"])    # (n,) array
"""

import numpy as np
from scipy.optimize import curve_fit


def fit(model, x, y, params, n_restarts=20, rng=None):
    """
    Fit `model` to data (x, y) by minimising the mean-squared error.

    Parameters
    ----------
    model : callable
        Signature ``model(x, params) -> array (n,)`` where ``x`` is the full
        input array of shape ``(n,)`` or ``(n, d)`` and ``params`` is a dict
        of named floats.  The model must support vectorised evaluation.
    x : array-like, shape (n,) or (n, d)
        Input features.  A 1-D array is treated as (n, 1).
    y : array-like, shape (n,) or (n, 1)
        Scalar targets.
    params : dict[str, tuple]
        Each entry ``name: (p0, lo, hi)`` specifies an initial guess and
        bounds.  Set ``p0=None`` to use the midpoint ``(lo + hi) / 2``.
    n_restarts : int
        Number of random restarts.  The best solution is returned.
    rng : np.random.Generator or None
        Random number generator for reproducibility.  Defaults to
        ``np.random.default_rng(0)``.

    Returns
    -------
    dict with keys:
        params     : dict[str, float]  — fitted parameter values
        stderr     : dict[str, float]  — 1-sigma uncertainty (from covariance)
        loss       : float             — MSE at the optimum
        residuals  : ndarray (n,)      — y - f(x)
    """
    rng = rng or np.random.default_rng(0)
    x = np.atleast_2d(np.asarray(x, dtype=float))
    if x.shape[0] == 1 and x.shape[1] > 1:   # (1, n) → (n, 1)
        x = x.T
    y = np.asarray(y, dtype=float).ravel()
    n = len(y)

    names  = list(params)
    p0_ref = np.array([(v[0] if v[0] is not None else (v[1] + v[2]) / 2.0)
                       for v in params.values()])
    lowers = np.array([v[1] for v in params.values()])
    uppers = np.array([v[2] for v in params.values()])

    def _f(_dummy, *p):
        return np.asarray(model(x, dict(zip(names, p)))).ravel()

    best_popt = best_pcov = None
    best_loss = np.inf

    starts = [p0_ref] + [rng.uniform(lowers, uppers) for _ in range(n_restarts - 1)]
    for p0 in starts:
        try:
            popt, pcov = curve_fit(_f, None, y, p0=p0, bounds=(lowers, uppers),
                                   full_output=False, check_finite=True)
        except (RuntimeError, ValueError):
            continue
        resid = y - _f(None, *popt)
        loss  = float(np.mean(resid ** 2))
        if loss < best_loss:
            best_loss, best_popt, best_pcov = loss, popt, pcov

    if best_popt is None:
        raise RuntimeError("curve_fit did not converge in any of the "
                           f"{n_restarts} restarts.")

    resid  = y - _f(None, *best_popt)
    stderr = np.sqrt(np.diag(best_pcov))

    return {
        "params":    dict(zip(names, best_popt)),
        "stderr":    dict(zip(names, stderr)),
        "loss":      float(np.mean(resid ** 2)),
        "residuals": resid,
    }
