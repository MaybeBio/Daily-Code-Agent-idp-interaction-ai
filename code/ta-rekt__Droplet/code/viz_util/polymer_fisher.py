"""
Smoothed polymer Fisher information rate from dwell-time distributions,
with run-level bootstrap error bars.

The inline notebook version estimates the score d/dphi log P(t|phi) by a 3-point
central difference between neighbouring concentrations.  Here the score is the
analytic derivative of a smoothing spline fitted along phi to the log-pdf
surface, which is less noisy and is defined at every concentration including the
two endpoints (the finite-difference version leaves them NaN).

Two things make this non-trivial and are handled explicitly:

1. I(phi) = INT P (dlogP/dphi)^2 dt is quadratic in the score, so score noise
   biases it *upward*.  A sharper derivative therefore does not mean a better
   I.  In particular, choosing the spline bandwidth by cross-validating the
   log-pdf itself is wrong -- that criterion optimises the function value, drives
   lambda to 0, and lands on an estimator worse than the finite difference.
   `select_lambda` instead minimises an estimated MSE of I: bootstrap variance
   plus squared deviation from a noisy-but-nearly-unbiased small-lambda pilot.

2. The 16 runs are the independent replicates, so the bootstrap resamples runs,
   not durations.  The same run selection is applied to the closed and open
   durations of a concentration, since both come from the same trajectory --
   otherwise the band on the summed rate is wrong.

Pre-binning each run's durations onto the shared log-time grid once means a
bootstrap replicate is a few matrix operations, so 400 replicates take seconds.

Self-test (lognormal dwell times with an analytically known Fisher information):
    python3 viz_util/polymer_fisher.py
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d


# --------------------------------------------------------------- time grids
def build_counts(times, n_grid=500):
    """Pre-bin every run onto one shared log-time grid.

    times : [conc][run] -> sequence of durations (the notebook's bound_times /
            unbound_times).
    Returns dict(grid (n_grid,), counts (K, R, n_grid) cumulative counts per run,
    nsamp (K, R)).  A pooled empirical CDF over any subset of runs is then
    counts[subset].sum(0)/nsamp[subset].sum(), which is what makes the bootstrap
    cheap -- the raw durations are never re-sorted.
    """
    K = len(times)
    R = max(len(times[k]) for k in range(K))
    runs = [[np.sort(np.asarray(times[k][r], float)) if r < len(times[k])
             else np.empty(0) for r in range(R)] for k in range(K)]
    runs = [[a[a > 0] for a in row] for row in runs]
    allt = np.concatenate([a for row in runs for a in row if a.size])
    grid = np.logspace(np.log10(allt.min()), np.log10(allt.max()), n_grid)
    counts = np.zeros((K, R, n_grid))
    nsamp = np.zeros((K, R))
    for k in range(K):
        for r in range(R):
            a = runs[k][r]
            nsamp[k, r] = a.size
            if a.size:
                counts[k, r] = np.searchsorted(a, grid, side="right")
    return dict(grid=grid, counts=counts, nsamp=nsamp)


def _pool(C, sel):
    counts, nsamp = C["counts"], C["nsamp"]
    if sel is None:
        return counts.sum(1), nsamp.sum(1)
    rows = np.arange(counts.shape[0])[:, None]
    return counts[rows, sel].sum(1), nsamp[rows, sel].sum(1)


def pdfs_from_counts(C, sel=None, sigma=10.0):
    """Pooled smoothed pdf f(t) per concentration -> (K, n_grid).

    Same recipe as the notebook cell: empirical CDF on the log grid, gaussian
    blur in log10 t, differentiate.  sel (K, R) selects runs; None pools all.
    """
    grid = C["grid"]
    tot, n = _pool(C, sel)
    F = tot / np.maximum(n, 1)[:, None]
    F = gaussian_filter1d(F, sigma=sigma, axis=1, mode="nearest")
    pdf = np.gradient(F, np.log10(grid), axis=1) / (grid * np.log(10.0))
    return np.clip(pdf, 0.0, None)


def mean_dwell(C, sel=None):
    """Pooled mean duration per concentration, E[T] = INT (1-F) dt."""
    grid = C["grid"]
    tot, n = _pool(C, sel)
    F = tot / np.maximum(n, 1)[:, None]
    return np.trapezoid(1.0 - F, grid, axis=1) + grid[0]


# ------------------------------------------------- smoothing spline along phi
def spline_ops(x, lam, n_sub=6):
    """Discrete cubic smoothing spline: minimise
    sum_k (u(x_k)-d_k)^2/K + lam * INT (u'')^2 dx  on a refined grid.
    Returns Sd (K, K), which maps K data values to du/dx at those same K points.
    """
    x = np.asarray(x, float)
    K = x.size
    xs = x[-1] - x[0]
    n = (K - 1) * n_sub + 1
    h = 1.0 / (n - 1)
    idx = np.arange(K) * n_sub
    E2 = (np.eye(n, k=1) * -2 + np.eye(n, k=2) + np.eye(n))[:-2]
    S = np.zeros((K, n))
    S[np.arange(K), idx] = 1.0
    A = S.T @ S / K + lam * (E2.T @ E2) / h**3
    A += 1e-12 * np.trace(A) / n * np.eye(n)
    H = np.linalg.solve(A, S.T / K)
    dx = xs * h
    G = np.zeros((n, n))                                  # np.gradient as a matrix
    for j in range(n):
        if j == 0:
            G[j, 0:3] = np.array([-3.0, 4.0, -1.0]) / (2 * dx)
        elif j == n - 1:
            G[j, -3:] = np.array([1.0, -4.0, 3.0]) / (2 * dx)
        else:
            G[j, j - 1], G[j, j + 1] = -1 / (2 * dx), 1 / (2 * dx)
    return (G @ H)[idx]


# ------------------------------------------------------- Fisher information
def _surface(P, floor_rel):
    fl = P.max(axis=1, keepdims=True) * floor_rel
    return np.log(np.maximum(P, fl)), np.all(P > fl, axis=0)


def fisher_from_pdfs(P, grid, Sd, floor_rel=1e-10, ok=None):
    """I(phi) = INT P(t|phi) (dlogP/dphi)^2 dt, one value per concentration.

    Times where any concentration's pdf sits on the numerical floor are dropped
    from the integral, as in the finite-difference version.  Pass `ok` to freeze
    that support mask across bootstrap replicates.
    """
    L, ok_ = _surface(P, floor_rel)
    if ok is None:
        ok = ok_
    return np.trapezoid(np.where(ok, P * (Sd @ L) ** 2, 0.0), grid, axis=1)


def cross_fisher(C, phi, lam, nsplit=24, seed=0, sigma=10.0, n_sub=6,
                 floor_rel=1e-10, ok=None):
    """Bias diagnostic: score taken from two disjoint halves of the runs.

    I = INT P score_A score_B dt.  The two scores carry independent noise, so
    the noise-squared term that inflates the ordinary estimator averages to
    zero.  On a channel with no real phi dependence this scatters around 0 and
    can go negative -- that negative excursion is the honest noise floor.  If a
    bump survives here, it is not a score-noise artefact.
    """
    rng = np.random.default_rng(seed)
    K, R, _ = C["counts"].shape
    P = pdfs_from_counts(C, sigma=sigma)
    L, ok_ = _surface(P, floor_rel)
    if ok is None:
        ok = ok_
    Sd = spline_ops(phi, lam, n_sub)
    acc = np.zeros(K)
    for _ in range(nsplit):
        p = rng.permutation(R)
        A, B = p[:R // 2], p[R // 2:]
        LA, _ = _surface(pdfs_from_counts(C, sel=np.tile(A, (K, 1)), sigma=sigma), floor_rel)
        LB, _ = _surface(pdfs_from_counts(C, sel=np.tile(B, (K, 1)), sigma=sigma), floor_rel)
        acc += np.trapezoid(np.where(ok, P * (Sd @ LA) * (Sd @ LB), 0.0), C["grid"], axis=1)
    return acc / nsplit


def select_lambda(C, phi, lams, nboot=40, seed=0, lam_pilot=None, sigma=10.0,
                  n_sub=6, floor_rel=1e-10, core=None):
    """Bandwidth by estimated MSE of I(phi): bootstrap variance + squared
    deviation from a small-lambda pilot (noisy, but nearly unbiased).

    Cross-validating the log-pdf instead would drive lambda to 0 and give an
    I(phi) worse than the finite-difference estimator -- see the module docstring.
    Returns (lam, mse curve).  `core` restricts the criterion to a slice of
    concentrations (default: drop the two endpoints, where the one-sided spline
    derivative is least reliable).
    """
    rng = np.random.default_rng(seed)
    grid = C["grid"]
    K, R, _ = C["counts"].shape
    core = slice(1, K - 1) if core is None else core
    P = pdfs_from_counts(C, sigma=sigma)
    L, ok = _surface(P, floor_rel)
    lam_pilot = lams[0] if lam_pilot is None else lam_pilot
    Sds = {lam: spline_ops(phi, lam, n_sub) for lam in lams}

    acc = {lam: [] for lam in lams}
    for _ in range(nboot):
        sel = rng.integers(0, R, size=(K, R))
        Pb = pdfs_from_counts(C, sel=sel, sigma=sigma)
        Lb, _ = _surface(Pb, floor_rel)
        for lam in lams:
            acc[lam].append(np.trapezoid(np.where(ok, Pb * (Sds[lam] @ Lb) ** 2, 0.0),
                                         grid, axis=1))
    pilot = fisher_from_pdfs(P, grid, spline_ops(phi, lam_pilot, n_sub), floor_rel, ok)
    mse = np.array([np.mean((fisher_from_pdfs(P, grid, Sds[lam], floor_rel, ok) - pilot)[core]**2
                            + np.var(np.array(acc[lam]), axis=0)[core]) for lam in lams])
    return lams[int(np.argmin(mse))], mse


# --------------------------------------------------------------- the estimator
def rate(Cb, Cu, phi, lam_b, lam_u, sel=None, sigma=10.0, n_sub=6,
         floor_rel=1e-10, norm=None, ok_b=None, ok_u=None):
    """Closed / open Fisher rate phi^2 I(phi)/norm for one (re)sample.

    norm=None reproduces the notebook's scalar normaliser: the pooled mean closed
    dwell time at the middle concentration.  Pass norm='tau' for the phi-resolved
    correlation time 1/(1/tau_closed + 1/tau_open), which changes the shape of
    the curve, not just its scale.
    """
    Pb = pdfs_from_counts(Cb, sel=sel, sigma=sigma)
    Pu = pdfs_from_counts(Cu, sel=sel, sigma=sigma)
    Ib = fisher_from_pdfs(Pb, Cb["grid"], spline_ops(phi, lam_b, n_sub), floor_rel, ok_b)
    Iu = fisher_from_pdfs(Pu, Cu["grid"], spline_ops(phi, lam_u, n_sub), floor_rel, ok_u)
    if norm is None:
        nrm = mean_dwell(Cb, sel)[len(phi) // 2]
    elif isinstance(norm, str) and norm == "tau":
        tb, tu = mean_dwell(Cb, sel), mean_dwell(Cu, sel)
        nrm = 1.0 / (1.0 / tb + 1.0 / tu)
    else:
        nrm = norm
    return phi**2 * Ib / nrm, phi**2 * Iu / nrm, nrm


def bootstrap_rate(Cb, Cu, phi, lam_b, lam_u, nboot=400, seed=0, **kw):
    """Resample runs with replacement; one selection shared by closed and open."""
    rng = np.random.default_rng(seed)
    K, R, _ = Cb["counts"].shape
    Fb = np.empty((nboot, K))
    Fu = np.empty((nboot, K))
    for b in range(nboot):
        sel = rng.integers(0, R, size=(K, R))
        Fb[b], Fu[b], _ = rate(Cb, Cu, phi, lam_b, lam_u, sel=sel, **kw)
    return Fb, Fu


def finite_difference_I(C, phi, sigma=10.0, floor_rel=1e-10):
    """The notebook's 3-point estimator, for comparison (endpoints are NaN)."""
    P = pdfs_from_counts(C, sigma=sigma)
    L, ok = _surface(P, floor_rel)
    dphi = phi[1] - phi[0]
    I = np.full(len(phi), np.nan)
    for k in range(1, len(phi) - 1):
        s = (L[k + 1] - L[k - 1]) / (2 * dphi)
        I[k] = np.trapezoid(np.where(ok, P[k] * s**2, 0.0), C["grid"])
    return I


# ------------------------------------------------------------------ self-test
if __name__ == "__main__":
    import time
    # Lognormal dwell times, log-mean sweeping with phi at fixed log-sd: the
    # Fisher information about phi carried by one observation is exactly
    # (mu'(phi)/s)^2, so the estimator has a ground truth to hit.
    phi = np.linspace(0.003, 0.018, 16)
    mu = 4.0 + 2.0 / (1.0 + np.exp(-(phi - 0.0095) / 0.002))
    s = 1.2
    I_true = np.gradient(mu, phi) ** 2 / s**2
    lams = np.logspace(-7, -1, 13)
    core = slice(3, 13)

    for nper in (100, 400):
        es, ef, ls = [], [], []
        for seed in range(5):
            rng = np.random.default_rng(seed)
            tb = [[np.maximum(rng.lognormal(mu[k], s, nper).round(), 1.0)
                   for _ in range(16)] for k in range(16)]
            C = build_counts(tb)
            lam, mse = select_lambda(C, phi, lams, seed=seed)
            I_sm = fisher_from_pdfs(pdfs_from_counts(C), C["grid"], spline_ops(phi, lam))
            I_fd = finite_difference_I(C, phi)
            ls.append(lam)
            es.append(np.median(np.abs(I_sm[core] / I_true[core] - 1)))
            ef.append(np.nanmedian(np.abs(I_fd[core] / I_true[core] - 1)))
        print(f"{nper:5d} durations/run: lam={np.median(ls):.0e}  "
              f"median |I/I_true - 1|: smoothed={np.mean(es):.2f}  "
              f"finite-diff={np.mean(ef):.2f}")

    rng = np.random.default_rng(0)
    tb = [[np.maximum(rng.lognormal(mu[k], s, 400).round(), 1.0) for _ in range(16)]
          for k in range(16)]
    tu = [[np.maximum(rng.lognormal(5.0, s, 400).round(), 1.0) for _ in range(16)]
          for k in range(16)]
    Cb, Cu = build_counts(tb), build_counts(tu)
    print("mean_dwell vs direct mean: %.1f vs %.1f"
          % (mean_dwell(Cb)[8], np.mean(np.concatenate(tb[8]))))
    lam_b, _ = select_lambda(Cb, phi, lams)
    lam_u, _ = select_lambda(Cu, phi, lams)
    t0 = time.time()
    fb, fu, nrm = rate(Cb, Cu, phi, lam_b, lam_u)
    Fb, Fu = bootstrap_rate(Cb, Cu, phi, lam_b, lam_u, nboot=200)
    lo, hi = np.percentile(Fb, [5, 95], axis=0)
    print("200 bootstrap replicates in %.1f s; band covers point estimate at %d/16"
          % (time.time() - t0, np.sum((fb >= lo) & (fb <= hi))))
    print("median relative half-width %.2f" % np.median((hi - lo) / (2 * np.abs(fb))))
    # the open channel was generated with no phi dependence: its Fisher rate is
    # the pure noise floor, and cross_fisher should put it near zero
    c_u = phi**2 * cross_fisher(Cu, phi, lam_u) / nrm
    print("null channel: ordinary %.2e, cross-split %+.2e, signal channel %.2e"
          % (np.median(fu), np.median(c_u), np.median(fb)))
