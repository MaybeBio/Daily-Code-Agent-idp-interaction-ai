"""The two-exponent distribution of the radius of gyration (Victor & Lhuillier, J. Chem. Phys. 92, 1362, 1990; d = 3):

    s = Rg / S_m,   p(s) ds  propto  exp(-A [ s^delta / delta + s^-g / g ]) ds,   delta = 1/(1 - nu),  g = 3/(3 nu - 1)

S_m is the mode of p(s), A sets the width, and the exponent nu in (1/3, 1) the asymmetry: the s^delta term is the
stretched-chain tail, the s^-g term the excluded-volume wall at small Rg. nu is parametrized as
nu = 1/3 + (2/3) expit(u). The normalisation and the moments are computed by quadrature on t = ln s in [-3, 2].

numpy functions take (ln A, ln S_m, nu) and broadcast over proteins; the jax functions take (ln A, ln S_m, u) and are
what the network training, the closed-form fit and the scoring use. Importing idprg enables float64 in JAX.
"""
import numpy as np
from scipy.special import expit, logit, logsumexp

import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp  # noqa: E402

D = 3
T = np.linspace(-3.0, 2.0, 4001)
_W = np.full(len(T), T[1] - T[0]); _W[[0, -1]] /= 2
LOG_W = np.log(_W)                                 # trapezoid weights on T


def nu_of(u):
    return 1 / 3 + (2 / 3) * expit(u)


def u_of(nu):
    return logit(1.5 * (np.asarray(nu, float) - 1 / 3))


def exponents(nu):
    """(g, delta): the exponents of the small-s wall and of the large-s tail."""
    return D / (nu * D - 1), 1.0 / (1.0 - nu)


def _f(t, g, dl, xp=np):
    return xp.exp(xp.minimum(dl * t, 700.0)) / dl + xp.exp(xp.minimum(-g * t, 700.0)) / g


def _log_mass(lnA, nu):
    """ln of the unnormalised mass per grid point of t = ln s, shape (..., len(T))."""
    lnA, nu = np.asarray(lnA, float)[..., None], np.asarray(nu, float)[..., None]
    g, dl = exponents(nu)
    return -np.exp(lnA) * _f(T, g, dl) + T + LOG_W


def logpdf(rg, lnA, lnSm, nu):
    """ln p(Rg) (density in Rg, nm^-1); rg broadcasts against the parameters' shape plus a trailing axis."""
    lnA, lnSm, nu = (np.asarray(a, float) for a in (lnA, lnSm, nu))
    g, dl = exponents(nu)
    logz = logsumexp(_log_mass(lnA, nu), axis=-1) + lnSm
    x = np.log(rg)
    if lnA.ndim:
        g, dl, lnA, lnSm, logz = (a[..., None] for a in (g, dl, lnA, lnSm, logz))
    return -np.exp(lnA) * _f(x - lnSm, g, dl) - logz


def pdf(rg, lnA, lnSm, nu):
    return np.exp(logpdf(rg, lnA, lnSm, nu))


def grid(lnA, lnSm, nu):
    """Rg (nm) on the quadrature grid and the normalised mass per unit ln Rg there, each (..., len(T))."""
    lm = _log_mass(lnA, nu)
    w = np.exp(lm - logsumexp(lm, axis=-1, keepdims=True)) / _W
    return np.exp(np.asarray(lnSm, float)[..., None] + T), w


def moments(lnA, lnSm, nu):
    """<Rg>, sigma(Rg) (nm), skewness and excess kurtosis."""
    r, w = grid(lnA, lnSm, nu)
    m = np.trapezoid(w * r, T, axis=-1)
    c = r - m[..., None]
    v = np.trapezoid(w * c ** 2, T, axis=-1)
    return (m, np.sqrt(v), np.trapezoid(w * c ** 3, T, axis=-1) / v ** 1.5,
            np.trapezoid(w * c ** 4, T, axis=-1) / v ** 2 - 3)


def cdf(rg, lnA, lnSm, nu):
    """P(Rg' <= rg), by interpolation of the quadrature CDF (scalar parameters)."""
    r, w = grid(lnA, lnSm, nu)
    c = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(T))])
    return np.interp(np.log(rg), np.log(r), c / c[-1])


def quantiles(q, lnA, lnSm, nu):
    """Rg (nm) at probabilities q (scalar parameters)."""
    r, w = grid(lnA, lnSm, nu)
    c = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(T))])
    return np.exp(np.interp(q, c / c[-1], np.log(r)))


def tail_below(ratio, lnA, nu, block=2000):
    """P(Rg < ratio x <Rg>) per chain, over arrays of chains. A function of the shape alone, so S_m plays no part."""
    ratio, lnA, nu = (np.asarray(a, float) for a in (ratio, lnA, nu))
    out, dT = np.empty(len(ratio)), T[1] - T[0]
    for b in range(0, len(ratio), block):
        i = slice(b, b + block)
        r, w = grid(lnA[i], np.zeros(len(lnA[i])), nu[i])
        mean = np.trapezoid(w * r, T, axis=-1)
        c = np.concatenate([np.zeros((w.shape[0], 1)), np.cumsum(0.5 * (w[:, 1:] + w[:, :-1]) * dT, axis=1)], axis=1)
        c /= c[:, -1:]
        t = np.log(ratio[i] * mean)
        k = np.clip(((t - T[0]) / dT).astype(int), 0, len(T) - 2)
        f, n = (t - T[k]) / dT, np.arange(len(k))
        out[i] = c[n, k] * (1 - f) + c[n, k + 1] * f
    return out


# Compact and extended populations. The bounds are the 5th and 95th percentiles of R_g / <R_g> of a typical chain of the
# same length: per training protein, the percentile over its MD frames; then the median per length bin, fitted by
# a + b / sqrt(N). So a typical chain has 5% of its conformations in each population at any N, and a chain's value reads
# directly against that. proteome/populations.py refits these constants and records why 5% (populations_thresholds.csv).
COMPACT_BOUND, EXTENDED_BOUND = (0.681, 0.311), (1.421, -0.571)


def population_bounds(N):
    """(t_c, t_e): a conformation is compact if Rg < t_c <Rg> and extended if Rg > t_e <Rg>."""
    s = np.asarray(N, float) ** -0.5
    return COMPACT_BOUND[0] + COMPACT_BOUND[1] * s, EXTENDED_BOUND[0] + EXTENDED_BOUND[1] * s


def populations(lnA, nu, N):
    """(P_c, P_e) per chain."""
    tc, te = population_bounds(N)
    return tail_below(tc, lnA, nu), 1 - tail_below(te, lnA, nu)


# ---------------------------------------------------------------- jax (training, fitting, scoring)
T_J, LOG_W_J = jnp.asarray(T), jnp.log(jnp.asarray(_W))


def jax_terms(lnA, lnSm, u):
    """(A, g, delta, ln Z) per protein, ln Z the normaliser of the density in Rg."""
    nu = 1 / 3 + (2 / 3) * jax.nn.sigmoid(u)
    g, dl = D / (nu * D - 1), 1.0 / (1.0 - nu)
    A = jnp.exp(lnA)
    lg = -A[:, None] * _f(T_J[None, :], g[:, None], dl[:, None], jnp) + T_J[None, :] + LOG_W_J[None, :]
    return A, g, dl, jax.scipy.special.logsumexp(lg, axis=1) + lnSm


def jax_logpdf(x, lnA, lnSm, u):
    """ln p(Rg) at x = ln Rg, x of shape (proteins, points)."""
    A, g, dl, logz = jax_terms(lnA, lnSm, u)
    return -A[:, None] * _f(x - lnSm[:, None], g[:, None], dl[:, None], jnp) - logz[:, None]


def jax_logpdf_flat(x, seg, lnA, lnSm, u):
    """ln p(Rg) at x = ln Rg for a flat array of points, seg[i] the protein of point i."""
    A, g, dl, logz = jax_terms(lnA, lnSm, u)
    return -A[seg] * _f(x - lnSm[seg], g[seg], dl[seg], jnp) - logz[seg]
