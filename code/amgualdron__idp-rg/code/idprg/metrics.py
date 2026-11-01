"""Scoring a predicted distribution against the MD frames, per protein:
    nll        -ln p(Rg) per frame over all 30,000 stored frames (density in Rg, nm^-1)
    mean, std, skew of the predicted distribution (quadrature on a 10,001-point ln Rg grid around the MD mean), and their
               errors against the full-trajectory MD moments: <Rg> and sigma error in %, skewness error
    w1_pct     Wasserstein-1 distance between the predicted distribution and the stored frames (exact, from the two
               CDFs), as a percentage of the MD <Rg>
    compact    MD <Rg> below 0.8 x the Flory baseline exp(LN_A_F + NU_F ln N)
"""
import jax.numpy as jnp
import numpy as np
import pandas as pd

from .data import LN_A_F, NU_F
from .distribution import jax_logpdf, nu_of

GRID = np.linspace(-3.0, 2.0, 10001)              # ln(Rg / MD <Rg>)


def logp(x, lnA, lnSm, u):
    return np.asarray(jax_logpdf(jnp.asarray(x), *(jnp.asarray(a) for a in (lnA, lnSm, u))))


def score(lnA, lnSm, u, proteins, F):
    """Per-protein scores; `proteins` are the scored rows of the protein table, the parameters aligned with them."""
    lnA, lnSm, u = (np.asarray(a, float) for a in (lnA, lnSm, u))
    rows, N = proteins.frame_row.to_numpy(), proteins.N.to_numpy(float)
    m_md, s_md, k_md = (proteins[c].to_numpy(float) for c in ("rg_mean", "rg_std", "rg_skew"))
    out = []
    for b in range(0, len(proteins), 64):
        i = np.arange(b, min(b + 64, len(proteins)))
        p = (lnA[i], lnSm[i], u[i])
        fr = np.sort(np.asarray(F[rows[i]], float), axis=1)
        nll = -logp(np.log(fr), *p).mean(1)
        x = np.log(m_md[i])[:, None] + GRID[None, :]
        r = np.exp(x)
        dm = np.exp(logp(x, *p)) * r                                   # d mass / d ln Rg
        mass = np.trapezoid(dm, x, axis=1)
        mean = np.trapezoid(dm * r, x, axis=1) / mass
        c = r - mean[:, None]
        var = np.trapezoid(dm * c ** 2, x, axis=1) / mass
        skew = np.trapezoid(dm * c ** 3, x, axis=1) / mass / var ** 1.5
        cdf = np.concatenate([np.zeros((len(i), 1)), np.cumsum(0.5 * (dm[:, 1:] + dm[:, :-1]) * np.diff(x), 1)], 1) / mass[:, None]
        emp = np.stack([np.searchsorted(fr[k], r[k], side="right") / fr.shape[1] for k in range(len(i))])
        w1 = np.trapezoid(np.abs(cdf - emp), r, axis=1)
        out.append(pd.DataFrame(dict(id=proteins.id.to_numpy()[i], N=N[i], lnA=lnA[i], lnSm=lnSm[i], nu=nu_of(u[i]),
                                     nll=nll, mass=mass, mean=mean, std=np.sqrt(var), skew=skew, md_mean=m_md[i],
                                     md_std=s_md[i], md_skew=k_md[i], w1_pct=100 * w1 / m_md[i])))
    d = pd.concat(out, ignore_index=True)
    d["mean_err_pct"] = 100 * (d["mean"] / d.md_mean - 1)
    d["std_err_pct"] = 100 * (d["std"] / d.md_std - 1)
    d["skew_err"] = d["skew"] - d.md_skew
    d["compact"] = d.md_mean / np.exp(LN_A_F + NU_F * np.log(d.N)) < 0.8
    return d


def summary(d):
    a = lambda c: float(np.abs(d[c]).mean())
    return dict(n=len(d), n_compact=int(d.compact.sum()), nll=float(d.nll.mean()), mean_err_pct=a("mean_err_pct"),
                bias_pct=float(d.mean_err_pct.mean()), std_err_pct=a("std_err_pct"), skew_mae=a("skew_err"),
                w1_pct_median=float(d.w1_pct.median()), w1_pct_mean=float(d.w1_pct.mean()),
                w1_pct_p90=float(d.w1_pct.quantile(0.9)),
                compact_mean_err_pct=float(np.abs(d.mean_err_pct[d.compact]).mean()) if d.compact.any() else np.nan,
                max_mass_dev=float(np.abs(d["mass"] - 1).max()))
