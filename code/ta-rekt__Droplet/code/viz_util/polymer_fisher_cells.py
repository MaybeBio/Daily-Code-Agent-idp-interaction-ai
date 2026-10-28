"""
Notebook chunks for the smoothed polymer Fisher rate with bootstrap error bars.

Paste each `# %% CELL n` block below into its own cell, after the cell that
defines `bound_times` / `unbound_times` (the compute_hysteresis_times cell).
Requires `fisher_rate_base` from the local-Fisher-rate cell for the dashed
reference curve in CELL 5; drop that line if you have not run it.

Convention taken from the notebook: bound_times = closed (contact) durations,
unbound_times = open durations.  Cell 108's swapped unpacking of
compute_hysteresis_times is what makes that true.
"""

# %% CELL 1 -- pre-bin the runs
import importlib
import numpy as np
import matplotlib.pyplot as plt
import viz_util.polymer_fisher as pf
importlib.reload(pf)

phi_inf = np.linspace(0.003, 0.018, 16)
n_grid, sigma = 500, 10.0          # as in the empirical-PDF cell

Cb = pf.build_counts(bound_times, n_grid=n_grid)     # closed durations
Cu = pf.build_counts(unbound_times, n_grid=n_grid)   # open durations
print("durations per concentration: closed %s" % Cb["nsamp"].sum(1).astype(int))
print("                             open   %s" % Cu["nsamp"].sum(1).astype(int))
print("runs per concentration: %d, time grid %.3g - %.3g"
      % (Cb["counts"].shape[1], Cb["grid"][0], Cb["grid"][-1]))

# %% CELL 2 -- pick the spline bandwidth along phi (one per channel)
lams = np.logspace(-7, -1, 13)
lam_b, mse_b = pf.select_lambda(Cb, phi_inf, lams, sigma=sigma, seed=0)
lam_u, mse_u = pf.select_lambda(Cu, phi_inf, lams, sigma=sigma, seed=0)
print("lambda closed = %.1e   open = %.1e" % (lam_b, lam_u))

fig, ax = plt.subplots(figsize=(4, 2.2))
ax.plot(lams, mse_b / mse_b.min(), "o-", ms=3, color="crimson", label="closed")
ax.plot(lams, mse_u / mse_u.min(), "o-", ms=3, color="green", label="open")
for lm, c in ((lam_b, "crimson"), (lam_u, "green")):
    ax.axvline(lm, color=c, ls="dotted", lw=1)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$\lambda$"); ax.set_ylabel("estimated MSE of $I(\\phi)$ (rel.)")
ax.legend(fontsize=8); plt.tight_layout()
# A minimum at either edge of the grid means the grid is too narrow -- widen it.

# %% CELL 3 -- point estimate + run-level bootstrap
rate_closed, rate_open, norm = pf.rate(Cb, Cu, phi_inf, lam_b, lam_u, sigma=sigma)
rate_sum = rate_closed + rate_open

Fb, Fu = pf.bootstrap_rate(Cb, Cu, phi_inf, lam_b, lam_u,
                           nboot=400, seed=1, sigma=sigma)
band_closed = np.percentile(Fb, [5, 95], axis=0)
band_open = np.percentile(Fu, [5, 95], axis=0)
band_sum = np.percentile(Fb + Fu, [5, 95], axis=0)    # same resample -> correlated
print("normaliser (mean closed dwell at mid phi) = %.4g" % norm)
print("median relative half-width: closed %.2f  open %.2f  sum %.2f"
      % tuple(np.median((b[1] - b[0]) / (2 * np.abs(r)))
              for b, r in ((band_closed, rate_closed), (band_open, rate_open),
                           (band_sum, rate_sum))))

# %% CELL 4 -- is the structure real, or the score-noise floor?
Ic = pf.cross_fisher(Cb, phi_inf, lam_b, sigma=sigma)
Iu_c = pf.cross_fisher(Cu, phi_inf, lam_u, sigma=sigma)
cross_closed = phi_inf**2 * Ic / norm
cross_open = phi_inf**2 * Iu_c / norm
print("closed: ordinary vs cross-split (ratio, bulk) %.2f"
      % np.median((rate_closed / np.where(cross_closed == 0, np.nan, cross_closed))[3:13]))
print("negative excursions of the cross-split estimate mark the noise floor:")
print("   closed min %.2e   open min %.2e   signal scale %.2e"
      % (cross_closed.min(), cross_open.min(), np.median(rate_sum)))

# %% CELL 5 -- figure
fig, ax = plt.subplots(figsize=(5, 2.6))
xticks = np.linspace(2.5, 17.5, 7)

for r, bnd, c, lab in ((rate_closed, band_closed, "crimson", "closed"),
                       (rate_open, band_open, "green", "open"),
                       (rate_sum, band_sum, "navy", "sum")):
    ax.fill_between(phi_inf, bnd[0], bnd[1], color=c, alpha=0.18, lw=0)
    ax.plot(phi_inf, r, "-", color=c, label=lab)

ax.plot(phi_inf, fisher_rate_base, ls="dashed", color="black",
        label=r"$\dot{\mathcal{F}}(c)$")
ax.plot(phi_inf, cross_closed, ":", color="crimson", lw=1, alpha=0.8)
ax.plot(phi_inf, cross_open, ":", color="green", lw=1, alpha=0.8)

ax.set_xticks(xticks / 1000); ax.set_xticklabels(xticks)
ax.set_xlabel(r"Concentration  $c\times 10^{-3}$", fontsize=12)
ax.set_ylabel(r"Polymer Fisher info. rate", fontsize=12)
ax.set_yscale("log"); ax.set_ylim(5e-4, 5)
ax.legend(fontsize=8)
plt.tight_layout()
plt.show()

# %% CELL 6 (optional) -- compare against the 3-point estimator you had
I_fd_b = pf.finite_difference_I(Cb, phi_inf, sigma=sigma)
I_fd_u = pf.finite_difference_I(Cu, phi_inf, sigma=sigma)
fig, ax = plt.subplots(figsize=(5, 2.6))
ax.fill_between(phi_inf, band_sum[0], band_sum[1], color="navy", alpha=0.18, lw=0)
ax.plot(phi_inf, rate_sum, color="navy", label="smoothed score")
ax.plot(phi_inf, phi_inf**2 * (I_fd_b + I_fd_u) / norm, "o--", ms=3,
        color="grey", label="3-point difference")
ax.set_xticks(xticks / 1000); ax.set_xticklabels(xticks)
ax.set_yscale("log"); ax.set_ylim(5e-4, 5)
ax.set_xlabel(r"Concentration  $c\times 10^{-3}$", fontsize=12)
ax.set_ylabel(r"Polymer Fisher info. rate", fontsize=12)
ax.legend(fontsize=8); plt.tight_layout(); plt.show()
