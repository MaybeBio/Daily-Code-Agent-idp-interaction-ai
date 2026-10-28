"""
Fisher information rate of the joint (T_open, T_closed) first-passage-time
distribution with respect to phi_inf, computed from Langevin trajectories of
the Rouse polymer in spec/rouse_mfpt_attractive.jl.

T_open and T_closed are conditionally independent once the inter-bead
distance crosses the potential ridge, so log P(T_open, T_closed | phi) splits
and the joint Fisher rate is the sum of the marginals.

Per the user's definition,
    dot{I}_X(phi) = -phi^2 / <T_X> * < d^2 log P(T_X|phi) / dphi^2 >_phi.
We use the equivalent first-derivative form
    -<d^2 log P / dphi^2> = <(d log P / dphi)^2> = integral P (d log P/dphi)^2 dT,
which is numerically more stable when P is estimated from histograms.

Pipeline:
  1. Load distance trajectories from data_dir, one record per phi_inf.
  2. For each record locate the barrier ridge on the recorded total_potential.
  3. Cut each trajectory into open/closed dwells and collect T_open, T_closed.
  4. Histogram in log T, smooth, interpolate to a common T grid, renormalise.
  5. Differentiate log P along phi_inf, integrate against P to get the
     per-event Fisher info; divide by <T_X> to get the rate.
"""

from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
import os
import glob
import json

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

matplotlib.use("Agg")


# -- Config ------------------------------------------------------------------
EXP_NAME = "rouse_mfpt_attractive"
DATA_DATE = "20260528"
OUT_DATE = "20260603"
DATA_DIR = f"data/{DATA_DATE}/{EXP_NAME}"
OUT_PNG = f"media/figures/{OUT_DATE}/{EXP_NAME}_fisher_hires.pdf"
OUT_PNG2 = f"media/figures/{OUT_DATE}/{EXP_NAME}_fisher_rates_hires.pdf"
OUT_JSON = f"media/figures/{OUT_DATE}/{EXP_NAME}_fisher_hires.json"

N_LOG_T_BINS = 60       # histogram bins in log T
SMOOTH_SIGMA = 1.5      # gaussian smoothing in log-T bins
N_T_INTERP = 400      # interpolation grid size in T
PDF_FLOOR = 1e-12    # protects log(P) and division at empty bins


# -- Load --------------------------------------------------------------------
files = sorted(glob.glob(os.path.join(DATA_DIR, "output_cpu*.json")))
if not files:
    sys.exit(f"no output_cpu*.json under {DATA_DIR}")

records = []
for fp in files:
    with open(fp) as f:
        blob = json.load(f)
    for _, by_beta in blob.items():
        for _, rec in by_beta.items():
            records.append(rec)
records.sort(key=lambda r: r["surface_params"]["phi_inf"])
phi_inf = np.array([r["surface_params"]["phi_inf"] for r in records])
print(f"loaded {len(records)} records from {len(files)} files")
print(f"phi_inf range: [{phi_inf.min():g}, {phi_inf.max():g}]")


# -- Per-record ridge -------------------------------------------------------
def potential_grid(rec):
    U = np.array([np.nan if x is None else x for x in rec["total_potential"]],
                 dtype=float)
    N = rec["rouse_params"]["N"]
    a = rec["rouse_params"]["kuhn_length"]
    d = np.linspace(0.0, 3.0 * np.sqrt(N) * a, len(U))
    return d, U


def locate_ridge(d, U):
    """First local maximum that lies between the inner (closed) and outer
    (open) wells. Falls back to argmax over the finite range."""
    m = np.isfinite(U)
    d_m, U_m = d[m], U[m]
    dU = np.diff(U_m)
    transitions = np.where(np.diff(np.sign(dU)))[0] + 1
    # walk: need min, then max, then min
    saw_min = False
    for i in transitions:
        is_max = dU[i - 1] > 0
        if is_max and saw_min:
            return float(d_m[i])
        if not is_max:
            saw_min = True
    return float(d_m[np.argmax(U_m)])


d_ridge_per_phi = np.array([
    locate_ridge(*potential_grid(r)) for r in records
])


# -- First-passage extraction -----------------------------------------------
def extract_fpts(traj, d_ridge, dt_per_sample):
    """Cut a trajectory at ridge crossings into alternating open/closed
    dwell times. Drops the leading and trailing (censored) partial dwells."""
    state = traj >= d_ridge                # True = open
    change = np.flatnonzero(np.diff(state.astype(np.int8)))
    if len(change) < 2:
        return np.empty(0), np.empty(0)
    seg_lens = np.diff(change) * dt_per_sample
    seg_states = state[change[:-1] + 1]
    return seg_lens[seg_states], seg_lens[~seg_states]


t_open_all = []
t_closed_all = []
for r, d_ridge in zip(records, d_ridge_per_phi):
    traj = np.asarray(r["distance_traj"], dtype=float)
    dt = r["rouse_params"]["dt"]
    nsteps = r["rouse_params"]["nsteps"]
    sample_skip = r.get("sample_skip", max(1, nsteps // len(traj)))
    dt_per_sample = dt * sample_skip
    to, tc = extract_fpts(traj, d_ridge, dt_per_sample)
    t_open_all.append(to)
    t_closed_all.append(tc)
    print(f"phi_inf={r['surface_params']['phi_inf']:.4f}  "
    f"d_ridge={d_ridge:.3f}  "
    f"open n={len(to):>5d} <T>={to.mean() if to.size else float('nan'):>8.3f}  "
    f"closed n={len(tc):>5d} <T>={tc.mean() if tc.size else float('nan'):>8.3f}")


# -- Smoothed PDFs and Fisher rate -----------------------------------------
def fisher_rate_family(samples_per_phi, phi_inf,
                       n_bins=N_LOG_T_BINS, sigma=SMOOTH_SIGMA,
                       n_interp=N_T_INTERP):
    """Fisher information rate (per unit time) for one FPT family,
    one value per phi_inf."""

    finite = np.concatenate([s for s in samples_per_phi if s.size > 0])
    finite = finite[finite > 0]
    log_edges = np.linspace(np.log(finite.min()),
                            np.log(finite.max()), n_bins + 1)
    log_t_centers = 0.5 * (log_edges[:-1] + log_edges[1:])
    T_grid = np.exp(np.linspace(log_t_centers[0], log_t_centers[-1], n_interp))

    pdfs = np.full((len(samples_per_phi), n_interp), np.nan)
    mean_T = np.full(len(samples_per_phi), np.nan)
    for i, s in enumerate(samples_per_phi):
        if s.size == 0:
            continue
        counts, _ = np.histogram(np.log(s), bins=log_edges)
        density_log = counts.astype(float) / (s.size * np.diff(log_edges))
        density_log = gaussian_filter1d(density_log, sigma=sigma)
        # density on log T -> density on T (Jacobian 1/T)
        f = interp1d(log_t_centers, density_log, kind="cubic",
                     bounds_error=False, fill_value=0.0)
        pT = np.clip(f(np.log(T_grid)), 0.0, None) / T_grid
        norm = np.trapz(pT, T_grid)
        if norm > 0:
            pT = pT / norm
        pdfs[i] = pT
        mean_T[i] = s.mean()

    log_pdfs = np.log(np.clip(pdfs, PDF_FLOOR, None))
    # d log P / d phi at each (phi, T) via central differences along phi axis
    dlogP_dphi = np.gradient(log_pdfs, phi_inf, axis=0)
    integrand = pdfs * dlogP_dphi ** 2
    fisher_per_event_phi = np.trapz(
        np.where(np.isfinite(integrand), integrand, 0.0), T_grid, axis=1)
    fisher_per_event_logphi = phi_inf ** 2 * fisher_per_event_phi
    with np.errstate(divide="ignore", invalid="ignore"):
        fisher_rate = fisher_per_event_logphi / mean_T

    return dict(
        T_grid=T_grid,
        pdfs=pdfs,
        mean_T=mean_T,
        fisher_per_event=fisher_per_event_logphi,  # Fisher info wrt log phi_inf
        fisher_rate=fisher_rate,                   # per unit simulation time
    )


open_res = fisher_rate_family(t_open_all,   phi_inf)
closed_res = fisher_rate_family(t_closed_all, phi_inf)
joint_rate = open_res["fisher_rate"] + closed_res["fisher_rate"]


# -- Dump JSON --------------------------------------------------------------
os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)


def _j(a):
    if isinstance(a, np.ndarray):
        return np.where(np.isfinite(a), a, None).tolist()
    return a


with open(OUT_JSON, "w") as f:
    json.dump(dict(
        blurb=records[0].get("blurb", ""),
        phi_inf=_j(phi_inf),
        d_ridge=_j(d_ridge_per_phi),
        open=dict(
            T_grid=_j(open_res["T_grid"]),
            pdfs=_j(open_res["pdfs"]),
            mean_T=_j(open_res["mean_T"]),
            fisher_per_event=_j(open_res["fisher_per_event"]),
            fisher_rate=_j(open_res["fisher_rate"]),
        ),
        closed=dict(
            T_grid=_j(closed_res["T_grid"]),
            pdfs=_j(closed_res["pdfs"]),
            mean_T=_j(closed_res["mean_T"]),
            fisher_per_event=_j(closed_res["fisher_per_event"]),
            fisher_rate=_j(closed_res["fisher_rate"]),
        ),
        joint_rate=_j(joint_rate),
    ), f)
print(f"wrote -> {OUT_JSON}")


# -- Plot: PDFs and MFPTs --------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
cmap = plt.get_cmap("viridis")

for col, (name, res) in enumerate([("T_open", open_res),
                                   ("T_closed", closed_res)]):
    ax = axes[0, col]
    for i in range(len(phi_inf)):
        if not np.all(np.isnan(res["pdfs"][i])):
            ax.loglog(res["T_grid"],
                      np.clip(res["pdfs"][i], PDF_FLOOR, None),
                      color=cmap(i / max(1, len(phi_inf) - 1)),
                      lw=1.1)
    ax.set_xlabel("T")
    ax.set_ylabel(f"P({name} | $\\phi_\\infty$)")
    ax.set_title(f"Smoothed PDF of {name}")
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[1, col]
    ax.plot(phi_inf, res["mean_T"], "o-", lw=1.5)
    ax.set_xlabel(r"$\phi_\infty$")
    ax.set_ylabel(rf"$\langle {name} \rangle$")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)

# colour bar for phi_inf
sm = plt.cm.ScalarMappable(cmap=cmap,
                           norm=plt.Normalize(vmin=phi_inf.min(), vmax=phi_inf.max()))
cb = fig.colorbar(sm, ax=axes[0, :], orientation="vertical",
                  fraction=0.025, pad=0.04)
cb.set_label(r"$\phi_\infty$")

fig.suptitle("Rouse-MFPT first-passage statistics", fontsize=12)
fig.tight_layout(rect=[0, 0, 0.95, 0.96])
os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
print(f"saved -> {OUT_PNG}")


# -- Plot: Fisher rate vs phi_inf ------------------------------------------
fig2, ax = plt.subplots(1, 1, figsize=(5.8, 4.0))
ax.plot(phi_inf, open_res["fisher_rate"],   "o-", color="C0",
        label=r"$\dot I_{T_{\rm open}}$")
ax.plot(phi_inf, closed_res["fisher_rate"], "s-", color="C1",
        label=r"$\dot I_{T_{\rm closed}}$")
ax.plot(phi_inf, joint_rate, "^-", color="k", lw=2,
        label=r"sum (joint, conditionally indep.)")
ax.set_xlabel(r"$\phi_\infty$")
ax.set_ylabel(r"Fisher information rate (per unit time)")
ax.set_yscale("log")
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)
fig2.tight_layout()
fig2.savefig(OUT_PNG2, dpi=150, bbox_inches="tight")
print(f"saved -> {OUT_PNG2}")
