"""
Full 3D Rouse chain vs the 1D Markovian reduction: exit and entry first-passage
times, their medians, and the ratio between the two models.

Left/middle columns: the two first-passage problems that the Fisher-rate pipeline
rests on -- EXIT (contact -> x_div) and ENTRY (x_div -> contact) -- with the Rouse
chain's sample mean and median (spec/rouse_mfpt_attractive.jl, one sojourn stream
per phi_inf) against the backward-FPE solver's MFPT and median
(spec/fisher_rate_from_fpt_solver.jl, the same U_tot, the same landmarks, the same
phi_inf grid, D_eff = 2 kBT/gamma).  Right column: the ratio between them vs the
barrier height, and the survival functions with each phi_inf's time rescaled by its
own solver MFPT, which separates a pure time-scale shift from a change of shape.

Rouse means are censoring-corrected: the sojourn still open when the trajectory ends
is counted as time-at-risk (MLE for the mean), so a phi_inf whose MFPT approaches the
run length is not biased low by discarding it.

    python viz_util/rouse_vs_1d_exit_entry.py [rouse_dir] [fisher_json] [out_dir]
"""

import glob
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

ROUSE_DIR = sys.argv[1] if len(sys.argv) > 1 else "data/20260912/rouse"
FISHER    = sys.argv[2] if len(sys.argv) > 2 else "data/20260912/test3/fisher_rate_from_fpt_solver.json"
OUT_DIR   = sys.argv[3] if len(sys.argv) > 3 else "media/figures/20260914"

C_ROUSE, C_1D, C_ENTRY = "#cc3311", "#1f77b4", "#2e7d5b"
LOW_BARRIER, HIGH_BARRIER = 4.5, (7.5, 10.7)     # the two survival groups

# ── load ─────────────────────────────────────────────────────────────────────
data = {}
for _p in sorted(glob.glob(f"{ROUSE_DIR}/output_cpu*.json")):
    for _keya, _runs in json.load(open(_p)).items():
        for _run in _runs.values():
            data[_run["alpha"]] = _run
fisher = json.load(open(FISHER))
solver = {r["phi_inf"]: r for r in fisher["per_phi"]}

phis = np.array(sorted(data))
recs = [solver[min(solver, key=lambda p: abs(p - a))] for a in phis]


def censored_mean(run, which):
    """MLE mean sojourn: the final, still-open interval counts as time-at-risk.

    The chain starts from an equilibrium free chain (<x> ~ a sqrt(s) >> x_contact),
    so it starts in B and the pushes alternate B,A,B,A...; the first interval of each
    list is pinned to step 0 and is dropped, exactly as the spec does."""
    A = np.asarray(run["A_intervals"], float)
    B = np.asarray(run["B_intervals"], float)
    dt, nst = run["rouse_params"]["dt"], run["rouse_params"]["nsteps"]
    open_steps = (nst - 1) - (A.sum() + B.sum())
    open_in_A = len(B) == len(A) + 1                  # last push was B => now in A
    ints, extra = (A, open_steps if open_in_A else 0.0) if which == "escape" \
                  else (B, 0.0 if open_in_A else open_steps)
    n = len(ints) - 1
    return (ints[1:].sum() + extra) / n * dt if n > 0 else np.nan


def solver_quantile(rec, key, q=0.5):
    """Quantile of the solver's FPT density, from its stored (t, pdf) log grid."""
    t, p = np.asarray(rec[f"{key}_t"]), np.asarray(rec[f"{key}_pdf"])
    c = np.concatenate([[0.0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(t))])
    return float(np.interp(q, c / c[-1], t))


def solver_survival(rec, key):
    t, p = np.asarray(rec[f"{key}_t"]), np.asarray(rec[f"{key}_pdf"])
    c = np.concatenate([[0.0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(t))])
    return t, np.clip(1.0 - c / c[-1], 0.0, 1.0)


barrier = np.array([d["barrier"] for d in (data[a] for a in phis)])
tab = {}
for key, nkey in (("escape", "n_escape"), ("return", "n_return")):
    t_r  = np.array([np.asarray(data[a][f"{key}_times"], float) for a in phis], dtype=object)
    n    = np.array([len(x) for x in t_r])
    mean = np.array([censored_mean(data[a], key) for a in phis])
    med  = np.array([np.median(x) if len(x) else np.nan for x in t_r])
    sem  = np.array([x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan for x in t_r])
    s_mean = np.array([r[f"{key}_mfpt"] for r in recs])
    s_med  = np.array([solver_quantile(r, key) for r in recs])
    tab[key] = dict(t=t_r, n=n, mean=mean, med=med, sem=sem, s_mean=s_mean, s_med=s_med)

# ── figure ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(2, 3, figsize=(16.5, 8.6))
p0 = fisher["params"]
fig.suptitle(
    f"Full 3D Rouse chain ({ROUSE_DIR}) vs the 1D Markovian reduction "
    f"({os.path.dirname(FISHER)})  --  "
    f"$s$={p0['rouse']['s']:.0f}, $a$={p0['rouse']['a']:.0f}, $R$={p0['surface']['R']:.0f}, "
    f"$x_{{\\rm div}}$={fisher['x_div']:.1f}, "
    f"$\\tau_{{\\rm Rouse}}$={data[phis[0]]['rouse_params']['tau_rouse']:.2e} s",
    fontsize=12)

for j, (key, title) in enumerate((("escape", "(a)  EXIT: contact $\\to x_{\\rm div}$"),
                                  ("return", "(b)  ENTRY: $x_{\\rm div} \\to$ contact"))):
    a_, T = ax[0][j], tab[key]
    ok = T["n"] > 0
    a_.plot(phis[ok], T["med"][ok], "^", mfc="none", ms=6, color=C_ROUSE, label="Rouse, median")
    a_.plot(phis, T["s_mean"], "-",  color=C_1D, lw=2, label="1D reduction, MFPT")
    a_.plot(phis, T["s_med"],  "--", color=C_1D, lw=2, label="1D reduction, median")
    a_.errorbar(phis[ok], T["mean"][ok], yerr=T["sem"][ok], fmt="o", ms=5,
                color=C_ROUSE, capsize=2, lw=1, label="Rouse 3D chain, mean")
    a_.set_yscale("log")
    if (~ok).any():
        a_.axvspan(phis[~ok].min(), phis.max(), color="0.85", zorder=0)
        a_.text(phis[~ok].min(), 0.02, " no\n events", transform=a_.get_xaxis_transform(),
                va="bottom", fontsize=9, color="0.4")
    a_.set_xlabel(r"$\phi_\infty$", fontsize=12)
    a_.set_ylabel(f"{'exit' if key == 'escape' else 'entry'} time  [s]", fontsize=11)
    a_.set_title(title, fontsize=12, loc="left")
    a_.legend(fontsize=9)
    a_.grid(alpha=0.25)

# (c) ratio vs barrier
a_ = ax[0][2]
for key, col, m1, m2, lab in (("escape", C_ROUSE, "o", "^", "exit"),
                              ("return", C_ENTRY, "s", "v", "entry")):
    T = tab[key]; ok = T["n"] > 0
    a_.plot(barrier[ok], (T["mean"] / T["s_mean"])[ok], m1 + "-", color=col, ms=5,
            lw=1.2, label=f"{lab}, mean ratio")
    a_.plot(barrier[ok], (T["med"] / T["s_med"])[ok], m2 + "--", color=col, ms=5,
            mfc="none", lw=1.0, alpha=0.8, label=f"{lab}, median ratio")
a_.axhline(1.0, color="k", lw=1.2)
lo = barrier <= LOW_BARRIER
hi = (barrier >= HIGH_BARRIER[0]) & (barrier <= HIGH_BARRIER[1])
for m, lab in ((lo, "low $c$"), (hi, "high $c$")):
    a_.axvspan(barrier[m].min(), barrier[m].max(), color="0.9", zorder=0)
    a_.text(barrier[m].min(), 0.97, f" {lab}", transform=a_.get_xaxis_transform(),
            fontsize=9, color="0.4", va="top")
good = (tab["escape"]["n"] >= 100)
gm = np.exp(np.log((tab["escape"]["mean"] / tab["escape"]["s_mean"])[good]).mean())
a_.axhline(gm, color=C_ROUSE, lw=0.9, ls=":")
a_.text(barrier.max(), gm, f"  $\\times${gm:.2f}", color=C_ROUSE, fontsize=9, va="center", ha="right")
a_.set_yscale("log")
a_.set_ylim(0.2, 35)
a_.set_xlabel(r"barrier $\Delta U / k_B T$", fontsize=12)
a_.set_ylabel("Rouse / 1D reduction", fontsize=11)
a_.set_title("(c)  reduction error vs barrier", fontsize=12, loc="left")
a_.legend(fontsize=8, loc="lower left")
a_.grid(alpha=0.25)

# (d)-(f) survival functions, each phi_inf rescaled by its own solver MFPT
def survival_panel(a_, key, groups, title):
    for m, col, lab in groups:
        pooled = np.concatenate([tab[key]["t"][i] / tab[key]["s_mean"][i]
                                 for i in np.flatnonzero(m) if tab[key]["n"][i]])
        pooled = np.sort(pooled[pooled > 0])
        a_.step(pooled, 1.0 - np.arange(pooled.size) / pooled.size, where="post",
                color=col, lw=1.6, label=f"{lab} (n={pooled.size})")
        for i in np.flatnonzero(m):
            t, S = solver_survival(recs[i], key)
            a_.plot(t / tab[key]["s_mean"][i], S, "--" if len(groups) > 1 else "-",
                    color=C_1D if len(groups) == 1 else col, lw=0.9, alpha=0.6,
                    zorder=0)
    a_.plot([], [], "--" if len(groups) > 1 else "-", color=C_1D if len(groups) == 1 else "0.5",
            lw=1.2, label="1D reduction (per $\\phi_\\infty$)")
    a_.set_xscale("log"); a_.set_yscale("log"); a_.set_ylim(8e-3, 1.6)
    a_.set_xlabel(r"$t\,/\,\langle t\rangle_{\rm 1D}$", fontsize=12)
    a_.set_ylabel(r"survival  $P(T>t)$", fontsize=11)
    a_.set_title(title, fontsize=12, loc="left")
    a_.legend(fontsize=8, loc="lower left")
    a_.grid(alpha=0.25)

b_lo, b_hi = barrier[lo], barrier[hi]
survival_panel(ax[1][0], "escape", [(lo, C_ROUSE, "Rouse")],
               f"(d)  EXIT, LOW $c$  ($\\Delta U$ = {b_lo.min():.1f}-{b_lo.max():.1f} $k_BT$)")
survival_panel(ax[1][1], "escape", [(hi, C_ROUSE, "Rouse")],
               f"(e)  EXIT, HIGH $c$  ($\\Delta U$ = {b_hi.min():.1f}-{b_hi.max():.1f} $k_BT$)")
survival_panel(ax[1][2], "return", [(lo, C_ENTRY, "Rouse, low $c$"),
                                    (hi, "#7b4fa8", "Rouse, high $c$")],
               "(f)  ENTRY, both groups")

fig.tight_layout(rect=[0, 0, 1, 0.965])
os.makedirs(OUT_DIR, exist_ok=True)
for ext in ("pdf", "png"):
    out = os.path.join(OUT_DIR, f"rouse_vs_1d_exit_entry.{ext}")
    fig.savefig(out, format=ext, bbox_inches="tight", dpi=150 if ext == "png" else None)
    print("saved", out)

# ── the numbers behind the figure ────────────────────────────────────────────
print(f"\nwell-sampled phi_inf (n >= 100): {good.sum()} of {len(phis)}")
for key in ("escape", "return"):
    T = tab[key]
    rm = np.log((T["mean"] / T["s_mean"])[good]); rq = np.log((T["med"] / T["s_med"])[good])
    print(f"  {key:6s}  mean ratio {np.exp(rm.mean()):5.2f} (x{np.exp(rm.std()):.2f})"
          f"   median ratio {np.exp(rq.mean()):5.2f} (x{np.exp(rq.std()):.2f})"
          f"   => D*/D_eff = {np.exp(-rm.mean()):.3f}")
K_r = tab["escape"]["mean"] / tab["return"]["mean"]
K_s = tab["escape"]["s_mean"] / tab["return"]["s_mean"]
print(f"  equilibrium constant K = T_exit/T_entry:  Rouse/1D = "
      f"{np.exp(np.log((K_r / K_s)[good]).mean()):.2f} "
      f"(x{np.exp(np.log((K_r / K_s)[good]).std()):.2f})  <- kinetics shift, thermodynamics does not")
