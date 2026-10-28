"""
Load output_cpu*.json files from a prewetting_diagram sweep and plot two
heatmaps over a chosen pair of axes from (R, omega_b, phi_inf):
  - left:  excess grand potential (surface free energy)
  - right: phi_0 (surface volume fraction)

The third axis is collapsed by selecting the slice closest to SLICE_VALUE
(or the unique value if there is only one).

Output: media/figures/<date>/prewetting_phase_diagram.png
"""

import sys, os, glob, json
from datetime import date

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Config ───────────────────────────────────────────────────────────────────
EXP_NAME = "prewetting_diagram_R_3"
Y_AXIS    = "omega_b"
X_AXIS    = "phi_inf"
SLICE_AT  = None                       # value for the third axis; None => median
# OUT_DATE  = date.today().strftime("%Y%m%d")
OUT_DATE = "20260706"
DATA_DIR  = f"data/{OUT_DATE}/{EXP_NAME}"
OUT_PATH  = f"media/figures/{OUT_DATE}/{EXP_NAME}.png"

assert X_AXIS in ("R", "omega_b", "phi_inf")
assert Y_AXIS in ("R", "omega_b", "phi_inf") and Y_AXIS != X_AXIS
SLICE_KEY = ({"R", "omega_b", "phi_inf"} - {X_AXIS, Y_AXIS}).pop()


# ── Load ─────────────────────────────────────────────────────────────────────
files = sorted(glob.glob(os.path.join(DATA_DIR, "output_cpu*.json")))
if not files:
    sys.exit(f"no output_cpu*.json under {DATA_DIR}")

records = []
for fp in files:
    with open(fp) as f:
        blob = json.load(f)
    for _, by_b in blob.items():
        for _, by_c in by_b.items():
            for _, rec in by_c.items():
                records.append(rec)

print(f"loaded {len(records)} records from {len(files)} files")


# ── Slice the third axis ─────────────────────────────────────────────────────
slice_vals = np.array(sorted({r[SLICE_KEY] for r in records}))
if SLICE_AT is None:
    slice_val = float(slice_vals[len(slice_vals) // 2])
else:
    slice_val = float(slice_vals[np.argmin(np.abs(slice_vals - SLICE_AT))])

eps = 1e-9 * max(1.0, abs(slice_val))
sliced = [r for r in records if abs(r[SLICE_KEY] - slice_val) <= eps]
print(f"sliced at {SLICE_KEY}={slice_val:g}  ->  {len(sliced)} records")


# ── Bin onto (x, y) grid ─────────────────────────────────────────────────────
xs = np.array(sorted({r[X_AXIS] for r in sliced}))
ys = np.array(sorted({r[Y_AXIS] for r in sliced}))

F   = np.full((len(ys), len(xs)), np.nan)
P0  = np.full((len(ys), len(xs)), np.nan)

x_idx = {v: i for i, v in enumerate(xs)}
y_idx = {v: i for i, v in enumerate(ys)}

for r in sliced:
    i = y_idx[r[Y_AXIS]]
    j = x_idx[r[X_AXIS]]
    if r.get("success", True):
        F[i, j]  = r["excess_grand_potential"]
        P0[i, j] = r["phi_0"]


# ── Plot ─────────────────────────────────────────────────────────────────────
def _extent(a):
    if len(a) == 1:
        return (a[0] - 0.5, a[0] + 0.5)
    d = (a[-1] - a[0]) / (len(a) - 1) / 2
    return (a[0] - d, a[-1] + d)

xlo, xhi = _extent(xs)
ylo, yhi = _extent(ys)

axis_label = {"R": r"$R$", "omega_b": r"$h_0$", "phi_inf": r"$\phi_\infty$"}

fig, axes = plt.subplots(1, 1, figsize=(5, 3.5))

im0 = axes.imshow(P0, origin="lower", aspect="auto",
                     extent=[xlo, xhi, ylo, yhi], cmap="viridis")
axes.set_xlabel(axis_label[X_AXIS], fontsize=12)
axes.set_ylabel(axis_label[Y_AXIS], fontsize=12)
axes.set_title(r"Surface composition $\phi_0$", fontsize=10)
fig.colorbar(im0, ax=axes, label=r"$\phi(0)$")
# axes.axvline(x=0.009, color="red")
# axes.scatter([0.009, 0.009], [-1.24, -1.28], color="red", s=50)

# im1 = axes[1].imshow(P0, origin="lower", aspect="auto",
#                      extent=[xlo, xhi, ylo, yhi], cmap="magma")
# axes[1].set_xlabel(axis_label[X_AXIS], fontsize=12)
# axes[1].set_ylabel(axis_label[Y_AXIS], fontsize=12)
# axes[1].set_title(r"Surface composition $\phi_0$", fontsize=13)
# fig.colorbar(im1, ax=axes[1], label=r"$\phi(0)$")

fig.suptitle(
    f"{axis_label[SLICE_KEY]} = {slice_val:g}, " 
        + r"$\chi = 4$, "
        + r"Binodal at $\phi_\infty=$" + str(round(r["phi_bin_low"], 3)),
    fontsize=13,
)
fig.tight_layout()

os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
fig.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
print(f"saved -> {OUT_PATH}")
