
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

# Input files
size_file = Path("size_statistics_20.csv")
selected_dir = Path("GA001/curve_1/pdbs")

# Read pool statistics
rg, ree, names = [], [], []

with open(size_file) as f:
    for line in f:
        p = line.split()
        if len(p) < 6:
            continue

        rg.append(float(p[0]))
        ree.append(float(p[2]))
        names.append(Path(p[5]).stem)

rg = np.array(rg)
ree = np.array(ree)

# Identify structures selected by GAJOE
selected_names = {
    p.stem for p in selected_dir.glob("*.cif")
}

mask = np.array([
    name in selected_names for name in names
])

print(f"Pool: {len(rg)}")
print(f"Selected by GAJOE: {mask.sum()}")
print("Selected structures:", sorted(selected_names))

assert mask.sum() == len(selected_names), \
    "Some selected structures were not found in the pool!"

# Figure
fig = plt.figure(figsize=(7, 6))

gs = GridSpec(
    2, 2,
    width_ratios=[4, 1.3],
    height_ratios=[1.3, 4],
    hspace=0.05,
    wspace=0.05
)

ax = fig.add_subplot(gs[1, 0])
ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
ax_right = fig.add_subplot(gs[1, 1], sharey=ax)

# 2D scatter
ax.scatter(
    ree, rg,
    s=20, color="coral", alpha=0.7,
    label="Full pool"
)

ax.scatter(
    ree[mask], rg[mask],
    s=65, facecolors="none",
    edgecolors="dodgerblue",
    linewidths=1.5,
    label="Selected by EOM"
)

# Histograms
bins_ree = np.linspace(0, max(ree) * 1.1, 12)
bins_rg = np.linspace(0, max(rg) * 1.1, 12)

ax_top.hist(
    ree, bins=bins_ree, density=True,
    color="coral", alpha=0.7
)
ax_top.hist(
    ree[mask], bins=bins_ree, density=True,
    color="dodgerblue", alpha=0.5
)

ax_right.hist(
    rg, bins=bins_rg, density=True,
    orientation="horizontal",
    color="coral", alpha=0.7
)
ax_right.hist(
    rg[mask], bins=bins_rg, density=True,
    orientation="horizontal",
    color="dodgerblue", alpha=0.5
)

# Formatting
ax.set_xlabel(r"$R_{ee}$ ($\AA$)", fontsize=14)
ax.set_ylabel(r"$R_g$ ($\AA$)", fontsize=14)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

ax_top.tick_params(labelbottom=False)
ax_right.tick_params(labelleft=False)

plt.savefig("EOM_2D.png", dpi=300, bbox_inches="tight")
plt.show()

