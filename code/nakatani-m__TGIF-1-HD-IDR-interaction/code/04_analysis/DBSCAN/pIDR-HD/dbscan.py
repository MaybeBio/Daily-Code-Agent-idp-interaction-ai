import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
from sklearn.cluster import DBSCAN
from collections import Counter
# =====================================
# parameters
# =====================================
mapping_file = "../../PCA/pca_traj_P.txt"
pc_cols = [1, 2]

eps = 0.013
min_samples = 650

out_png = f"dbscan_eps{eps:.3f}_ms{min_samples}.png"
# =====================================
# load data
# =====================================
rows = []
with open(mapping_file) as f:
    for line in f:
        if not line.strip():
             continue
        if line.lstrip().startswith(("#", "@")):
             continue

        rows.append(line.split())

data = np.array(rows, dtype=float)
frames = data[:, 0].astype(int)

X_raw = data[:, pc_cols]

# z-score normalization
X = (X_raw - X_raw.mean(axis=0)) / ( X_raw.std(axis=0) + 1e-12 )
# =====================================
# DBSCAN
# =====================================
db = DBSCAN(
    eps=eps,
    min_samples=min_samples,
    n_jobs=-1
)

labels = db.fit_predict(X)

# =====================================
# relabel by population
# largest cluster -> C1
# second largest -> C2
# ...
# noise remains -1
# =====================================

cnt = Counter(labels)

# noise removal
cluster_sizes = [
    (lab, size)
    for lab, size in cnt.items()
    if lab != -1
]

# sort by size
cluster_sizes.sort(
    key=lambda x: x[1],
    reverse=True
)

# old label -> new label
label_map = {
    old_lab: new_lab
    for new_lab, (old_lab, _) in enumerate(cluster_sizes, start=1)
}

# relabel
labels_re = np.array([
    -1 if lab == -1 else label_map[lab]
    for lab in labels
])

# =====================================
# statistics
# =====================================

n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
noise_ratio = np.mean(labels_re == -1)

print()
print(f"eps = {eps}")
print(f"min_samples = {min_samples}")
print(f"clusters = {n_clusters}")
print(f"noise ratio = {noise_ratio:.2%}")
print()


cnt = Counter(labels_re)

for k, v in sorted(
    cnt.items(),
    key=lambda x: x[1],
    reverse=True):

    print(
        f"label={k:3d} "
        f"count={v:8d} "
        f"fraction={v/len(labels):8.3%}"
    )

# =====================================
# plot
# =====================================

fm.fontManager.addfont("/home/biostr1/SOTU23/nak/fonts/Arial.ttf")
mpl.rcParams["font.family"] = "Arial"

plt.figure(figsize=(8, 8))
plot_noise = True

# noise first

if plot_noise:

    mask_noise = labels_re == -1

    plt.scatter(
        X_raw[mask_noise, 0],
        X_raw[mask_noise, 1],
        s=1,
        c="lightgray",
        alpha=0.3,
        rasterized=True,
        label="noise"
    )

# clusters
unique_labels = sorted(
    [x for x in set(labels_re) if x != -1]
    )

cmap = plt.cm.get_cmap("tab20", len(unique_labels))

plt.tick_params(direction="in", labelsize=18, pad=10, length=10, width=2)

plt.xlabel("PC1", fontsize=20, labelpad=12)
plt.ylabel("PC2", fontsize=20, labelpad=12)

plt.xlim(-2500, 500)
plt.ylim(-600, 600)

ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2)

total = len(labels)

for i, lab in enumerate(unique_labels):

    mask = labels_re == lab
    n = np.sum(mask)
    frac = 100 * n / total

    plt.scatter(
        X_raw[mask, 0],
        X_raw[mask, 1],
        s=2,
        color=cmap(i),
        alpha=0.8,
        rasterized=True,
        label=f"C{lab} ({frac:.1f}%)"
    )

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=14,
    markerscale=4
)

plt.tight_layout()

plt.savefig(
    out_png,
    dpi=300,
    bbox_inches="tight"
)

print()
print(f"saved -> {out_png}")
