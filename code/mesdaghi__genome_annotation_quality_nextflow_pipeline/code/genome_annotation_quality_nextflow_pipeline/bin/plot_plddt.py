import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import sys
import glob
import pandas as pd
import os

from sklearn.mixture import GaussianMixture
from adjustText import adjust_text

USAGE = """Usage:
  Single input:
    python plot_plddt.py <dataset_name> [metapredict_csv]
      (reads plddt_all_values_<dataset_name>*.pkl from the current directory)

  Several inputs on shared plots:
    python plot_plddt.py --combined <output_prefix> \\
        --query NAME PKL [METAPREDICT_CSV] \\
        --query NAME PKL [METAPREDICT_CSV] ...
"""

# ---------- ARGUMENTS ----------
# combined_mode: every --query is drawn on the same figures, labelled by
# NAME, and the outputs are named <output_prefix>_<plot>.png.
# Otherwise the original single-input behaviour is kept unchanged.
combined_mode = len(sys.argv) > 1 and sys.argv[1] == "--combined"

if combined_mode:
    argv = sys.argv[2:]
    if not argv or argv[0].startswith("--"):
        print(USAGE)
        sys.exit(1)
    output_prefix = argv.pop(0)
    queries = []  # list of (name, pkl_path, metapredict_csv_or_None)
    while argv:
        flag = argv.pop(0)
        if flag != "--query":
            print(f"ERROR: unexpected argument {flag!r}\n{USAGE}")
            sys.exit(1)
        vals = []
        while argv and argv[0] != "--query":
            vals.append(argv.pop(0))
        if len(vals) not in (2, 3):
            print(f"ERROR: --query takes NAME PKL [METAPREDICT_CSV], got {vals}\n{USAGE}")
            sys.exit(1)
        queries.append((vals[0], vals[1], vals[2] if len(vals) == 3 else None))
    if not queries:
        print(f"ERROR: --combined needs at least one --query\n{USAGE}")
        sys.exit(1)
    names = [q[0] for q in queries]
    if len(set(names)) != len(names):
        print(f"ERROR: --query names must be unique, got {names}")
        sys.exit(1)
    dataset_name = output_prefix
    # Disorder for the coloured scatter is only drawn if every query has it
    query_metapredict = {name: csv for name, _, csv in queries}
    metapredict_csv = "per-query" if all(query_metapredict.values()) else None
    plot_title_suffix = f"{len(queries)} inputs" if len(queries) > 1 else queries[0][0]
else:
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print(USAGE)
        sys.exit(1)
    dataset_name = sys.argv[1]
    metapredict_csv = sys.argv[2] if len(sys.argv) == 3 else None
    plot_title_suffix = dataset_name


def output_png(kind):
    """kind: e.g. 'plddt_density_scipy'."""
    if combined_mode:
        return f"{output_prefix}_{kind}.png"
    return f"{kind}_{dataset_name}.png"

# ---------- CSV FILE ----------
# Resolve CSV path robustly. Search order:
#   1. $PLDDT_CSV environment variable (if set)
#   2. ./plddt_model_organisms.csv               (current working dir)
#   3. ./reference/plddt_model_organisms.csv     (common project layout)
#   4. <script_dir>/plddt_model_organisms.csv    (legacy behaviour)
#   5. <script_dir>/reference/plddt_model_organisms.csv
_script_dir = os.path.dirname(os.path.abspath(__file__))
_csv_candidates = [
    os.environ.get("PLDDT_CSV"),
    "plddt_model_organisms.csv",
    os.path.join("reference", "plddt_model_organisms.csv"),
    os.path.join(_script_dir, "plddt_model_organisms.csv"),
    os.path.join(_script_dir, "reference", "plddt_model_organisms.csv"),
    # Nextflow-style layout: script in bin/, CSV in sibling reference/
    os.path.join(_script_dir, "..", "reference", "plddt_model_organisms.csv"),
]
csv_file = next((p for p in _csv_candidates if p and os.path.exists(p)), None)
if csv_file is None:
    print("ERROR: plddt_model_organisms.csv not found. Looked in:")
    for p in _csv_candidates:
        if p:
            print(f"  - {p}")
    print("Set $PLDDT_CSV to override, or place the file next to the script.")
    sys.exit(1)
print(f"Using CSV: {csv_file}")

# ---------- DISORDER PKL (for disorder-coloured scatter) ----------
# Same search logic as the pLDDT CSV — works whether run standalone,
# from the project root, or as a Nextflow process.
_disorder_candidates = [
    os.environ.get("DISORDER_PKL"),
    "combined_proteome_disorder.pkl",
    os.path.join("reference", "combined_proteome_disorder.pkl"),
    os.path.join(_script_dir, "combined_proteome_disorder.pkl"),
    os.path.join(_script_dir, "reference", "combined_proteome_disorder.pkl"),
    os.path.join(_script_dir, "..", "reference", "combined_proteome_disorder.pkl"),
]
disorder_pkl_path = next((p for p in _disorder_candidates if p and os.path.exists(p)), None)
if disorder_pkl_path is None:
    print("NOTE: combined_proteome_disorder.pkl not found in any of:")
    for p in _disorder_candidates:
        if p:
            print(f"  - {p}")
    print("Disorder-coloured scatter will be skipped.")
else:
    print(f"Using disorder PKL: {disorder_pkl_path}")

# Cache file — bumped to v2 because schema changed (Skew column removed).
# Lives in cwd; in Nextflow each task gets a fresh work dir so the cache
# is effectively single-run, but the CSV GMM fits are cheap so this is fine.
cache_file = "cached_csv_gmm_metrics_model_organisms_v2.csv"

# Species to exclude from CSV
species_to_exclude = ["Cauris6684", "Homo_sapiens_2k", "afdb_Homo_sapiens_2k", "CaurisB8441", "Pan_troglodytes", "Rattus_norvegicus"]

# ---------- TAXONOMY ----------
# Grouping: kingdom-level, but protists are split by supergroup since
# Apicomplexa (Toxoplasma, Plasmodium) and Euglenozoa (Trypanosoma) are
# not meaningfully related and tend to behave differently.
SPECIES_TAXONOMY = {
    "Homo_sapiens":             "Metazoa",
    "Mus_musculus":             "Metazoa",
    "Drosophila_melanogaster":  "Metazoa",
    "Saccharomyces_cerevisiae": "Fungi",
    "Arabidopsis_thaliana":     "Plants",
    "Toxoplasma_gondii":        "Apicomplexa",
    "Plasmodium_falciparum":    "Apicomplexa",
    "Trypanosoma_brucei":       "Euglenozoa",
}

# Fixed colour mapping so colours stay stable across runs and datasets
TAXON_COLOURS = {
    "Metazoa":     "#1f77b4",  # blue
    "Fungi":       "#ff7f0e",  # orange
    "Plants":      "#2ca02c",  # green
    "Apicomplexa": "#d62728",  # red
    "Euglenozoa":  "#9467bd",  # purple
    "Unknown":     "#7f7f7f",  # grey fallback
}

def get_taxon(species):
    return SPECIES_TAXONOMY.get(species, "Unknown")

# ---------- LOAD CSV (all species except excluded, HS renamed) ----------
df_csv = pd.read_csv(csv_file)
df_csv["Species"] = df_csv["Species"].replace({"HS": "Homo_sapiens"})
df_csv = df_csv[~df_csv["Species"].isin(species_to_exclude)]
csv_species_values = (
    df_csv.groupby("Species")["Mean_pLDDT"]
    .apply(list)
    .to_dict()
)
print("Loaded CSV species:", list(csv_species_values.keys()))
if not csv_species_values:
    print("ERROR: CSV loaded but contains no species after filtering. "
          "Check species_to_exclude and CSV contents.")
    sys.exit(1)

# ---------- LOAD PKL ----------
if combined_mode:
    # One entry per query, keyed by its display name. A PKL normally holds
    # a single dataset; if it holds several they are pooled under NAME.
    species_plddt = {}
    for name, pkl_path, _ in queries:
        if not os.path.exists(pkl_path):
            print(f"ERROR: file not found: {pkl_path}")
            sys.exit(1)
        with open(pkl_path, "rb") as f:
            loaded = pickle.load(f)
        pooled = {}
        for sub in loaded.values():
            pooled.update(sub)
        species_plddt[name] = pooled
else:
    pkl_files = glob.glob(f"plddt_all_values_{dataset_name}*.pkl")

    if len(pkl_files) == 0:
        print(f"ERROR: file not found: plddt_all_values_{dataset_name}*.pkl")
        sys.exit(1)

    pkl_file = pkl_files[0]

    with open(pkl_file, "rb") as f:
        species_plddt = pickle.load(f)

print("Loaded PKL species:", list(species_plddt.keys()))

def scipy_kde_curve(values, grid):
    kde = gaussian_kde(values)
    return kde.evaluate(grid)

plt.rcParams["figure.figsize"] = (10, 6)

# ---------- CONVERT PKL DICT ----------
species_values = {}
for species, protein_dict in species_plddt.items():
    species_values[species] = list(protein_dict.values())

# =========================================================
# GUIDED GMM INITIALISATION from reference species
# =========================================================
reference_species = [
    "Drosophila_melanogaster",
    "Homo_sapiens",
    "Saccharomyces_cerevisiae",
    "Mus_musculus",
    "Arabidopsis_thaliana"
]

print("Fitting reference species to compute initial GMM parameters...")

gmm_means_list = []
gmm_cov_list = []
gmm_weights_list = []

for sp in reference_species:
    if sp in species_values:
        values = np.array(species_values[sp])
    elif sp in csv_species_values:
        values = np.array(csv_species_values[sp])
    else:
        print(f"  Reference species not found, skipping: {sp}")
        continue

    if len(values) < 20:
        continue

    X = values.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, random_state=0).fit(X)

    means = gmm.means_.flatten()
    order = np.argsort(means)

    gmm_means_list.append(means[order])
    gmm_cov_list.append(gmm.covariances_.flatten()[order])
    gmm_weights_list.append(gmm.weights_[order])

if len(gmm_means_list) == 0:
    print("WARNING: No reference species found — falling back to default GMM init.")
    reference_means = None
    reference_weights = None
else:
    reference_means = np.mean(gmm_means_list, axis=0)
    reference_weights = np.mean(gmm_weights_list, axis=0)
    print(f"  Reference means:   {reference_means}")
    print(f"  Reference weights: {reference_weights}")

print("Reference GMM parameters computed.")

def fit_gmm(values):
    """Fit 2-component GMM with guided init if available, else default."""
    X = values.reshape(-1, 1)
    if reference_means is not None:
        gmm = GaussianMixture(
            n_components=2,
            means_init=reference_means.reshape(-1, 1),
            weights_init=reference_weights,
            random_state=0
        ).fit(X)
    else:
        gmm = GaussianMixture(n_components=2, random_state=0).fit(X)
    return gmm

# ---------- CREATE KDE GRID ----------
all_values_list = [vals for vals in species_values.values() if len(vals) > 0]
if len(csv_species_values) > 0:
    all_values_list.extend([vals for vals in csv_species_values.values() if len(vals) > 0])

all_values = np.concatenate(all_values_list)
grid = np.linspace(min(all_values), max(all_values), 400)

# =========================================================
# FIGURE — SciPy KDE density plot
# =========================================================
plt.figure()

for species, values in species_values.items():
    if len(values) < 2:
        continue
    y = scipy_kde_curve(values, grid)
    plt.plot(grid, y, label=species)

for species, values in csv_species_values.items():
    if len(values) < 2:
        continue
    y = scipy_kde_curve(values, grid)
    plt.plot(grid, y, linestyle="--", alpha=0.3, label=f"{species} (CSV)")

plt.title(f"pLDDT Density (SciPy KDE) - {plot_title_suffix}")
plt.xlabel("Mean pLDDT")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig(output_png("plddt_density_scipy"), dpi=300)
plt.close()

# =========================================================
# GMM SCATTER
# =========================================================
metrics = []

# ---------------------------
# LOAD OR COMPUTE CSV METRICS
# ---------------------------
if os.path.exists(cache_file):
    print("Loading cached CSV GMM metrics...")
    csv_metrics_df = pd.read_csv(cache_file)
    csv_metrics_df["Species"] = csv_metrics_df["Species"].replace({"HS": "Homo_sapiens"})
    csv_metrics_df = csv_metrics_df[~csv_metrics_df["Species"].isin(species_to_exclude)]
else:
    print("Computing CSV GMM metrics (first run)...")
    csv_metrics = []

    for sp, values in csv_species_values.items():
        values = np.array(values)
        if len(values) < 20:
            continue
        gmm = fit_gmm(values)
        means = gmm.means_.flatten()
        order = np.argsort(means)
        upper_prop = gmm.weights_[order][1]
        csv_metrics.append({
            "Species": sp,
            "GMM_upper_prop": upper_prop,
            "Prop_ge_70": np.mean(values >= 70)
        })

    csv_metrics_df = pd.DataFrame(csv_metrics)
    csv_metrics_df.to_csv(cache_file, index=False)
    print("CSV metrics cached.")

for _, row in csv_metrics_df.iterrows():
    metrics.append({
        "Species": row["Species"],
        "GMM_upper_prop": row["GMM_upper_prop"],
        "Prop_ge_70": row["Prop_ge_70"],
        "Is_PKL": False
    })

# ---------------------------
# PKL METRICS (always fresh)
# ---------------------------
for sp, values in species_values.items():
    values = np.array(values)
    if len(values) < 20:
        continue
    gmm = fit_gmm(values)
    means = gmm.means_.flatten()
    order = np.argsort(means)
    upper_prop = gmm.weights_[order][1]
    metrics.append({
        "Species": sp,
        "GMM_upper_prop": upper_prop,
        "Prop_ge_70": np.mean(values >= 70),
        "Is_PKL": True
    })

df_metrics = pd.DataFrame(metrics)
df_metrics["Taxon"] = df_metrics["Species"].map(get_taxon)
df_metrics["Colour"] = df_metrics["Taxon"].map(TAXON_COLOURS)

# ---------------------------
# SCATTER — coloured by taxon, PKL query highlighted
# ---------------------------
fig, ax = plt.subplots(figsize=(9, 7))

# Plot one scatter per taxon so we get a clean categorical legend
for taxon, group in df_metrics.groupby("Taxon"):
    ax.scatter(
        group["GMM_upper_prop"],
        group["Prop_ge_70"],
        c=TAXON_COLOURS.get(taxon, TAXON_COLOURS["Unknown"]),
        label=taxon,
        alpha=0.85,
        edgecolors="black",
        linewidths=0.3,
        s=60,
        zorder=2
    )

ax.legend(title="Taxon", loc="best", frameon=True)

# Expand axis limits by 10% on each side to prevent clipping
x_vals = df_metrics["GMM_upper_prop"]
y_vals = df_metrics["Prop_ge_70"]
x_margin = (x_vals.max() - x_vals.min()) * 0.10
y_margin = (y_vals.max() - y_vals.min()) * 0.10
ax.set_xlim(x_vals.min() - x_margin, x_vals.max() + x_margin)
ax.set_ylim(y_vals.min() - y_margin, y_vals.max() + y_margin)

# Build text objects for adjustText — CSV species
texts = []
for _, row in df_metrics[~df_metrics["Is_PKL"]].iterrows():
    t = ax.text(
        row["GMM_upper_prop"],
        row["Prop_ge_70"],
        row["Species"],
        fontsize=8
    )
    texts.append(t)

# Overlay PKL query species with red ring — no legend entry
pkl_points = df_metrics[df_metrics["Is_PKL"]]

if len(pkl_points) > 0:
    ax.scatter(
        pkl_points["GMM_upper_prop"],
        pkl_points["Prop_ge_70"],
        facecolors="none",
        edgecolors="red",
        s=200,
        linewidths=2,
        zorder=3
    )
    for _, row in pkl_points.iterrows():
        t = ax.text(
            row["GMM_upper_prop"],
            row["Prop_ge_70"],
            row["Species"],
            fontsize=9,
            fontweight="bold",
            color="red"
        )
        texts.append(t)

# Auto-adjust all labels to avoid overlap
adjust_text(
    texts,
    ax=ax,
    expand=(1.3, 1.5),
    arrowprops=dict(arrowstyle="-", color="grey", lw=0.5)
)

ax.set_xlabel("GMM Upper Component Proportion")
ax.set_ylabel("Proportion ≥70 pLDDT")
ax.set_title(f"GMM Upper Peak vs ≥70 Proportion ({plot_title_suffix})")
plt.tight_layout()
plt.savefig(output_png("plddt_gmm_scatter"), dpi=300)
plt.close()

print("GMM scatter diagnostics generated (CSV cached, PKL highlighted).")


# =========================================================
# DISORDER-COLOURED SCATTER — same axes/points as the taxon scatter,
# but coloured by proteome-wide mean disorder.
# Skipped (with a warning) if either the disorder reference PKL or
# the query metapredict CSV is unavailable.
# =========================================================

def _compute_query_mean_disorder(csv_path):
    """Per-protein mean disorder from a metapredict CSV (fields 2+ per row),
    then averaged across proteins. Matches plot_metapredict.py's parsing."""
    per_protein = []
    with open(csv_path) as fh:
        for line in fh:
            parts = line.strip().split(",")
            try:
                scores = [float(x) for x in parts[2:] if x.strip()]
            except ValueError:
                continue
            if scores:
                per_protein.append(float(np.mean(scores)))
    if not per_protein:
        return None
    return float(np.mean(per_protein))


can_make_disorder_plot = disorder_pkl_path is not None and metapredict_csv is not None
if not can_make_disorder_plot:
    reasons = []
    if disorder_pkl_path is None:
        reasons.append("no disorder reference PKL")
    if metapredict_csv is None:
        reasons.append("no metapredict CSV passed for the query")
    print(f"Skipping disorder-coloured scatter ({', '.join(reasons)}).")
else:
    print("Building disorder-coloured scatter...")

    # Reference disorder: per-species mean of mean_disorder
    df_disorder = pd.read_pickle(disorder_pkl_path)
    df_disorder["species"] = df_disorder["species"].replace({"HS": "Homo_sapiens"})
    ref_disorder = (
        df_disorder.groupby("species", observed=True)["mean_disorder"]
        .mean()
        .to_dict()
    )

    # Query disorder: one value per query species (keyed as in species_values)
    def _load_query_disorder(csv_path):
        if not os.path.exists(csv_path):
            print(f"  WARNING: metapredict CSV not found: {csv_path}")
            return None
        value = _compute_query_mean_disorder(csv_path)
        if value is None:
            print(f"  WARNING: no usable rows in {csv_path}")
        return value

    if combined_mode:
        query_disorder = {name: _load_query_disorder(csv)
                          for name, csv in query_metapredict.items()}
    else:
        single = _load_query_disorder(metapredict_csv)
        query_disorder = {sp: single for sp in species_values}
    for name, value in query_disorder.items():
        if value is not None:
            print(f"  Query mean disorder ({name}): {value:.4f}")

    # Attach a MeanDisorder value to each row of df_metrics
    def _disorder_for_row(row):
        if row["Is_PKL"]:
            value = query_disorder.get(row["Species"])
            return value if value is not None else np.nan
        return ref_disorder.get(row["Species"], np.nan)

    df_metrics["MeanDisorder"] = df_metrics.apply(_disorder_for_row, axis=1)

    valid_mask = df_metrics["MeanDisorder"].notna()
    n_missing = int((~valid_mask).sum())
    if n_missing > 0:
        missing_species = df_metrics.loc[~valid_mask, "Species"].tolist()
        print(f"  {n_missing} point(s) without disorder data (plotted in grey): {missing_species}")

    fig, ax = plt.subplots(figsize=(9, 7))

    # Points with disorder data — coloured by viridis
    pts_valid = df_metrics[valid_mask]
    sc = ax.scatter(
        pts_valid["GMM_upper_prop"],
        pts_valid["Prop_ge_70"],
        c=pts_valid["MeanDisorder"],
        cmap="viridis",
        alpha=0.9,
        edgecolors="black",
        linewidths=0.3,
        s=70,
        zorder=2,
    )
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Mean disorder (proteome-wide)")

    # Points without disorder data — grey, with a small legend entry
    pts_missing = df_metrics[~valid_mask]
    if len(pts_missing) > 0:
        ax.scatter(
            pts_missing["GMM_upper_prop"],
            pts_missing["Prop_ge_70"],
            c="lightgrey",
            alpha=0.85,
            edgecolors="black",
            linewidths=0.3,
            s=70,
            zorder=2,
            label="No disorder data",
        )
        ax.legend(loc="best", frameon=True, fontsize=9)

    # Same axis limits as the taxon scatter (recompute, since axis() not stored)
    x_vals = df_metrics["GMM_upper_prop"]
    y_vals = df_metrics["Prop_ge_70"]
    x_margin = (x_vals.max() - x_vals.min()) * 0.10
    y_margin = (y_vals.max() - y_vals.min()) * 0.10
    ax.set_xlim(x_vals.min() - x_margin, x_vals.max() + x_margin)
    ax.set_ylim(y_vals.min() - y_margin, y_vals.max() + y_margin)

    # CSV species labels
    disorder_texts = []
    for _, row in df_metrics[~df_metrics["Is_PKL"]].iterrows():
        disorder_texts.append(ax.text(
            row["GMM_upper_prop"], row["Prop_ge_70"], row["Species"],
            fontsize=8
        ))

    # Query species — red ring + bold label
    pkl_points = df_metrics[df_metrics["Is_PKL"]]
    if len(pkl_points) > 0:
        ax.scatter(
            pkl_points["GMM_upper_prop"],
            pkl_points["Prop_ge_70"],
            facecolors="none",
            edgecolors="red",
            s=200,
            linewidths=2,
            zorder=3,
        )
        for _, row in pkl_points.iterrows():
            disorder_texts.append(ax.text(
                row["GMM_upper_prop"], row["Prop_ge_70"], row["Species"],
                fontsize=9, fontweight="bold", color="red"
            ))

    adjust_text(
        disorder_texts,
        ax=ax,
        expand=(1.3, 1.5),
        arrowprops=dict(arrowstyle="-", color="grey", lw=0.5),
    )

    ax.set_xlabel("GMM Upper Component Proportion")
    ax.set_ylabel("Proportion ≥70 pLDDT")
    ax.set_title(f"GMM Scatter — coloured by Mean Disorder ({plot_title_suffix})")
    plt.tight_layout()
    plt.savefig(output_png("plddt_gmm_scatter_disorder"), dpi=300)
    plt.close()
    print("Disorder-coloured scatter generated.")

