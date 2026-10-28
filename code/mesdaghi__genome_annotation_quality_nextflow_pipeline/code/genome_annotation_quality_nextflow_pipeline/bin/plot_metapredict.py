import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pandas as pd
import sys
import os

# =========================================================
# CONFIG
# =========================================================
REFERENCE_PKL = os.path.join(os.path.dirname(__file__), "..", "reference", "combined_proteome_disorder.pkl")

# Species set must match the pLDDT plot (plot_plddt.py).
# These are the species we *want* to show.
SPECIES_TO_INCLUDE = {
    "Homo_sapiens",
    "Mus_musculus",
    "Drosophila_melanogaster",
    "Saccharomyces_cerevisiae",
    "Arabidopsis_thaliana",
    "Toxoplasma_gondii",
    "Plasmodium_falciparum",
    "Trypanosoma_brucei",
}

# Species explicitly dropped (kept in sync with plot_plddt.py for clarity)
SPECIES_TO_EXCLUDE = {
    "Cauris6684", "Homo_sapiens_2k", "afdb_Homo_sapiens_2k",
    "CaurisB8441", "Pan_troglodytes", "Rattus_norvegicus",
}

# =========================================================
# ARGUMENTS
# =========================================================
# Single input (original behaviour):
#     python plot_metapredict.py <metapredict_csv>
#       -> <csv basename>_mean_disorder_density_scipy.png
# Several inputs on one plot:
#     python plot_metapredict.py A.csv B.csv ... [--names A B ...] [--out-prefix combined]
#       -> <out-prefix>_metapredict_mean_disorder_density_scipy.png
import argparse

ap = argparse.ArgumentParser(description="Overlay query disorder distributions on the reference proteomes.")
ap.add_argument("csvs", nargs="+", help="Query metapredict CSV(s)")
ap.add_argument("--names", nargs="+", default=None,
                help="Display names, one per CSV (default: CSV basenames)")
ap.add_argument("--out-prefix", default=None,
                help="Write <prefix>_metapredict_mean_disorder_density_scipy.png "
                     "(default for a single CSV: <csv basename>_mean_disorder_density_scipy.png)")
args = ap.parse_args()

if args.names is not None and len(args.names) != len(args.csvs):
    sys.exit(f"ERROR: got {len(args.names)} --names for {len(args.csvs)} CSVs")

query_labels = args.names or [os.path.basename(c).replace(".csv", "") for c in args.csvs]
if len(set(query_labels)) != len(query_labels):
    sys.exit(f"ERROR: query names must be unique, got {query_labels}")

if args.out_prefix:
    output_png = f"{args.out_prefix}_metapredict_mean_disorder_density_scipy.png"
elif len(args.csvs) == 1:
    output_png = f"{query_labels[0]}_mean_disorder_density_scipy.png"
else:
    sys.exit("ERROR: --out-prefix is required with more than one CSV")

print("Datasets:", query_labels)

# =========================================================
# LOAD REFERENCE PKL
# =========================================================
if not os.path.exists(REFERENCE_PKL):
    raise FileNotFoundError(f"Reference PKL not found: {REFERENCE_PKL}")

df = pd.read_pickle(REFERENCE_PKL)
print("Loaded reference PKL")

# Normalise any aliases (mirrors HS -> Homo_sapiens rename in plot_plddt.py)
df["species"] = df["species"].replace({"HS": "Homo_sapiens"})

available_species = set(df["species"].unique())
print("Species found in PKL:", sorted(available_species))

# Apply exclusion list first (drops the unwanted extras like Pan_troglodytes etc.)
df = df[~df["species"].isin(SPECIES_TO_EXCLUDE)]

# Then restrict to the target set
df = df[df["species"].isin(SPECIES_TO_INCLUDE)]

kept_species = sorted(df["species"].unique())
print("Species kept after filtering:", kept_species)

# Warn about anything in the target set with no data in the PKL
missing = sorted(SPECIES_TO_INCLUDE - set(kept_species))
if missing:
    print("WARNING: the following target species are NOT in the disorder PKL "
          "and will be absent from the plot:")
    for sp in missing:
        print(f"  - {sp}")

# =========================================================
# LOAD HIGHLIGHT CSV(S)
# =========================================================
def load_mean_disorder(csv_path):
    values = []
    if not os.path.exists(csv_path):
        print("WARNING: highlight CSV not found:", csv_path)
        return values
    print("Loading highlight CSV:", csv_path)
    with open(csv_path) as f:
        for line in f:
            parts = line.strip().split(",")
            try:
                scores = [float(x) for x in parts[2:] if x]
            except:
                continue
            if scores:
                values.append(np.mean(scores))
    return values

highlights = [(label, load_mean_disorder(csv)) for label, csv in zip(query_labels, args.csvs)]

# =========================================================
# PREP DATA
# =========================================================
species_mean_values = {}
for species, subdf in df.groupby("species"):
    species_mean_values[species] = subdf["mean_disorder"].to_numpy()

# =========================================================
# STYLE FUNCTIONS
# =========================================================
def style_background():
    return dict(linestyle="--", alpha=0.3, linewidth=1.5)

def style_highlight():
    return dict(linestyle="-", alpha=1.0, linewidth=2.5)

# =========================================================
# SCIPY KDE
# =========================================================
plt.figure()
grid = np.linspace(0, 1, 400)

for species, values in species_mean_values.items():
    if len(values) < 2:
        continue
    y = gaussian_kde(values)(grid)
    plt.plot(grid, y, label=species, **style_background())

for label, values in highlights:
    if len(values) > 1:
        y = gaussian_kde(values)(grid)
        plt.plot(grid, y, label=label, **style_highlight())

plt.title("Mean Disorder Density (SciPy KDE)")
plt.xlabel("Mean Disorder")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig(output_png, dpi=300)
plt.close()

print("Done — overlay plots generated")
