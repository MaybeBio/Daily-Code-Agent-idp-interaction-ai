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
if len(sys.argv) < 2:
    print("Usage: python plot_metapredict.py <metapredict_csv>")
    sys.exit(1)

highlight_csv = sys.argv[1]
dataset_label = os.path.basename(highlight_csv).replace(".csv", "")

print("Dataset:", dataset_label)

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
# LOAD HIGHLIGHT CSV
# =========================================================
highlight_mean_values = []

if os.path.exists(highlight_csv):
    print("Loading highlight CSV:", highlight_csv)

    with open(highlight_csv) as f:
        for line in f:
            parts = line.strip().split(",")
            try:
                scores = [float(x) for x in parts[2:] if x]
            except:
                continue
            if scores:
                highlight_mean_values.append(np.mean(scores))
else:
    print("WARNING: highlight CSV not found")

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

if highlight_mean_values and len(highlight_mean_values) > 1:
    y = gaussian_kde(highlight_mean_values)(grid)
    plt.plot(grid, y, label=dataset_label, **style_highlight())

plt.title("Mean Disorder Density (SciPy KDE)")
plt.xlabel("Mean Disorder")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig(f"{dataset_label}_mean_disorder_density_scipy.png", dpi=300)
plt.close()

print("Done — overlay plots generated")
