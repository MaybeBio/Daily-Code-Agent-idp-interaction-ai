#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# STYLE
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 11

# =========================================================
# SPECIES CONFIG — kept in sync with plot_plddt.py and plot_metapredict.py
# =========================================================
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

SPECIES_TO_EXCLUDE = {
    "Cauris6684", "Homo_sapiens_2k", "afdb_Homo_sapiens_2k",
    "CaurisB8441", "Pan_troglodytes", "Rattus_norvegicus",
}

# The psauron CSV labels the three protists with database-style strings.
# Map them onto the canonical species names used by the other plots.
SPECIES_RENAMES = {
    "HS": "Homo_sapiens",
    "TgondiiME49_proteins_psauron_score":     "Toxoplasma_gondii",
    "Pfalciparum3D7_proteins_psauron_score":  "Plasmodium_falciparum",
    "TbruceiTREU927_proteins_psauron_score":  "Trypanosoma_brucei",
}

# ARGUMENTS
if len(sys.argv) < 2:
    print("Usage: python plot_psauron_distribution.py <combined_csv> [highlight_csv]")
    sys.exit(1)

csv_file = sys.argv[1]
highlight_csv = sys.argv[2] if len(sys.argv) >= 3 else None

print("=" * 60)
print("PSAURON CATEGORICAL DISTRIBUTION BY SPECIES")
print("=" * 60)

# LOAD CSV
df = pd.read_csv(csv_file)
print(f"✓ Loaded {len(df):,} rows from {csv_file}")

if "in-frame_score" not in df.columns:
    print("Column 'in-frame_score' not found.")
    sys.exit(1)

# Normalise species labels
df["species"] = df["species"].replace(SPECIES_RENAMES)

# Apply exclusions, then restrict to target set
available_species = set(df["species"].unique())
df = df[~df["species"].isin(SPECIES_TO_EXCLUDE)]
df = df[df["species"].isin(SPECIES_TO_INCLUDE)]

kept_species = sorted(df["species"].unique())
print(f"✓ Species kept after filtering: {kept_species}")

missing = sorted(SPECIES_TO_INCLUDE - set(kept_species))
if missing:
    print("⚠ WARNING: target species NOT present in the psauron CSV:")
    for sp in missing:
        print(f"    - {sp}")

print(f"✓ {len(df):,} rows remaining")

# OPTIONAL HIGHLIGHT DATASET
highlight_values = None
highlight_label = None

if highlight_csv and os.path.exists(highlight_csv):
    highlight_label = os.path.basename(highlight_csv).replace(".csv", "")
    tmp_df = pd.read_csv(highlight_csv, skiprows=2)

    if "in-frame_score" in tmp_df.columns:
        highlight_values = tmp_df["in-frame_score"].dropna().to_numpy()
        print(f"✓ Loaded highlight dataset: {highlight_label} (n={len(highlight_values):,})")
    else:
        print("⚠ Highlight CSV missing 'in-frame_score' column")


# CATEGORY CALCULATION
def calculate_category_percentages(scores):
    total = len(scores)
    if total == 0:
        return [0, 0, 0]

    low = (scores < 0.1).sum()
    moderate = ((scores >= 0.1) & (scores < 0.9)).sum()
    high = (scores >= 0.9).sum()

    return [
        (low / total) * 100,
        (moderate / total) * 100,
        (high / total) * 100,
    ]


# COMPUTE STATS PER SPECIES (preserve a stable order — same as SPECIES_TO_INCLUDE)
species_order = [sp for sp in [
    "Homo_sapiens",
    "Mus_musculus",
    "Drosophila_melanogaster",
    "Saccharomyces_cerevisiae",
    "Arabidopsis_thaliana",
    "Toxoplasma_gondii",
    "Plasmodium_falciparum",
    "Trypanosoma_brucei",
] if sp in kept_species]

species_stats = {}
species_counts = {}

for sp in species_order:
    values = df.loc[df["species"] == sp, "in-frame_score"].dropna().to_numpy()
    species_stats[sp] = calculate_category_percentages(values)
    species_counts[sp] = len(values)

# PLOTTING
categories = ['Low\n(<0.1)', 'Moderate\n(0.1-0.9)', 'High\n(≥0.9)']
x = np.arange(len(species_stats))
width = 0.25

fig, ax = plt.subplots(figsize=(14, 6))

colors = ['#d62728', '#ff7f0e', '#2ca02c']  # red, orange, green

# Plot bars per category
for i, category in enumerate(categories):
    values = [species_stats[sp][i] for sp in species_stats.keys()]
    ax.bar(x + (i - 1) * width,
           values,
           width,
           label=category,
           color=colors[i],
           alpha=0.8)

# Optional highlight overlay
if highlight_values is not None:
    highlight_stats = calculate_category_percentages(highlight_values)
    for i in range(3):
        ax.bar(len(x) + (i - 1) * width,
               highlight_stats[i],
               width,
               color=colors[i],
               edgecolor='black',
               linewidth=2.5,
               alpha=1.0)

    species_labels = list(species_stats.keys()) + [highlight_label]
    x = np.append(x, len(x))
else:
    species_labels = list(species_stats.keys())

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(species_labels, rotation=45, ha='right')
ax.set_ylabel("Percentage of Proteins (%)", fontweight='bold')
ax.set_title("Psauron In-frame Score Quality Categories by Species",
             fontweight='bold', fontsize=13)
ax.set_ylim(0, 100)
ax.legend(title="Quality Level",
          loc="upper center",
          bbox_to_anchor=(0.5, -0.35),
          ncol=3,
          frameon=False)
ax.grid(True, axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig("psauron_categorical_distribution.png", dpi=300, bbox_inches='tight')
plt.close()

print("Figure saved: psauron_categorical_distribution.png")


# PRINT SUMMARY TABLE
print("\nSUMMARY STATISTICS")
print("─" * 70)
print(f"{'Species':<28} {'Low (%)':<10} {'Moderate (%)':<14} {'High (%)':<10} {'n':<8}")
print("─" * 70)

for sp in species_stats:
    low, mod, high = species_stats[sp]
    n = species_counts[sp]
    print(f"{sp:<28} {low:<10.1f} {mod:<14.1f} {high:<10.1f} {n:<8}")

if highlight_values is not None:
    low, mod, high = highlight_stats
    print("─" * 70)
    print(f"{highlight_label:<28} {low:<10.1f} {mod:<14.1f} {high:<10.1f} {len(highlight_values):<8}")

print("\nDone.")
