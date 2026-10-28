"""Figures (PNG, 200 dpi, white background) from the tables written by 08_compare.py."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *

OURS, REF, SP = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e3e3e0"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
                     "legend.frameon": False, "font.family": "DejaVu Sans"})
S = json.loads((TAB / "summary.json").read_text())
N = S["n_clean"]


def save(fig, name):
    fig.savefig(FIG / name, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ 1 class / architecture composition
cls = pd.read_csv(TAB / "class_composition.csv")
arch = pd.read_csv(TAB / "architecture_composition.csv", dtype={"arch": str})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1, 1.35]})
x = np.arange(len(cls)); w = 0.38
a1.bar(x - w / 2 - 0.01, cls.ours_share, w, color=OURS, label=f"SAXS set ({S['composition']['n_units']} protein-superfamily units)")
a1.errorbar(x - w / 2 - 0.01, cls.ours_share, yerr=[cls.ours_share - cls.ours_lo, cls.ours_hi - cls.ours_share],
            fmt="none", ecolor=INK, capsize=3, lw=1)
a1.bar(x + w / 2 + 0.01, cls.cath_share, w, color=REF, label=f"CATH {S['cath_reference']['version']} S35 clusters (n={S['composition']['cath_s35_clusters_classes_1_4']:,})")
for i, r in cls.iterrows():
    a1.text(i, max(r.ours_hi, r.cath_share) + 0.015, f"×{r.enrichment:.2f}", ha="center", color=INK2, fontsize=9)
a1.set_xticks(x, ["Mainly\nalpha", "Mainly\nbeta", "Alpha\nbeta", "Few secondary\nstructures"])
a1.set_ylabel("Share of domains"); a1.set_ylim(0, 0.62); a1.grid(axis="x", visible=False)
a1.legend(loc="upper left", fontsize=8.5)
a1.set_title("a  CATH class composition (enrichment = ours / CATH)", loc="left", fontsize=11)
ar = arch.sort_values("enrichment")
y = np.arange(len(ar))
a2.axvline(1, color=INK2, lw=1, ls="--")
lo = np.clip(ar.enr_lo, 0.08, None)
a2.errorbar(ar.enrichment.clip(lower=0.08), y, xerr=[ar.enrichment.clip(lower=0.08) - lo, ar.enr_hi - ar.enrichment.clip(lower=0.08)],
            fmt="o", color=OURS, ecolor=OURS, ms=6, capsize=2, lw=1.2)
a2.set_xscale("log"); a2.set_xlim(0.08, 8)
a2.set_xticks([0.1, 0.25, 0.5, 1, 2, 4, 8], ["0.1", "0.25", "0.5", "1", "2", "4", "8"])
a2.set_yticks(y, [f"{r.arch} {r['name']}  (ours {r.ours_units}, CATH {r.cath_s35:,})" for _, r in ar.iterrows()], fontsize=8.5)
a2.set_xlabel("Enrichment ratio, ours / CATH S35 (log scale, 95% bootstrap CI)")
a2.set_title("b  CATH architectures (CATH share ≥1% or ≥3 of ours)", loc="left", fontsize=11)
fig.tight_layout()
save(fig, "fig1_class_composition.png")

# ------------------------------------------------------------------ 2 rank-abundance
rk = pd.read_csv(TAB / "cath_superfamily_rank.csv")
fig, (b1, b2) = plt.subplots(2, 1, figsize=(10, 7.5), gridspec_kw={"height_ratios": [2.2, 1]}, sharex=True)
b1.plot(rk["rank"], rk.s35_clusters, color="#9a9994", lw=1.5, label="CATH superfamilies (all)")
pres = rk[rk.in_ours_gene3d]
b1.scatter(pres["rank"], pres.s35_clusters, s=22, color=OURS, edgecolor="white", linewidth=0.6, zorder=3,
           label=f"present in SAXS set (CATH-Gene3D, {len(pres)} superfamilies)")
b1.set_xscale("log"); b1.set_yscale("log")
b1.set_ylabel("S35 sequence clusters per superfamily")
for k in (50, 200):
    b1.axvline(k, color=INK2, lw=0.8, ls=":")
    c = S["coverage"][f"top{k}_superfamilies_present_gene3d"]
    b1.text(k * 1.06, 12 if k == 50 else 4, f"top {k}:\n{c['k']}/{k} present", fontsize=8.5, color=INK2, va="top")
for (_, r), dy in zip(pres.sort_values("rank").head(3).iterrows(), (2, 8, -14)):
    b1.annotate(f"{r['name'][:40]} ({r.sf})", (r["rank"], r.s35_clusters), xytext=(8, dy), textcoords="offset points",
                fontsize=7.5, color=INK2)
b1.legend(loc="lower left", fontsize=9)
b1.set_title(f"a  CATH {S['cath_reference']['version']} superfamilies ranked by S35 clusters (n={len(rk):,})", loc="left", fontsize=11)
edges = [0, 10, 50, 200, 500, 1000, 2000, len(rk)]
shares, labs, mids = [], [], []
for lo_, hi_ in zip(edges[:-1], edges[1:]):
    sub = rk[(rk["rank"] > lo_) & (rk["rank"] <= hi_)]
    shares.append(sub.in_ours_gene3d.mean()); labs.append(f"{lo_ + 1}–{hi_}"); mids.append((lo_ + 1, hi_))
for (l_, h_), s_ in zip(mids, shares):
    b2.fill_between([l_, h_ + 1], 0, s_, step="pre", color=OURS, alpha=0.85, edgecolor="white", linewidth=1.5)
    b2.text(np.sqrt(l_ * (h_ + 1)), s_ + 0.03, f"{s_:.0%}", ha="center", fontsize=8.5, color=INK)
b2.set_ylim(0, 1.05); b2.set_ylabel("Share present")
b2.set_xlabel("Superfamily rank by S35 clusters (log scale)")
b2.set_title("b  Share of CATH superfamilies present, by rank bin", loc="left", fontsize=11)
fig.tight_layout()
save(fig, "fig2_superfamily_rank_abundance.png")

# ------------------------------------------------------------------ 3 organisms
org = pd.read_csv(TAB / "organism_comparison.csv", index_col=0)
fig, ax = plt.subplots(figsize=(10, 5.2))
x = np.arange(len(org)); w = 0.27
ax.bar(x - w - 0.02, org.ours_share, w, color=OURS, label=f"SAXS set (n={N})")
ax.errorbar(x - w - 0.02, org.ours_share, yerr=[org.ours_share - org.ours_lo, org.ours_hi - org.ours_share], fmt="none", ecolor=INK, capsize=2.5, lw=1)
ax.bar(x, org.pdb_nr30_share, w, color=REF, label=f"PDB, 30% identity clusters, ≥30 aa (n={S['organism']['n_pdb_nr30_ge30']:,})")
ax.bar(x + w + 0.02, org.swissprot_share, w, color=SP, label=f"UniProtKB/Swiss-Prot (n={S['organism']['n_swissprot']:,})")
for i, r in enumerate(org.itertuples()):
    ax.text(i - w - 0.02, r.ours_hi + 0.01, f"{r.ours}", ha="center", fontsize=8, color=INK2)
ax.set_xticks(x, org.index); ax.set_ylabel("Share of proteins"); ax.grid(axis="x", visible=False)
ax.legend(loc="upper right", fontsize=9)
ax.set_title("Source organism: SAXS set vs non-redundant PDB vs Swiss-Prot (numbers = SAXS entries; bars ±95% Wilson CI)", loc="left", fontsize=10.5)
fig.tight_layout()
save(fig, "fig3_organisms.png")

# ------------------------------------------------------------------ 4 lengths and domains per protein
df = pd.read_csv(TAB / "per_protein_mapping.csv")
C = df[df.clean.astype(str) == "True"]
pdb = pd.read_csv(TAB / "rcsb_nr30_representatives.csv.gz")
sp = pd.read_csv(TAB / "swissprot_summary.csv.gz")
fig, (c1, c2) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1.4, 1]})
bins = np.logspace(np.log10(10), np.log10(5000), 34)
for vals, col, lab in [(C.construct_length, OURS, f"SAXS constructs (n={N}, median {int(C.construct_length.median())})"),
                       (pdb.length[pdb.length >= 30], REF, f"PDB 30% clusters ≥30 aa (median {int(pdb.length[pdb.length >= 30].median())})"),
                       (sp.length, SP, f"Swiss-Prot (median {int(sp.length.median())})")]:
    h, e = np.histogram(vals, bins=bins)
    c1.step(e[:-1], h / h.sum(), where="post", color=col, lw=2, label=lab)
c1.set_xscale("log"); c1.set_xlabel("Chain length (residues, log scale)"); c1.set_ylabel("Fraction per log-length bin")
c1.legend(fontsize=8.5, loc="upper left")
c1.set_title("a  Chain length distribution (SAXS constructs are all ≤500 aa)", loc="left", fontsize=11)
dpp = pd.read_csv(TAB / "domains_per_protein.csv", index_col=0)
x = np.arange(len(dpp)); w = 0.27
c2.bar(x - w - 0.02, dpp.ours_construct, w, color=OURS, label="SAXS construct region")
c2.bar(x, dpp.ours_parent_fulllength, w, color="#9cc3ee", label="SAXS parent protein, full length")
c2.bar(x + w + 0.02, dpp.swissprot, w, color=SP, label="Swiss-Prot, full length")
c2.set_xticks(x, dpp.index.astype(str)); c2.set_xlabel("CATH-Gene3D domains per protein"); c2.set_ylabel("Share of proteins")
c2.grid(axis="x", visible=False); c2.legend(fontsize=8.5)
c2.set_title("b  Domains per protein", loc="left", fontsize=11)
fig.tight_layout()
save(fig, "fig4_lengths_domains.png")

# ------------------------------------------------------------------ 5 disorder and structure availability
dis = pd.read_csv(TAB / "disorder_comparison.csv", index_col=0)
sh = dis / dis.sum()
cols = ["ours_construct", "ours_parent_fulllength", "swissprot_all", "swissprot_with_pdb"]
labels = ["SAXS\nconstructs", "SAXS parent\nproteins", "Swiss-Prot\nall", "Swiss-Prot with\nPDB entry"]
fig, (d1, d2) = plt.subplots(1, 2, figsize=(13, 5))
x = np.arange(len(cols)); w = 0.36
d1.bar(x - w / 2 - 0.01, sh.loc["mostly disordered (>=50% residues)", cols], w, color="#1c5aa6", label="≥50% of residues disordered")
d1.bar(x + w / 2 + 0.01, sh.loc["long IDR (>=30 aa), <50%", cols], w, color="#8fb9ea", label="IDR ≥30 aa, <50% disordered")
for i, c in enumerate(cols):
    for off, k in [(-w / 2 - 0.01, "mostly disordered (>=50% residues)"), (w / 2 + 0.01, "long IDR (>=30 aa), <50%")]:
        d1.text(i + off, sh.loc[k, c] + 0.008, f"{sh.loc[k, c]:.1%}", ha="center", fontsize=8, color=INK2)
d1.set_xticks(x, labels); d1.set_ylabel("Share of proteins"); d1.grid(axis="x", visible=False); d1.legend(fontsize=9)
d1.set_title("a  UniProt 'Disordered' regions (mainly MobiDB-lite)", loc="left", fontsize=11)
sa = S["structure_availability"]["by_disorder_class_no_pdb_ge50"]
order = ["folded", "partly disordered", "disordered"]
v = [sa[k]["share"] for k in order]
lo_ = [sa[k]["share"] - sa[k]["ci95"][0] for k in order]; hi_ = [sa[k]["ci95"][1] - sa[k]["share"] for k in order]
d2.bar(range(3), v, 0.55, color=OURS)
d2.errorbar(range(3), v, yerr=[lo_, hi_], fmt="none", ecolor=INK, capsize=3, lw=1)
for i, k in enumerate(order):
    d2.text(i, sa[k]["ci95"][1] + 0.02, f"{sa[k]['k']}/{sa[k]['n']}", ha="center", fontsize=9, color=INK2)
d2.set_xticks(range(3), ["Folded", "Partly disordered", "Disordered"]); d2.set_ylim(0, 0.9)
d2.set_xlabel("Disorder class of the SAXS entry (ADOPT2)")
d2.set_ylabel("Share without PDB coverage of ≥50% of construct")
d2.grid(axis="x", visible=False)
d2.set_title("b  Constructs lacking an experimental PDB structure", loc="left", fontsize=11)
fig.tight_layout()
save(fig, "fig5_disorder_structure.png")
print("figures written to", FIG)
