from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = Path("data/depmap")
RAW = Path("data/raw")
SUPP = Path("results/depmap_essentiality/supp")
SUPP.mkdir(parents=True, exist_ok=True)
TMP = Path(tempfile.gettempdir()) / "idrfeat_mmseqs_tmp"
SEED = 1729
MIN_ID = 0.30
THREADS = "3"
rng = np.random.default_rng(SEED)


def strip_entrez(name: str) -> str:
    return re.sub(r"\s*\(\d+\)$", "", str(name)).strip()


def load_sequences(symbols: set[str]) -> dict[str, str]:
    best: dict[str, str] = {}
    sym = None
    seq: list[str] = []
    with gzip.open(RAW / "UP000005640_reviewed.fasta.gz", "rt") as f:
        for line in f:
            if line.startswith(">"):
                if sym and sym in symbols:
                    s = "".join(seq)
                    if sym not in best or len(s) > len(best[sym]):
                        best[sym] = s
                m = re.search(r"GN=(\S+)", line)
                sym = m.group(1) if m else None
                seq = []
            else:
                seq.append(line.strip())
    if sym and sym in symbols:
        s = "".join(seq)
        if sym not in best or len(s) > len(best[sym]):
            best[sym] = s
    return best


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 50:
        return np.nan
    x, y = a[m], b[m]
    if x.std() == 0 or y.std() == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def main() -> None:
    gl = pd.read_csv(DATA / "gene_level_idr_features.csv", index_col=0)
    eff = pd.read_csv(DATA / "mean_gene_effect.csv", index_col=0)["mean_gene_effect"]
    genes = set(gl.index.astype(str))

    print("Loading proteome sequences...")
    seqs = load_sequences(genes)
    print(f"  mapped sequences for {len(seqs)}/{len(genes)} feature genes")

    fasta = SUPP / "genes.fasta"
    with open(fasta, "w") as f:
        for g, s in seqs.items():
            f.write(f">{g}\n{s}\n")

    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)

    print(f"Clustering with mmseqs (min-seq-id {MIN_ID})...")
    cl_prefix = str(SUPP / "clust")
    subprocess.run(["mmseqs", "easy-cluster", str(fasta), cl_prefix, str(TMP / "c"),
                    "--min-seq-id", str(MIN_ID), "-c", "0.5", "--threads", THREADS, "-v", "1"],
                   check=True)
    cl = pd.read_csv(cl_prefix + "_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
    gene2cluster = dict(zip(cl["member"], cl["rep"]))
    pd.DataFrame({"gene": list(gene2cluster), "cluster": list(gene2cluster.values())}).to_csv(
        SUPP / "gene_clusters.csv", index=False)
    sizes = cl.groupby("rep").size().sort_values(ascending=False)
    n_singleton = int((sizes == 1).sum())
    print(f"  clusters {len(sizes)} | multi-gene {int((sizes > 1).sum())} | singletons {n_singleton}")

    ras = ["KRAS", "HRAS", "NRAS"]
    ras_rows = []
    for g in ras:
        ras_rows.append({"gene": g, "in_feature_set": g in genes,
                         "mean_gene_effect": round(float(eff[g]), 3) if g in eff.index else np.nan})
    pd.DataFrame(ras_rows).to_csv(SUPP / "ras_check.csv", index=False)
    print("\nRAS family:")
    for r in ras_rows:
        print(f"  {r['gene']}: in IDR feature set = {r['in_feature_set']}, "
              f"pan-cancer mean effect = {r['mean_gene_effect']}")
    print("\nLargest sequence-similarity clusters among feature genes:")
    for rep, n in sizes.head(8).items():
        members = cl.loc[cl["rep"] == rep, "member"].tolist()
        print(f"  {rep} (n={n}): {', '.join(members[:8])}{' ...' if n > 8 else ''}")

    print("\nAll-vs-all homology search...")
    m8 = SUPP / "allvall.m8"
    subprocess.run(["mmseqs", "easy-search", str(fasta), str(fasta), str(m8), str(TMP / "s"),
                    "--threads", THREADS, "-s", "5.7", "--max-seqs", "300",
                    "--format-output", "query,target,fident,alnlen,evalue,bits", "-v", "1"],
                   check=True)
    hits = pd.read_csv(m8, sep="\t", header=None,
                       names=["q", "t", "fident", "alnlen", "evalue", "bits"])
    hits = hits[(hits["q"] != hits["t"]) & (hits["alnlen"] >= 50) & (hits["fident"] >= MIN_ID)]
    a = np.minimum(hits["q"], hits["t"])
    b = np.maximum(hits["q"], hits["t"])
    hits = hits.assign(ga=a, gb=b).sort_values("fident", ascending=False).drop_duplicates(["ga", "gb"])
    print(f"  homologous gene pairs (id>={MIN_ID}, alnlen>=50): {len(hits)}")
    if len(hits) > 150000:
        hits = hits.sample(150000, random_state=SEED)
        print(f"  sampled {len(hits)} pairs for the diagnostic")

    with open(DATA / "CRISPRGeneEffect.csv") as f:
        header = f.readline().rstrip("\n").split(",")
    colmap: dict[str, str] = {}
    name2idx: dict[str, int] = {}
    for i, col in enumerate(header):
        name2idx[col] = i
    sym2idx: dict[str, int] = {}
    for col in header[1:]:
        sym2idx.setdefault(strip_entrez(col), name2idx[col])
    colmap = {g: g for g in sym2idx}

    involved = set(hits["ga"]) | set(hits["gb"])
    involved = {g for g in involved if g in sym2idx}
    usecols = [0] + sorted(sym2idx[g] for g in involved)
    prof = pd.read_csv(DATA / "CRISPRGeneEffect.csv", index_col=0, usecols=usecols)
    prof.columns = [strip_entrez(c) for c in prof.columns]
    P = {g: prof[g].to_numpy(dtype=float) for g in prof.columns}

    pair_rows = []
    for ga, gb, fid in zip(hits["ga"], hits["gb"], hits["fident"]):
        if ga in P and gb in P:
            r = pearson(P[ga], P[gb])
            if r == r:
                pair_rows.append((ga, gb, float(fid), r))
    pairs_df = pd.DataFrame(pair_rows, columns=["gene_a", "gene_b", "identity", "essentiality_corr"])
    pairs_df.to_csv(SUPP / "homology_vs_essentiality.csv", index=False)

    glist = list(P.keys())
    homo_set = {(ga, gb) for ga, gb in zip(hits["ga"], hits["gb"])}
    bg = []
    tries = 0
    while len(bg) < min(5000, len(pairs_df)) and tries < 200000:
        i, j = rng.integers(0, len(glist), 2)
        tries += 1
        if i == j:
            continue
        ga, gb = sorted((glist[i], glist[j]))
        if (ga, gb) in homo_set:
            continue
        r = pearson(P[ga], P[gb])
        if r == r:
            bg.append(r)
    bg = np.array(bg)

    from scipy.stats import spearmanr
    rho = spearmanr(pairs_df["identity"], pairs_df["essentiality_corr"]).statistic
    bins = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01]
    pairs_df["bin"] = pd.cut(pairs_df["identity"], bins, right=False)
    binned = pairs_df.groupby("bin", observed=True)["essentiality_corr"].agg(["count", "mean", "median"])
    binned.to_csv(SUPP / "homology_vs_essentiality_binned.csv")

    print("\nIdentity vs essentiality-profile diagnostic:")
    print(f"  usable homologous pairs with profiles: {len(pairs_df)}")
    print(f"  Spearman(sequence identity, essentiality-profile correlation) = {rho:.3f}")
    print(f"  mean essentiality corr, homologous pairs  = {pairs_df['essentiality_corr'].mean():.3f}")
    print(f"  mean essentiality corr, random background = {bg.mean():.3f} (n={len(bg)})")
    print("  by identity bin:")
    print(binned.round(3).to_string())

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].scatter(pairs_df["identity"], pairs_df["essentiality_corr"], s=6, alpha=0.25, color="#0f766e")
    ax[0].axhline(bg.mean(), color="#c0392b", lw=1.4, ls="--", label=f"random pairs mean ({bg.mean():.2f})")
    ax[0].set_xlabel("pairwise sequence identity")
    ax[0].set_ylabel("essentiality-profile correlation")
    ax[0].set_title(f"Homologous pairs (Spearman {rho:.2f})")
    ax[0].legend(fontsize=8)
    centers = [(bins[i] + bins[i + 1]) / 2 for i in range(len(binned))]
    ax[1].bar(centers, binned["mean"].to_numpy(), width=0.08, color="#1f6feb")
    ax[1].axhline(bg.mean(), color="#c0392b", lw=1.4, ls="--")
    ax[1].set_xlabel("sequence identity bin")
    ax[1].set_ylabel("mean essentiality correlation")
    ax[1].set_title("Mean essentiality correlation by identity")
    fig.tight_layout()
    fig.savefig(SUPP / "homology_vs_essentiality.png", dpi=170)
    fig.savefig(SUPP / "homology_vs_essentiality.pdf")

    shutil.rmtree(TMP, ignore_errors=True)
    print(f"\nwrote gene_clusters.csv, ras_check.csv, homology_vs_essentiality.* to {SUPP}")


if __name__ == "__main__":
    main()
