#!/usr/bin/env python3
"""
Cluster the SAE decoder vectors into feature families.

The compositional MobiDB triplets are too sparse to resolve families by F1, so the
grouping is taken from decoder geometry instead. Uses the SAE checkpoint, the feature
atlas, and the descriptions; no retraining and no GPU.

Usage:
    python scripts/09_feature_clusters.py paths=colab model=esm2_8m sae=topk_32x
    python scripts/09_feature_clusters.py paths=colab model=esm2_8m sae=topk_32x +cluster_min_f1=0.3
"""
import logging
import sys
from pathlib import Path

import hydra
import numpy as np
import pandas as pd
import torch
from omegaconf import DictConfig

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from idpinterp.sae.architecture import build_sae

logger = logging.getLogger(__name__)


def _run_name(cfg):
    layer = cfg.model.get("layer") or cfg.model.layers_to_extract[-1]
    tag = cfg.sae.get("tag", None)
    return f"{cfg.model.name}_{cfg.sae.type}_{cfg.sae.expansion_factor}x{'_'+tag if tag else ''}_layer{layer}"


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s")
    run = _run_name(cfg)
    atlas_dir = Path(cfg.paths.output_dir) / "atlas" / run
    out_dir = Path(cfg.paths.output_dir) / "diagnostics" / run
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Decoder vectors from the checkpoint (rows of W_dec = feature dictionaries)
    ckpt = Path(cfg.paths.checkpoint_dir) / f"{run}_final.pt"
    if not ckpt.exists():
        logger.error(f"No checkpoint: {ckpt}"); return
    sae = build_sae(cfg.sae, cfg.model)
    sae.load_state_dict(torch.load(ckpt, map_location="cpu", weights_only=True)["model_state_dict"])
    W_dec = sae.W_dec.detach().cpu().numpy()                      # (d_latent, d_in)
    W_dec = W_dec / (np.linalg.norm(W_dec, axis=1, keepdims=True) + 1e-8)

    # ── Restrict to interpretable features (aligned in the atlas) ──────────────
    atlas = pd.read_parquet(atlas_dir / "feature_atlas.parquet")
    min_f1 = float(cfg.get("cluster_min_f1", 0.3))
    feats = atlas[atlas["f1"] >= min_f1].copy()
    idx = feats["feature_idx"].to_numpy()
    X = W_dec[idx]
    logger.info(f"Clustering {len(idx)} features (F1>={min_f1}) by decoder geometry.")

    # ── Embeddings: moderate-dim for CLUSTERING, 2D for PLOTTING ───────────────
    # (Clustering on the 2D projection collapses near-duplicate disorder features
    #  into one blob → degenerate 1–2 clusters; cluster on higher-dim instead.)
    seed = int(cfg.seed)
    n_clu_dim = min(int(cfg.get("cluster_dims", 10)), len(idx) - 2, X.shape[1] - 1)
    try:
        import umap
        Xc = umap.UMAP(n_components=n_clu_dim, metric="cosine", n_neighbors=15,
                       min_dist=0.0, random_state=seed).fit_transform(X)
        XY = umap.UMAP(n_components=2, metric="cosine", n_neighbors=15,
                       min_dist=0.1, random_state=seed).fit_transform(X)
        emb = "UMAP"
    except Exception as e:
        from sklearn.decomposition import PCA
        Xc = PCA(n_components=min(n_clu_dim, X.shape[1]), random_state=seed).fit_transform(X)
        XY = PCA(n_components=2, random_state=seed).fit_transform(X)
        emb = f"PCA (umap unavailable: {e})"
    logger.info(f"Embedding: {emb} | clustering in {Xc.shape[1]}d, plotting in 2d")

    # ── Cluster on the moderate-dim embedding; KMeans fallback if degenerate ───
    min_size = int(cfg.get("cluster_min_size", 6))
    try:
        from sklearn.cluster import HDBSCAN
        labels = HDBSCAN(min_cluster_size=min_size, min_samples=1,
                         metric="euclidean").fit_predict(Xc)
        nclu = len(set(labels)) - (1 if -1 in labels else 0)
        clu = f"HDBSCAN"
        if nclu < 5:
            raise ValueError(f"HDBSCAN degenerate ({nclu} clusters)")
    except Exception as e:
        from sklearn.cluster import KMeans
        k = min(int(cfg.get("cluster_k", 25)), len(idx))
        labels = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(Xc)
        clu = f"KMeans({k}) [{e}]"
    feats["cluster"] = labels
    feats["umap_x"], feats["umap_y"] = XY[:, 0], XY[:, 1]
    nclu = len(set(labels)) - (1 if -1 in labels else 0)
    logger.info(f"Clustering: {clu} → {nclu} families")

    # ── Attach LLM motif snippets if present ───────────────────────────────────
    desc_path = atlas_dir / "feature_descriptions.parquet"
    if desc_path.exists():
        dd = pd.read_parquet(desc_path)[["feature_idx", "description"]]
        dd["motif"] = dd["description"].str.split("MOTIF").str[-1].str[:80].str.replace("\n", " ")
        feats = feats.merge(dd[["feature_idx", "motif"]], on="feature_idx", how="left")

    feats.to_csv(out_dir / "feature_clusters.csv", index=False)

    # ── Per-cluster summary: dominant triplet + a representative motif ─────────
    summ = []
    for c, g in feats.groupby("cluster"):
        if c == -1:
            continue
        top_trip = g["triplet"].mode().iat[0] if len(g) else "?"
        motif = g.sort_values("f1", ascending=False)["motif"].dropna().iat[0] if "motif" in g and g["motif"].notna().any() else ""
        summ.append({"cluster": int(c), "n": len(g), "dominant_triplet": top_trip,
                     "max_f1": round(g["f1"].max(), 3), "example_motif": motif})
    summ_df = pd.DataFrame(summ).sort_values("n", ascending=False)
    summ_df.to_csv(out_dir / "feature_cluster_summary.csv", index=False)
    logger.info("Cluster summary:\n" + summ_df.to_string(index=False))

    # ── Plot ───────────────────────────────────────────────────────────────────
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(11, 9))
        sc = ax.scatter(XY[:, 0], XY[:, 1], c=labels, cmap="tab20", s=14, alpha=0.8)
        ax.set_title(f"SAE decoder-vector clusters — disorder feature families ({run})")
        ax.set_xlabel("dim 1"); ax.set_ylabel("dim 2")
        fig.colorbar(sc, ax=ax, label="cluster")
        fig.tight_layout(); fig.savefig(out_dir / "feature_clusters.png", dpi=150)
        logger.info(f"Saved → {out_dir / 'feature_clusters.png'}")
    except Exception as e:
        logger.warning(f"Plot skipped: {e}")


if __name__ == "__main__":
    main()
