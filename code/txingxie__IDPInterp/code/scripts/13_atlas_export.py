#!/usr/bin/env python3
"""
Consolidate one run's atlas into a single publication-ready table + summary.

Joins, per feature:
  feature_atlas.parquet        (best axis, F1, precision, region_recall)
  feature_descriptions.parquet (LLM description, held-out pearson_r)   [optional]
  diagnostics/<run>/feature_clusters.csv (decoder family id)            [optional]

Outputs → outputs/atlas/<run>/:
  atlas_table.csv / .parquet   — one row per feature, sorted by F1
  atlas_summary.json           — SAE-quality-independent atlas stats:
                                 #features F1>0.5 / >0.3, median held-out r,
                                 #families, per-axis best F1 (from triplet CSVs)

Usage:
  python scripts/13_atlas_export.py paths=colab model=esm2_8m  model.layer=4  sae=topk_32x
  python scripts/13_atlas_export.py paths=colab model=seqdance model.layer=10 sae=topk_32x
"""
import glob
import json
import logging
import os
import sys
from pathlib import Path

import hydra
import numpy as np
import pandas as pd
from omegaconf import DictConfig

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
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
    atlas_path = atlas_dir / "feature_atlas.parquet"
    if not atlas_path.exists():
        logger.error(f"No atlas at {atlas_path}; run scripts/04 first.")
        return
    atlas = pd.read_parquet(atlas_path)

    # ── merge descriptions (optional) ──────────────────────────────────────
    desc_path = atlas_dir / "feature_descriptions.parquet"
    if desc_path.exists():
        d = pd.read_parquet(desc_path)
        keep = [c for c in ("feature_idx", "description", "pearson_r", "top_proteins")
                if c in d.columns]
        atlas = atlas.merge(d[keep], on="feature_idx", how="left")
    else:
        logger.warning("No feature_descriptions.parquet — table will lack descriptions.")

    # ── merge decoder family id (optional) ─────────────────────────────────
    clust_path = Path(cfg.paths.output_dir) / "diagnostics" / run / "feature_clusters.csv"
    if clust_path.exists():
        c = pd.read_csv(clust_path)
        cl_col = next((x for x in ("cluster", "family", "label") if x in c.columns), None)
        if cl_col:
            atlas = atlas.merge(
                c[["feature_idx", cl_col]].rename(columns={cl_col: "family"}),
                on="feature_idx", how="left")
    else:
        logger.warning("No feature_clusters.csv — run scripts/09 for family structure.")

    atlas = atlas.sort_values("f1", ascending=False)
    atlas.to_parquet(atlas_dir / "atlas_table.parquet", index=False)
    atlas.to_csv(atlas_dir / "atlas_table.csv", index=False)

    # ── per-axis best F1 from the triplet CSVs ─────────────────────────────
    per_axis = {}
    for f in glob.glob(str(atlas_dir / "triplets" / "*.csv")):
        axis = os.path.basename(f)[:-4]
        df = pd.read_csv(f)
        col = next((x for x in ("f1", "test_f1", "best_f1") if x in df.columns), None)
        if col:
            per_axis[axis] = round(float(df[col].fillna(0).max()), 4)

    r = atlas["pearson_r"].dropna() if "pearson_r" in atlas.columns else pd.Series([], dtype=float)
    r = r[np.isfinite(r)]
    summary = {
        "run": run,
        "n_features_total": int(len(atlas)),
        "n_F1_gt_0.5": int((atlas["f1"] > 0.5).sum()),
        "n_F1_gt_0.3": int((atlas["f1"] > 0.3).sum()),
        "n_described": int(atlas["description"].apply(
            lambda x: isinstance(x, str) and len(x) > 0 and not x.startswith("[API error")
        ).sum()) if "description" in atlas.columns else 0,
        "median_heldout_pearson_r": round(float(r.median()), 4) if len(r) else None,
        "n_features_validated": int(len(r)),
        "n_families": int(atlas["family"].nunique()) if "family" in atlas.columns else None,
        "per_axis_best_f1": dict(sorted(per_axis.items(), key=lambda kv: -kv[1])),
    }
    json.dump(summary, open(atlas_dir / "atlas_summary.json", "w"), indent=2)
    logger.info("Atlas summary:\n" + json.dumps(summary, indent=2))
    logger.info(f"Saved atlas_table + atlas_summary → {atlas_dir}")


if __name__ == "__main__":
    main()
