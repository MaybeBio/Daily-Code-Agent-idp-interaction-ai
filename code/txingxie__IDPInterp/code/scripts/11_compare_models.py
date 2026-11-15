#!/usr/bin/env python3
"""
Cross-model SAE feature comparison by activation pattern.

The two decoders live in different residual spaces (ESM-2 320-d, SeqDance 480-d) and
share no basis, so features are compared by where they fire on the shared evaluation
corpus rather than geometrically.

Method:
  1. Align residues of shared proteins by (uniprot_id, position) through the per-model
     metadata; truncation may differ, so take the minimum length per protein.
  2. Sample S residues and z-score each feature's activations across them.
  3. Cross-correlation C = Zr^T Zc / S (ref_features x cmp_features).
  4. Best |corr| match per feature, split into shared (> thr) and model-unique (<= thr).
  5. Linear CKA between the two activation matrices as overall similarity.
  6. Annotate cmp features with the disorder atlas to identify disorder and sub-state
     features unique to the cmp model.

Outputs -> outputs/compare/<ref>__vs__<cmp>/:
  summary.json, cmp_feature_match.parquet, cmp_disorder_unique.csv

Usage:
  python scripts/11_compare_models.py paths=colab sae=topk_32x \
     +compare.ref_model=esm2_8m  +compare.ref_layer=4 \
     +compare.cmp_model=seqdance +compare.cmp_layer=10 \
     paths.embedding_dir=/mnt/local-scratch/embeddings
"""
import json
import logging
import sys
from pathlib import Path

import hydra
import numpy as np
import pandas as pd
from omegaconf import DictConfig

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from idpinterp.models.memmap_store import MetadataIndex

logger = logging.getLogger(__name__)

try:
    import torch
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False


def _run_name(model_name, sae_cfg, layer):
    tag = sae_cfg.get("tag", None)
    return f"{model_name}_{sae_cfg.type}_{sae_cfg.expansion_factor}x{'_'+tag if tag else ''}_layer{layer}"


def _load_side(cfg, model_name, layer, dataset):
    """Load one model's eval-corpus activation cache + metadata + atlas."""
    run = _run_name(model_name, cfg.sae, layer)
    scratch = cfg.paths.get("scratch_dir", None) or cfg.paths.embedding_dir
    act_path = Path(scratch) / "act_cache" / f"{run}.mmap"
    cfg_path = act_path.with_suffix(".json")
    if not (act_path.exists() and cfg_path.exists()):
        raise FileNotFoundError(
            f"Missing activation cache for {run}:\n  {act_path}\n"
            f"Rebuild it by running scripts/04_evaluate_atlas.py for "
            f"model={model_name} model.layer={layer} (the encoding pass writes this cache)."
        )
    ac = json.load(open(cfg_path))
    acts = np.memmap(act_path, dtype=np.float16, mode="r",
                     shape=(ac["total_residues"], ac["d_latent"]), order="C")
    meta = MetadataIndex.load(
        Path(cfg.paths.embedding_dir) / model_name / dataset / f"layer{layer}_metadata.parquet"
    )
    atlas_path = Path(cfg.paths.output_dir) / "atlas" / run / "feature_atlas.parquet"
    atlas = pd.read_parquet(atlas_path) if atlas_path.exists() else None
    return {"run": run, "acts": acts, "meta": meta, "d_latent": ac["d_latent"], "atlas": atlas}


def _aligned_indices(meta_ref, meta_cmp):
    """Row indices into each act-cache for residues of shared proteins (aligned)."""
    m = meta_ref.merge(meta_cmp, on="uniprot_id", suffixes=("_r", "_c"))
    ref_idx, cmp_idx = [], []
    for row in m.itertuples():
        n = int(min(row.stored_length_r, row.stored_length_c))
        ref_idx.append(np.arange(row.mmap_start_idx_r, row.mmap_start_idx_r + n))
        cmp_idx.append(np.arange(row.mmap_start_idx_c, row.mmap_start_idx_c + n))
    if not ref_idx:
        raise RuntimeError("No shared proteins between the two metadata indices.")
    return np.concatenate(ref_idx), np.concatenate(cmp_idx)


def _gather(acts, idx):
    """Read mmap rows for idx (sorted reads are far faster), return float32."""
    order = np.argsort(idx)
    out = np.empty((len(idx), acts.shape[1]), dtype=np.float32)
    out[order] = np.asarray(acts[idx[order]], dtype=np.float32)
    return out


def _standardize(mat):
    """z-score columns; zero-variance (dead/constant) features → all-zero column."""
    sd = mat.std(0, keepdims=True)
    dead = (sd < 1e-8).ravel()
    z = (mat - mat.mean(0, keepdims=True)) / np.where(sd < 1e-8, 1.0, sd)
    z[:, dead] = 0.0
    return z, dead


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s")
    if "compare" not in cfg:
        raise SystemExit("Provide +compare.ref_model/+compare.ref_layer/"
                         "+compare.cmp_model/+compare.cmp_layer overrides.")
    c = cfg.compare
    dataset = c.get("dataset", "gold_curated")

    ref = _load_side(cfg, c.ref_model, c.ref_layer, dataset)
    cmp = _load_side(cfg, c.cmp_model, c.cmp_layer, dataset)
    logger.info(f"REF {ref['run']} (F={ref['d_latent']})  vs  CMP {cmp['run']} (F={cmp['d_latent']})")

    ref_idx, cmp_idx = _aligned_indices(ref["meta"], cmp["meta"])
    logger.info(f"Aligned residues (shared proteins): {len(ref_idx):,}")

    n_sample = int(c.get("n_sample", 100_000))
    if len(ref_idx) > n_sample:
        sel = np.random.default_rng(0).choice(len(ref_idx), n_sample, replace=False)
        ref_idx, cmp_idx = ref_idx[sel], cmp_idx[sel]

    Zr, _ = _standardize(_gather(ref["acts"], ref_idx))
    Zc, _ = _standardize(_gather(cmp["acts"], cmp_idx))
    S = Zr.shape[0]
    logger.info(f"Computing {Zr.shape[1]}x{Zc.shape[1]} cross-correlation over {S:,} residues...")

    # Cross-correlation + self-Gram norms for linear CKA (GPU if available).
    if HAS_TORCH and torch.cuda.is_available():
        tr = torch.from_numpy(Zr).cuda()
        tc = torch.from_numpy(Zc).cuda()
        C = (tr.T @ tc) / S
        cross_fro2 = float((C * C).sum()) * S * S            # ||Zr^T Zc||_F^2
        denom = float((tr.T @ tr).norm() * (tc.T @ tc).norm())
        C = C.cpu().numpy()
        del tr, tc
        torch.cuda.empty_cache()
    else:
        C = (Zr.T @ Zc) / S
        cross_fro2 = float((C * C).sum()) * S * S
        denom = float(np.linalg.norm(Zr.T @ Zr) * np.linalg.norm(Zc.T @ Zc))
    cka = cross_fro2 / denom if denom > 0 else 0.0

    absC = np.abs(C)
    ref_best = absC.max(1)
    cmp_best = absC.max(0)
    cmp_match = absC.argmax(0)
    thr = float(c.get("match_thr", 0.5))

    summary = {
        "ref_run": ref["run"], "cmp_run": cmp["run"],
        "aligned_residues_sampled": S,
        "linear_cka": round(cka, 4),
        "ref_features": int(C.shape[0]), "cmp_features": int(C.shape[1]),
        "match_thr": thr,
        "ref_shared": int((ref_best > thr).sum()), "ref_unique": int((ref_best <= thr).sum()),
        "cmp_shared": int((cmp_best > thr).sum()), "cmp_unique": int((cmp_best <= thr).sum()),
        "ref_best_corr_median": round(float(np.median(ref_best)), 3),
        "cmp_best_corr_median": round(float(np.median(cmp_best)), 3),
    }
    out_dir = Path(cfg.paths.output_dir) / "compare" / f"{ref['run']}__vs__{cmp['run']}"
    out_dir.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(out_dir / "summary.json", "w"), indent=2)
    logger.info("Summary:\n" + json.dumps(summary, indent=2))

    # Per cmp-feature best cross-model match, annotated with the disorder atlas.
    cmp_tbl = pd.DataFrame({
        "cmp_feature": np.arange(C.shape[1]),
        "best_ref_feature": cmp_match,
        "best_abs_corr": np.round(cmp_best, 3),
        "shared": cmp_best > thr,
    })
    if cmp["atlas"] is not None:
        a = cmp["atlas"].rename(columns={"feature_idx": "cmp_feature"})
        keep = [col for col in ("cmp_feature", "triplet", "f1") if col in a.columns]
        cmp_tbl = cmp_tbl.merge(a[keep], on="cmp_feature", how="left")
    cmp_tbl.to_parquet(out_dir / "cmp_feature_match.parquet", index=False)

    # Disorder/sub-state-aligned cmp features that are UNIQUE to the cmp model.
    if "f1" in cmp_tbl.columns:
        uniq = cmp_tbl[(cmp_tbl["f1"] > float(c.get("disorder_f1", 0.2))) & (~cmp_tbl["shared"])]
        uniq = uniq.sort_values("f1", ascending=False)
        uniq.to_csv(out_dir / "cmp_disorder_unique.csv", index=False)
        logger.info(f"{c.cmp_model}-UNIQUE disorder/sub-state features (best_abs_corr<={thr}):\n"
                    + uniq.head(20).to_string(index=False))

    logger.info(f"Saved → {out_dir}")


if __name__ == "__main__":
    main()
