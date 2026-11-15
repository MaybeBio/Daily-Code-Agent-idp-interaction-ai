#!/usr/bin/env python3
"""
Diagnostics computed from an existing feature atlas; no SAE retraining.

Sub-analyses (run independently; each guards its own dependencies):

  1. specificity  - Feature x triplet F1 matrix from the per-triplet CSVs, measuring
                    how triplet-exclusive each feature is. Needs only pandas.
  2. thresholds   - For features whose best threshold is below 1e-3, compares activation
                    distributions on disordered versus ordered residues, testing whether
                    region-aware recall is inflated by always-on features. Needs the
                    activation cache and a disorder triplet.
  3. deadfeatures - Forward-passes the SAE over the eval corpus and reports how many
                    training-time dead features activate on this corpus. Needs torch,
                    the SAE checkpoint, and embeddings.

Usage:
    python scripts/06_diagnostics.py paths=colab model=esm2_8m sae=l1_32x +diag=specificity
    python scripts/06_diagnostics.py paths=colab model=esm2_8m sae=l1_32x +diag=thresholds
    python scripts/06_diagnostics.py paths=colab model=esm2_8m sae=l1_32x +diag=deadfeatures
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

from idpinterp.constants import ALL_TARGET_TRIPLETS

logger = logging.getLogger(__name__)

# Generic-disorder reference triplets (a sub-state-specific feature should beat
# its own generic-disorder F1 by a clear margin).
GENERIC_DISORDER = [
    "prediction-disorder-th_50",
    "curated-disorder-merge",
    "curated-disorder-priority",
    "homology-disorder-merge",
    "prediction-disorder-alphafold",
]
# Triplets that name a specific disorder sub-state.
SUBSTATE_TRIPLETS = [
    "derived-binding_mode_disorder_to_order-mobi",      # MoRF
    "derived-binding_mode_disorder_to_disorder-mobi",    # fuzzy complex
    "curated-binding_mode_disorder_to_disorder-fuzdb",
    "prediction-compact-mobidb_lite_sub",               # molten-globule proxy
    "prediction-extended-mobidb_lite_sub",              # random-coil proxy
    "curated-phase_separation-merge",                    # LLPS
    "curated-lip-mfib",                                  # amyloid
]


def _scratch_act_cache(cfg, run_name: str) -> Path:
    root = cfg.paths.get("scratch_dir", None) or cfg.paths.output_dir
    return Path(root) / "act_cache" / f"{run_name}.mmap"


def _run_name(cfg: DictConfig) -> tuple[str, int]:
    layer = cfg.model.get("layer")
    if layer is None:
        layer = cfg.model.layers_to_extract[-1]
    tag = cfg.sae.get("tag", None)
    suffix = f"_{tag}" if tag else ""
    return (
        f"{cfg.model.name}_{cfg.sae.type}_{cfg.sae.expansion_factor}x{suffix}_layer{layer}",
        layer,
    )


# ── 1. Feature x triplet specificity matrix ─────────────────────────────────

def diag_specificity(cfg: DictConfig, run_name: str, out_dir: Path) -> None:
    triplets_dir = Path(cfg.paths.output_dir) / "atlas" / run_name / "triplets"
    if not triplets_dir.exists():
        logger.error(f"No per-triplet CSVs at {triplets_dir}. Run scripts/04_evaluate_atlas.py first.")
        return

    frames = []
    for csv in sorted(triplets_dir.glob("*.csv")):
        df = pd.read_csv(csv)
        if "triplet" not in df.columns:
            df["triplet"] = csv.stem
        frames.append(df[["feature_idx", "triplet", "f1"]])
    if not frames:
        logger.error("No triplet CSVs found.")
        return

    long = pd.concat(frames, ignore_index=True)
    mat = long.pivot_table(
        index="feature_idx", columns="triplet", values="f1", fill_value=0.0
    )
    mat.to_parquet(out_dir / "feature_triplet_f1_matrix.parquet")

    # Specificity stats at a meaningful-alignment threshold
    thr = 0.3
    aligned = mat >= thr
    n_aligned_triplets = aligned.sum(axis=1)           # per feature
    features_with_any = (n_aligned_triplets > 0).sum()
    exclusive = (n_aligned_triplets == 1).sum()
    promiscuous = (n_aligned_triplets >= 3).sum()

    summary = {
        "n_features_total": int(mat.shape[0]),
        "n_triplets": int(mat.shape[1]),
        "alignment_threshold": thr,
        "features_aligned_to_any": int(features_with_any),
        "features_triplet_exclusive": int(exclusive),
        "features_promiscuous_ge3": int(promiscuous),
        "exclusive_fraction_of_aligned": (
            float(exclusive / features_with_any) if features_with_any else 0.0
        ),
        "mean_triplets_per_aligned_feature": (
            float(n_aligned_triplets[n_aligned_triplets > 0].mean())
            if features_with_any
            else 0.0
        ),
    }
    with open(out_dir / "specificity_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    logger.info("Specificity summary:\n" + json.dumps(summary, indent=2))
    logger.info(f"Saved matrix → {out_dir / 'feature_triplet_f1_matrix.parquet'}")

    _substate_analysis(mat, out_dir)


def _substate_analysis(mat: pd.DataFrame, out_dir: Path) -> None:
    """
    For each feature, contrast its best disorder
    SUB-STATE F1 against its best GENERIC-disorder F1. Features with a large
    positive margin are candidates that resolve a specific sub-state (MoRF /
    fuzzy / molten-globule / random-coil / LLPS / amyloid) beyond what a generic
    disorder label captures — exactly what current scalar predictors cannot do.
    """
    sub_cols = [c for c in SUBSTATE_TRIPLETS if c in mat.columns]
    gen_cols = [c for c in GENERIC_DISORDER if c in mat.columns]
    if not sub_cols:
        logger.warning("No sub-state triplets present in the matrix; "
                       "re-run scripts/04_evaluate_atlas.py with the expanded ALL_TARGET_TRIPLETS.")
        return

    best_sub = mat[sub_cols].max(axis=1)
    best_sub_name = mat[sub_cols].idxmax(axis=1)
    best_gen = mat[gen_cols].max(axis=1) if gen_cols else pd.Series(0.0, index=mat.index)

    ranking = pd.DataFrame({
        "feature_idx": mat.index,
        "best_substate": best_sub_name.values,
        "substate_f1": best_sub.values,
        "generic_disorder_f1": best_gen.values,
        "substate_margin": (best_sub - best_gen).values,
    }).sort_values("substate_margin", ascending=False)
    ranking.to_csv(out_dir / "substate_specific_features.csv", index=False)
    logger.info(
        "Top sub-state-specific features (substate_f1 high, generic low):\n"
        + ranking.head(20).to_string(index=False)
    )

    # Heatmap of the top sub-state-specific features × all triplets
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        top = ranking.head(40)["feature_idx"].tolist()
        cols = gen_cols + sub_cols
        hm = mat.loc[top, cols]
        fig, ax = plt.subplots(figsize=(max(8, len(cols) * 0.6), max(6, len(top) * 0.25)))
        im = ax.imshow(hm.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([c.replace("-", "\n") for c in cols], fontsize=6, rotation=90)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top, fontsize=5)
        ax.axvline(len(gen_cols) - 0.5, color="red", lw=1)  # generic | sub-state divider
        ax.set_title("Top sub-state-specific features  (left: generic disorder | right: sub-states)")
        fig.colorbar(im, ax=ax, label="test F1")
        fig.tight_layout()
        fig.savefig(out_dir / "substate_heatmap.png", dpi=150)
        logger.info(f"Saved heatmap → {out_dir / 'substate_heatmap.png'}")
    except Exception as e:
        logger.warning(f"Heatmap skipped ({e}); CSV ranking still written.")


# ── 2. Near-zero-threshold diagnostic ───────────────────────────────────────

def diag_thresholds(cfg: DictConfig, run_name: str, out_dir: Path) -> None:
    atlas_path = Path(cfg.paths.output_dir) / "atlas" / run_name / "feature_atlas.parquet"
    if not atlas_path.exists():
        # fall back to the bundled notebooks copy if present
        atlas_path = Path("notebooks") / run_name / "feature_atlas.parquet"
    if not atlas_path.exists():
        logger.error(f"No feature_atlas.parquet found for {run_name}.")
        return

    atlas = pd.read_parquet(atlas_path)
    flagged = atlas[atlas["best_threshold"] < 1e-3].sort_values("f1", ascending=False)
    logger.info(
        f"{len(flagged)} features with best_threshold < 1e-3 (always-on suspects):\n"
        + flagged[["feature_idx", "triplet", "f1", "precision",
                   "region_recall", "best_threshold"]].head(20).to_string()
    )

    # Activation-distribution comparison needs the cache + a disorder triplet.
    act_cache = _scratch_act_cache(cfg, run_name)
    cfg_json = act_cache.with_suffix(".json")
    if not (act_cache.exists() and cfg_json.exists()):
        logger.warning(
            "Activation cache not found — reporting flagged feature list only. "
            "Run scripts/04_evaluate_atlas.py to enable the disordered-vs-ordered distribution test."
        )
        flagged.to_csv(out_dir / "near_zero_threshold_features.csv", index=False)
        return

    with open(cfg_json) as f:
        ac = json.load(f)
    acts = np.memmap(
        act_cache, dtype=np.float16, mode="r",
        shape=(ac["total_residues"], ac["d_latent"]), order="C",
    )

    # Build a disordered-residue mask from a curated triplet on the eval corpus.
    eval_df = pd.read_parquet(Path(cfg.paths.dataset_dir) / "gold_curated.parquet")
    embed_dir = Path(cfg.paths.embedding_dir) / cfg.model.name / cfg.dataset.name
    _, layer = _run_name(cfg)
    from idpinterp.models.memmap_store import MetadataIndex
    meta = MetadataIndex.load(embed_dir / f"layer{layer}_metadata.parquet").set_index("uniprot_id")

    dis_mask = np.zeros(ac["total_residues"], dtype=bool)
    col = "mask_curated_disorder_merge"
    for _, p in eval_df.iterrows():
        uid, m = p["uniprot_id"], p.get(col)
        if not isinstance(m, str) or uid not in meta.index:
            continue
        s = int(meta.loc[uid, "mmap_start_idx"]); L = int(meta.loc[uid, "stored_length"])
        arr = (np.frombuffer(m[:L].encode(), dtype=np.uint8) - ord("0")).astype(bool)
        dis_mask[s:s + L] = arr

    rows = []
    for fidx in flagged["feature_idx"].head(10):
        col_act = np.asarray(acts[:, int(fidx)], dtype=np.float32)
        rows.append({
            "feature_idx": int(fidx),
            "mean_disordered": float(col_act[dis_mask].mean()),
            "mean_ordered": float(col_act[~dis_mask].mean()),
            "frac_active_disordered": float((col_act[dis_mask] > 1e-6).mean()),
            "frac_active_ordered": float((col_act[~dis_mask] > 1e-6).mean()),
        })
    res = pd.DataFrame(rows)
    res["separation_ratio"] = res["mean_disordered"] / (res["mean_ordered"] + 1e-9)
    res.to_csv(out_dir / "near_zero_threshold_distributions.csv", index=False)
    logger.info("Disordered-vs-ordered activation for near-zero features:\n"
                + res.to_string())


# ── 3. Dead-feature revival probe ───────────────────────────────────────────

def diag_deadfeatures(cfg: DictConfig, run_name: str, out_dir: Path) -> None:
    import torch
    from idpinterp.models.memmap_store import load_mmap_config
    from idpinterp.sae.architecture import build_sae
    from idpinterp.sae.dataset import MemmapDataset

    device = "cuda" if torch.cuda.is_available() else "cpu"
    _, layer = _run_name(cfg)
    embed_dir = Path(cfg.paths.embedding_dir) / cfg.model.name / cfg.dataset.name
    mmap_path = embed_dir / f"layer{layer}.mmap"
    mcfg = load_mmap_config(mmap_path)

    ckpt_dir = Path(cfg.paths.checkpoint_dir)
    ckpt = ckpt_dir / f"{run_name}_final.pt"
    if not ckpt.exists():
        cands = sorted(ckpt_dir.glob(f"{run_name}_step*.pt"))
        if not cands:
            logger.error(f"No checkpoint for {run_name}.")
            return
        ckpt = cands[-1]

    sae = build_sae(cfg.sae, cfg.model)
    sae.load_state_dict(torch.load(ckpt, map_location=device, weights_only=True)["model_state_dict"])
    sae = sae.to(device).eval()

    ds = MemmapDataset(mmap_path, mcfg["total_residues"], mcfg["hidden_dim"])
    loader = torch.utils.data.DataLoader(ds, batch_size=4096, shuffle=False, num_workers=2)
    ever_active = torch.zeros(sae.d_latent, dtype=torch.bool, device=device)
    with torch.no_grad():
        for batch in loader:
            z = sae.encode(batch.to(device))
            ever_active |= (z.max(dim=0).values > 1e-5)

    n_active = int(ever_active.sum().item())
    n_dead = sae.d_latent - n_active
    logger.info(
        f"On eval corpus: {n_active}/{sae.d_latent} features active, "
        f"{n_dead} never fire ({100 * n_dead / sae.d_latent:.1f}% dead on eval)."
    )
    np.save(out_dir / "eval_active_features.npy", ever_active.cpu().numpy())


# ── 0. Cheap sparsity report for tier-1 SAE selection ───────────────────────

def diag_sparsity(cfg: DictConfig, run_name: str, out_dir: Path) -> None:
    """
    """
    import torch
    from idpinterp.models.memmap_store import load_mmap_config
    from idpinterp.sae.architecture import build_sae
    from idpinterp.sae.dataset import MemmapDataset

    device = "cuda" if torch.cuda.is_available() else "cpu"
    _, layer = _run_name(cfg)
    embed_dir = Path(cfg.paths.embedding_dir) / cfg.model.name / cfg.dataset.name
    mmap_path = embed_dir / f"layer{layer}.mmap"
    mcfg = load_mmap_config(mmap_path)

    ckpt_dir = Path(cfg.paths.checkpoint_dir)
    ckpt = ckpt_dir / f"{run_name}_final.pt"
    if not ckpt.exists():
        cands = sorted(ckpt_dir.glob(f"{run_name}_step*.pt"))
        if not cands:
            logger.error(f"No checkpoint for {run_name}.")
            return
        ckpt = cands[-1]
    sae = build_sae(cfg.sae, cfg.model)
    sae.load_state_dict(torch.load(ckpt, map_location=device, weights_only=True)["model_state_dict"])
    sae = sae.to(device).eval()

    cap = int(cfg.get("diag_max_residues", 300_000))
    ds = MemmapDataset(mmap_path, mcfg["total_residues"], mcfg["hidden_dim"])
    loader = torch.utils.data.DataLoader(ds, batch_size=4096, shuffle=False, num_workers=2)

    seen = 0
    l0_sum = 0.0
    ever_active = torch.zeros(sae.d_latent, dtype=torch.bool, device=device)
    sse = 0.0  # residual sum of squares
    tss = 0.0  # total sum of squares (var around mean)
    with torch.no_grad():
        for batch in loader:
            x = batch.to(device)
            z = sae.encode(x)
            x_hat = sae.decode(z)
            l0_sum += (z > 1e-6).float().sum(dim=1).sum().item()
            ever_active |= (z.max(dim=0).values > 1e-6)
            sse += ((x - x_hat) ** 2).sum().item()
            tss += ((x - x.mean(dim=0, keepdim=True)) ** 2).sum().item()
            seen += x.shape[0]
            if seen >= cap:
                break

    n_dead = int((~ever_active).sum().item())
    report = {
        "run_name": run_name,
        "sae_type": cfg.sae.type,
        "k": cfg.sae.get("k", None),
        "l1_coeff": cfg.sae.get("l1_coeff", None),
        "residues_scanned": seen,
        "L0_mean": l0_sum / max(1, seen),
        "dead_features": n_dead,
        "dead_pct": 100.0 * n_dead / sae.d_latent,
        "fve": 1.0 - sse / max(1e-9, tss),
    }
    with open(out_dir / "sparsity_report.json", "w") as f:
        json.dump(report, f, indent=2)
    logger.info("Sparsity report:\n" + json.dumps(report, indent=2))


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(name)s] %(levelname)s: %(message)s"
    )
    run_name, _ = _run_name(cfg)
    diag = cfg.get("diag", "specificity")
    out_dir = Path(cfg.paths.output_dir) / "diagnostics" / run_name
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"=== Diagnostic: {diag} [{run_name}] ===")

    {
        "sparsity": diag_sparsity,
        "specificity": diag_specificity,
        "thresholds": diag_thresholds,
        "deadfeatures": diag_deadfeatures,
    }[diag](cfg, run_name, out_dir)


if __name__ == "__main__":
    main()
