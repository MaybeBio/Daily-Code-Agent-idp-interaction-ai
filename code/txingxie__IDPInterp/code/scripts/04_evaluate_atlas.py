#!/usr/bin/env python3
"""
Evaluate trained SAE features against MobiDB annotations.

Usage:
    python scripts/04_evaluate_atlas.py paths=colab model=esm2_8m sae=l1_32x
    python scripts/04_evaluate_atlas.py paths=colab model=esm2_8m sae=l1_32x eval.force_reencode=true
"""

import json
import logging
import sys
from pathlib import Path

import hydra
import numpy as np
import pandas as pd
import torch
from omegaconf import DictConfig
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from idpinterp.constants import ALL_TARGET_TRIPLETS, ESM2_MAX_SEQ_LEN
from idpinterp.eval.atlas import AtlasBuilder
from idpinterp.eval.metrics import RegionAwareEvaluator
from idpinterp.eval.normalize import FeatureNormalizer
from idpinterp.models.memmap_store import MemmapReader, MetadataIndex, load_mmap_config
from idpinterp.sae.architecture import build_sae
from idpinterp.sae.dataset import MemmapDataset

logger = logging.getLogger(__name__)


# ── Helpers ────────────────────────────────────────────────────────────────────

def _build_region_lookup(region_df: pd.DataFrame) -> dict:
    """
    Precompute {(uniprot_id, triplet): [(start, end), ...]} once so triplet/label
    construction is O(regions) instead of O(proteins × triplets × |region_df|)
    pandas filters (which becomes a multi-minute host stall at 20 triplets).
    """
    lookup: dict = {}
    if region_df is None or len(region_df) == 0:
        return lookup
    for (uid, trip), grp in region_df.groupby(["uniprot_id", "triplet"], sort=False):
        lookup[(uid, trip)] = list(
            zip(grp["region_start"].astype(int), grp["region_end"].astype(int))
        )
    return lookup


def _build_triplet_data(
    eval_df: pd.DataFrame,
    region_lookup: dict,
    meta_df: pd.DataFrame,
    total_residues: int,
    triplet: str,
) -> tuple[np.ndarray, list[tuple[int, int]], np.ndarray]:
    """
    Build per-residue label array and region list for one triplet.

    Returns:
        labels: (total_residues,) float32 binary array
        regions: list of (start, end) in global memmap coordinates
        region_val_mask: (len(regions),) bool, True where the region's protein
            is in the validation split (used for leakage-free evaluation).
    """
    col = f"mask_{triplet.replace('-', '_')}"
    if col not in eval_df.columns:
        return None, None, None

    labels  = np.zeros(total_residues, dtype=np.float32)
    regions = []
    region_val = []

    # Build a fast lookup: uniprot_id → meta row
    meta_lookup = meta_df.set_index("uniprot_id")
    has_split = "split" in eval_df.columns

    for _, protein in eval_df.iterrows():
        uid      = protein["uniprot_id"]
        mask_str = protein.get(col)
        if mask_str is None or uid not in meta_lookup.index:
            continue

        meta_row   = meta_lookup.loc[uid]
        start_idx  = int(meta_row["mmap_start_idx"])
        stored_len = int(meta_row["stored_length"])
        is_val     = bool(protein.get("split") == "val") if has_split else False

        # Per-residue labels
        mask = np.frombuffer(
            mask_str[:stored_len].encode(), dtype=np.uint8
        ) - ord("0")                                        # "010..." → [0,1,0,...]
        labels[start_idx : start_idx + stored_len] = mask.astype(np.float32)

        # Region boundaries (O(1) dict lookup, precomputed)
        for r_start, r_end in region_lookup.get((uid, triplet), []):
            r_end = min(int(r_end), stored_len)
            r_start = int(r_start)
            if r_start < stored_len:
                regions.append((start_idx + r_start, start_idx + r_end))
                region_val.append(is_val)

    return labels, regions, np.array(region_val, dtype=bool)


def _build_residue_split_mask(
    eval_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    total_residues: int,
) -> np.ndarray | None:
    """
    Build a (total_residues,) bool mask, True where the residue's protein is in
    the validation split. Returns None if no 'split' column exists.
    """
    if "split" not in eval_df.columns:
        return None
    val_mask = np.zeros(total_residues, dtype=bool)
    meta_lookup = meta_df.set_index("uniprot_id")
    for _, protein in eval_df.iterrows():
        uid = protein["uniprot_id"]
        if protein.get("split") != "val" or uid not in meta_lookup.index:
            continue
        meta_row   = meta_lookup.loc[uid]
        start_idx  = int(meta_row["mmap_start_idx"])
        stored_len = int(meta_row["stored_length"])
        val_mask[start_idx : start_idx + stored_len] = True
    return val_mask


def _cache_paths(act_cache_dir: Path, run_name: str, layer: int) -> dict:
    """Return a dict of all cache file paths for this run."""
    base = act_cache_dir / f"{run_name}"
    return {
        "mmap":      base.with_suffix(".mmap"),
        "cfg":       base.with_suffix(".json"),
        "max_acts":  base.parent / f"{run_name}_max_acts.npy",
        "reg_maxes": base.parent / f"{run_name}_region_maxes.npz",
    }


# ── Main ───────────────────────────────────────────────────────────────────────

@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    )

    layer = cfg.model.get("layer")
    if layer is None:
        layer = cfg.model.layers_to_extract[-1]
    _tag = cfg.sae.get("tag", None)
    _tag_suffix = f"_{_tag}" if _tag else ""
    run_name = f"{cfg.model.name}_{cfg.sae.type}_{cfg.sae.expansion_factor}x{_tag_suffix}_layer{layer}"
    logger.info(f"=== Atlas evaluation [{run_name}] ===")

    device = "cuda" if torch.cuda.is_available() else "cpu"

    # ── Load eval dataset first (needed for region maxes during encoding) ─────
    eval_dataset_path = Path(cfg.paths.dataset_dir) / "gold_curated.parquet"
    if not eval_dataset_path.exists():
        logger.error(f"Evaluation dataset not found: {eval_dataset_path}")
        return

    logger.info("Loading evaluation dataset and annotation regions...")
    eval_df     = pd.read_parquet(eval_dataset_path)
    regions_path = Path(cfg.paths.dataset_dir) / "gold_curated_regions.parquet"
    region_df   = pd.read_parquet(regions_path) if regions_path.exists() else pd.DataFrame()

    # ── Embedding mmap (source for SAE encoding) ───────────────────────────────
    embed_dir   = Path(cfg.paths.embedding_dir) / cfg.model.name / cfg.dataset.name
    mmap_path   = embed_dir / f"layer{layer}.mmap"
    mmap_cfg    = load_mmap_config(mmap_path)
    total_residues = mmap_cfg["total_residues"]
    hidden_dim     = mmap_cfg["hidden_dim"]
    meta_df        = MetadataIndex.load(embed_dir / f"layer{layer}_metadata.parquet")

    # ── Residue-level validation mask (for leakage-free evaluation) ────────────
    val_residue_mask = _build_residue_split_mask(eval_df, meta_df, total_residues)
    if val_residue_mask is None:
        logger.warning(
            "No 'split' column in eval dataset — falling back to FULL-CORPUS "
            "evaluation (threshold selection leakage). Re-fetch with create_splits."
        )
    else:
        logger.info(
            f"Val/test split enforced: {int(val_residue_mask.sum())} val residues, "
            f"{int((~val_residue_mask).sum())} test residues."
        )

    # ── Build triplet data (labels + regions) up front ─────────────────────────
    # Optionally restrict to a subset (e.g. for a fast layer sweep): pass
    # +eval_triplets=[prediction-disorder-alphafold,curated-phase_separation-merge,...]
    target_triplets = list(cfg.get("eval_triplets", []) or ALL_TARGET_TRIPLETS)
    logger.info(
        f"Building region lookup + per-triplet labels for {len(target_triplets)} triplets..."
    )
    region_lookup = _build_region_lookup(region_df)
    triplet_data: dict[str, tuple[np.ndarray, list, np.ndarray]] = {}
    for triplet in target_triplets:
        labels, regions, region_val_mask = _build_triplet_data(
            eval_df, region_lookup, meta_df, total_residues, triplet
        )
        if labels is not None:
            triplet_data[triplet] = (labels, regions, region_val_mask)
            logger.info(
                f"  {triplet}: {int(labels.sum())} labeled residues, "
                f"{len(regions)} regions"
            )

    if not triplet_data:
        logger.error("No valid triplets found in eval dataset.")
        return

    # ── Activation cache ────────────────────────────────────────────────────────
    scratch_root = cfg.paths.get("scratch_dir", None) or cfg.paths.output_dir
    act_cache_dir = Path(scratch_root) / "act_cache"
    act_cache_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Activation cache dir: {act_cache_dir}")
    paths = _cache_paths(act_cache_dir, run_name, layer)

    force_reencode = cfg.get("eval", {}).get("force_reencode", False)
    cache_valid = (
        paths["mmap"].exists()
        and paths["cfg"].exists()
        and paths["max_acts"].exists()
        and paths["reg_maxes"].exists()
        and not force_reencode
    )

    # ── Load SAE ───────────────────────────────────────────────
    ckpt_dir  = Path(cfg.paths.checkpoint_dir)
    final_ckpt = ckpt_dir / f"{run_name}_final.pt"
    if not final_ckpt.exists():
        ckpts = sorted(ckpt_dir.glob(f"{run_name}_step*.pt"))
        if not ckpts:
            logger.error(f"No checkpoint found for {run_name}.")
            return
        final_ckpt = ckpts[-1]

    logger.info(f"Loading SAE from: {final_ckpt}")
    sae = build_sae(cfg.sae, cfg.model)
    ckpt = torch.load(final_ckpt, map_location=device, weights_only=True)
    sae.load_state_dict(ckpt["model_state_dict"])
    sae = sae.to(device).eval()

    if cache_valid:
        # ── Cache hit: skip encoding entirely ─────────────────────────────────
        logger.info(f"✓ Using cached activations from: {paths['mmap']}")
        with open(paths["cfg"]) as f:
            act_cfg = json.load(f)

        all_activations = np.memmap(
            paths["mmap"],
            dtype=np.float16,
            mode="r",
            shape=(act_cfg["total_residues"], act_cfg["d_latent"]),
            order="C",
        )

        # Restore SAE normalization so atlas thresholds are meaningful
        max_acts_np = np.load(paths["max_acts"])
        max_acts    = torch.from_numpy(max_acts_np).to(device)
        FeatureNormalizer.normalize_sae(sae, max_acts)

        # Load cached region maxes
        reg_maxes_cache = np.load(paths["reg_maxes"])
        region_maxes_by_triplet: dict[str, np.ndarray] = {
            k: reg_maxes_cache[k] for k in reg_maxes_cache.files
        }
        logger.info(f"  Loaded region maxes for {len(region_maxes_by_triplet)} triplets.")

    else:
        # ── Cache miss: encode and cache ──────────────────────────────────────
        logger.info("Cache miss — running encoding pass (will be cached for next run).")

        dataset = MemmapDataset(mmap_path, total_residues, hidden_dim)

        # Normalize first (required before encoding so activations are in [0,1])
        logger.info("Normalizing SAE features...")
        max_acts = FeatureNormalizer.compute_max_activations(
            sae, dataset, batch_size=4096, device=device
        )
        FeatureNormalizer.normalize_sae(sae, max_acts)
        np.save(paths["max_acts"], max_acts.cpu().numpy())
        logger.info(f"  Saved max_acts → {paths['max_acts']}")

        # Preallocate activation cache as C-order memmap (fast sequential row writes during encoding).
        # Where RAM exceeds the mmap size the OS page cache holds the whole file
        # after the first region-max scan, so column reads during evaluation are then RAM-speed.
        d_latent = sae.d_latent
        logger.info(
            f"Allocating activation cache: ({total_residues}, {d_latent}) float16 C-order"
            f" ≈ {total_residues * d_latent * 2 / 1e9:.1f} GB"
        )
        all_activations = np.memmap(
            paths["mmap"],
            dtype=np.float16,
            mode="w+",
            shape=(total_residues, d_latent),
            order="C",
        )

        # Precompute batch→region overlap lookup to avoid inner Python loops
        batch_size_enc = 4096
        n_batches      = (total_residues + batch_size_enc - 1) // batch_size_enc
        triplet_names  = list(triplet_data.keys())
        n_triplets     = len(triplet_names)

        # region_maxes_by_triplet: dict[triplet_name → (R, F) float32]
        region_maxes_by_triplet = {
            t: np.zeros((len(triplet_data[t][1]), d_latent), dtype=np.float32)
            for t in triplet_names
        }

        # batch_overlaps[b] = list of (triplet_name, region_idx, local_start, local_end)
        logger.info("Building batch→region overlap lookup...")
        batch_overlaps: list[list] = [[] for _ in range(n_batches)]
        for trip_name in triplet_names:
            _, regions, _ = triplet_data[trip_name]
            for r_idx, (rstart, rend) in enumerate(regions):
                first_b = rstart // batch_size_enc
                last_b  = (rend - 1) // batch_size_enc
                for b in range(first_b, min(last_b + 1, n_batches)):
                    b_start = b * batch_size_enc
                    b_end   = min(b_start + batch_size_enc, total_residues)
                    ls = max(rstart, b_start) - b_start
                    le = min(rend,   b_end)   - b_start
                    if le > ls:
                        batch_overlaps[b].append((trip_name, r_idx, ls, le))

        # Encoding loop: encode → write mmap + update region maxes simultaneously
        logger.info("Encoding residues and computing region maxes...")
        loader = torch.utils.data.DataLoader(
            dataset, batch_size=batch_size_enc, shuffle=False, num_workers=2
        )
        row_idx = 0
        with torch.no_grad():
            for b_idx, batch in enumerate(tqdm(loader, desc="Encoding", mininterval=30)):
                z = sae.encode(batch.to(device)).half()    # (bz, F)
                bz = z.shape[0]

                # Write to F-order mmap
                all_activations[row_idx : row_idx + bz] = z.cpu().numpy()

                # Update region maxes for overlapping regions
                overlaps = batch_overlaps[b_idx]
                if overlaps:
                    z_f32 = z.float().cpu().numpy()        # (bz, F) temporary
                    for trip_name, r_idx, ls, le in overlaps:
                        region_maxes_by_triplet[trip_name][r_idx] = np.maximum(
                            region_maxes_by_triplet[trip_name][r_idx],
                            z_f32[ls:le].max(axis=0),
                        )

                row_idx += bz

        # Flush mmap
        all_activations.flush()
        logger.info(f"  Saved activation cache → {paths['mmap']}")

        # Save config and region maxes
        with open(paths["cfg"], "w") as f:
            json.dump({"total_residues": total_residues, "d_latent": d_latent}, f)
        np.savez(paths["reg_maxes"], **region_maxes_by_triplet)
        logger.info(f"  Saved region maxes → {paths['reg_maxes']}")

        # Re-open as read-only so subsequent code uses the mmap, not the w+ handle
        del all_activations
        all_activations = np.memmap(
            paths["mmap"],
            dtype=np.float16,
            mode="r",
            shape=(total_residues, d_latent),
            order="C",
        )

    # ── Evaluate all triplets ─────────────────────────────────────────────────
    # ── Pre-warm OS page cache ─────────────────────────────────────────
    # A sequential sum scan loads every mmap page into RAM before the chunked
    # column reads begin, trading one sequential pass for repeated page-fault
    # stalls across triplets. Matters most on network-backed storage.
    logger.info("Pre-warming mmap into OS page cache (sequential scan)...")
    _warm = all_activations.sum(axis=1)   # row-wise sum, sequential C-order read
    del _warm

    evaluator      = RegionAwareEvaluator(chunk_size=int(cfg.get("eval_chunk_size", 256)))
    triplet_results: dict[str, pd.DataFrame] = {}

    atlas_dir = Path(cfg.paths.output_dir) / "atlas" / run_name
    triplets_dir = atlas_dir / "triplets"
    triplets_dir.mkdir(parents=True, exist_ok=True)

    for triplet, (labels, regions, region_val_mask) in triplet_data.items():
        triplet_csv = triplets_dir / f"{triplet}.csv"

        # Checkpoint resume
        if triplet_csv.exists() and not force_reencode:
            logger.info(f"✓ Loaded cached evaluation for triplet: {triplet}")
            triplet_results[triplet] = pd.read_csv(triplet_csv)
            continue

        logger.info(
            f"Evaluating triplet: {triplet} "
            f"({int(labels.sum())} labeled residues, {len(regions)} regions)"
        )

        # Use cached region maxes (precomputed during encoding)
        reg_mx = region_maxes_by_triplet.get(triplet)

        result = evaluator.evaluate_all_features(
            all_activations,
            labels,
            regions,
            triplet,
            region_maxes_matrix=reg_mx,
            device=device,
            val_residue_mask=val_residue_mask,
            val_region_mask=region_val_mask if val_residue_mask is not None else None,
        )

        if result:
            df = pd.DataFrame(result)
            triplet_results[triplet] = df
            df.to_csv(triplet_csv, index=False)
            logger.info(f"  Saved {triplet} checkpoint → {triplet_csv}")

    # ── Build and save atlas ───────────────────────────────────────────────────
    atlas_dir = Path(cfg.paths.output_dir) / "atlas" / run_name
    builder   = AtlasBuilder(atlas_dir)
    atlas     = builder.build_atlas(triplet_results)

    # Flag candidate missing annotations
    # Use last triplet's labels for the missing-annotation scan (any labels work)
    last_labels = list(triplet_data.values())[-1][0]
    missing = builder.flag_missing_annotations(
        atlas, all_activations, last_labels, meta_df
    )
    builder.save(atlas, missing)

    # ── Print top features ─────────────────────────────────────────────────────
    logger.info(f"\n{'='*60}")
    logger.info("TOP 20 FEATURES BY F1:")
    logger.info(f"{'='*60}")
    if len(atlas) > 0:
        for _, row in atlas.head(20).iterrows():
            logger.info(
                f"  Feature {row['feature_idx']:5d} | "
                f"{row['triplet']:45s} | "
                f"F1={row['f1']:.3f} | "
                f"P={row['precision']:.3f} | "
                f"R={row['region_recall']:.3f} | "
                f"T={row['best_threshold']:.3f}"
            )

    logger.info(f"\n=== Evaluation complete: {run_name} ===")


if __name__ == "__main__":
    main()
