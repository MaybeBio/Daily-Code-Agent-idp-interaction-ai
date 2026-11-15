#!/usr/bin/env python3
"""
Fetch protein records from MobiDB and build the Parquet datasets.

Usage:
    python scripts/01_fetch_mobidb.py paths=local dataset=gold_curated
    python scripts/01_fetch_mobidb.py paths=colab dataset=full_swissprot
"""

import logging
import sys
from pathlib import Path

import hydra
from omegaconf import DictConfig

# Add project root to path so idpinterp is importable without pip install
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from idpinterp.data.fetch import MobiDBClient
from idpinterp.data.process import (
    build_region_index,
    deduplicate_against,
    parse_mobidb_entries,
    save_dataset,
)
from idpinterp.data.splits import create_splits

logger = logging.getLogger(__name__)


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    )

    logger.info(f"=== MobiDB data fetch ===")
    logger.info(f"Dataset config: {cfg.dataset.name}")
    logger.info(f"Cache dir: {cfg.paths.cache_dir}")
    logger.info(f"Output dir: {cfg.paths.dataset_dir}")

    # Initialize client
    client = MobiDBClient(
        cache_dir=Path(cfg.paths.cache_dir),
        page_size=100,  # Conservative for reliability
    )

    # Fetch data based on dataset config
    max_length = cfg.dataset.get("max_length", None)
    max_entries = cfg.dataset.get("max_proteins", None)
    feature_filter = cfg.dataset.get("query_filter", None)
    reviewed_only = cfg.dataset.get("reviewed_only", False)

    entries = client.fetch_by_filter(
        feature_filter=feature_filter,
        reviewed_only=reviewed_only,
        max_length=max_length,
        max_entries=max_entries,
    )

    # If there are additional filters, fetch those too (union)
    additional_filters = cfg.dataset.get("additional_filters", [])
    if additional_filters:
        seen_accs = {e["acc"] for e in entries}
        for filt in additional_filters:
            logger.info(f"Fetching additional filter: {filt}")
            extra = client.fetch_by_filter(
                feature_filter=filt,
                reviewed_only=reviewed_only,
                max_length=max_length,
            )
            for e in extra:
                if e["acc"] not in seen_accs:
                    entries.append(e)
                    seen_accs.add(e["acc"])
            logger.info(f"  -> Total unique proteins so far: {len(entries)}")

    if not entries:
        logger.error("No entries fetched. Check your filters and network connection.")
        return

    # Parse into DataFrame with annotation masks
    logger.info("Parsing MobiDB entries into structured DataFrame...")
    df = parse_mobidb_entries(entries)

    # Train/eval deduplication: drop eval proteins present in the training corpus
    dedup_against = cfg.dataset.get("dedup_against", None)
    if dedup_against:
        logger.info(f"Deduplicating against training set: {dedup_against}")
        df = deduplicate_against(
            df, Path(cfg.paths.dataset_dir) / f"{dedup_against}.parquet"
        )

    # Build region index for evaluation 
    split_mode = cfg.dataset.get("split_mode", "eval_only")
    if split_mode == "train_only":
        logger.info("train_only dataset → skipping region index build.")
        region_idx = None
    else:
        logger.info("Building region index for region-aware recall...")
        region_idx = build_region_index(df)

    # Apply splits
    test_fraction = cfg.dataset.get("test_fraction", 0.5)
    df = create_splits(
        df,
        split_mode=split_mode,
        test_fraction=test_fraction,
        random_seed=cfg.seed,
    )

    # Save outputs
    dataset_dir = Path(cfg.paths.dataset_dir)
    dataset_name = cfg.dataset.name

    save_dataset(
        df=df,
        output_path=dataset_dir / f"{dataset_name}.parquet",
        region_index=region_idx,
        region_index_path=dataset_dir / f"{dataset_name}_regions.parquet",
    )

    # Print summary
    logger.info(f"\n{'='*60}")
    logger.info(f"DATASET SUMMARY: {dataset_name}")
    logger.info(f"{'='*60}")
    logger.info(f"  Total proteins: {len(df)}")
    logger.info(f"  Total residues: {df['length'].sum():,}")
    logger.info(f"  Mean length: {df['length'].mean():.0f}")
    logger.info(f"  Max length: {df['length'].max()}")
    if "split" in df.columns:
        logger.info(f"  Splits: {df['split'].value_counts().to_dict()}")
    logger.info(f"  Total annotation regions: {len(region_idx) if region_idx is not None else 0}")
    logger.info(f"  Output: {dataset_dir / f'{dataset_name}.parquet'}")
    logger.info(f"{'='*60}")


if __name__ == "__main__":
    main()
