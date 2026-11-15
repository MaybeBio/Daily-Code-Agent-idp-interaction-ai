#!/usr/bin/env python3
"""
Standalone dataset fetch, independent of Hydra and the training environment.

Run wherever mobidb.org is reachable; the built Parquet files can then be moved to
whichever machine runs the pipeline.

Usage:
    python scripts/local_fetch.py --dataset gold_curated
    python scripts/local_fetch.py --dataset full_swissprot
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from idpinterp.constants import ALL_TARGET_TRIPLETS
from idpinterp.data.process import build_binary_mask, build_region_index

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(name)s] %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# ── Dataset configs (mirrors the YAML but avoids Hydra dependency) ──────────
DATASET_CONFIGS = {
    "gold_curated": {
        "query_filter": "curated-disorder-merge",
        "additional_filters": [
            "curated-lip-merge",
            "curated-phase_separation-merge",
            "curated-binding_mode_disorder_to_disorder-fuzdb",
        ],
        "max_length": 1022,
        "reviewed_only": False,
        "max_proteins": None,
        "split_mode": "eval_only",
        "test_fraction": 0.5,
    },
    "full_swissprot": {
        "query_filter": None,
        "additional_filters": [],
        "max_length": 1022,
        "reviewed_only": True,
        "max_proteins": 50000,
        "split_mode": "train_only",
        "test_fraction": 0.0,
    },
}


def create_session() -> requests.Session:
    """Create a requests session with retries and a browser User-Agent."""
    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 IDPInterp/0.1",
        "Accept": "application/json",
    })
    retries = Retry(total=5, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session


def fetch_from_mobidb(
    session: requests.Session,
    feature_filter: str | None,
    max_length: int | None,
    reviewed_only: bool,
    max_entries: int | None,
) -> list[dict]:
    """Fetch entries from MobiDB v6 API using cursor pagination."""
    BASE_URL = "https://mobidb.org/api"
    entries = []
    last_id = None

    while True:
        params = {
            "limit": 1000,
            "sort_asc": "_id",
            "format": "json",
        }
        if feature_filter:
            params[feature_filter] = "exists"
        if max_length:
            params["length_max"] = max_length
        if last_id:
            params["_id_gt"] = last_id

        try:
            resp = session.get(f"{BASE_URL}/download_page", params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            logger.error(f"Request failed: {e}")
            time.sleep(10)
            continue

        batch = data.get("data", [])
        if not batch:
            break

        if reviewed_only:
            batch = [e for e in batch if e.get("reviewed") is True]

        entries.extend(batch)
        logger.info(f"  Fetched {len(entries)} entries so far...")

        if max_entries and len(entries) >= max_entries:
            entries = entries[:max_entries]
            break

        last_id = data.get("metadata", {}).get("last_id")
        if not last_id:
            break

        # Be polite to the server
        time.sleep(0.5)

    return entries


def parse_entries(entries: list[dict]) -> pd.DataFrame:
    """Parse MobiDB v6 entries into a DataFrame with binary masks."""
    rows = []
    for entry in entries:
        acc = entry.get("acc")
        seq = entry.get("sequence", "")
        length = len(seq)

        if not acc or not seq:
            continue

        row = {
            "uniprot_id": acc,
            "sequence": seq,
            "length": length,
        }

        for triplet in ALL_TARGET_TRIPLETS:
            col_name = f"mask_{triplet.replace('-', '_')}"
            feature_data = entry.get(triplet)

            if feature_data and "regions" in feature_data:
                row[col_name] = build_binary_mask(length, feature_data["regions"])
            else:
                row[col_name] = None

        rows.append(row)

    return pd.DataFrame(rows)


def create_splits(df: pd.DataFrame, split_mode: str, test_fraction: float, seed: int = 42) -> pd.DataFrame:
    """Assign splits to the dataframe."""
    np.random.seed(seed)
    df = df.copy()

    if split_mode == "train_only":
        df["split"] = "train"
    elif split_mode == "eval_only":
        mask = np.random.rand(len(df)) < test_fraction
        df["split"] = np.where(mask, "test", "val")
    else:
        rand = np.random.rand(len(df))
        df["split"] = np.where(rand < 0.8, "train", np.where(rand < 0.9, "val", "test"))

    return df


def main():
    parser = argparse.ArgumentParser(description="Fetch MobiDB data locally")
    parser.add_argument("--dataset", required=True, choices=list(DATASET_CONFIGS.keys()))
    parser.add_argument("--output-dir", default="data/datasets")
    args = parser.parse_args()

    config = DATASET_CONFIGS[args.dataset]
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    session = create_session()

    # ── Main fetch ──
    logger.info(f"=== Fetching dataset: {args.dataset} ===")
    logger.info(f"  Filter: {config['query_filter']}")
    logger.info(f"  Max length: {config['max_length']}")
    logger.info(f"  Reviewed only: {config['reviewed_only']}")

    entries = fetch_from_mobidb(
        session,
        feature_filter=config["query_filter"],
        max_length=config["max_length"],
        reviewed_only=config["reviewed_only"],
        max_entries=config["max_proteins"],
    )

    # ── Additional filters (union) ──
    if config["additional_filters"]:
        seen_accs = {e["acc"] for e in entries}
        for filt in config["additional_filters"]:
            logger.info(f"Fetching additional filter: {filt}")
            extra = fetch_from_mobidb(
                session,
                feature_filter=filt,
                max_length=config["max_length"],
                reviewed_only=config["reviewed_only"],
                max_entries=None,
            )
            for e in extra:
                if e["acc"] not in seen_accs:
                    entries.append(e)
                    seen_accs.add(e["acc"])
            logger.info(f"  -> Total unique proteins: {len(entries)}")

    if not entries:
        logger.error("No entries fetched! Check your network connection.")
        return

    # ── Parse ──
    logger.info("Parsing entries into structured DataFrame...")
    df = parse_entries(entries)

    # ── Region index ──
    logger.info("Building region index...")
    region_idx = build_region_index(df)

    # ── Splits ──
    df = create_splits(df, config["split_mode"], config["test_fraction"])

    # ── Save ──
    df.to_parquet(output_dir / f"{args.dataset}.parquet", engine="pyarrow", index=False)
    region_idx.to_parquet(output_dir / f"{args.dataset}_regions.parquet", engine="pyarrow", index=False)

    # ── Summary ──
    logger.info(f"\n{'='*60}")
    logger.info(f"DATASET SUMMARY: {args.dataset}")
    logger.info(f"{'='*60}")
    logger.info(f"  Total proteins: {len(df)}")
    logger.info(f"  Total residues: {df['length'].sum():,}")
    logger.info(f"  Mean length: {df['length'].mean():.0f}")
    logger.info(f"  Max length: {df['length'].max()}")
    if "split" in df.columns:
        logger.info(f"  Splits: {df['split'].value_counts().to_dict()}")
    logger.info(f"  Total annotation regions: {len(region_idx)}")
    logger.info(f"  Saved to: {output_dir / f'{args.dataset}.parquet'}")
    logger.info(f"{'='*60}")


if __name__ == "__main__":
    main()
