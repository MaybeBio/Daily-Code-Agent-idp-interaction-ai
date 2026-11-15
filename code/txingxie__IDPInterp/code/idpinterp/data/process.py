import logging
from pathlib import Path

import numpy as np
import pandas as pd

from idpinterp.constants import ALL_TARGET_TRIPLETS

logger = logging.getLogger(__name__)


def build_binary_mask(seq_len: int, regions: list[list[int]]) -> str:
    """Convert 1-indexed [start, end] region pairs to a binary string mask."""
    mask = np.zeros(seq_len, dtype=int)
    for start, end in regions:
        if start >= 1 and end <= seq_len:
            mask[start - 1 : end] = 1
    return "".join(map(str, mask))


def parse_mobidb_entries(entries: list[dict]) -> pd.DataFrame:
    """
    Parse a list of MobiDB JSON entries into a DataFrame with binary masks.
    
    MobiDB v6 stores features as top-level keys in each entry, e.g.:
        entry["curated-disorder-merge"] = {"regions": [[233, 246]], "content_fraction": 0.057, ...}
    """
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

        # In MobiDB v6, each triplet is a top-level key in the entry dict
        for triplet in ALL_TARGET_TRIPLETS:
            col_name = f"mask_{triplet.replace('-', '_')}"
            feature_data = entry.get(triplet)
            
            if feature_data and "regions" in feature_data:
                row[col_name] = build_binary_mask(length, feature_data["regions"])
            else:
                row[col_name] = None

        rows.append(row)

    return pd.DataFrame(rows)


def build_region_index(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a flat DataFrame of continuous regions for evaluation.
    Each row is a single contiguous region of 1s in a mask.
    """
    regions = []
    
    for _, row in df.iterrows():
        uid = row["uniprot_id"]
        
        for triplet in ALL_TARGET_TRIPLETS:
            col_name = f"mask_{triplet.replace('-', '_')}"
            mask_str = row.get(col_name)
            
            if not isinstance(mask_str, str):
                continue
                
            mask = np.array([int(c) for c in mask_str])
            if mask.sum() == 0:
                continue
                
            # Find contiguous regions of 1s
            diffs = np.diff(np.concatenate([[0], mask, [0]]))
            starts = np.where(diffs == 1)[0]
            ends = np.where(diffs == -1)[0]
            
            for s, e in zip(starts, ends):
                regions.append({
                    "uniprot_id": uid,
                    "triplet": triplet,
                    "region_start": int(s),
                    "region_end": int(e),
                    "region_length": int(e - s),
                })
                
    return pd.DataFrame(regions)


def deduplicate_against(
    df: pd.DataFrame,
    other_parquet: Path,
) -> pd.DataFrame:
    """
    Drop rows whose ``uniprot_id`` appears in another dataset (e.g. the SAE
    training corpus). Prevents train/eval leakage where the SAE may have
    partially memorised an evaluation protein's embedding.

    Returns the deduplicated DataFrame. If ``other_parquet`` does not exist,
    returns ``df`` unchanged with a warning (caller should fetch the training
    set first for dedup to take effect).
    """
    other_parquet = Path(other_parquet)
    if not other_parquet.exists():
        logger.warning(
            f"  Dedup target not found: {other_parquet}. "
            f"Skipping dedup (fetch the training set first to enable it)."
        )
        return df

    other_ids = set(pd.read_parquet(other_parquet, columns=["uniprot_id"])["uniprot_id"])
    before = len(df)
    df = df[~df["uniprot_id"].isin(other_ids)].reset_index(drop=True)
    logger.info(
        f"  Deduplicated against {other_parquet.name}: "
        f"removed {before - len(df)} overlapping proteins ({before} → {len(df)})."
    )
    return df


def save_dataset(
    df: pd.DataFrame, 
    output_path: Path, 
    region_index: pd.DataFrame | None = None,
    region_index_path: Path | None = None
) -> None:
    """Save the dataset and region index as Parquet files."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path, engine="pyarrow", index=False)
    
    if region_index is not None and region_index_path is not None:
        region_index.to_parquet(region_index_path, engine="pyarrow", index=False)
