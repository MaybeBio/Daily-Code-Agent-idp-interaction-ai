#!/usr/bin/env python3
"""
Confirm that every triplet in ALL_TARGET_TRIPLETS produced a usable per-residue mask
in the evaluation dataset.

The check runs against the built parquet rather than the MobiDB API: the fetch pulls
full protein records through the curated filters, and process.py parses every triplet
key out of each record, so the parquet is authoritative and the check is offline. A
triplet with zero proteins is either a bad key or absent from the corpus; evaluation
skips empty triplets.

Usage:
    python scripts/07_probe_triplets.py paths=colab
    python scripts/07_probe_triplets.py paths=colab probe_dataset=gold_curated
"""
import json
import logging
import sys
from pathlib import Path

import hydra
import pandas as pd
from omegaconf import DictConfig

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from idpinterp.constants import ALL_TARGET_TRIPLETS

logger = logging.getLogger(__name__)


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(name)s] %(levelname)s: %(message)s"
    )
    name = cfg.get("probe_dataset", "gold_curated")
    path = Path(cfg.paths.dataset_dir) / f"{name}.parquet"
    if not path.exists():
        logger.error(f"Dataset not found: {path}. Run scripts/01_fetch_mobidb.py first.")
        return

    df = pd.read_parquet(path)
    logger.info(f"Validating {len(df)} proteins in {path}")

    rows = []
    for triplet in ALL_TARGET_TRIPLETS:
        col = f"mask_{triplet.replace('-', '_')}"
        if col not in df.columns:
            rows.append({"triplet": triplet, "present": False,
                         "proteins": 0, "annotated_residues": 0})
            continue
        masks = df[col].dropna()
        n_res = int(sum(s.count("1") for s in masks))
        rows.append({
            "triplet": triplet,
            "present": True,
            "proteins": int(len(masks)),
            "annotated_residues": n_res,
        })

    empty = [r["triplet"] for r in rows if r["proteins"] == 0]
    logger.info(f"\n{'='*78}\nTRIPLET COVERAGE in {name} "
                f"({len(rows) - len(empty)}/{len(rows)} populated)\n{'='*78}")
    for r in rows:
        flag = "OK " if r["proteins"] > 0 else "!! "
        logger.info(f"  {flag}{r['triplet']:55s} proteins={r['proteins']:<6} "
                    f"residues={r['annotated_residues']}")
    if empty:
        logger.warning(
            "Empty triplets (bad key name or absent in this corpus; the evaluator will "
            "skip them):\n  " + "\n  ".join(empty)
        )

    out = Path(cfg.paths.output_dir) / "triplet_coverage.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(rows, f, indent=2)
    logger.info(f"Saved → {out}")


if __name__ == "__main__":
    main()
