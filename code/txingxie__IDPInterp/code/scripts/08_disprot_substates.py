#!/usr/bin/env python3
"""
Add DisProt experimental disorder sub-states (molten globule, pre-molten globule,
MoRF / disorder-to-order) as evaluation axes.

Two modes:
  probe (default) - fetch DisProt and report per-state protein and region counts, both
                    database-wide and within the gold_curated corpus.
  build           - build per-residue masks for the viable states and merge them into
                    gold_curated.parquet (and _regions.parquet) as additive triplets
                    named `disprot-<state>`. Add those names to ALL_TARGET_TRIPLETS and
                    re-run the atlas evaluation to score them.

Usage:
    python scripts/08_disprot_substates.py paths=colab +disprot_mode=probe
    python scripts/08_disprot_substates.py paths=colab +disprot_mode=build
"""
import json
import logging
import sys
from pathlib import Path

import hydra
import pandas as pd
from omegaconf import DictConfig

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from idpinterp.data.disprot import (
    DEFAULT_TARGET_STATES,
    fetch_released,
    probe_counts,
    region_table,
)
from idpinterp.data.process import build_binary_mask

logger = logging.getLogger(__name__)

# Sparse conformational states worth pulling extra proteins for (MoRF/d2o is
# already well-powered DB-wide, so we don't bloat the corpus with its 500+).
EXPAND_STATES = ["molten_globule", "pre_molten_globule", "order_to_disorder"]


def _union_disprot_proteins(eval_df, entries, regions, states, train_ids, seed):
    """
    Add DisProt proteins carrying any of `states` that are NOT already in the
    eval corpus (and not in the training set, to stay leakage-free) as new rows.
    New rows carry the DisProt sequence; MobiDB mask columns are left null (they
    contribute only true-negatives to the MobiDB axes). DisProt masks are built
    later by the normal build loop once these accs are in eval_accs.
    Splits are assigned to the NEW rows only (existing splits untouched).
    """
    needles = [n for s in states for n in DEFAULT_TARGET_STATES[s]]
    hit = regions[regions["term_name"].apply(lambda t: any(n in t for n in needles))]
    have = set(eval_df["uniprot_id"])
    candidate = set(hit["acc"]) - have - set(train_ids)
    seq_by_acc = {e.get("acc"): e.get("sequence", "") for e in entries}

    new_rows = []
    for acc in sorted(candidate):
        seq = seq_by_acc.get(acc)
        if not isinstance(seq, str) or not seq:
            continue
        new_rows.append({"uniprot_id": acc, "sequence": seq, "length": len(seq)})
    if not new_rows:
        return eval_df, 0

    new_df = pd.DataFrame(new_rows)
    for col in eval_df.columns:
        if col not in new_df.columns:
            new_df[col] = None
    rng = np.random.default_rng(seed)
    new_df["split"] = np.where(rng.random(len(new_df)) < 0.5, "test", "val")
    out = pd.concat([eval_df, new_df[eval_df.columns]], ignore_index=True)
    return out, len(new_df)


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(name)s] %(levelname)s: %(message)s"
    )
    mode = cfg.get("disprot_mode", "probe")
    dataset_dir = Path(cfg.paths.dataset_dir)
    gold_path = dataset_dir / "gold_curated.parquet"
    regions_path = dataset_dir / "gold_curated_regions.parquet"

    logger.info("Fetching DisProt released set...")
    entries = fetch_released(cfg.get("disprot_release", "current"))
    regions = region_table(entries)

    eval_df = pd.read_parquet(gold_path) if gold_path.exists() else None
    eval_accs = set(eval_df["uniprot_id"]) if eval_df is not None else None

    if mode == "probe":
        counts = probe_counts(regions, eval_accs)
        logger.info("DisProt sub-state counts:\n" + json.dumps(counts, indent=2))
        out = Path(cfg.paths.output_dir) / "disprot_substate_counts.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as f:
            json.dump(counts, f, indent=2)
        logger.info(f"Saved → {out}")
        if eval_accs is not None:
            viable = [s for s, r in counts.items() if r.get("viable")]
            logger.info(f"Viable in-corpus states (>=50 regions): {viable or 'NONE'}")
        return

    # ── build ────────────────────────────────────────────────────────────────
    if eval_df is None:
        logger.error("gold_curated.parquet not found; fetch the eval corpus first.")
        return

    # Opt-in: union out-of-corpus DisProt conformational-state proteins so the
    # molten-globule/pre-molten axes aren't capped at whatever overlaps MobiDB
    # curated disorder. Dedup against the SAE training corpus to stay leakage-free.
    if cfg.get("disprot_expand_corpus", False):
        train_name = cfg.get("disprot_train_dedup", "full_swissprot")
        train_path = dataset_dir / f"{train_name}.parquet"
        train_ids = (
            set(pd.read_parquet(train_path, columns=["uniprot_id"])["uniprot_id"])
            if train_path.exists() else set()
        )
        if not train_ids:
            logger.warning(
                f"  Training set {train_path} not found — corpus union will NOT "
                f"dedup against training (leakage risk). Fetch it first."
            )
        before = len(eval_df)
        eval_df, n_added = _union_disprot_proteins(
            eval_df, entries, regions, EXPAND_STATES, train_ids, cfg.seed
        )
        eval_accs = set(eval_df["uniprot_id"])
        logger.info(
            f"  Corpus union: added {n_added} out-of-corpus DisProt proteins "
            f"({before} → {len(eval_df)}). NOTE: regenerate gold_curated embeddings "
            f"+ act caches before re-evaluating."
        )

    target_states = DEFAULT_TARGET_STATES
    built_triplets = []
    new_region_rows = []
    seq_lookup = eval_df.set_index("uniprot_id")["sequence"].to_dict()

    for state, needles in target_states.items():
        triplet = f"disprot-{state}"
        col = f"mask_{triplet.replace('-', '_')}"
        state_regs = regions[regions["term_name"].apply(
            lambda t: any(n in t for n in needles))]
        state_regs = state_regs[state_regs["acc"].isin(eval_accs)]
        if state_regs.empty:
            logger.info(f"  {state}: no in-corpus regions, skipping.")
            continue

        masks = {}
        for acc, grp in state_regs.groupby("acc"):
            seq = seq_lookup.get(acc)
            if not isinstance(seq, str):
                continue
            region_pairs = [[int(r.start), int(r.end)] for r in grp.itertuples()]
            masks[acc] = build_binary_mask(len(seq), region_pairs)
            # region index rows (0-indexed contiguous spans within the mask)
            mask_arr = masks[acc]
            import numpy as np
            m = np.array([int(c) for c in mask_arr])
            diffs = np.diff(np.concatenate([[0], m, [0]]))
            for s, e in zip(np.where(diffs == 1)[0], np.where(diffs == -1)[0]):
                new_region_rows.append({
                    "uniprot_id": acc, "triplet": triplet,
                    "region_start": int(s), "region_end": int(e),
                    "region_length": int(e - s),
                })

        eval_df[col] = eval_df["uniprot_id"].map(masks)
        built_triplets.append(triplet)
        logger.info(f"  {state}: built masks for {len(masks)} proteins → {col}")

    if not built_triplets:
        logger.warning("No DisProt states had in-corpus data; nothing built.")
        return

    eval_df.to_parquet(gold_path, index=False)
    if regions_path.exists() and new_region_rows:
        old = pd.read_parquet(regions_path)
        old = old[~old["triplet"].isin(built_triplets)]  # idempotent re-build
        pd.concat([old, pd.DataFrame(new_region_rows)], ignore_index=True).to_parquet(
            regions_path, index=False)

    triplets_file = dataset_dir / "disprot_triplets.txt"
    triplets_file.write_text("\n".join(built_triplets) + "\n")
    logger.info(
        f"Built DisProt triplets: {built_triplets}\n"
        f"→ Add these to ALL_TARGET_TRIPLETS (idpinterp/constants.py) and re-run "
        f"the atlas evaluation to score them. List saved to {triplets_file}."
    )


if __name__ == "__main__":
    main()
