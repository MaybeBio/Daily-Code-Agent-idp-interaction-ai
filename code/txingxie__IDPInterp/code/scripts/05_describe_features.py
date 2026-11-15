#!/usr/bin/env python3
"""
Generate and validate natural-language descriptions for SAE features.

Examples span the full activation range (deciles plus zero-activation negatives). For
validation, the model predicts activation on a disjoint set of proteins from its own
description, and we report the Pearson correlation against measured peak activation.

Usage:
    python scripts/05_describe_features.py paths=colab model=esm2_8m sae=l1_32x \
        +describe_model=claude-sonnet-4-6 +describe_max_calls=100
"""

import json
import logging
import os
import sys
import time
from pathlib import Path

import hydra
import numpy as np
import pandas as pd
import torch
from omegaconf import DictConfig
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from idpinterp.models.memmap_store import MetadataIndex

logger = logging.getLogger(__name__)

try:
    import anthropic
    HAS_ANTHROPIC = True
except ImportError:
    HAS_ANTHROPIC = False

try:
    from google import genai
    HAS_GEMINI = True
except ImportError:
    HAS_GEMINI = False


# ── Example selection (decile-spanning) ──────────────────────────────────────

def _protein_peaks(feat_acts: np.ndarray, starts: np.ndarray) -> np.ndarray:
    """Per-protein max activation for one feature (segments defined by starts)."""
    return np.maximum.reduceat(feat_acts, starts)


def _make_example(
    prot_row: pd.Series,
    feat_acts: np.ndarray,
    uid_to_row: dict,
    triplet: str,
    window_radius: int,
) -> dict | None:
    """Build a marked sequence window around a protein's peak residue."""
    uid = prot_row["uniprot_id"]
    if uid not in uid_to_row:
        return None
    start = int(prot_row["mmap_start_idx"])
    stored_len = int(prot_row["stored_length"])
    seg = feat_acts[start : start + stored_len]
    if seg.size == 0:
        return None
    local_idx = int(np.argmax(seg))
    peak = float(seg[local_idx])

    seq = uid_to_row[uid]["sequence"]
    ws = max(0, local_idx - window_radius)
    we = min(stored_len, local_idx + window_radius + 1)
    window = seq[ws:we]
    p = local_idx - ws
    marked = window[:p] + f"[{window[p]}]" + window[p + 1:]

    mask_col = f"mask_{triplet.replace('-', '_')}"
    mask_str = uid_to_row[uid].get(mask_col)
    annotation = "unknown"
    if isinstance(mask_str, str) and local_idx < len(mask_str):
        annotation = "True" if mask_str[local_idx] == "1" else "False"

    return {
        "uniprot_id": uid,
        "sequence_window": marked,
        "peak_activation": peak,
        "annotation": annotation,
    }


def _select_decile_examples(
    feat_idx: int,
    all_activations: np.ndarray,
    meta_df: pd.DataFrame,
    starts: np.ndarray,
    uid_to_row: dict,
    triplet: str,
    n_per_bin: int = 2,
    top_bin_n: int = 8,
    n_negatives: int = 8,
    window_radius: int = 15,
    seed: int = 0,
) -> tuple[list[dict], list[dict], int]:
    """
    Returns (generation_examples, validation_examples, n_high) where examples
    span activation deciles plus zero-activation negatives. The two lists are
    disjoint sets of proteins (validation is held out for Pearson r).
    n_high = number of proteins reaching the top-3 activation bins (>=0.7).
    """
    rng = np.random.default_rng(seed + feat_idx)
    feat_acts = all_activations[:, feat_idx].astype(np.float32)
    peaks = _protein_peaks(feat_acts, starts)            # (n_proteins,) aligned to meta_df

    gen, val = [], []
    n_high = int((peaks >= 0.7).sum())

    def add_from(idxs: np.ndarray, k: int):
        rng.shuffle(idxs)
        added = 0
        for pi in idxs:
            ex = _make_example(meta_df.iloc[int(pi)], feat_acts, uid_to_row,
                               triplet, window_radius)
            if ex is None:
                continue
            (gen if (added % 2 == 0) else val).append(ex)
            added += 1
            if added >= k:
                break

    for b in range(10):
        lo, hi = b * 0.1, (b + 1) * 0.1 + (1e-6 if b == 9 else 0.0)
        bin_idx = np.nonzero((peaks >= lo) & (peaks < hi))[0]
        if bin_idx.size == 0:
            continue
        k = top_bin_n if b == 9 else n_per_bin
        add_from(bin_idx, k)

    # Zero-activation negatives
    zero_idx = np.nonzero(peaks < 1e-6)[0]
    if zero_idx.size:
        add_from(zero_idx, n_negatives)

    return gen, val, n_high


# ── Prompts ──────────────────────────────────────────────────────────────────

def _build_description_prompt(feat_idx, triplet, f1, examples) -> str:
    lines = []
    for i, p in enumerate(sorted(examples, key=lambda x: -x["peak_activation"])):
        lines.append(
            f"  P{i+1} ({p['uniprot_id']}): activation={p['peak_activation']:.2f} "
            f"| MobiDB[{triplet}]={p['annotation']} | window: {p['sequence_window']}"
        )
    body = "\n".join(lines)
    return f"""You are analyzing latent Feature #{feat_idx} from a Sparse Autoencoder trained on a protein language model (ESM-2), studied for how it represents intrinsically disordered protein regions. It best matches the MobiDB annotation "{triplet}" (region-aware F1={f1:.3f}).

Below are proteins spanning the full activation range (0=never fires, 1=max), with a ±15-residue window around each peak; the peak residue is in [brackets]:

{body}

Using the contrast between high- and low-activation examples:
1. MOTIF: What specific amino-acid composition, charge pattern, or sequence motif does this feature detect? Be precise.
2. SUMMARY: One sentence that would let someone predict this feature's activation (0-1) on a new sequence window.
3. SIGNIFICANCE: Likely structural/functional role in disorder biology (e.g. linker, LLPS, MoRF, amyloid, NLS)."""


def _build_prediction_prompt(description: str, val_examples) -> tuple[str, list[dict]]:
    shown = [{"uniprot_id": p["uniprot_id"], "sequence_window": p["sequence_window"]}
             for p in val_examples]
    lines = [f"  P{i+1}: {p['sequence_window']}" for i, p in enumerate(shown)]
    body = "\n".join(lines)
    prompt = f"""A protein-language-model feature was described as follows:

\"\"\"{description}\"\"\"

For each sequence window below (peak residue in [brackets]), predict this feature's activation as a number in [0, 1] (0=will not fire, 1=fires maximally). Use ONLY the description.

{body}

Respond with ONLY a JSON array of {len(shown)} numbers in order P1..P{len(shown)}, e.g. [0.9, 0.1, ...]."""
    return prompt, shown


def _parse_predictions(text: str, n: int) -> list[float] | None:
    try:
        s = text[text.index("["): text.rindex("]") + 1]
        arr = json.loads(s)
        if isinstance(arr, list) and len(arr) == n:
            return [float(x) for x in arr]
    except Exception:
        return None
    return None


def _pearson(a: list[float], b: list[float]) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.size < 3 or np.std(a) < 1e-9 or np.std(b) < 1e-9:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


# ── LLM call wrappers ────────────────────────────────────────────────────────

def _make_caller(cfg, anthropic_key, gemini_key):
    """Return (call_fn, model_name). Prefers Anthropic (capable default)."""
    model = cfg.get("describe_model", "claude-sonnet-4-6")
    max_tokens = int(cfg.get("describe_max_tokens", 600))
    sleep = float(cfg.get("describe_sleep", 0.0))

    if anthropic_key and HAS_ANTHROPIC:
        client = anthropic.Anthropic(api_key=anthropic_key)

        def call(prompt: str) -> str:
            for attempt in range(4):
                try:
                    if sleep:
                        time.sleep(sleep)
                    r = client.messages.create(
                        model=model, max_tokens=max_tokens,
                        messages=[{"role": "user", "content": prompt}],
                    )
                    return r.content[0].text
                except Exception as e:
                    if attempt == 3:
                        return f"[API error: {e}]"
                    time.sleep(2 ** attempt)
        return call, model

    if gemini_key and HAS_GEMINI:
        gclient = genai.Client(api_key=gemini_key)
        gmodel = cfg.get("describe_model", "gemini-3.1-flash-lite")
        # Guard: if describe_model is a non-Gemini id (e.g. the Claude default),
        # don't send it to the Gemini endpoint — fall back to a Gemini model.
        if "gemini" not in gmodel.lower():
            gmodel = "gemini-3.1-flash-lite"

        def call(prompt: str) -> str:
            time.sleep(max(sleep, 4.1))  # free-tier 15 RPM
            try:
                return gclient.models.generate_content(model=gmodel, contents=prompt).text
            except Exception as e:
                return f"[API error: {e}]"
        return call, gmodel

    return None, None


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(name)s] %(levelname)s: %(message)s"
    )

    layer = cfg.model.get("layer") or cfg.model.layers_to_extract[-1]
    _tag = cfg.sae.get("tag", None)
    _tag_suffix = f"_{_tag}" if _tag else ""
    run_name = f"{cfg.model.name}_{cfg.sae.type}_{cfg.sae.expansion_factor}x{_tag_suffix}_layer{layer}"
    logger.info(f"=== Feature description [{run_name}] ===")

    anthropic_key = os.environ.get("ANTHROPIC_API_KEY")
    gemini_key = os.environ.get("GEMINI_API_KEY")
    call, model = _make_caller(cfg, anthropic_key, gemini_key)
    if call is None:
        logger.warning("No API key/client available — will write prompts only.")
    else:
        logger.info(f"LLM: {model}")

    max_calls = int(cfg.get("describe_max_calls", 50))
    validate = bool(cfg.get("describe_validate", True))
    min_high = int(cfg.get("describe_min_high", 6))
    logger.info(f"API budget: {max_calls} calls | validation={validate}")

    atlas_dir = Path(cfg.paths.output_dir) / "atlas" / run_name
    atlas_path = atlas_dir / "feature_atlas.parquet"
    if not atlas_path.exists():
        logger.error(f"Atlas not found: {atlas_path}. Run scripts/04_evaluate_atlas.py first.")
        return
    atlas = pd.read_parquet(atlas_path)

    min_f1 = float(cfg.get("describe_min_f1", 0.3))
    high = atlas[atlas["f1"] >= min_f1].sort_values("f1", ascending=False)
    # Optional: restrict to features whose best axis is a specific triplet
    # (e.g. +describe_triplet=disprot-molten_globule to describe just the MG family).
    only_triplet = cfg.get("describe_triplet", None)
    if only_triplet:
        high = high[high["triplet"] == only_triplet]
    # Override: describe an explicit feature-id list regardless of f1/triplet
    # (e.g. +describe_features=[12669,5229,3966,14136] for the MG-unique family).
    feat_list = cfg.get("describe_features", None)
    if feat_list:
        ids = [int(x) for x in feat_list]
        high = atlas[atlas["feature_idx"].isin(ids)].sort_values("f1", ascending=False)
    logger.info(f"Features to describe: {len(high)}"
                f"{f' (explicit id list)' if feat_list else f' with F1>={min_f1}'}"
                f"{f' on {only_triplet}' if only_triplet and not feat_list else ''}")
    if high.empty:
        return

    scratch_root = cfg.paths.get("scratch_dir", None) or cfg.paths.output_dir
    act_cache_dir = Path(scratch_root) / "act_cache"
    act_mmap_path = act_cache_dir / f"{run_name}.mmap"
    act_cfg_path = act_cache_dir / f"{run_name}.json"
    if not act_mmap_path.exists():
        logger.error(f"Activation cache not found: {act_mmap_path}. Run scripts/04_evaluate_atlas.py first.")
        return
    with open(act_cfg_path) as f:
        act_cfg = json.load(f)
    all_activations = np.memmap(
        act_mmap_path, dtype=np.float16, mode="r",
        shape=(act_cfg["total_residues"], act_cfg["d_latent"]), order="C",
    )

    embed_dir = Path(cfg.paths.embedding_dir) / cfg.model.name / cfg.dataset.name
    meta_df = MetadataIndex.load(embed_dir / f"layer{layer}_metadata.parquet").reset_index(drop=True)
    starts = meta_df["mmap_start_idx"].values.astype(np.int64)
    eval_path = Path(cfg.paths.dataset_dir) / "gold_curated.parquet"
    eval_df = pd.read_parquet(eval_path) if eval_path.exists() else pd.DataFrame()
    uid_to_row = (
        eval_df.drop_duplicates("uniprot_id").set_index("uniprot_id").to_dict(orient="index")
        if not eval_df.empty else {}
    )

    logger.info("Pre-warming mmap into OS page cache (sequential scan)...")
    _warm = all_activations.sum(axis=1)
    del _warm

    descriptions, prompts_for_review = [], []
    calls_made = 0

    for _, row in tqdm(high.iterrows(), total=len(high), desc="Describing", mininterval=30):
        feat_idx, triplet, f1 = int(row["feature_idx"]), row["triplet"], float(row["f1"])
        gen, val, n_high = _select_decile_examples(
            feat_idx, all_activations, meta_df, starts, uid_to_row, triplet,
            n_per_bin=int(cfg.get("describe_n_per_bin", 2)),
            top_bin_n=int(cfg.get("describe_top_bin_n", 8)),
            n_negatives=int(cfg.get("describe_n_negatives", 8)),
        )
        if len(gen) < 3:
            continue

        desc_prompt = _build_description_prompt(feat_idx, triplet, f1, gen)
        # Save the held-out validation examples (with measured peak activations)
        # alongside the prompt so the Pearson-r validation can run later on a
        # CPU-only runtime (scripts/10) without the activation cache.
        prompts_for_review.append({
            "feature_idx": feat_idx, "prompt": desc_prompt,
            "val": val, "n_high": int(n_high),
        })

        description, pearson_r = "", float("nan")
        if call is not None and calls_made < max_calls:
            description = call(desc_prompt)
            calls_made += 1

            if validate and n_high >= min_high and len(val) >= 3 and calls_made < max_calls:
                pred_prompt, shown = _build_prediction_prompt(description, val)
                pred_text = call(pred_prompt)
                calls_made += 1
                preds = _parse_predictions(pred_text, len(shown))
                if preds is not None:
                    measured = [p["peak_activation"] for p in val]
                    pearson_r = _pearson(preds, measured)

        descriptions.append({
            "feature_idx": feat_idx, "triplet": triplet, "f1": f1,
            "description": description, "pearson_r": pearson_r,
            "n_examples_gen": len(gen), "n_examples_val": len(val), "n_high": n_high,
            "top_proteins": json.dumps([p["uniprot_id"] for p in gen[:5]]),
        })

    desc_df = pd.DataFrame(descriptions)
    atlas_dir.mkdir(parents=True, exist_ok=True)
    # Always (re)write the prompts — that's the emit artifact.
    with open(atlas_dir / "description_prompts.json", "w") as f:
        json.dump(prompts_for_review, f, indent=2)
    # Only (over)write the descriptions file if we ACTUALLY described something.
    # An emit-only pass (max_calls=0) must NOT clobber an existing, populated
    # feature_descriptions.parquet — scripts/10 owns incremental description writes.
    if calls_made > 0:
        desc_df.to_parquet(atlas_dir / "feature_descriptions.parquet", index=False)
        desc_df.to_csv(atlas_dir / "feature_descriptions.csv", index=False)
    else:
        logger.info("Emit-only (0 calls): wrote description_prompts.json; left "
                    "feature_descriptions.parquet untouched.")

    valid_r = desc_df["pearson_r"].dropna()
    valid_r = valid_r[np.isfinite(valid_r)]
    if len(valid_r):
        logger.info(
            f"Description validation: median Pearson r = {valid_r.median():.3f} "
            f"over {len(valid_r)} features (mean={valid_r.mean():.3f})."
        )
    logger.info(f"=== Description complete: {len(descriptions)} features, {calls_made} API calls ===")


if __name__ == "__main__":
    main()
