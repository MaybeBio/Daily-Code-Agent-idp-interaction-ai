#!/usr/bin/env python3
"""
Describe features from the prompts saved by scripts/05
(outputs/atlas/<run>/description_prompts.json).

Held-out Pearson-r validation runs here as well: scripts/05 stores the held-out
validation examples with their measured peak activations inside the prompt file, so the
describe-predict-correlate loop needs no activation cache. Disable with
+describe_validate=false.

Usage:
    python scripts/10_describe_from_prompts.py paths=colab model=esm2_8m sae=topk_32x \
        +describe_max_calls=2000
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
from omegaconf import DictConfig

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
logger = logging.getLogger(__name__)


# ── Held-out validation helpers (mirror scripts/05; pure functions, CPU-only) ──
def _build_prediction_prompt(description: str, val_examples):
    shown = [{"uniprot_id": p["uniprot_id"], "sequence_window": p["sequence_window"]}
             for p in val_examples]
    body = "\n".join(f"  P{i+1}: {p['sequence_window']}" for i, p in enumerate(shown))
    prompt = (
        "A protein-language-model feature was described as follows:\n\n"
        f'"""{description}"""\n\n'
        "For each sequence window below (peak residue in [brackets]), predict this "
        "feature's activation as a number in [0, 1] (0=will not fire, 1=fires "
        "maximally). Use ONLY the description.\n\n"
        f"{body}\n\n"
        f"Respond with ONLY a JSON array of {len(shown)} numbers in order "
        f"P1..P{len(shown)}, e.g. [0.9, 0.1, ...]."
    )
    return prompt, shown


def _parse_predictions(text: str, n: int):
    try:
        s = text[text.index("["): text.rindex("]") + 1]
        arr = json.loads(s)
        if isinstance(arr, list) and len(arr) == n:
            return [float(x) for x in arr]
    except Exception:
        return None
    return None


def _pearson(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.size < 3 or np.std(a) < 1e-9 or np.std(b) < 1e-9:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])

try:
    from google import genai
    HAS_GEMINI = True
except ImportError:
    HAS_GEMINI = False
try:
    import anthropic
    HAS_ANTHROPIC = True
except ImportError:
    HAS_ANTHROPIC = False


def _run_name(cfg):
    layer = cfg.model.get("layer") or cfg.model.layers_to_extract[-1]
    tag = cfg.sae.get("tag", None)
    return f"{cfg.model.name}_{cfg.sae.type}_{cfg.sae.expansion_factor}x{'_'+tag if tag else ''}_layer{layer}"


def _make_caller(cfg):
    gkey, akey = os.environ.get("GEMINI_API_KEY"), os.environ.get("ANTHROPIC_API_KEY")
    model = cfg.get("describe_model", "gemini-3.1-flash-lite")
    sleep = float(cfg.get("describe_sleep", 4.1))
    if gkey and HAS_GEMINI:
        gc = genai.Client(api_key=gkey)
        if "gemini" not in model.lower():
            model = "gemini-3.1-flash-lite"
        def call(p):
            time.sleep(sleep)
            try:
                return gc.models.generate_content(model=model, contents=p).text
            except Exception as e:
                return f"[API error: {e}]"
        return call, model
    if akey and HAS_ANTHROPIC:
        ac = anthropic.Anthropic(api_key=akey)
        if "claude" not in model.lower():
            model = "claude-sonnet-4-6"
        mt = int(cfg.get("describe_max_tokens", 600))
        def call(p):
            try:
                return ac.messages.create(model=model, max_tokens=mt,
                    messages=[{"role": "user", "content": p}]).content[0].text
            except Exception as e:
                return f"[API error: {e}]"
        return call, model
    return None, None


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s")
    run = _run_name(cfg)
    atlas_dir = Path(cfg.paths.output_dir) / "atlas" / run
    prompts_path = atlas_dir / "description_prompts.json"
    if not prompts_path.exists():
        logger.error(f"No prompts at {prompts_path}. Run scripts/05 once (GPU) to emit them.")
        return
    prompts = json.load(open(prompts_path))
    feat_list = cfg.get("describe_features", None)
    if feat_list:
        ids = {int(x) for x in feat_list}
        prompts = [p for p in prompts if int(p["feature_idx"]) in ids]
        logger.info(f"Restricted to {len(prompts)} explicitly-requested features.")
    atlas = pd.read_parquet(atlas_dir / "feature_atlas.parquet").set_index("feature_idx")

    call, model = _make_caller(cfg)
    if call is None:
        logger.error("No API key (set GEMINI_API_KEY or ANTHROPIC_API_KEY).")
        return
    logger.info(f"LLM: {model} | {len(prompts)} prompts available")

    out_path = atlas_dir / "feature_descriptions.parquet"
    cur_ids = {int(p["feature_idx"]) for p in prompts}
    done, prev_keep = {}, None
    if out_path.exists():
        prev = pd.read_parquet(out_path)
        for r in prev.itertuples():
            d = getattr(r, "description", None)
            if isinstance(d, str) and d and not d.startswith("[API error"):
                done[int(r.feature_idx)] = (d, float(getattr(r, "pearson_r", float("nan"))))
        # ADDITIVE: rows for features outside this prompt set are preserved as-is,
        # so a targeted pass never clobbers earlier descriptions.
        prev_keep = prev[~prev["feature_idx"].isin(cur_ids)].copy()
        logger.info(f"Resuming — {len(done)} described; preserving {len(prev_keep)} "
                    f"prior rows outside this prompt set.")

    def _write(rows):
        df = pd.DataFrame(rows)
        if prev_keep is not None and len(prev_keep):
            df = pd.concat([prev_keep, df], ignore_index=True)
        if "f1" in df.columns:
            df = df.sort_values("f1", ascending=False)
        df.to_parquet(out_path, index=False)
        return df

    validate = bool(cfg.get("describe_validate", True))
    min_high = int(cfg.get("describe_min_high", 6))
    max_calls = int(cfg.get("describe_max_calls", 2000))
    rows, n_new = [], 0
    for item in prompts:
        fi = int(item["feature_idx"])
        pearson_r = float("nan")
        if fi in done:
            desc, pearson_r = done[fi]
        elif n_new < max_calls:
            desc = call(item["prompt"]); n_new += 1
            # Held-out validation: predict activations from the description alone
            # on disjoint val examples, correlate with the saved measured peaks.
            val = item.get("val") or []
            n_high = int(item.get("n_high", 0))
            if (validate and isinstance(desc, str)
                    and not desc.startswith("[API error")
                    and n_high >= min_high and len(val) >= 3 and n_new < max_calls):
                pred_prompt, shown = _build_prediction_prompt(desc, val)
                preds = _parse_predictions(call(pred_prompt), len(shown))
                n_new += 1
                if preds is not None:
                    pearson_r = _pearson(preds, [p["peak_activation"] for p in val])
        else:
            desc = ""
        rows.append({
            "feature_idx": fi,
            "triplet": atlas.loc[fi, "triplet"] if fi in atlas.index else "",
            "f1": float(atlas.loc[fi, "f1"]) if fi in atlas.index else 0.0,
            "description": desc,
            "pearson_r": pearson_r,
        })
        if n_new and n_new % 25 == 0:
            _write(rows)  # incremental checkpoint (merged with preserved prior rows)

    df = _write(rows)
    df.to_csv(atlas_dir / "feature_descriptions.csv", index=False)
    valid_r = df["pearson_r"].dropna()
    valid_r = valid_r[np.isfinite(valid_r)]
    med = float(valid_r.median()) if len(valid_r) else float("nan")
    logger.info(f"Described {n_new} new this run; {len(df)} total in table → {out_path}")
    logger.info(f"Held-out validation: median Pearson r = {med:.3f} over {len(valid_r)} "
                f"features (atlas-wide)")


if __name__ == "__main__":
    main()
