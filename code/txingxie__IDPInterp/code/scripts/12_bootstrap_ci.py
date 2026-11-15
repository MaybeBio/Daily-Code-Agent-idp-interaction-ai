#!/usr/bin/env python3
"""
Region-level bootstrap CI for a sub-state axis (default: molten globule).

Resamples the test-split annotated regions with replacement and recomputes the best
feature's F1 as the evaluator defines it: per-residue precision on the test split with
region-level recall on the test split, at the feature's validation-selected threshold.
Reports a percentile CI for each of the two runs (ref and cmp).

Usage:
  python scripts/12_bootstrap_ci.py paths=colab sae=topk_32x \
     +boot.ref_model=esm2_8m  +boot.ref_layer=4 \
     +boot.cmp_model=seqdance +boot.cmp_layer=10 \
     +boot.triplet=disprot-molten_globule \
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


def _run_name(model_name, sae_cfg, layer):
    tag = sae_cfg.get("tag", None)
    return f"{model_name}_{sae_cfg.type}_{sae_cfg.expansion_factor}x{'_'+tag if tag else ''}_layer{layer}"


def _load(cfg, model, layer, triplet, dataset, feat_override=None):
    """Return (run, act-cache mmap, feature_idx, threshold, reported_f1) for one
    model. If feat_override is given, use that feature instead of the best-F1 one."""
    run = _run_name(model, cfg.sae, layer)
    scratch = cfg.paths.get("scratch_dir", None) or cfg.paths.embedding_dir
    act_path = Path(scratch) / "act_cache" / f"{run}.mmap"
    if not act_path.with_suffix(".json").exists():
        raise FileNotFoundError(f"No activation cache for {run}: {act_path} "
                                f"(rebuild via scripts/04_evaluate_atlas.py).")
    ac = json.load(open(act_path.with_suffix(".json")))
    acts = np.memmap(act_path, dtype=np.float16, mode="r",
                     shape=(ac["total_residues"], ac["d_latent"]), order="C")
    meta = MetadataIndex.load(
        Path(cfg.paths.embedding_dir) / model / dataset / f"layer{layer}_metadata.parquet")
    csv = Path(cfg.paths.output_dir) / "atlas" / run / "triplets" / f"{triplet}.csv"
    tdf = pd.read_csv(csv)
    feat_col = next(c for c in ("f1", "test_f1", "best_f1") if c in tdf.columns)
    if feat_override is not None:
        row = tdf[tdf["feature_idx"] == int(feat_override)]
        if row.empty:
            raise ValueError(f"feature {feat_override} not in {csv}")
        best = row.iloc[0]
    else:
        best = tdf.sort_values(feat_col, ascending=False).iloc[0]
    feat = int(best["feature_idx"])
    thr = float(best["best_threshold"])
    return run, acts, meta, feat, thr, float(best[feat_col])


def _regions_and_fp(cfg, meta, acts, feat, thr, triplet, dataset):
    """Return {(uid, local_start, local_end): (region_max_act, region_TP_count)}
    over TEST regions for the feature (keyed by biological region so the two
    models can be paired), plus the fixed background FP count (firing test
    residues outside any region)."""
    gold = pd.read_parquet(Path(cfg.paths.dataset_dir) / f"{dataset}.parquet")
    col = f"mask_{triplet.replace('-', '_')}"
    if col not in gold.columns:
        raise KeyError(f"{col} not in {dataset}.parquet — build the axis first (08).")
    split_by = dict(zip(gold["uniprot_id"], gold["split"]))
    mask_by = dict(zip(gold["uniprot_id"], gold[col]))

    acts_col = np.asarray(acts[:, feat], dtype=np.float32)
    fire = acts_col >= thr

    regions = {}
    test_fire = 0      # all test residues firing
    test_pos_fire = 0  # test residues inside a region AND firing (= TP residues)
    for r in meta.itertuples():
        if split_by.get(r.uniprot_id) != "test":
            continue
        s0, L = int(r.mmap_start_idx), int(r.stored_length)
        seg_fire = fire[s0:s0 + L]
        test_fire += int(seg_fire.sum())
        m = mask_by.get(r.uniprot_id)
        if not isinstance(m, str):
            continue
        a = (np.frombuffer(m[:L].encode(), dtype=np.uint8) - ord("0")).astype(bool)
        if a.size < L:  # mask shorter than stored length: pad
            a = np.concatenate([a, np.zeros(L - a.size, bool)])
        seg_act = acts_col[s0:s0 + L]
        d = np.diff(np.concatenate([[0], a.astype(int), [0]]))
        for st, en in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
            sub = seg_act[st:en]
            regions[(r.uniprot_id, int(st), int(en))] = (
                float(sub.max()) if sub.size else 0.0, int((sub >= thr).sum()))
        test_pos_fire += int((seg_fire & a).sum())
    fp = test_fire - test_pos_fire  # firing test residues outside any region
    return regions, int(fp)


def _f1(reg_hit, reg_tp, fp):
    R = reg_hit.mean() if reg_hit.size else 0.0
    TP = reg_tp.sum()
    P = TP / (TP + fp + 1e-10)
    return 2 * P * R / (P + R + 1e-10), P, R


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s")
    b = cfg.boot
    dataset = b.get("dataset", "gold_curated")
    triplet = b.get("triplet", "disprot-molten_globule")
    n_boot = int(b.get("n_boot", 5000))
    rng = np.random.default_rng(int(b.get("seed", 0)))

    sides = {}
    for tag, model, layer, ov in (
        ("ref", b.ref_model, b.ref_layer, b.get("ref_feature", None)),
        ("cmp", b.cmp_model, b.cmp_layer, b.get("cmp_feature", None)),
    ):
        run, acts, meta, feat, thr, reported = _load(cfg, model, layer, triplet, dataset, ov)
        regions, fp = _regions_and_fp(cfg, meta, acts, feat, thr, triplet, dataset)
        sides[tag] = {"run": run, "feat": feat, "thr": thr, "reported": reported,
                      "regions": regions, "fp": fp}

    # Pair on the SAME biological regions (uid, start, end) present for both models.
    common = sorted(set(sides["ref"]["regions"]) & set(sides["cmp"]["regions"]))
    if len(common) < 5:
        raise RuntimeError(f"Only {len(common)} shared regions — too few to bootstrap.")

    def arrays(side):
        rm = np.array([side["regions"][k][0] for k in common])
        tp = np.array([side["regions"][k][1] for k in common], dtype=np.float64)
        return (rm >= side["thr"]).astype(float), tp, side["fp"]

    rh_r, tp_r, fp_r = arrays(sides["ref"])
    rh_c, tp_c, fp_c = arrays(sides["cmp"])
    n = len(common)

    f1_r, P_r, R_r = _f1(rh_r, tp_r, fp_r)
    f1_c, P_c, R_c = _f1(rh_c, tp_c, fp_c)

    # PAIRED bootstrap: resample region indices ONCE, apply to both models.
    gap_b = np.empty(n_boot)
    fr_b = np.empty(n_boot)
    fc_b = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        fr_b[i] = _f1(rh_r[idx], tp_r[idx], fp_r)[0]
        fc_b[i] = _f1(rh_c[idx], tp_c[idx], fp_c)[0]
        gap_b[i] = fc_b[i] - fr_b[i]
    glo, ghi = np.percentile(gap_b, [2.5, 97.5])
    p_gap_le0 = float((gap_b <= 0).mean())  # bootstrap one-sided p that cmp<=ref

    summary = {
        "triplet": triplet, "n_paired_regions": int(n),
        "ref": {"run": sides["ref"]["run"], "feature": sides["ref"]["feat"],
                "point_f1": round(float(f1_r), 4),
                "ci95": [round(float(x), 4) for x in np.percentile(fr_b, [2.5, 97.5])]},
        "cmp": {"run": sides["cmp"]["run"], "feature": sides["cmp"]["feat"],
                "point_f1": round(float(f1_c), 4),
                "ci95": [round(float(x), 4) for x in np.percentile(fc_b, [2.5, 97.5])]},
        "gap_cmp_minus_ref": round(float(f1_c - f1_r), 4),
        "gap_ci95": [round(float(glo), 4), round(float(ghi), 4)],
        "gap_excludes_0": bool(glo > 0),
        "boot_p_cmp_le_ref": round(p_gap_le0, 4),
    }
    logger.info(f"PAIRED gap (cmp-ref) = {f1_c - f1_r:+.3f}  95%CI[{glo:+.3f},{ghi:+.3f}]  "
                f"excludes 0: {glo > 0}  (boot p[cmp<=ref]={p_gap_le0:.4f}, n={n} paired regions)")

    out = Path(cfg.paths.output_dir) / "bootstrap" / f"{triplet}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(out, "w"), indent=2)
    logger.info("Summary:\n" + json.dumps(summary, indent=2))
    logger.info(f"Saved → {out}")


if __name__ == "__main__":
    main()
