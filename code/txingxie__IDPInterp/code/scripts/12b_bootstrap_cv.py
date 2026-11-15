#!/usr/bin/env python3
"""
K-fold cross-validated region-level bootstrap CI for a sub-state axis.

Usage:
  python scripts/12b_bootstrap_cv.py paths=colab sae=topk_32x \
     +boot.ref_model=esm2_8m  +boot.ref_layer=4  +boot.ref_feature=2594 \
     +boot.cmp_model=seqdance +boot.cmp_layer=10 +boot.cmp_feature=188 \
     +boot.triplet=disprot-molten_globule +boot.k_folds=5 \
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


def _load(cfg, model, layer, triplet, dataset, feat_override):
    """Return (run, act mmap, meta, feature_idx) for one model. Feature must be pinned."""
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
    if feat_override is None:
        raise ValueError(f"{model}: a pinned feature is required "
                         f"(+boot.{'ref' if 'ref' else 'cmp'}_feature=...). CV validates a "
                         f"fixed feature; it must not re-select.")
    return run, acts, meta, int(feat_override)


def _collect(cfg, meta, acts, feat, triplet, dataset):
    """"""
    gold = pd.read_parquet(Path(cfg.paths.dataset_dir) / f"{dataset}.parquet")
    col = f"mask_{triplet.replace('-', '_')}"
    if col not in gold.columns:
        raise KeyError(f"{col} not in {dataset}.parquet — build the axis first (08).")
    mask_by = dict(zip(gold["uniprot_id"], gold[col]))

    acts_col = np.asarray(acts[:, feat], dtype=np.float32)

    regions = {}
    bg_by_uid = {}
    for r in meta.itertuples():
        uid = r.uniprot_id
        s0, L = int(r.mmap_start_idx), int(r.stored_length)
        seg_act = acts_col[s0:s0 + L]
        m = mask_by.get(uid)
        if not isinstance(m, str):
            # no annotation string → entire protein is background
            bg_by_uid[uid] = np.concatenate([bg_by_uid.get(uid, np.empty(0, np.float32)), seg_act])
            continue
        a = (np.frombuffer(m[:L].encode(), dtype=np.uint8) - ord("0")).astype(bool)
        if a.size < L:
            a = np.concatenate([a, np.zeros(L - a.size, bool)])
        # region residues → per-region arrays
        d = np.diff(np.concatenate([[0], a.astype(int), [0]]))
        for st, en in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
            regions[(uid, int(st), int(en))] = seg_act[st:en].astype(np.float32)
        # non-region residues of this protein → background
        bg = seg_act[~a]
        bg_by_uid[uid] = np.concatenate([bg_by_uid.get(uid, np.empty(0, np.float32)), bg])
    return regions, bg_by_uid


def _sorted(a):
    return np.sort(np.asarray(a, dtype=np.float32))


def _ge_count(sorted_arr, thr):
    """# elements >= thr in an ascending-sorted array."""
    return int(sorted_arr.size - np.searchsorted(sorted_arr, thr, side="left"))


def _fit_threshold(train_regions, train_bg_sorted):
    """Pick the threshold that maximises region-F1 on the training folds — the same
    objective the split-aware evaluator uses on its val split. Candidates = distinct
    region-max activations (a threshold between region maxes never changes recall)."""
    maxes = np.array([r.max() if r.size else 0.0 for r in train_regions], dtype=np.float32)
    n_reg = len(train_regions)
    cands = np.unique(maxes[maxes > 1e-6])
    if cands.size == 0:
        return None
    fp_total = train_bg_sorted  # sorted
    best_thr, best_f1 = None, -1.0
    for thr in cands:
        recall = float((maxes >= thr).mean())
        tp = sum(_ge_count(_sorted(r), thr) for r in train_regions)
        fp = _ge_count(fp_total, thr)
        prec = tp / (tp + fp + 1e-10)
        f1 = 2 * prec * recall / (prec + recall + 1e-10)
        if f1 > best_f1:
            best_thr, best_f1 = float(thr), f1
    return best_thr


def _oof_score(cfg, side, common, fold_of, k, rng):
    """K-fold out-of-fold scoring for one model. Returns per-region (hit, tp) arrays
    aligned to `common`, plus the pooled out-of-fold background FP count."""
    regions, bg_by_uid = side["regions"], side["bg_by_uid"]
    reg_sorted = {key: _sorted(regions[key]) for key in common}
    reg_max = {key: (regions[key].max() if regions[key].size else 0.0) for key in common}

    hit = np.zeros(len(common))
    tp = np.zeros(len(common))
    fp_total = 0
    for f in range(k):
        heldout_uids = {u for u, fo in fold_of.items() if fo == f}
        # training regions = those whose protein is NOT held out
        train_keys = [key for key in common if key[0] not in heldout_uids]
        train_regions = [regions[key] for key in train_keys]
        train_bg = np.concatenate(
            [bg_by_uid[u] for u, fo in fold_of.items() if fo != f and u in bg_by_uid]
            or [np.empty(0, np.float32)])
        thr = _fit_threshold(train_regions, _sorted(train_bg))
        if thr is None:
            thr = np.inf  # degenerate fold: nothing fires → no detections
        # score held-out regions
        for i, key in enumerate(common):
            if key[0] in heldout_uids:
                hit[i] = 1.0 if reg_max[key] >= thr else 0.0
                tp[i] = _ge_count(reg_sorted[key], thr)
        # held-out background FP at this fold's threshold
        held_bg = np.concatenate(
            [bg_by_uid[u] for u in heldout_uids if u in bg_by_uid]
            or [np.empty(0, np.float32)])
        fp_total += _ge_count(_sorted(held_bg), thr)
    return hit, tp, int(fp_total)


def _f1_from(hit, tp, fp):
    R = hit.mean() if hit.size else 0.0
    TP = tp.sum()
    P = TP / (TP + fp + 1e-10)
    return 2 * P * R / (P + R + 1e-10)


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s")
    b = cfg.boot
    dataset = b.get("dataset", "gold_curated")
    triplet = b.get("triplet", "disprot-molten_globule")
    n_boot = int(b.get("n_boot", 5000))
    k = int(b.get("k_folds", 5))
    rng = np.random.default_rng(int(b.get("seed", 0)))

    sides = {}
    for tag, model, layer, ov in (
        ("ref", b.ref_model, b.ref_layer, b.get("ref_feature", None)),
        ("cmp", b.cmp_model, b.cmp_layer, b.get("cmp_feature", None)),
    ):
        run, acts, meta, feat = _load(cfg, model, layer, triplet, dataset, ov)
        regions, bg = _collect(cfg, meta, acts, feat, triplet, dataset)
        sides[tag] = {"run": run, "feat": feat, "regions": regions, "bg_by_uid": bg}
        logger.info(f"{tag}: {run} feat {feat} — {len(regions)} regions (all splits), "
                    f"{len(bg)} proteins with background.")

    # Pair on biological regions present for both models.
    common = sorted(set(sides["ref"]["regions"]) & set(sides["cmp"]["regions"]))
    if len(common) < k:
        raise RuntimeError(f"Only {len(common)} shared regions — need >= k_folds={k}.")

    # Fold assignment is by PROTEIN (uid), over the union of both models' proteins, so a
    # protein's residues never straddle folds. Deterministic given seed.
    all_uids = sorted(set(sides["ref"]["bg_by_uid"]) | set(sides["cmp"]["bg_by_uid"])
                      | {key[0] for key in common})
    perm = rng.permutation(len(all_uids))
    fold_of = {all_uids[i]: int(perm[i] % k) for i in range(len(all_uids))}

    hit_r, tp_r, fp_r = _oof_score(cfg, sides["ref"], common, fold_of, k, rng)
    hit_c, tp_c, fp_c = _oof_score(cfg, sides["cmp"], common, fold_of, k, rng)
    n = len(common)

    f1_r = _f1_from(hit_r, tp_r, fp_r)
    f1_c = _f1_from(hit_c, tp_c, fp_c)

    # PAIRED bootstrap over regions: resample indices once, apply to both models.
    # Background FP is a fixed nuisance constant per model (as in script 12).
    gap_b = np.empty(n_boot)
    fr_b = np.empty(n_boot)
    fc_b = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        fr_b[i] = _f1_from(hit_r[idx], tp_r[idx], fp_r)
        fc_b[i] = _f1_from(hit_c[idx], tp_c[idx], fp_c)
        gap_b[i] = fc_b[i] - fr_b[i]
    glo, ghi = np.percentile(gap_b, [2.5, 97.5])
    p_gap_le0 = float((gap_b <= 0).mean())

    summary = {
        "method": f"{k}-fold CV (feature pinned, threshold re-fit per fold), paired region bootstrap",
        "triplet": triplet, "n_regions_all_splits": int(n), "k_folds": k,
        "ref": {"run": sides["ref"]["run"], "feature": sides["ref"]["feat"],
                "point_f1": round(float(f1_r), 4), "bg_fp": int(fp_r),
                "ci95": [round(float(x), 4) for x in np.percentile(fr_b, [2.5, 97.5])]},
        "cmp": {"run": sides["cmp"]["run"], "feature": sides["cmp"]["feat"],
                "point_f1": round(float(f1_c), 4), "bg_fp": int(fp_c),
                "ci95": [round(float(x), 4) for x in np.percentile(fc_b, [2.5, 97.5])]},
        "gap_cmp_minus_ref": round(float(f1_c - f1_r), 4),
        "gap_ci95": [round(float(glo), 4), round(float(ghi), 4)],
        "gap_excludes_0": bool(glo > 0),
        "boot_p_cmp_le_ref": round(p_gap_le0, 4),
    }
    logger.info(f"CV gap (cmp-ref) = {f1_c - f1_r:+.3f}  95%CI[{glo:+.3f},{ghi:+.3f}]  "
                f"excludes 0: {glo > 0}  (boot p[cmp<=ref]={p_gap_le0:.4f}, "
                f"n={n} regions over {k} folds)")

    out = Path(cfg.paths.output_dir) / "bootstrap" / f"{triplet}_cv{k}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(out, "w"), indent=2)
    logger.info("Summary:\n" + json.dumps(summary, indent=2))
    logger.info(f"Saved → {out}")


if __name__ == "__main__":
    main()
