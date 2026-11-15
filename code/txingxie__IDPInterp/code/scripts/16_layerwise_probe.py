#!/usr/bin/env python3
"""
Layer-wise linear probes on raw residue embeddings.

Usage (one model at a time; embeddings must exist for the layers requested)
  python scripts/16_layerwise_probe.py paths=colab \
     +probe.model=esm2_8m +probe.layers=[1,2,3,4,5,6] \
     +probe.axes=[disprot-molten_globule,disprot-pre_molten_globule] \
     paths.embedding_dir=/mnt/local-scratch/embeddings

  # composition baseline (no embeddings, no layer sweep)
  python scripts/16_layerwise_probe.py paths=colab +probe.features=composition \
     +probe.axes=[disprot-molten_globule,disprot-pre_molten_globule]
"""
import json
import logging
import sys
import time
from pathlib import Path

import hydra
import numpy as np
import pandas as pd
from omegaconf import DictConfig

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from idpinterp.models.memmap_store import MetadataIndex, load_mmap_config

logger = logging.getLogger(__name__)

# Kyte-Doolittle hydropathy
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
AAS = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AAS)}
CHARGE = {"K": 1.0, "R": 1.0, "D": -1.0, "E": -1.0}
AROMATIC = set("FWY")
DISORDER_PROMOTING = set("PGSQEKRA")


# Region-aware F1, mirrored from 12b_bootstrap_cv.py so probe and SAE-feature scores
# are directly comparable. Keep the two in sync.
def _sorted(a):
    return np.sort(np.asarray(a, dtype=np.float32))


def _ge_count(sorted_arr, thr):
    """# elements >= thr in an ascending-sorted array."""
    return int(sorted_arr.size - np.searchsorted(sorted_arr, thr, side="left"))


def _fit_threshold(train_regions, train_bg_sorted, max_cands=2000):
    """Threshold maximising region-F1 on the threshold-fitting proteins. Candidates are
    distinct region-max scores (a threshold between region maxes never changes recall)."""
    maxes = np.array([r.max() if r.size else -np.inf for r in train_regions], dtype=np.float32)
    finite = maxes[np.isfinite(maxes)]
    if finite.size == 0:
        return None
    cands = np.unique(finite)
    if cands.size > max_cands:  # probe scores are continuous; subsample the grid
        cands = np.quantile(cands, np.linspace(0, 1, max_cands)).astype(np.float32)

    # Vectorised over candidates: TP at a threshold is the number of in-region residues
    # at or above it, so every region concatenates into one sorted array and a candidate
    # costs a searchsorted instead of a re-sort. Probe scores are continuous, which makes
    # the per-region loop the bottleneck. Identical objective, computed without it.
    n_reg = maxes.size
    maxes_sorted = np.sort(finite)
    all_res = (_sorted(np.concatenate([r for r in train_regions if r.size]))
               if any(r.size for r in train_regions) else np.empty(0, np.float32))

    rec = 1.0 - np.searchsorted(maxes_sorted, cands, side="left") / max(n_reg, 1)
    tp = all_res.size - np.searchsorted(all_res, cands, side="left")
    fp = train_bg_sorted.size - np.searchsorted(train_bg_sorted, cands, side="left")
    prec = tp / (tp + fp + 1e-10)
    f1 = 2 * prec * rec / (prec + rec + 1e-10)
    return float(cands[int(np.argmax(f1))]) if f1.size else None


def _f1_from(hit, tp, fp):
    R = hit.mean() if hit.size else 0.0
    TP = tp.sum()
    P = TP / (TP + fp + 1e-10)
    return 2 * P * R / (P + R + 1e-10)


_GOLD_CACHE = {}


def _gold(cfg, dataset):
    """Read the annotation parquet once per process; re-reading it per (layer, axis)
    is slow on network-backed storage."""
    if dataset not in _GOLD_CACHE:
        _GOLD_CACHE[dataset] = pd.read_parquet(
            Path(cfg.paths.dataset_dir) / f"{dataset}.parquet")
    return _GOLD_CACHE[dataset]


def _prepare(cfg, dataset, axis, meta):
    """Precompute, once per axis and independent of any score:
        slices  {uid: (mmap_start, stored_length)}
        regions {uid: [(start, end), ...]}   local coordinates
        mask    {uid: bool array of length stored_length}
    A protein with no annotation string contributes background only."""
    gold = _gold(cfg, dataset)
    col = f"mask_{axis.replace('-', '_')}"
    if col not in gold.columns:
        raise KeyError(f"{col} not in {dataset}.parquet — build the axis first (08).")
    mask_by = dict(zip(gold["uniprot_id"], gold[col]))
    seq_by = dict(zip(gold["uniprot_id"], gold["sequence"]))

    slices, regions, masks = {}, {}, {}
    for r in meta.itertuples():
        uid = r.uniprot_id
        s0, L = int(r.mmap_start_idx), int(r.stored_length)
        slices[uid] = (s0, L)
        m = mask_by.get(uid)
        if not isinstance(m, str):
            a = np.zeros(L, bool)
        else:
            a = (np.frombuffer(m[:L].encode(), dtype=np.uint8) - ord("0")).astype(bool)
            if a.size < L:
                a = np.concatenate([a, np.zeros(L - a.size, bool)])
        masks[uid] = a
        d = np.diff(np.concatenate([[0], a.astype(int), [0]]))
        regions[uid] = [(int(st), int(en))
                        for st, en in zip(np.where(d == 1)[0], np.where(d == -1)[0])]
    return slices, regions, masks, seq_by


def _shuffle_regions(regions, masks, rng):
    """Null control: relocate every region to a random position inside the SAME protein,
    preserving region count and length. The corpus geometry (how many regions, how long,
    how much background) is untouched, so probing these labels measures the F1 this
    pipeline reports when there is no signal to find.

    This matters because the region-aware metric is lenient by construction — one firing
    residue recovers a whole region — so its chance level is well above zero and must be
    quoted alongside any small F1."""
    r_out, m_out = {}, {}
    for uid, regs in regions.items():
        L = masks[uid].size
        a = np.zeros(L, bool)
        placed = []
        for (st, en) in regs:
            w = en - st
            if w <= 0 or w >= L:
                continue
            for _ in range(50):                       # rejection-sample a free slot
                s = int(rng.integers(0, L - w))
                if not a[s:s + w].any():
                    a[s:s + w] = True
                    placed.append((s, s + w))
                    break
        m_out[uid] = a
        d = np.diff(np.concatenate([[0], a.astype(int), [0]]))
        r_out[uid] = [(int(s), int(e))
                      for s, e in zip(np.where(d == 1)[0], np.where(d == -1)[0])]
    return r_out, m_out


def _score_sets(uids, scores_by_uid, regions, masks):
    """Turn per-protein score arrays into the region / background structures the
    metric consumes."""
    reg, bg = {}, {}
    for uid in uids:
        sc = scores_by_uid[uid]
        for (st, en) in regions[uid]:
            reg[(uid, st, en)] = sc[st:en]
        bg[uid] = sc[~masks[uid]]
    return reg, bg


def _composition_matrix(seq, L, window):
    """Per-residue handcrafted features from sequence alone: windowed AA frequency,
    hydropathy, net charge, FCR, aromatic and disorder-promoting fraction.
    Uses cumulative sums so it is O(L) per protein."""
    s = (seq or "")[:L]
    if len(s) < L:
        s = s + "A" * (L - len(s))
    idx = np.array([AA_IDX.get(c, -1) for c in s])
    onehot = np.zeros((L, len(AAS)), np.float32)
    ok = idx >= 0
    onehot[np.arange(L)[ok], idx[ok]] = 1.0

    scal = np.stack([
        np.array([KD.get(c, 0.0) for c in s], np.float32),
        np.array([CHARGE.get(c, 0.0) for c in s], np.float32),
        np.array([abs(CHARGE.get(c, 0.0)) for c in s], np.float32),      # FCR
        np.array([1.0 if c in AROMATIC else 0.0 for c in s], np.float32),
        np.array([1.0 if c in DISORDER_PROMOTING else 0.0 for c in s], np.float32),
    ], axis=1)

    raw = np.concatenate([onehot, scal], axis=1)                          # (L, 25)
    cs = np.concatenate([np.zeros((1, raw.shape[1]), np.float32), np.cumsum(raw, 0)], 0)
    lo = np.clip(np.arange(L) - window, 0, L)
    hi = np.clip(np.arange(L) + window + 1, 0, L)
    win = (cs[hi] - cs[lo]) / np.maximum(hi - lo, 1)[:, None]
    return np.concatenate([onehot, win, scal[:, :2]], axis=1).astype(np.float32)


class FeatureBank:
    """Uniform accessor: .rows(uid) -> (L, d) float32 for one protein."""

    def __init__(self, kind, slices, seq_by, mmap=None, window=10):
        self.kind, self.slices, self.seq_by, self.mmap, self.window = \
            kind, slices, seq_by, mmap, window
        any_uid = next(iter(slices))
        self.d = self.rows(any_uid).shape[1]

    def rows(self, uid):
        s0, L = self.slices[uid]
        if self.kind == "embedding":
            return np.asarray(self.mmap[s0:s0 + L], dtype=np.float32)
        return _composition_matrix(self.seq_by.get(uid, ""), L, self.window)


def _fit_probe(bank, uids, masks, rng, neg_ratio, max_rows, max_iter):
    """Balanced logistic probe on residues of `uids`. Keeps all positives, subsamples
    negatives to neg_ratio:1, caps total rows."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    Xs, ys = [], []
    for uid in uids:
        a = masks[uid]
        pos = np.where(a)[0]
        neg = np.where(~a)[0]
        if pos.size == 0 and neg.size == 0:
            continue
        n_neg = min(neg.size, max(int(pos.size * neg_ratio), 1 if pos.size else 0))
        if pos.size == 0:      # background-only protein: contribute a few negatives
            n_neg = min(neg.size, 64)
        take_neg = rng.choice(neg, n_neg, replace=False) if n_neg else np.empty(0, int)
        take = np.concatenate([pos, take_neg]).astype(int)
        if take.size == 0:
            continue
        R = bank.rows(uid)
        Xs.append(R[take])
        ys.append(a[take].astype(np.int8))
    if not Xs:
        return None
    X = np.concatenate(Xs); y = np.concatenate(ys)
    if y.sum() == 0 or y.sum() == y.size:
        return None
    if X.shape[0] > max_rows:
        keep = rng.choice(X.shape[0], max_rows, replace=False)
        X, y = X[keep], y[keep]
    sc = StandardScaler().fit(X)
    clf = LogisticRegression(max_iter=max_iter, class_weight="balanced", n_jobs=-1)
    clf.fit(sc.transform(X), y)
    return sc, clf


def _apply(model, bank, uids):
    sc, clf = model
    out = {}
    for uid in uids:
        R = bank.rows(uid)
        out[uid] = clf.decision_function(sc.transform(R)).astype(np.float32)
    return out


def run_axis(cfg, bank, axis, slices, regions, masks, k, rng, n_boot,
             neg_ratio, max_rows, max_iter, dump=None):
    """K-fold, protein-held-out probe; returns pooled out-of-fold region F1 + CI."""
    all_uids = sorted(slices)
    keys = sorted([(u, st, en) for u in all_uids for (st, en) in regions[u]])
    n = len(keys)
    if n < k:
        logger.warning(f"  {axis}: only {n} regions (< k={k}) — skipping.")
        return None

    perm = rng.permutation(len(all_uids))
    fold_of = {all_uids[i]: int(perm[i] % k) for i in range(len(all_uids))}
    pos_of = {key: i for i, key in enumerate(keys)}

    hit = np.zeros(n); tp = np.zeros(n); fp_total = 0
    for f in range(k):
        held = [u for u in all_uids if fold_of[u] == f]
        train = [u for u in all_uids if fold_of[u] != f]
        if not held or not train:
            continue
        tperm = rng.permutation(len(train))
        cut = max(int(0.8 * len(train)), 1)
        fit_uids = [train[i] for i in tperm[:cut]]
        thr_uids = [train[i] for i in tperm[cut:]] or fit_uids

        model = _fit_probe(bank, fit_uids, masks, rng, neg_ratio, max_rows, max_iter)
        if model is None:
            logger.warning(f"  {axis}: fold {f} had no usable training signal.")
            continue

        thr_scores = _apply(model, bank, thr_uids)
        treg, tbg = _score_sets(thr_uids, thr_scores, regions, masks)
        thr = _fit_threshold(list(treg.values()),
                             _sorted(np.concatenate(list(tbg.values()) or [np.empty(0, np.float32)])))
        if thr is None:
            thr = np.inf

        held_scores = _apply(model, bank, held)
        hreg, hbg = _score_sets(held, held_scores, regions, masks)
        for key, arr in hreg.items():
            i = pos_of[key]
            hit[i] = 1.0 if (arr.size and arr.max() >= thr) else 0.0
            tp[i] = _ge_count(_sorted(arr), thr)
        fp_total += _ge_count(
            _sorted(np.concatenate(list(hbg.values()) or [np.empty(0, np.float32)])), thr)

    f1 = _f1_from(hit, tp, int(fp_total))

    if dump is not None:
        # Per-region out-of-fold outcomes keyed by (uid, start, end), so two models can
        # be paired on identical regions and the gap bootstrapped. Marginal CIs on each
        # model separately cannot support a claim about the difference.
        Path(dump).parent.mkdir(parents=True, exist_ok=True)
        np.savez(dump,
                 uid=np.array([k_[0] for k_ in keys]),
                 start=np.array([k_[1] for k_ in keys]),
                 end=np.array([k_[2] for k_ in keys]),
                 hit=hit, tp=tp, fp=np.array([int(fp_total)]))

    boot = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boot[i] = _f1_from(hit[idx], tp[idx], int(fp_total))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return {"n_regions": int(n), "k_folds": int(k),
            "region_recall": round(float(hit.mean()), 4),
            "f1": round(float(f1), 4),
            "ci95": [round(float(lo), 4), round(float(hi), 4)],
            "bg_fp": int(fp_total)}


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s")
    p = cfg.get("probe", {})
    dataset = p.get("dataset", "gold_curated")
    model_name = p.get("model", "esm2_8m")
    kind = p.get("features", "embedding")          # embedding | composition
    layers = list(p.get("layers", [])) or None
    axes = list(p.get("axes", ["disprot-molten_globule", "disprot-pre_molten_globule"]))
    k = int(p.get("k_folds", 5))
    n_boot = int(p.get("n_boot", 2000))
    neg_ratio = float(p.get("neg_ratio", 20))
    max_rows = int(p.get("max_train_rows", 300_000))
    max_iter = int(p.get("max_iter", 300))
    window = int(p.get("window", 10))
    seed = int(p.get("seed", 0))
    null_control = bool(p.get("null_control", True))
    dump_regions = bool(p.get("dump_regions", False))

    embed_dir = Path(cfg.paths.embedding_dir) / model_name / dataset
    if kind == "composition":
        layers = [layers[0]] if layers else [None]   # one pass; layer only locates metadata
    elif not layers:
        raise ValueError("+probe.layers=[...] is required for features=embedding.")

    out = Path(cfg.paths.output_dir) / "probe" / f"{model_name}_{kind}.json"
    out.parent.mkdir(parents=True, exist_ok=True)

    # Resume from a partial sweep: completed layers are reused rather than recomputed,
    # since embeddings may live on scratch that does not survive a restart.
    # Pass +probe.overwrite=true to force a clean sweep.
    results = {"model": model_name, "features": kind, "dataset": dataset,
               "k_folds": k, "neg_ratio": neg_ratio, "layers": {}}
    if out.exists() and not bool(p.get("overwrite", False)):
        try:
            prev = json.load(open(out))
            if prev.get("k_folds") == k and prev.get("features") == kind:
                results = prev
                if results["layers"]:
                    logger.info(f"Resuming from {out}: "
                                f"{sorted(results['layers'])} already done.")
        except (json.JSONDecodeError, OSError) as e:
            logger.warning(f"Could not reuse {out} ({e}); starting fresh.")

    for layer in layers:
        # metadata defines residue order; for composition any available layer works
        meta_path = embed_dir / f"layer{layer}_metadata.parquet" if layer is not None else None
        if meta_path is None or not meta_path.exists():
            cand = sorted(embed_dir.glob("layer*_metadata.parquet"))
            if not cand:
                raise FileNotFoundError(
                    f"No metadata under {embed_dir}. Generate embeddings first (02).")
            if kind == "embedding":
                logger.warning(f"  layer {layer}: metadata missing — skipping.")
                continue
            meta_path = cand[0]
        meta = MetadataIndex.load(meta_path)

        mm = None
        if kind == "embedding":
            mpath = embed_dir / f"layer{layer}.mmap"
            if not mpath.exists():
                logger.warning(f"  layer {layer}: {mpath} missing — skipping.")
                continue
            mc = load_mmap_config(mpath)
            mm = np.memmap(mpath, dtype=np.float16, mode="r",
                           shape=(mc["total_residues"], mc["hidden_dim"]), order="C")

        tag = f"layer{layer}" if kind == "embedding" else "composition"
        done = results["layers"].get(tag, {})
        def _have(ax):
            npz = out.parent / "regions" / f"{model_name}_{tag}_{ax}.npz"
            return (ax in done
                    and (not null_control or (ax + "__null") in done)
                    and (not dump_regions or npz.exists()))
        if all(_have(ax) for ax in axes):
            logger.info(f"=== {model_name} / {tag} — already complete, skipping ===")
            if mm is not None:
                del mm
            continue
        logger.info(f"\n=== {model_name} / {tag} ===")
        for axis in axes:
            # Skip only if the null is present too when it was asked for — this is what
            # lets a fast null-free sweep be topped up later by a null-only pass.
            if _have(axis):
                logger.info(f"  {axis}: already done (F1={done[axis]['f1']:.3f}) — skipping.")
                continue
            slices, regions, masks, seq_by = _prepare(cfg, dataset, axis, meta)
            bank = FeatureBank(kind, slices, seq_by, mmap=mm, window=window)

            npz_ok = (not dump_regions or
                      (out.parent / "regions" / f"{model_name}_{tag}_{axis}.npz").exists())
            if axis in done and npz_ok:   # real result cached; only the null is missing
                res = done[axis]
                logger.info(f"  {axis}: reusing F1={res['f1']:.3f}; adding null only.")
            else:
                t0 = time.time()
                logger.info(f"  {axis}: probing {len(slices)} proteins, d={bank.d} ...")
                dump = (out.parent / "regions" /
                        f"{model_name}_{tag}_{axis}.npz") if dump_regions else None
                res = run_axis(cfg, bank, axis, slices, regions, masks, k,
                               np.random.default_rng(seed), n_boot,
                               neg_ratio, max_rows, max_iter, dump=dump)
                if res is None:
                    continue
                res["d_features"] = int(bank.d)
                res["seconds"] = round(time.time() - t0, 1)
                results["layers"].setdefault(tag, {})[axis] = res
                logger.info(f"  {axis}: F1={res['f1']:.3f} "
                            f"95%CI[{res['ci95'][0]:.3f},{res['ci95'][1]:.3f}] "
                            f"recall={res['region_recall']:.3f} n={res['n_regions']} "
                            f"[{res['seconds']:.0f}s]")

            if null_control and (axis + "__null") not in done:
                t0 = time.time()
                nreg, nmsk = _shuffle_regions(regions, masks, np.random.default_rng(seed + 1))
                nres = run_axis(cfg, bank, axis + " [null]", slices, nreg, nmsk, k,
                                np.random.default_rng(seed), n_boot, neg_ratio,
                                max_rows, max_iter)
                if nres is not None:
                    nres["seconds"] = round(time.time() - t0, 1)
                    results["layers"].setdefault(tag, {})[axis + "__null"] = nres
                    logger.info(f"  {axis} [shuffled null]: F1={nres['f1']:.3f} "
                                f"95%CI[{nres['ci95'][0]:.3f},{nres['ci95'][1]:.3f}] "
                                f"[{nres['seconds']:.0f}s] "
                                f"— chance level for this corpus geometry")
        if mm is not None:
            del mm
        json.dump(results, open(out, "w"), indent=2)   # checkpoint after every layer
        logger.info(f"  checkpointed → {out}")

    logger.info("\n" + json.dumps(results, indent=2))
    logger.info(f"Saved → {out}")


if __name__ == "__main__":
    main()
