"""Protein-level paired bootstrap significance test (Section "Experimental Settings").

In each of N (=1000) iterations, test proteins are sampled with replacement (same sample for both
methods → paired), residues of the sampled proteins are pooled, and micro-averaged metrics are
computed for both methods. The one-sided p-value is the fraction of iterations in which method A
(FlexBind) does NOT outperform method B:  p = #{metric_A - metric_B <= 0} / N.
Threshold-dependent metrics (Precision, Recall, MCC) use each method's own pre-specified threshold.

Prediction files are pickles in the format written by scripts/train.py / scripts/evaluate.py:
    {"tasks": [...], "proteins": {id: {"prob": (L, T), "label": (L, T)}}, "thresholds": {task: thr}}
Baselines (web servers, other tools) must be converted to this format; tasks a baseline does not
predict can be filled with NaN probabilities and are skipped.

Example:
    python scripts/bootstrap_test.py --pred_a checkpoints/DP81/seed42/test_predictions.pkl \
        --pred_b baselines/DP81/DisoFLAG.pkl --out results/DP81_vs_DisoFLAG.json
"""
import argparse
import json
import pickle

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.evaluate import METRIC_NAMES, best_f1_threshold, binary_metrics  # noqa: E402


def _metrics(y, p, thr):
    out = {m: np.nan for m in METRIC_NAMES}
    if len(y) == 0 or y.min() == y.max():
        return out
    out["AUPR"] = average_precision_score(y, p)
    out["AUROC"] = roc_auc_score(y, p)
    out["F1max"] = best_f1_threshold(y, p)[0]
    out.update({k: v for k, v in binary_metrics(y, p, thr).items() if k in METRIC_NAMES})
    return out


def _per_protein(preds, task, ids):
    j = preds["tasks"].index(task)
    ys, ps = [], []
    for pid in ids:
        rec = preds["proteins"][pid]
        y = np.asarray(rec["label"])[:, j]
        pr = np.asarray(rec["prob"], dtype=np.float64)[:, j]
        keep = y >= 0
        ys.append(y[keep].astype(np.int64))
        ps.append(pr[keep])
    return ys, ps


def paired_bootstrap(pred_a, pred_b, task, n_iter=1000, seed=0):
    ids = sorted(set(pred_a["proteins"]) & set(pred_b["proteins"]))
    ya, pa = _per_protein(pred_a, task, ids)
    yb, pb = _per_protein(pred_b, task, ids)
    if any(np.isnan(x).any() for x in pb) or any(np.isnan(x).any() for x in pa):
        return None
    for a, b in zip(ya, yb):
        if not np.array_equal(a, b):
            raise ValueError(f"Label mismatch between prediction files for task {task}")
    ta, tb = pred_a["thresholds"][task], pred_b["thresholds"][task]

    cat = lambda xs, idx: np.concatenate([xs[i] for i in idx])  # noqa: E731
    all_idx = np.arange(len(ids))
    obs_a = _metrics(cat(ya, all_idx), cat(pa, all_idx), ta)
    obs_b = _metrics(cat(yb, all_idx), cat(pb, all_idx), tb)

    rng = np.random.RandomState(seed)
    worse = {m: 0 for m in METRIC_NAMES}
    valid = {m: 0 for m in METRIC_NAMES}
    for _ in range(n_iter):
        idx = rng.randint(0, len(ids), size=len(ids))
        y = cat(ya, idx)
        ma = _metrics(y, cat(pa, idx), ta)
        mb = _metrics(y, cat(pb, idx), tb)
        for m in METRIC_NAMES:
            if np.isfinite(ma[m]) and np.isfinite(mb[m]):
                valid[m] += 1
                worse[m] += int(ma[m] - mb[m] <= 0)
    return {m: {"A": obs_a[m], "B": obs_b[m], "diff": obs_a[m] - obs_b[m],
                "p_value": worse[m] / valid[m] if valid[m] else np.nan,
                "n_valid_iter": valid[m]} for m in METRIC_NAMES}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pred_a", required=True, help="FlexBind predictions (.pkl)")
    p.add_argument("--pred_b", required=True, help="second-best method / configuration (.pkl)")
    p.add_argument("--n_iter", type=int, default=1000)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default=None)
    args = p.parse_args()

    with open(args.pred_a, "rb") as f:
        A = pickle.load(f)
    with open(args.pred_b, "rb") as f:
        B = pickle.load(f)

    report = {}
    for task in A["tasks"]:
        if task not in B["tasks"]:
            continue
        res = paired_bootstrap(A, B, task, args.n_iter, args.seed)
        if res is None:
            print(f"[{task}] skipped (missing predictions)")
            continue
        report[task] = res
        print(f"[{task}]")
        for m, r in res.items():
            star = "*" if np.isfinite(r["p_value"]) and r["p_value"] < args.alpha else ""
            print(f"  {m:<10} A={r['A']:.4f}  B={r['B']:.4f}  diff={r['diff']:+.4f}  p={r['p_value']:.4f} {star}")
    if args.out:
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)


if __name__ == "__main__":
    main()
