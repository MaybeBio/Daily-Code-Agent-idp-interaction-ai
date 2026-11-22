"""Evaluation utilities and CLI.

Metrics (all micro-averaged, i.e. residues pooled over all test proteins):
    AUPR, F1max, AUROC  — threshold-free (threshold sweep)
    Precision, Recall, MCC — at a pre-specified operating point: the threshold that maximised
                             F1 on the *validation* set (no test data used for threshold tuning).

CLI example (evaluate a trained checkpoint on another test set, e.g. DP94 with the DP93 model):
    python scripts/evaluate.py --checkpoint checkpoints/DP93/seed42/best.pt \
        --npz_path data/processed/DP94/test.npz \
        --disorder_txt ... --protein_txt ... --rna_txt ... --dna_txt ... \
        --out_dir results/DP94_seed42
"""
import argparse
import json
import os
import pickle
import sys
from typing import Dict, Optional, Sequence

import numpy as np
import torch
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.flexbind import FlexBind  # noqa: E402
from models.losses import asymmetric_focal_loss  # noqa: E402
from scripts.dataset import TASKS, build_loader, log  # noqa: E402

NAN = float("nan")
METRIC_NAMES = ["AUPR", "F1max", "AUROC", "Precision", "Recall", "MCC"]


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def _has_both_classes(y_true: np.ndarray) -> bool:
    return len(y_true) > 0 and y_true.min() != y_true.max()


def best_f1_threshold(y_true: np.ndarray, y_prob: np.ndarray):
    """Returns (F1max, threshold) from a full threshold sweep."""
    if not _has_both_classes(y_true):
        return NAN, 0.5
    P, R, TH = precision_recall_curve(y_true, y_prob)
    f1 = 2 * P[:-1] * R[:-1] / (P[:-1] + R[:-1] + 1e-12)   # last point (P=1, R=0) has no threshold
    i = int(np.nanargmax(f1))
    return float(f1[i]), float(TH[i])


def binary_metrics(y_true: np.ndarray, y_prob: np.ndarray, threshold: float) -> Dict[str, float]:
    pred = (y_prob >= threshold).astype(np.int64)
    y = y_true.astype(np.int64)
    tp = float(np.sum((pred == 1) & (y == 1)))
    fp = float(np.sum((pred == 1) & (y == 0)))
    fn = float(np.sum((pred == 0) & (y == 1)))
    tn = float(np.sum((pred == 0) & (y == 0)))
    precision = tp / (tp + fp) if tp + fp > 0 else 0.0
    recall = tp / (tp + fn) if tp + fn > 0 else 0.0
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0
    return {"Precision": precision, "Recall": recall, "MCC": float(mcc),
            "TP": tp, "FP": fp, "FN": fn, "TN": tn}


def compute_metrics(y_true: np.ndarray, y_prob: np.ndarray, threshold: float) -> Dict[str, float]:
    out = {"n_residues": int(len(y_true)), "n_positive": int(np.sum(y_true)), "threshold": float(threshold)}
    if _has_both_classes(y_true):
        out["AUPR"] = float(average_precision_score(y_true, y_prob))
        out["AUROC"] = float(roc_auc_score(y_true, y_prob))
        out["F1max"] = best_f1_threshold(y_true, y_prob)[0]
    else:
        out["AUPR"] = out["AUROC"] = out["F1max"] = NAN
    out.update(binary_metrics(y_true, y_prob, threshold))
    return out


# --------------------------------------------------------------------------- #
# Predictions
# --------------------------------------------------------------------------- #
@torch.no_grad()
def collect_predictions(model, loader, device, gammas: Optional[Dict[str, Dict[str, float]]] = None):
    """Runs the model and returns per-protein predictions (+ mean per-task ASL if `gammas`).

    preds = {"tasks": [...], "proteins": {id: {"prob": (L, T) float32, "label": (L, T) int8}}}
    """
    model.eval()
    proteins = {}
    loss_sums = {t: 0.0 for t in TASKS}
    n_batches = 0
    for batch in loader:
        H = batch["emb"].to(device, non_blocking=True)
        Y = batch["label"].to(device, non_blocking=True)
        M = batch["mask"].to(device, non_blocking=True)
        logits = model(H, M)
        if gammas is not None:
            for j, t in enumerate(TASKS):
                valid = M & (Y[:, :, j] > -50)
                loss_sums[t] += float(asymmetric_focal_loss(
                    logits[t], Y[:, :, j], valid, gammas["gamma_pos"][t], gammas["gamma_neg"][t]))
            n_batches += 1
        probs = torch.stack([torch.sigmoid(logits[t]) for t in TASKS], dim=-1).float().cpu().numpy()
        labels = Y.cpu().numpy()
        lengths = M.sum(1).cpu().numpy()
        for b, pid in enumerate(batch["ids"]):
            L = int(lengths[b])
            proteins[pid] = {"prob": probs[b, :L].astype(np.float32), "label": labels[b, :L].astype(np.int8)}
    preds = {"tasks": list(TASKS), "proteins": proteins}
    losses = {t: v / max(n_batches, 1) for t, v in loss_sums.items()} if gammas is not None else None
    return preds, losses


def pooled(preds: dict, task: str, ids: Optional[Sequence[str]] = None):
    """Pools residues of the given proteins (micro-averaging) for one task."""
    j = preds["tasks"].index(task)
    ids = list(preds["proteins"].keys()) if ids is None else ids
    ys, ps = [], []
    for pid in ids:
        rec = preds["proteins"][pid]
        y = rec["label"][:, j]
        keep = y >= 0
        ys.append(y[keep])
        ps.append(rec["prob"][keep, j])
    if not ys:
        return np.array([]), np.array([])
    return np.concatenate(ys).astype(np.int64), np.concatenate(ps).astype(np.float64)


def select_thresholds(preds: dict) -> Dict[str, float]:
    """Threshold maximising F1 for each task (to be computed on the validation set)."""
    return {t: best_f1_threshold(*pooled(preds, t))[1] for t in preds["tasks"]}


def evaluate_predictions(preds: dict, thresholds: Dict[str, float]) -> Dict[str, Dict[str, float]]:
    return {t: compute_metrics(*pooled(preds, t), thresholds.get(t, 0.5)) for t in preds["tasks"]}


def mean_aupr(preds: dict) -> float:
    vals = []
    for t in preds["tasks"]:
        y, p = pooled(preds, t)
        if _has_both_classes(y):
            vals.append(average_precision_score(y, p))
    return float(np.mean(vals)) if vals else NAN


def save_predictions(preds: dict, thresholds: Dict[str, float], path: str, method: str = "FlexBind"):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump({**preds, "thresholds": thresholds, "method": method}, f)


def format_table(metrics: Dict[str, Dict[str, float]]) -> str:
    lines = ["task      " + "".join(f"{m:>10}" for m in METRIC_NAMES)]
    for t, m in metrics.items():
        lines.append(f"{t:<10}" + "".join(f"{m[k]:>10.4f}" for k in METRIC_NAMES))
    return "\n".join(lines)


def plot_pr_curves(preds: dict, out_dir: str, suffix: str = ""):
    """Optional PR-curve plots (requires matplotlib)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log("matplotlib not installed; skipping PR plots")
        return
    os.makedirs(out_dir, exist_ok=True)
    for t in preds["tasks"]:
        y, p = pooled(preds, t)
        if not _has_both_classes(y):
            continue
        P, R, _ = precision_recall_curve(y, p)
        plt.figure(figsize=(5, 5))
        plt.plot(R, P, label=f"AUPR={average_precision_score(y, p):.4f}")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(f"PR curve - {t}")
        plt.legend(loc="upper right")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"pr_{t}{suffix}.png"), dpi=200)
        plt.close()


# --------------------------------------------------------------------------- #
# Checkpoints
# --------------------------------------------------------------------------- #
def load_flexbind(checkpoint_path: str, device, model_cfg: Optional[dict] = None):
    """Loads a FlexBind checkpoint written by scripts/train.py (or a bare state_dict)."""
    ckpt = torch.load(checkpoint_path, map_location=device)
    if isinstance(ckpt, dict) and "model_state" in ckpt:
        cfg = ckpt.get("model_config", model_cfg or {})
        state, thresholds = ckpt["model_state"], ckpt.get("thresholds")
    else:
        cfg, state, thresholds = (model_cfg or {}), ckpt, None
    model = FlexBind(**cfg).to(device)
    model.load_state_dict(state)
    model.eval()
    return model, thresholds


def main():
    p = argparse.ArgumentParser(description="Evaluate a trained FlexBind checkpoint.")
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--npz_path", required=True)
    p.add_argument("--disorder_txt", required=True)
    p.add_argument("--protein_txt", required=True)
    p.add_argument("--rna_txt", required=True)
    p.add_argument("--dna_txt", required=True)
    p.add_argument("--out_dir", default="./results")
    p.add_argument("--batch_size", type=int, default=32)
    p.add_argument("--num_workers", type=int, default=2)
    p.add_argument("--plot", action="store_true")
    args = p.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model, thresholds = load_flexbind(args.checkpoint, device)
    if thresholds is None:
        log("Checkpoint has no validation thresholds; falling back to 0.5")
        thresholds = {t: 0.5 for t in TASKS}

    loader = build_loader(args.npz_path, {"disorder": args.disorder_txt, "protein": args.protein_txt,
                                          "rna": args.rna_txt, "dna": args.dna_txt},
                          batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers)
    preds, _ = collect_predictions(model, loader, device)
    metrics = evaluate_predictions(preds, thresholds)
    os.makedirs(args.out_dir, exist_ok=True)
    save_predictions(preds, thresholds, os.path.join(args.out_dir, "predictions.pkl"))
    with open(os.path.join(args.out_dir, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)
    if args.plot:
        plot_pr_curves(preds, args.out_dir)
    print(format_table(metrics))


if __name__ == "__main__":
    main()
