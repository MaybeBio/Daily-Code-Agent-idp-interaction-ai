"""Training / evaluation of FlexBind (Sections "Training Objective" and "Experimental Settings").

* frozen ProtT5 embeddings (precomputed .npz), full-length sequences, dynamic padding
* AdamW (lr 1e-4, weight decay 0.01), batch size 32
* task-specific asymmetric focal loss + uncertainty-based task weighting (sigma_t init 1)
* early stopping / checkpoint selection on validation AUPR (mean over the four tasks), patience 10
* operating thresholds = F1-maximising thresholds on the validation set
* repeated over independent random seeds (10 by default) → mean ± SD

Example:
    python scripts/train.py --config configs/default_config.yaml --name DP81 \
        --npz_path data/processed/DP81/train_val.npz --test_npz_path data/processed/DP81/test.npz \
        --disorder_txt ... --protein_txt ... --rna_txt ... --dna_txt ... \
        --test_disorder_txt ... --test_protein_txt ... --test_rna_txt ... --test_dna_txt ...
"""
import argparse
import copy
import json
import os
import random
import sys
import time

import numpy as np
import torch
import yaml
from tqdm.auto import tqdm

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.flexbind import FlexBind  # noqa: E402
from models.layers import HAS_MAMBA_SSM  # noqa: E402
from models.losses import UncertaintyWeightedLoss, asymmetric_focal_loss  # noqa: E402
from scripts.dataset import TASKS, build_loader, load_npz_ids, log, split_ids  # noqa: E402
from scripts.evaluate import (METRIC_NAMES, collect_predictions, evaluate_predictions,  # noqa: E402
                              format_table, mean_aupr, plot_pr_curves, save_predictions,
                              select_thresholds)

DEFAULT_CONFIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              "configs", "default_config.yaml")


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def load_config(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def train_one_epoch(model, criterion, loader, optimizer, cfg, device, epoch):
    model.train()
    criterion.train()
    gp, gn = cfg["loss"]["gamma_pos"], cfg["loss"]["gamma_neg"]
    total, n = 0.0, 0
    pbar = tqdm(loader, desc=f"Train epoch {epoch}", leave=False)
    for batch in pbar:
        H = batch["emb"].to(device, non_blocking=True)
        Y = batch["label"].to(device, non_blocking=True)
        M = batch["mask"].to(device, non_blocking=True)

        logits = model(H, M)
        task_losses = {}
        for j, t in enumerate(TASKS):
            valid = M & (Y[:, :, j] > -50)        # padded / unlabelled residues excluded
            task_losses[t] = asymmetric_focal_loss(logits[t], Y[:, :, j], valid, gp[t], gn[t])
        loss = criterion(task_losses)

        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        if cfg["optim"].get("grad_clip"):
            torch.nn.utils.clip_grad_norm_(list(model.parameters()) + list(criterion.parameters()),
                                           cfg["optim"]["grad_clip"])
        optimizer.step()

        total += loss.item()
        n += 1
        pbar.set_postfix(loss=f"{total / n:.4f}")
    return total / max(n, 1)


def run_seed(seed, cfg, args, device, loaders):
    train_loader, val_loader, test_loader = loaders
    set_seed(seed)
    out_dir = os.path.join(args.out_dir, args.name, f"seed{seed}")
    os.makedirs(out_dir, exist_ok=True)

    model = FlexBind(**cfg["model"]).to(device)
    criterion = UncertaintyWeightedLoss(TASKS, init_sigma=cfg["loss"]["init_sigma"]).to(device)
    # Uncertainty parameters are excluded from weight decay (decay would bias sigma_t toward 1).
    optimizer = torch.optim.AdamW([
        {"params": model.parameters(), "weight_decay": cfg["optim"]["weight_decay"]},
        {"params": criterion.parameters(), "weight_decay": 0.0},
    ], lr=cfg["optim"]["lr"])

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    log(f"[seed {seed}] trainable parameters: {n_params / 1e6:.2f}M | fused Mamba kernel: "
        f"{HAS_MAMBA_SSM and device == 'cuda'}")

    best_path = os.path.join(out_dir, "best.pt")
    best_aupr, best_epoch, patience_left = -np.inf, 0, cfg["optim"]["patience"]
    history = []
    for epoch in range(1, cfg["optim"]["max_epochs"] + 1):
        t0 = time.time()
        train_loss = train_one_epoch(model, criterion, train_loader, optimizer, cfg, device, epoch)
        val_preds, val_losses = collect_predictions(model, val_loader, device, gammas=cfg["loss"])
        val_aupr = mean_aupr(val_preds)
        history.append({"epoch": epoch, "train_loss": train_loss, "val_mean_aupr": val_aupr,
                        "val_task_loss": val_losses, "sigma": criterion.sigmas()})
        log(f"[seed {seed}] epoch {epoch:3d} | train_loss {train_loss:.4f} | val mean AUPR {val_aupr:.4f} "
            f"| sigma {criterion.sigmas()} | {time.time() - t0:.0f}s")

        if np.isfinite(val_aupr) and val_aupr > best_aupr:
            best_aupr, best_epoch, patience_left = val_aupr, epoch, cfg["optim"]["patience"]
            torch.save({"model_state": model.state_dict(), "loss_state": criterion.state_dict(),
                        "model_config": model.config, "epoch": epoch, "seed": seed,
                        "val_mean_aupr": val_aupr}, best_path)
        else:
            patience_left -= 1
            if patience_left <= 0:
                log(f"[seed {seed}] early stopping at epoch {epoch} (best epoch {best_epoch})")
                break

    # ---- reload best checkpoint, fix thresholds on validation, evaluate on test ----
    ckpt = torch.load(best_path, map_location=device)
    model.load_state_dict(ckpt["model_state"])
    val_preds, _ = collect_predictions(model, val_loader, device)
    thresholds = select_thresholds(val_preds)
    ckpt["thresholds"] = thresholds
    torch.save(ckpt, best_path)

    test_preds, _ = collect_predictions(model, test_loader, device)
    test_metrics = evaluate_predictions(test_preds, thresholds)
    save_predictions(test_preds, thresholds, os.path.join(out_dir, "test_predictions.pkl"))
    if args.plot:
        plot_pr_curves(test_preds, out_dir, suffix="_test")

    result = {"seed": seed, "best_epoch": best_epoch, "val_mean_aupr": best_aupr,
              "thresholds": thresholds, "test": test_metrics, "history": history}
    with open(os.path.join(out_dir, "results.json"), "w") as f:
        json.dump(result, f, indent=2)
    log(f"[seed {seed}] TEST (thresholds from validation: {thresholds})\n{format_table(test_metrics)}")
    return result


def summarize(results, path):
    summary = {}
    for t in TASKS:
        summary[t] = {}
        for m in METRIC_NAMES:
            vals = np.array([r["test"][t][m] for r in results], dtype=np.float64)
            summary[t][m] = {"mean": float(np.nanmean(vals)),
                             "sd": float(np.nanstd(vals, ddof=1)) if len(vals) > 1 else 0.0,
                             "values": vals.tolist()}
    with open(path, "w") as f:
        json.dump({"seeds": [r["seed"] for r in results], "summary": summary}, f, indent=2)
    lines = ["task      " + "".join(f"{m:>18}" for m in METRIC_NAMES)]
    for t in TASKS:
        lines.append(f"{t:<10}" + "".join(
            f"{summary[t][m]['mean']:>10.4f}±{summary[t][m]['sd']:<7.4f}" for m in METRIC_NAMES))
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=DEFAULT_CONFIG)
    p.add_argument("--name", default="run", help="experiment name, e.g. DP81 / DP93")
    p.add_argument("--npz_path", required=True, help="training (or train+validation) embeddings")
    p.add_argument("--test_npz_path", required=True)
    p.add_argument("--disorder_txt", required=True)
    p.add_argument("--protein_txt", required=True)
    p.add_argument("--rna_txt", required=True)
    p.add_argument("--dna_txt", required=True)
    p.add_argument("--test_disorder_txt", required=True)
    p.add_argument("--test_protein_txt", required=True)
    p.add_argument("--test_rna_txt", required=True)
    p.add_argument("--test_dna_txt", required=True)
    # Optional explicit validation set (benchmark-defined split). If omitted, a fixed random
    # split of --npz_path (val_ratio, split_seed) is used, identical for every model seed.
    p.add_argument("--val_npz_path", default=None)
    p.add_argument("--val_disorder_txt", default=None)
    p.add_argument("--val_protein_txt", default=None)
    p.add_argument("--val_rna_txt", default=None)
    p.add_argument("--val_dna_txt", default=None)
    p.add_argument("--out_dir", default="./checkpoints")
    p.add_argument("--seeds", type=int, nargs="+", default=None, help="override the seeds in the config")
    p.add_argument("--lr", type=float, default=None)
    p.add_argument("--weight_decay", type=float, default=None)
    p.add_argument("--batch_size", type=int, default=None)
    p.add_argument("--max_epochs", type=int, default=None)
    p.add_argument("--patience", type=int, default=None)
    p.add_argument("--plot", action="store_true")
    args = p.parse_args()

    cfg = copy.deepcopy(load_config(args.config))
    for k in ["lr", "weight_decay", "batch_size", "max_epochs", "patience"]:
        if getattr(args, k) is not None:
            cfg["optim"][k] = getattr(args, k)
    seeds = args.seeds or cfg["random_seeds"]
    device = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"device={device} | seeds={seeds} | optim={cfg['optim']}")

    os.makedirs(os.path.join(args.out_dir, args.name), exist_ok=True)
    with open(os.path.join(args.out_dir, args.name, "config_used.yaml"), "w") as f:
        yaml.safe_dump({**cfg, "random_seeds": seeds, "args": vars(args)}, f)

    labels = {"disorder": args.disorder_txt, "protein": args.protein_txt,
              "rna": args.rna_txt, "dna": args.dna_txt}
    test_labels = {"disorder": args.test_disorder_txt, "protein": args.test_protein_txt,
                   "rna": args.test_rna_txt, "dna": args.test_dna_txt}
    bs, nw = cfg["optim"]["batch_size"], cfg.get("num_workers", 2)

    if args.val_npz_path:
        val_labels = {"disorder": args.val_disorder_txt, "protein": args.val_protein_txt,
                      "rna": args.val_rna_txt, "dna": args.val_dna_txt}
        if any(v is None for v in val_labels.values()):
            p.error("--val_npz_path requires all four --val_*_txt label files")
        train_ids, val_npz, val_ids = None, args.val_npz_path, None
    else:
        train_ids, val_ids = split_ids(load_npz_ids(args.npz_path), cfg["val_ratio"], cfg["split_seed"])
        val_npz, val_labels = args.npz_path, labels
        log(f"Random train/val split (split_seed={cfg['split_seed']}): train={len(train_ids)}, val={len(val_ids)}")

    results = []
    for seed in seeds:
        loaders = (
            build_loader(args.npz_path, labels, bs, shuffle=True, include_ids=train_ids, num_workers=nw, seed=seed),
            build_loader(val_npz, val_labels, bs, shuffle=False, include_ids=val_ids, num_workers=nw),
            build_loader(args.test_npz_path, test_labels, bs, shuffle=False, num_workers=nw),
        )
        results.append(run_seed(seed, cfg, args, device, loaders))

    table = summarize(results, os.path.join(args.out_dir, args.name, "summary.json"))
    log(f"TEST over {len(results)} seed(s), mean ± SD:\n{table}")


if __name__ == "__main__":
    main()
