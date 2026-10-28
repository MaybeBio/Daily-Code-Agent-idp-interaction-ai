from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import random
import sys
import time
from pathlib import Path

import numpy as np
import torch
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from utils import Status, atomic_json, sha256_file, utc_now
from model import PhaseLACE, BaselineHead, residual_regularizer, paper_regularizer
from protocol import CoverageCycledSampler, EarlyStopping, assert_finite_async, build_inner_split, make_length_batches, minimum_epochs_for_coverage, objective_parts, select_trial_ids


def digest_json(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def evaluate(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float | int]:
    if len(labels) == 0 or set(np.unique(labels)) != {0, 1}:
        raise ValueError("Evaluation requires both classes")
    if not np.isfinite(probabilities).all() or np.any((probabilities < 0) | (probabilities > 1)):
        raise ValueError("Invalid prediction probabilities")
    return {
        "n": int(len(labels)),
        "positive": int(labels.sum()),
        "ap": float(average_precision_score(labels, probabilities)),
        "auroc": float(roc_auc_score(labels, probabilities)),
        "brier": float(brier_score_loss(labels, probabilities)),
    }


def save_checkpoint(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("wb") as handle:
        torch.save(payload, handle)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def append_metric(path: Path, record: dict[str, object]) -> None:
    existing_epochs: set[int] = set()
    if path.exists():
        with path.open(encoding="utf-8") as handle:
            existing_epochs = {int(json.loads(line)["epoch"]) for line in handle if line.strip()}
    if int(record["epoch"]) not in existing_epochs:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
            handle.flush()
            os.fsync(handle.fileno())


def write_predictions(path: Path, rows: list[dict[str, str]], indices: np.ndarray, probabilities: np.ndarray) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row", "original_row", "accession", "sequence_sha256", "label", "fold", "probability"])
        for index, probability in zip(indices, probabilities):
            row = rows[int(index)]
            writer.writerow([int(index), row["original_row"], row["accession"], row["sequence_sha256"], int(row["label"]), int(row["fold"]), float(probability)])
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Train the frozen base head or the final PhaseLACE residual members.")
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stage", choices=["baseline", "phaselace"], required=True)
    parser.add_argument("--baseline-root", type=Path)
    parser.add_argument("--trial-ids-file", type=Path)
    parser.add_argument("--seeds", nargs="+", type=int, default=[711016, 711017, 711018])
    parser.add_argument("--folds", nargs="+", type=int, default=list(range(5)))
    parser.add_argument("--inner-seed", type=int, default=20260904)
    parser.add_argument("--validation-fraction", type=float, default=0.2)
    parser.add_argument("--max-epochs", type=int, default=80)
    parser.add_argument("--min-epochs", type=int, default=10)
    parser.add_argument("--patience", type=int, default=8)
    parser.add_argument("--min-delta", type=float, default=0.002)
    parser.add_argument("--learning-rate", type=float, default=0.003)
    parser.add_argument("--residual-scale", type=float, default=0.3)
    parser.add_argument("--residual-regularizer-target", type=float, default=0.001374)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--gpu-cache", action="store_true")
    parser.add_argument("--eval-residue-budget", type=int, default=65536)
    parser.add_argument("--eval-max-batch-size", type=int, default=128)
    parser.add_argument("--max-hours", type=float, default=36.0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--skip-cache-content-hash", action="store_true", help="Only for synthetic smoke fixtures.")
    args = parser.parse_args()
    args.variants = [args.stage]
    if len(set(args.variants)) != len(args.variants) or not set(args.folds) <= set(range(5)):
        raise ValueError("Duplicate variants or invalid folds")
    if args.max_epochs < args.min_epochs or args.max_hours <= 0 or args.learning_rate <= 0 or args.eval_residue_budget < 1 or args.eval_max_batch_size < 1:
        raise ValueError("Invalid training limits")
    residual_requested = any(variant != "baseline" for variant in args.variants)
    if residual_requested and args.baseline_root is None:
        raise ValueError("Residual variants require --baseline-root")
    if args.output.exists() and not args.resume:
        raise FileExistsError("Output exists; use --resume only for the same frozen run")
    requested_trials = None
    if args.trial_ids_file is not None:
        requested_trials = [line.strip() for line in args.trial_ids_file.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
    selected_trials = select_trial_ids(args.variants, args.seeds, args.folds, requested_trials)
    if not selected_trials:
        raise ValueError("No trials selected")
    selected_trial_set = set(selected_trials)
    args.output.mkdir(parents=True, exist_ok=True)
    from utils import flock, LOCK_EX, LOCK_NB, LOCK_UN

    lock_handle = (args.output / "runner.lock").open("a")
    flock(lock_handle, LOCK_EX | LOCK_NB)
    status = Status(args.output / "status.json", "validate_inputs")
    try:
        torch.set_num_threads(4)
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.benchmark = False
        torch.use_deterministic_algorithms(True)
        manifest_path = args.cache / "manifest.json"
        index_path = args.cache / "train_index.csv"
        manifest = json.loads(manifest_path.read_text())
        with index_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if sha256_file(index_path) != manifest["train_index_sha256"]:
            raise ValueError("Training index checksum mismatch")
        cache_name = "pooled_fp32.npy" if set(args.variants) == {"baseline"} else "raw_fp32.npy"
        if "baseline" in args.variants and residual_requested:
            raise ValueError("Run baseline and residual stages in separate output roots/cache types")
        cache_path = args.cache / cache_name
        if cache_path.stat().st_size != manifest["cache_bytes"]:
            raise ValueError("Cache size mismatch")
        expected_cache_hash = manifest.get("cache_sha256", manifest.get("pooled_sha256"))
        if expected_cache_hash is None:
            raise ValueError("Cache manifest has no content checksum")
        if not args.skip_cache_content_hash and sha256_file(cache_path) != expected_cache_hash:
            raise ValueError("Cache content checksum mismatch")
        cache = np.load(cache_path, mmap_mode="r", allow_pickle=False)
        if list(cache.shape) != manifest["cache_shape"] or cache.dtype != np.float32:
            raise ValueError("Cache schema mismatch")
        labels = np.asarray([int(row["label"]) for row in rows], dtype=np.int64)
        offsets = np.asarray([int(row["offset"]) for row in rows], dtype=np.int64)
        lengths = np.asarray([int(row["length"]) for row in rows], dtype=np.int64)
        source_hashes = {path.name: sha256_file(path) for path in sorted(Path(__file__).parent.glob("*.py")) if not path.name.startswith("test_")}
        scientific_config = {
            "schema": "pstp_protocol_v2",
            "variants": args.variants,
            "seeds": args.seeds,
            "folds": args.folds,
            "trial_ids": selected_trials,
            "inner_seed": args.inner_seed,
            "validation_fraction": args.validation_fraction,
            "max_epochs": args.max_epochs,
            "min_epochs": args.min_epochs,
            "early_stopping_coverage_gate": "at least one complete inner-training negative cycle",
            "patience": args.patience,
            "min_delta": args.min_delta,
            "learning_rate": args.learning_rate,
            "scheduler": {"name": "ReduceLROnPlateau", "mode": "max", "factor": 0.5, "patience": 3},
            "negative_sampling": "coverage-cycled 1:1 within inner train",
            "selection": "inner-validation AP; outer fold evaluated once after best checkpoint reload",
            "loss_baseline": "2 * BCE + 1e-4 * original paper regularizer",
            "residual_scale": args.residual_scale,
            "residual_regularizer_target": args.residual_regularizer_target,
            "evaluation_batching": {"length_bucketed": True, "residue_budget": args.eval_residue_budget, "max_batch_size": args.eval_max_batch_size},
            "cache_manifest_sha256": sha256_file(manifest_path),
            "source_sha256": source_hashes,
            "runtime": {"python": sys.version, "torch": torch.__version__, "numpy": np.__version__, "platform": platform.platform()},
        }
        config_hash = digest_json(scientific_config)
        config_path = args.output / "config.json"
        if config_path.exists():
            if json.loads(config_path.read_text()) != scientific_config:
                raise ValueError("Resume requires unchanged code, input, and scientific configuration")
        else:
            atomic_json(config_path, scientific_config)
        split_dir = args.output / "splits"
        split_dir.mkdir(exist_ok=True)
        splits: dict[int, dict[str, object]] = {}
        for fold in args.folds:
            split = build_inner_split(rows, fold, args.inner_seed, args.validation_fraction)
            path = split_dir / f"outer_fold{fold}.json"
            if path.exists() and json.loads(path.read_text()) != split:
                raise ValueError("Frozen inner split changed")
            if not path.exists():
                atomic_json(path, split)
            splits[fold] = split
        cache_gpu = None
        offsets_gpu = None
        lengths_gpu = None
        labels_gpu = None
        if args.gpu_cache:
            free, _ = torch.cuda.mem_get_info(args.device)
            if free < cache.nbytes + 3 * 1024**3:
                raise RuntimeError("GPU cache exceeds free memory with 3 GiB reserve")
            status.update(phase="load_gpu_cache", cache_gib=cache.nbytes / 1024**3)
            # The mmap is read-only, but this operation only reads it and creates an
            # independent GPU tensor. Avoiding a second 66 GiB host copy shortens startup.
            cache_gpu = torch.from_numpy(cache).to(args.device)
            offsets_gpu = torch.as_tensor(offsets, device=args.device)
            lengths_gpu = torch.as_tensor(lengths, device=args.device)
            labels_gpu = torch.as_tensor(labels, dtype=torch.float32, device=args.device)

        def features(index: int) -> torch.Tensor:
            start, end = int(offsets[index]), int(offsets[index] + lengths[index])
            if cache_gpu is not None:
                values = cache_gpu[start:end]
            else:
                values = torch.from_numpy(np.array(cache[start:end], copy=True)).to(args.device)
            return values.permute(1, 0, 2) if cache_name == "pooled_fp32.npy" else values

        def evaluate_indices(model: torch.nn.Module, indices: np.ndarray, phase: str) -> tuple[np.ndarray, dict[str, float | int]]:
            model.eval()
            probabilities = np.empty(len(indices), dtype=np.float64)
            positions = {int(index): position for position, index in enumerate(indices)}
            with torch.inference_mode():
                if cache_gpu is not None and offsets_gpu is not None and lengths_gpu is not None and hasattr(model, "forward_batch"):
                    batches = make_length_batches(indices, lengths, args.eval_residue_budget, args.eval_max_batch_size)
                    completed_count = 0
                    for batch in batches:
                        batch_indices = torch.as_tensor(batch, device=args.device)
                        batch_lengths = lengths_gpu[batch_indices]
                        padded_length = int(batch_lengths.max())
                        residue_positions = torch.arange(padded_length, device=args.device)
                        mask = residue_positions[None, :] < batch_lengths[:, None]
                        gather = offsets_gpu[batch_indices, None] + torch.minimum(residue_positions[None, :], batch_lengths[:, None] - 1)
                        gathered = cache_gpu[gather]
                        batch_features = gathered.permute(0, 2, 1, 3) if cache_name == "pooled_fp32.npy" else gathered
                        batch_probabilities = model.forward_batch(batch_features, mask).cpu().numpy()
                        for index, probability in zip(batch, batch_probabilities):
                            probabilities[positions[index]] = probability
                        completed_count += len(batch)
                        status.update(phase=phase, evaluation_completed=completed_count, evaluation_total=len(indices))
                else:
                    for position, index in enumerate(indices):
                        probabilities[position] = float(model(features(int(index))))
                        if position % 250 == 0:
                            status.update(phase=phase, evaluation_completed=position + 1, evaluation_total=len(indices))
            return probabilities, evaluate(labels[indices], probabilities)

        deadline = time.monotonic() + args.max_hours * 3600
        completed: list[dict[str, object]] = []
        results_path = args.output / "results.json"
        if results_path.exists():
            completed = json.loads(results_path.read_text())["trials"]
        total_trials = len(selected_trials)
        for seed in args.seeds:
            for fold in args.folds:
                split = splits[fold]
                inner_train = np.asarray(split["inner_train"], dtype=np.int64)
                inner_validation = np.asarray(split["inner_validation"], dtype=np.int64)
                outer_evaluation = np.asarray(split["outer_evaluation"], dtype=np.int64)
                positive = inner_train[labels[inner_train] == 1]
                negative = inner_train[labels[inner_train] == 0]
                effective_min_epochs = minimum_epochs_for_coverage(len(positive), len(negative), args.min_epochs, args.max_epochs)
                for variant in args.variants:
                    trial_id = f"{variant}_seed{seed}_fold{fold}"
                    if trial_id not in selected_trial_set:
                        continue
                    trial = args.output / trial_id
                    trial.mkdir(exist_ok=True)
                    completed_path = trial / "completed.json"
                    if completed_path.exists():
                        result = json.loads(completed_path.read_text())
                        if result["config_hash"] != config_hash:
                            raise ValueError("Completed trial config mismatch")
                        if not any(item["trial"] == trial_id for item in completed):
                            completed.append(result)
                        continue
                    random.seed(seed)
                    np.random.seed(seed)
                    torch.manual_seed(seed)
                    torch.cuda.manual_seed_all(seed)
                    baseline_checkpoint_sha256 = None
                    if variant == "baseline":
                        model: torch.nn.Module = BaselineHead().to(args.device)
                        regularizer_weight = 1e-4
                    else:
                        baseline_path = args.baseline_root / trial_id.replace(variant, "baseline", 1) / "best.pt"
                        baseline_checkpoint_sha256 = sha256_file(baseline_path)
                        baseline_checkpoint = torch.load(baseline_path, map_location=args.device, weights_only=False)
                        baseline = BaselineHead().to(args.device)
                        baseline.load_state_dict(baseline_checkpoint["model"])
                        model = PhaseLACE(baseline, args.residual_scale).to(args.device)
                        initial_regularizer = float(residual_regularizer(model).detach())
                        regularizer_weight = args.residual_regularizer_target / initial_regularizer
                    optimizer = torch.optim.Adam((parameter for parameter in model.parameters() if parameter.requires_grad), lr=args.learning_rate)
                    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="max", factor=0.5, patience=3)
                    stopping = EarlyStopping(effective_min_epochs, args.patience, args.min_delta)
                    sampler = CoverageCycledSampler(positive, negative, seed + 1009 * fold)
                    start_epoch = 0
                    last_path = trial / "last.pt"
                    best_path = trial / "best.pt"
                    if last_path.exists():
                        checkpoint = torch.load(last_path, map_location=args.device, weights_only=False)
                        if checkpoint["config_hash"] != config_hash or checkpoint["split_sha256"] != split["sha256"]:
                            raise ValueError("Checkpoint config/split mismatch")
                        if checkpoint.get("baseline_checkpoint_sha256") != baseline_checkpoint_sha256:
                            raise ValueError("Frozen baseline checkpoint changed")
                        model.load_state_dict(checkpoint["model"])
                        optimizer.load_state_dict(checkpoint["optimizer"])
                        scheduler.load_state_dict(checkpoint["scheduler"])
                        stopping.load_state_dict(checkpoint["early_stopping"])
                        sampler.load_state_dict(checkpoint["sampler"])
                        regularizer_weight = float(checkpoint["regularizer_weight"])
                        start_epoch = int(checkpoint["completed_epochs"])
                        if checkpoint.get("last_record"):
                            append_metric(trial / "metrics.jsonl", checkpoint["last_record"])
                    trial_started = time.monotonic()
                    for epoch_index in range(start_epoch, args.max_epochs):
                        epoch = epoch_index + 1
                        order = sampler.next_epoch()
                        model.train()
                        total_names = ("objective", "bce", "regularizer", "scaled_regularizer")
                        totals_gpu = torch.zeros(len(total_names), device=args.device)
                        epoch_started = time.monotonic()
                        status.update(phase="training", trial=trial_id, epoch=epoch, max_epochs=args.max_epochs, completed_trials=len(completed), total_trials=total_trials, samples=len(order))
                        order_hash = hashlib.sha256(np.asarray(order, dtype=np.int64).tobytes()).hexdigest()
                        for step, index in enumerate(order):
                            optimizer.zero_grad(set_to_none=True)
                            probability = model(features(index))
                            target = labels_gpu[index].to(probability.dtype) if labels_gpu is not None else probability.new_tensor(float(labels[index]))
                            regularizer = paper_regularizer(model) if variant == "baseline" else residual_regularizer(model)
                            parts = objective_parts(probability, target, regularizer, regularizer_weight)
                            loss = parts["objective"]
                            assert_finite_async(loss, f"loss in {trial_id}, epoch {epoch}, step {step}")
                            loss.backward()
                            optimizer.step()
                            with torch.no_grad():
                                totals_gpu += torch.stack([parts[name].detach() for name in total_names])
                            if step % 100 == 0:
                                status.update(step=step + 1)
                        validation_probabilities, validation_metrics = evaluate_indices(model, inner_validation, "inner_validation")
                        if sampler.selection_ready:
                            decision = stopping.step(epoch, float(validation_metrics["ap"]))
                            scheduler.step(float(validation_metrics["ap"]))
                        else:
                            decision = {"improved": False, "should_stop": False}
                        total_values = totals_gpu.cpu().tolist()
                        record: dict[str, object] = {
                            "trial": trial_id,
                            "epoch": epoch,
                            "train": {name: value / len(order) for name, value in zip(total_names, total_values)},
                            "inner_validation": validation_metrics,
                            "learning_rate": optimizer.param_groups[0]["lr"],
                            "regularizer_weight": regularizer_weight,
                            "negative_coverage": sampler.coverage_fraction,
                            "negative_unique": sampler.unique_negative_count,
                            "negative_total": len(negative),
                            "effective_min_epochs": effective_min_epochs,
                            "sampling_cycle": sampler.cycle,
                            "selection_active": sampler.selection_ready,
                            "order_sha256": order_hash,
                            "improved": decision["improved"],
                            "bad_epochs": stopping.bad_epochs,
                            "seconds": time.monotonic() - epoch_started,
                            "timestamp": utc_now(),
                        }
                        checkpoint_payload: dict[str, object] = {
                            "model": model.state_dict(),
                            "optimizer": optimizer.state_dict(),
                            "scheduler": scheduler.state_dict(),
                            "early_stopping": stopping.state_dict(),
                            "sampler": sampler.state_dict(),
                            "completed_epochs": epoch,
                            "config_hash": config_hash,
                            "split_sha256": split["sha256"],
                            "variant": variant,
                            "seed": seed,
                            "fold": fold,
                            "regularizer_weight": regularizer_weight,
                            "baseline_checkpoint_sha256": baseline_checkpoint_sha256,
                            "last_record": record,
                        }
                        if decision["improved"]:
                            save_checkpoint(best_path, checkpoint_payload)
                            write_predictions(trial / "best_inner_predictions.csv", rows, inner_validation, validation_probabilities)
                        save_checkpoint(last_path, checkpoint_payload)
                        append_metric(trial / "metrics.jsonl", record)
                        print(json.dumps(record, allow_nan=False), flush=True)
                        if stopping.should_stop:
                            break
                        if time.monotonic() > deadline:
                            status.close("PAUSED_BUDGET", trial=trial_id, completed_epochs=epoch, resume="same command with --resume")
                            return
                    if not best_path.exists():
                        raise RuntimeError("Best checkpoint was not created")
                    best = torch.load(best_path, map_location=args.device, weights_only=False)
                    model.load_state_dict(best["model"])
                    save_checkpoint(trial / "inference.pt", {"model": model.state_dict(), "stage": variant, "seed": seed, "fold": fold, "residual_scale": args.residual_scale})
                    status.update(phase="outer_evaluation_once", trial=trial_id, best_epoch=best["completed_epochs"])
                    outer_probabilities, outer_metrics = evaluate_indices(model, outer_evaluation, "outer_evaluation_once")
                    write_predictions(trial / "outer_predictions.csv", rows, outer_evaluation, outer_probabilities)
                    result = {
                        "trial": trial_id,
                        "variant": variant,
                        "seed": seed,
                        "fold": fold,
                        "epochs_completed": int(torch.load(last_path, map_location="cpu", weights_only=False)["completed_epochs"]),
                        "best_epoch": int(best["completed_epochs"]),
                        "best_inner_ap": float(best["early_stopping"]["best_score"]),
                        "outer_metrics": outer_metrics,
                        "trainable_parameters": sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad),
                        "regularizer_weight": regularizer_weight,
                        "effective_min_epochs": effective_min_epochs,
                        "split_sha256": split["sha256"],
                        "config_hash": config_hash,
                        "best_checkpoint_sha256": sha256_file(best_path),
                        "baseline_checkpoint_sha256": baseline_checkpoint_sha256,
                        "seconds_this_attempt": time.monotonic() - trial_started,
                        "finished_at": utc_now(),
                    }
                    atomic_json(completed_path, result)
                    completed.append(result)
                    atomic_json(results_path, {"trials": completed, "expected_trials": total_trials})
                    print(json.dumps(result, allow_nan=False), flush=True)
                    if time.monotonic() > deadline and len(completed) < total_trials:
                        status.close("PAUSED_BUDGET", completed_trials=len(completed), total_trials=total_trials, resume="same command with --resume")
                        return
        status.close("COMPLETED", phase="complete", completed_trials=len(completed), total_trials=total_trials)
    except BaseException as exc:
        status.close("FAILED", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        flock(lock_handle, LOCK_UN)
        lock_handle.close()


if __name__ == "__main__":
    main()
