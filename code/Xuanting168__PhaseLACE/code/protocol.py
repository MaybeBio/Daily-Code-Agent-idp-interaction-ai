from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict

import numpy as np
import torch
from torch.nn import functional as F


def _digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def minimum_epochs_for_coverage(positive_count: int, negative_count: int, requested_minimum: int, maximum: int) -> int:
    if min(positive_count, negative_count, requested_minimum, maximum) < 1:
        raise ValueError("Counts and epoch limits must be positive")
    required = math.ceil(negative_count / positive_count)
    effective = max(requested_minimum, required)
    if effective > maximum:
        raise ValueError("Maximum epochs cannot complete one negative-coverage cycle")
    return effective


def make_length_batches(
    indices: np.ndarray,
    lengths: np.ndarray,
    residue_budget: int,
    max_batch_size: int,
) -> list[list[int]]:
    if residue_budget < 1 or max_batch_size < 1:
        raise ValueError("Batch limits must be positive")
    ordered = sorted((int(index) for index in indices), key=lambda index: (int(lengths[index]), index))
    batches: list[list[int]] = []
    current: list[int] = []
    current_max = 0
    for index in ordered:
        candidate_max = max(current_max, int(lengths[index]))
        if current and (len(current) >= max_batch_size or candidate_max * (len(current) + 1) > residue_budget):
            batches.append(current)
            current = []
            current_max = 0
        current.append(index)
        current_max = max(current_max, int(lengths[index]))
    if current:
        batches.append(current)
    return batches


def assert_finite_async(value: torch.Tensor, context: str) -> None:
    finite = torch.isfinite(value)
    if value.is_cuda:
        torch._assert_async(finite, f"Non-finite value: {context}")
    elif not bool(finite):
        raise FloatingPointError(f"Non-finite value: {context}")


def select_trial_ids(variants: list[str], seeds: list[int], folds: list[int], requested: list[str] | None = None) -> list[str]:
    available = [f"{variant}_seed{seed}_fold{fold}" for seed in seeds for fold in folds for variant in variants]
    if requested is None:
        return available
    if len(requested) != len(set(requested)):
        raise ValueError("Duplicate trial assignment")
    unknown = set(requested) - set(available)
    if unknown:
        raise ValueError(f"Unknown assigned trials: {sorted(unknown)}")
    requested_set = set(requested)
    return [trial_id for trial_id in available if trial_id in requested_set]


def build_inner_split(
    rows: list[dict[str, object]],
    outer_fold: int,
    seed: int,
    validation_fraction: float = 0.2,
) -> dict[str, object]:
    """Create one frozen inner split without touching the outer fold."""
    if outer_fold not in range(5) or not 0 < validation_fraction < 0.5:
        raise ValueError("Invalid outer fold or validation fraction")
    candidates = [index for index, row in enumerate(rows) if int(row["fold"]) != outer_fold]
    lengths = np.asarray([int(rows[index]["length"]) for index in candidates])
    boundaries = np.quantile(lengths, [0.25, 0.5, 0.75])
    strata: dict[str, list[int]] = defaultdict(list)
    for index in candidates:
        row = rows[index]
        length_bin = int(np.searchsorted(boundaries, int(row["length"]), side="right"))
        key = f'{int(row["label"])}|{row.get("category", "unknown")}|q{length_bin}'
        strata[key].append(index)
    rng = np.random.RandomState(seed + 1009 * outer_fold)
    validation: list[int] = []
    for key in sorted(strata):
        members = np.asarray(sorted(strata[key]), dtype=np.int64)
        rng.shuffle(members)
        if len(members) < 2:
            continue
        count = min(len(members) - 1, max(1, int(round(len(members) * validation_fraction))))
        validation.extend(int(value) for value in members[:count])
    validation_set = set(validation)
    training = sorted(index for index in candidates if index not in validation_set)
    validation = sorted(validation)
    if not training or not validation:
        raise ValueError("Inner split is empty")
    for name, indices in (("inner train", training), ("inner validation", validation)):
        if {int(rows[index]["label"]) for index in indices} != {0, 1}:
            raise ValueError(f"{name} does not contain both labels")
    split = {
        "schema": "pstp_inner_split_v2",
        "outer_fold": outer_fold,
        "seed": seed,
        "validation_fraction": validation_fraction,
        "length_quantile_boundaries": [float(value) for value in boundaries],
        "inner_train": training,
        "inner_validation": validation,
        "outer_evaluation": [index for index, row in enumerate(rows) if int(row["fold"]) == outer_fold],
    }
    split["sha256"] = _digest(split)
    return split


class CoverageCycledSampler:
    """Draw 1:1 negatives and finish each background cycle before repeating it."""

    def __init__(self, positive: list[int] | np.ndarray, negative: list[int] | np.ndarray, seed: int) -> None:
        self.positive = np.asarray(positive, dtype=np.int64)
        self.negative = np.asarray(negative, dtype=np.int64)
        if not len(self.positive) or len(self.negative) < len(self.positive):
            raise ValueError("Sampler requires positives and at least as many negatives")
        if set(self.positive) & set(self.negative):
            raise ValueError("Positive and negative pools overlap")
        self.rng = np.random.RandomState(seed)
        self.order = self.rng.permutation(self.negative)
        self.cursor = 0
        self.epoch = 0
        self.cycle = 0
        self.seen: set[int] = set()

    @property
    def unique_negative_count(self) -> int:
        return len(self.seen)

    @property
    def coverage_fraction(self) -> float:
        return len(self.seen) / len(self.negative)

    @property
    def selection_ready(self) -> bool:
        return len(self.seen) == len(self.negative)

    def _take_negatives(self, count: int) -> list[int]:
        result: list[int] = []
        while len(result) < count:
            available = len(self.order) - self.cursor
            take = min(count - len(result), available)
            result.extend(int(value) for value in self.order[self.cursor : self.cursor + take])
            self.cursor += take
            if self.cursor == len(self.order):
                self.cycle += 1
                self.order = self.rng.permutation(self.negative)
                self.cursor = 0
        self.seen.update(result)
        return result

    def next_epoch(self) -> list[int]:
        negatives = self._take_negatives(len(self.positive))
        values = np.concatenate((self.positive, np.asarray(negatives, dtype=np.int64)))
        self.rng.shuffle(values)
        self.epoch += 1
        return [int(value) for value in values]

    def state_dict(self) -> dict[str, object]:
        return {
            "positive": self.positive.copy(),
            "negative": self.negative.copy(),
            "rng_state": self.rng.get_state(),
            "order": self.order.copy(),
            "cursor": self.cursor,
            "epoch": self.epoch,
            "cycle": self.cycle,
            "seen": sorted(self.seen),
        }

    def load_state_dict(self, state: dict[str, object]) -> None:
        if not np.array_equal(self.positive, state["positive"]) or not np.array_equal(self.negative, state["negative"]):
            raise ValueError("Sampler pools changed")
        self.rng.set_state(state["rng_state"])
        self.order = np.asarray(state["order"], dtype=np.int64)
        self.cursor = int(state["cursor"])
        self.epoch = int(state["epoch"])
        self.cycle = int(state["cycle"])
        self.seen = {int(value) for value in state["seen"]}


class EarlyStopping:
    def __init__(self, min_epochs: int, patience: int, min_delta: float) -> None:
        if min_epochs < 1 or patience < 1 or min_delta < 0:
            raise ValueError("Invalid early-stopping settings")
        self.min_epochs = min_epochs
        self.patience = patience
        self.min_delta = min_delta
        self.best_score = -float("inf")
        self.best_epoch = 0
        self.bad_epochs = 0
        self.should_stop = False

    def step(self, epoch: int, score: float) -> dict[str, object]:
        if not np.isfinite(score):
            raise ValueError("Early-stopping score must be finite")
        improved = score >= self.best_score + self.min_delta
        if improved:
            self.best_score = float(score)
            self.best_epoch = epoch
            self.bad_epochs = 0
        else:
            self.bad_epochs += 1
        self.should_stop = epoch >= self.min_epochs and self.bad_epochs >= self.patience
        return {"improved": improved, "should_stop": self.should_stop}

    def state_dict(self) -> dict[str, object]:
        return dict(self.__dict__)

    def load_state_dict(self, state: dict[str, object]) -> None:
        expected = {"min_epochs", "patience", "min_delta", "best_score", "best_epoch", "bad_epochs", "should_stop"}
        if set(state) != expected:
            raise ValueError("Invalid early-stopping state")
        self.__dict__.update(state)


def objective_parts(
    probability: torch.Tensor,
    target: torch.Tensor,
    regularizer: torch.Tensor,
    regularizer_weight: float,
) -> dict[str, torch.Tensor]:
    bce = F.binary_cross_entropy(probability, target)
    scaled_regularizer = regularizer_weight * regularizer
    return {"bce": bce, "regularizer": regularizer, "scaled_regularizer": scaled_regularizer, "objective": 2 * bce + scaled_regularizer}
