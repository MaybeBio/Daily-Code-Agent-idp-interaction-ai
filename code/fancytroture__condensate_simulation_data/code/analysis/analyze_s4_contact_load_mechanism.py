#!/usr/bin/env python3
"""Analyze per-partner bead-contact loads for supplementary Figure S4."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE_DIR = ROOT / "source_data/supplementary/s4_contact_load_cache"
DEFAULT_OUTDIR = ROOT / "source_data/supplementary/s4_contact_load_analysis"
SYSTEMS = ("FUS", "A1", "TDP43")
ROLES = ("Tl", "Th")
TEMPS = {
    "FUS": {"Tl": "280K", "Th": "325K"},
    "A1": {"Tl": "230K", "Th": "260K"},
    "TDP43": {"Tl": "260K", "Th": "300K"},
}
LAGS_NS = np.asarray([0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0])
EPSILONS = (0.05, 0.10)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260723)
    return parser.parse_args()


def find_cache(cache_dir: Path, system: str, role: str) -> Path:
    temp = TEMPS[system][role]
    matches = sorted(
        cache_dir.glob(f"partner_contact_load_{system}_{temp}_parts*_pstride005.npz")
    )
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one S4 cache for {system} {role}, found {len(matches)}"
        )
    return matches[0]


def percentile_interval(values: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan"), float("nan")
    return tuple(float(value) for value in np.percentile(finite, [2.5, 97.5]))


def bootstrap_equal_mean(
    values: np.ndarray,
    indices: np.ndarray,
) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    output = np.full(indices.shape[0], np.nan, dtype=np.float64)
    for boot_index, selected in enumerate(indices):
        sample = arr[selected]
        if np.any(np.isfinite(sample)):
            output[boot_index] = float(np.nanmean(sample))
    return output


def bootstrap_ratio_of_sums(
    numerator: np.ndarray,
    denominator: np.ndarray,
    indices: np.ndarray,
    *,
    complement: bool = False,
) -> np.ndarray:
    numerator = np.asarray(numerator, dtype=np.float64)
    denominator = np.asarray(denominator, dtype=np.float64)
    output = np.full(indices.shape[0], np.nan, dtype=np.float64)
    for boot_index, selected in enumerate(indices):
        den = float(np.nansum(denominator[selected]))
        if den <= 0:
            continue
        value = float(np.nansum(numerator[selected])) / den
        output[boot_index] = 1.0 - value if complement else value
    return output


def add_metric(
    rows: list[dict[str, object]],
    system: str,
    role: str,
    lag_ns: float,
    metric: str,
    value: float,
    bootstrap: np.ndarray,
    n_segments: int,
) -> None:
    low, high = percentile_interval(bootstrap)
    rows.append(
        {
            "system": system,
            "role": role,
            "lag_ns": lag_ns,
            "metric": metric,
            "value": value,
            "ci_low": low,
            "ci_high": high,
            "n_segments": n_segments,
        }
    )


def analyze_condition(
    cache_path: Path,
    system: str,
    role: str,
    args: argparse.Namespace,
) -> tuple[
    list[dict[str, object]],
    dict[str, np.ndarray],
    dict[tuple[float, str], np.ndarray],
]:
    with np.load(cache_path, allow_pickle=False) as data:
        contact_load = np.asarray(data["contact_load"], dtype=np.int32)
        valid = np.asarray(data["valid"], dtype=bool)
        dt_ns = float(np.asarray(data["dt_ns"]).reshape(-1)[0])
        cutoff_nm = float(np.asarray(data["cutoff_nm"]).reshape(-1)[0])

    n_segments = contact_load.shape[0]
    rng = np.random.default_rng(
        int(args.seed) + sum(ord(char) for char in f"{system}-{role}")
    )
    bootstrap_indices = rng.integers(
        0,
        n_segments,
        size=(int(args.n_bootstrap), n_segments),
    )
    rows: list[dict[str, object]] = []
    distribution: dict[str, np.ndarray] = {}
    segment_metrics: dict[tuple[float, str], np.ndarray] = {}

    for lag_ns in LAGS_NS:
        lag_frames = int(round(float(lag_ns) / dt_ns))
        if not np.isclose(lag_frames * dt_ns, lag_ns):
            raise ValueError(f"Lag {lag_ns} ns is not representable at dt={dt_ns}")

        comp_abs_delta = np.zeros(n_segments)
        comp_turnover = np.zeros(n_segments)
        comp_window_mean = np.full(n_segments, np.nan)
        q_load = np.full(n_segments, np.nan)
        partner_survival = np.full(n_segments, np.nan)
        pair_count = np.zeros(n_segments)
        conditional_sum = {
            epsilon: np.zeros(n_segments)
            for epsilon in EPSILONS
        }
        conditional_count = {
            epsilon: np.zeros(n_segments)
            for epsilon in EPSILONS
        }
        distribution_lost: list[np.ndarray] = []
        distribution_gained: list[np.ndarray] = []

        exit_loss = np.zeros(n_segments)
        enter_gain = np.zeros(n_segments)
        retained_loss = np.zeros(n_segments)
        retained_gain = np.zeros(n_segments)

        for segment_index in range(n_segments):
            segment_valid = valid[segment_index]
            if lag_frames >= segment_valid.size:
                continue
            pair_mask = (
                segment_valid[:-lag_frames]
                & segment_valid[lag_frames:]
            )
            if not np.any(pair_mask):
                continue

            start = contact_load[segment_index, :-lag_frames][pair_mask]
            future = contact_load[segment_index, lag_frames:][pair_mask]
            difference = future - start
            lost = np.maximum(-difference, 0).sum(axis=1, dtype=np.float64)
            gained = np.maximum(difference, 0).sum(axis=1, dtype=np.float64)
            total_start = start.sum(axis=1, dtype=np.float64)
            total_future = future.sum(axis=1, dtype=np.float64)
            turnover = lost + gained
            absolute_delta = np.abs(total_future - total_start)

            start_present = start > 0
            future_present = future > 0
            exiting = start_present & ~future_present
            entering = ~start_present & future_present
            retained = start_present & future_present
            exit_loss[segment_index] = float(np.sum(start[exiting]))
            enter_gain[segment_index] = float(np.sum(future[entering]))
            retained_loss[segment_index] = float(
                np.sum(np.maximum(-difference, 0)[retained])
            )
            retained_gain[segment_index] = float(
                np.sum(np.maximum(difference, 0)[retained])
            )

            comp_abs_delta[segment_index] = float(np.sum(absolute_delta))
            comp_turnover[segment_index] = float(np.sum(turnover))
            positive_turnover = turnover > 0
            if np.any(positive_turnover):
                comp_window_mean[segment_index] = float(
                    np.mean(
                        1.0
                        - absolute_delta[positive_turnover]
                        / turnover[positive_turnover]
                    )
                )

            intersection = np.sum(start_present & future_present, axis=1)
            union = np.sum(start_present | future_present, axis=1)
            initial_partners = np.sum(start_present, axis=1)
            replacement = np.full(union.size, np.nan, dtype=np.float64)
            survival = np.full(union.size, np.nan, dtype=np.float64)
            union_positive = union > 0
            initial_positive = initial_partners > 0
            replacement[union_positive] = (
                1.0 - intersection[union_positive] / union[union_positive]
            )
            survival[initial_positive] = (
                intersection[initial_positive] / initial_partners[initial_positive]
            )
            if np.any(np.isfinite(survival)):
                partner_survival[segment_index] = float(np.nanmean(survival))

            denominator = total_start + total_future
            symmetric_change = np.full(denominator.size, np.nan)
            denominator_positive = denominator > 0
            symmetric_change[denominator_positive] = (
                2.0
                * absolute_delta[denominator_positive]
                / denominator[denominator_positive]
            )
            pair_count[segment_index] = float(np.sum(np.isfinite(replacement)))
            for epsilon in EPSILONS:
                selected = (
                    np.isfinite(replacement)
                    & np.isfinite(symmetric_change)
                    & (symmetric_change <= epsilon)
                )
                conditional_count[epsilon][segment_index] = float(np.sum(selected))
                conditional_sum[epsilon][segment_index] = float(
                    np.nansum(replacement[selected])
                )

            all_totals = contact_load[segment_index].sum(axis=1, dtype=np.float64)
            finite_totals = all_totals[segment_valid]
            standard_deviation = float(np.std(finite_totals))
            if standard_deviation > 1e-12:
                z_values = (
                    all_totals - float(np.mean(finite_totals))
                ) / standard_deviation
                states = np.full(z_values.size, -1, dtype=np.int8)
                states[z_values < -0.5] = 0
                states[(z_values >= -0.5) & (z_values <= 0.5)] = 1
                states[z_values > 0.5] = 2
                finite_states = states[segment_valid]
                populations = np.asarray(
                    [np.mean(finite_states == state) for state in range(3)]
                )
                baseline = float(np.sum(populations**2))
                same = states[:-lag_frames][pair_mask] == states[lag_frames:][pair_mask]
                if baseline < 1.0:
                    q_load[segment_index] = (
                        float(np.mean(same)) - baseline
                    ) / (1.0 - baseline)

            if np.isclose(lag_ns, 2.0):
                mean_load = float(np.mean(finite_totals))
                if mean_load > 0:
                    distribution_lost.append(lost / mean_load)
                    distribution_gained.append(gained / mean_load)

        comp_value = (
            1.0 - np.sum(comp_abs_delta) / np.sum(comp_turnover)
            if np.sum(comp_turnover) > 0
            else float("nan")
        )
        comp_bootstrap = bootstrap_ratio_of_sums(
            comp_abs_delta,
            comp_turnover,
            bootstrap_indices,
            complement=True,
        )
        add_metric(
            rows,
            system,
            role,
            float(lag_ns),
            "contact_load_compensation_pooled",
            comp_value,
            comp_bootstrap,
            n_segments,
        )
        segment_metrics[
            (float(lag_ns), "load_state_persistence_q")
        ] = q_load.copy()
        segment_metrics[
            (float(lag_ns), "partner_survival")
        ] = partner_survival.copy()

        add_metric(
            rows,
            system,
            role,
            float(lag_ns),
            "contact_load_compensation_window_mean",
            float(np.nanmean(comp_window_mean)),
            bootstrap_equal_mean(comp_window_mean, bootstrap_indices),
            n_segments,
        )
        add_metric(
            rows,
            system,
            role,
            float(lag_ns),
            "load_state_persistence_q",
            float(np.nanmean(q_load)),
            bootstrap_equal_mean(q_load, bootstrap_indices),
            n_segments,
        )
        add_metric(
            rows,
            system,
            role,
            float(lag_ns),
            "partner_survival",
            float(np.nanmean(partner_survival)),
            bootstrap_equal_mean(partner_survival, bootstrap_indices),
            n_segments,
        )

        for epsilon in EPSILONS:
            count = conditional_count[epsilon]
            replacement_value = (
                np.sum(conditional_sum[epsilon]) / np.sum(count)
                if np.sum(count) > 0
                else float("nan")
            )
            replacement_bootstrap = bootstrap_ratio_of_sums(
                conditional_sum[epsilon],
                count,
                bootstrap_indices,
            )
            add_metric(
                rows,
                system,
                role,
                float(lag_ns),
                f"conditional_replacement_eps{epsilon:.2f}",
                replacement_value,
                replacement_bootstrap,
                n_segments,
            )
            selected_fraction = (
                np.sum(count) / np.sum(pair_count)
                if np.sum(pair_count) > 0
                else float("nan")
            )
            fraction_bootstrap = bootstrap_ratio_of_sums(
                count,
                pair_count,
                bootstrap_indices,
            )
            add_metric(
                rows,
                system,
                role,
                float(lag_ns),
                f"condition_fraction_eps{epsilon:.2f}",
                selected_fraction,
                fraction_bootstrap,
                n_segments,
            )

        if np.isclose(lag_ns, 2.0):
            components = {
                "exit_loss": exit_loss,
                "enter_gain": enter_gain,
                "retained_loss": retained_loss,
                "retained_gain": retained_gain,
            }
            component_total = sum(np.sum(values) for values in components.values())
            for metric, values in components.items():
                value = (
                    float(np.sum(values)) / component_total
                    if component_total > 0
                    else float("nan")
                )
                component_bootstrap = bootstrap_ratio_of_sums(
                    values,
                    sum(components.values()),
                    bootstrap_indices,
                )
                add_metric(
                    rows,
                    system,
                    role,
                    float(lag_ns),
                    f"turnover_fraction_{metric}",
                    value,
                    component_bootstrap,
                    n_segments,
                )
            distribution["lost_normalized"] = np.concatenate(distribution_lost)
            distribution["gained_normalized"] = np.concatenate(distribution_gained)

    distribution["cutoff_nm"] = np.asarray(cutoff_nm)
    distribution["dt_ns"] = np.asarray(dt_ns)
    return rows, distribution, segment_metrics


def temperature_difference_rows(
    segment_metrics: dict[
        tuple[str, str],
        dict[tuple[float, str], np.ndarray],
    ],
    args: argparse.Namespace,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for system in SYSTEMS:
        for metric in ("load_state_persistence_q", "partner_survival"):
            for lag_ns in LAGS_NS:
                low_values = segment_metrics[(system, "Tl")][
                    (float(lag_ns), metric)
                ]
                high_values = segment_metrics[(system, "Th")][
                    (float(lag_ns), metric)
                ]
                rng = np.random.default_rng(
                    int(args.seed)
                    + sum(
                        ord(char)
                        for char in f"temperature-difference-{system}-{metric}-{lag_ns}"
                    )
                )
                bootstrap = np.full(int(args.n_bootstrap), np.nan)
                for index in range(int(args.n_bootstrap)):
                    low_sample = low_values[
                        rng.integers(0, low_values.size, low_values.size)
                    ]
                    high_sample = high_values[
                        rng.integers(0, high_values.size, high_values.size)
                    ]
                    bootstrap[index] = (
                        float(np.nanmean(high_sample))
                        - float(np.nanmean(low_sample))
                    )
                value = (
                    float(np.nanmean(high_values))
                    - float(np.nanmean(low_values))
                )
                add_metric(
                    rows,
                    system,
                    "Th-Tl",
                    float(lag_ns),
                    f"temperature_difference_{metric}",
                    value,
                    bootstrap,
                    int(low_values.size + high_values.size),
                )
    return rows


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    distributions: dict[str, np.ndarray] = {}
    segment_metrics: dict[
        tuple[str, str],
        dict[tuple[float, str], np.ndarray],
    ] = {}

    for system in SYSTEMS:
        for role in ROLES:
            (
                condition_rows,
                condition_distribution,
                condition_segment_metrics,
            ) = analyze_condition(
                find_cache(args.cache_dir, system, role),
                system,
                role,
                args,
            )
            rows.extend(condition_rows)
            for key, values in condition_distribution.items():
                distributions[f"{system}_{role}_{key}"] = values
            segment_metrics[(system, role)] = condition_segment_metrics

    rows.extend(temperature_difference_rows(segment_metrics, args))

    summary_path = args.outdir / "FigureS4_contact_load_mechanism_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    distribution_path = args.outdir / "FigureS4_contact_load_distributions.npz"
    np.savez_compressed(distribution_path, **distributions)
    print(f"[save] {summary_path}")
    print(f"[save] {distribution_path}")


if __name__ == "__main__":
    main()
