#!/usr/bin/env python3
"""Bootstrap the interchain contact-number t0.2 temperature-response ratio."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE_DIR = PROJECT_ROOT / "source_data/main/interchain_ratio_bootstrap"
DEFAULT_OUTPUT = DEFAULT_CACHE_DIR / "interchain_contact_relaxation_ratio_bootstrap.csv"

CONDITIONS = {
    "FUS": {
        "Tl": "interchain_contact_relaxation_FUS_280K_parts081-150.npz",
        "Th": "interchain_contact_relaxation_FUS_325K_parts081-150.npz",
    },
    "A1": {
        "Tl": "interchain_contact_relaxation_A1_230K_parts057-126.npz",
        "Th": "interchain_contact_relaxation_A1_260K_parts033-102.npz",
    },
    "TDP43": {
        "Tl": "interchain_contact_relaxation_TDP43_260K_parts081-150.npz",
        "Th": "interchain_contact_relaxation_TDP43_300K_parts021-090.npz",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--n-bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260717)
    return parser.parse_args()


def crossing_time(time_ns: np.ndarray, acf: np.ndarray, threshold: float) -> float:
    indices = np.flatnonzero(acf <= threshold)
    if indices.size == 0:
        return float("nan")
    idx = int(indices[0])
    if idx == 0:
        return float(time_ns[0])
    x0, x1 = float(time_ns[idx - 1]), float(time_ns[idx])
    y0, y1 = float(acf[idx - 1]), float(acf[idx])
    if np.isclose(y0, y1):
        return x1
    return x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)


def load_cache(path: Path) -> tuple[np.ndarray, np.ndarray, float]:
    with np.load(path, allow_pickle=False) as data:
        if "segment_acfs" not in data.files:
            raise KeyError(f"{path} has no segment_acfs; rerun the trajectory calculation")
        time_ns = np.asarray(data["time_ns"], dtype=np.float64)
        acfs = np.asarray(data["segment_acfs"], dtype=np.float64)
        threshold = float(np.asarray(data["threshold"]).reshape(-1)[0])
    if acfs.ndim != 2 or acfs.shape[1] != time_ns.size:
        raise ValueError(f"Inconsistent segment ACF shape in {path}: {acfs.shape}")
    return time_ns, acfs, threshold


def bootstrap_ratio(
    tl_time: np.ndarray,
    tl_acfs: np.ndarray,
    th_time: np.ndarray,
    th_acfs: np.ndarray,
    threshold: float,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> np.ndarray:
    ratios = np.full(n_bootstrap, np.nan, dtype=np.float64)
    for boot_idx in range(n_bootstrap):
        tl_idx = rng.integers(0, tl_acfs.shape[0], size=tl_acfs.shape[0])
        th_idx = rng.integers(0, th_acfs.shape[0], size=th_acfs.shape[0])
        tl_tau = crossing_time(tl_time, np.nanmean(tl_acfs[tl_idx], axis=0), threshold)
        th_tau = crossing_time(th_time, np.nanmean(th_acfs[th_idx], axis=0), threshold)
        if np.isfinite(tl_tau) and np.isfinite(th_tau) and tl_tau > 0:
            ratios[boot_idx] = th_tau / tl_tau
    return ratios[np.isfinite(ratios)]


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    rows: list[dict[str, str]] = []

    for system, paths in CONDITIONS.items():
        tl_time, tl_acfs, tl_threshold = load_cache(args.cache_dir / paths["Tl"])
        th_time, th_acfs, th_threshold = load_cache(args.cache_dir / paths["Th"])
        if not np.isclose(tl_threshold, th_threshold):
            raise ValueError(f"Threshold mismatch for {system}: {tl_threshold} vs {th_threshold}")

        tl_tau = crossing_time(tl_time, np.nanmean(tl_acfs, axis=0), tl_threshold)
        th_tau = crossing_time(th_time, np.nanmean(th_acfs, axis=0), th_threshold)
        ratios = bootstrap_ratio(
            tl_time,
            tl_acfs,
            th_time,
            th_acfs,
            tl_threshold,
            args.n_bootstrap,
            rng,
        )
        if ratios.size < 0.95 * args.n_bootstrap:
            raise RuntimeError(
                f"Only {ratios.size}/{args.n_bootstrap} finite bootstrap ratios for {system}"
            )
        low, high = np.percentile(ratios, [2.5, 97.5])
        rows.append(
            {
                "system": system,
                "tau_tl_ns": f"{tl_tau:.8f}",
                "tau_th_ns": f"{th_tau:.8f}",
                "ratio_th_over_tl": f"{th_tau / tl_tau:.8f}",
                "ratio_ci_low": f"{low:.8f}",
                "ratio_ci_high": f"{high:.8f}",
                "threshold": f"{tl_threshold:.6f}",
                "n_segments_tl": str(tl_acfs.shape[0]),
                "n_segments_th": str(th_acfs.shape[0]),
                "n_bootstrap": str(args.n_bootstrap),
                "n_finite_bootstrap": str(ratios.size),
                "method": "independent_segment_acf_bootstrap_of_tH_over_tL_ratio",
            }
        )
        print(
            f"{system}: ratio={th_tau / tl_tau:.4f}, "
            f"95% CI=[{low:.4f}, {high:.4f}] "
            f"(segments {tl_acfs.shape[0]}/{th_acfs.shape[0]})"
        )

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
