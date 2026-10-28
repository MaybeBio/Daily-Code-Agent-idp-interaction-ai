#!/usr/bin/env python3
"""Compute condensate intra-chain contact-number relaxation.

The observable is the number of nonlocal intra-chain residue pairs with
distance < 0.66 nm. Pairs with |i-j| <= 2 are absent from the existing
combined_features.npy arrays. Each trajectory is clipped to the same bulk-chain
time window and selected part range used by the corrected inter-chain analysis.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INTERCHAIN_DIR = ROOT / "storyline_figures_for_ppt_corrected_v3"
DEFAULT_OUTDIR = INTERCHAIN_DIR / "intra_inter_contact_relaxation"
SYSTEMS = ("FUS", "A1", "TDP43")
ROLES = ("Tl", "Th")
TEMP_COLORS = {"Tl": "#2f7f73", "Th": "#b54b5f"}
BAR_COLORS = {"intrachain": "#2f7f73", "interchain": "#c96d28"}


@dataclass(frozen=True)
class Condition:
    system: str
    temp: str
    role: str
    feature_path: Path
    lengths_path: Path
    trajectory_dir: Path
    metadata_path: Path
    interchain_cache: Path


def condition(
    system: str,
    temp: str,
    role: str,
    tica_dir: str,
    intermediate_name: str,
    metadata_path: str,
    interchain_cache: str,
) -> Condition:
    tica = ROOT / tica_dir
    intermediate = tica / intermediate_name
    return Condition(
        system=system,
        temp=temp,
        role=role,
        feature_path=intermediate / "combined_features.npy",
        lengths_path=intermediate / "combined_tica_lengths.npy",
        trajectory_dir=tica / "extracted_500ns_bulk_chains",
        metadata_path=ROOT / metadata_path,
        interchain_cache=ROOT / interchain_cache,
    )


CONDITIONS = (
    condition(
        "FUS",
        "280K",
        "Tl",
        "fus/fus_protein_single_simulation/condensate/tica_100dimension_10frames_280K",
        "intermediate_fus",
        "fus/fus_condensate_simulation/condensate_280K_200chains/extracted_500ns_bulk_chains/segments_meta.json",
        "ppt_5_11_figures/results/interchain_contact_relaxation_FUS_280K_parts081-150.npz",
    ),
    condition(
        "FUS",
        "325K",
        "Th",
        "fus/fus_protein_single_simulation/condensate/tica_100dimension_10frames_325K",
        "intermediate_fus",
        "fus/fus_condensate_simulation/condensate_325K_200chains/extracted_500ns_bulk_chains/segments_meta.json",
        "ppt_5_11_figures/results/interchain_contact_relaxation_FUS_325K_parts081-150.npz",
    ),
    condition(
        "A1",
        "230K",
        "Tl",
        "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_230K",
        "intermediate_a1",
        "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_230K/extracted_500ns_bulk_chains/segments_meta.json",
        "ppt_5_11_figures/results/interchain_contact_relaxation_A1_230K_parts057-126.npz",
    ),
    condition(
        "A1",
        "260K",
        "Th",
        "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_260K",
        "intermediate_a1",
        "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_260K/extracted_500ns_bulk_chains/segments_meta.json",
        "ppt_5_11_figures/results/interchain_contact_relaxation_A1_260K_parts033-102.npz",
    ),
    condition(
        "TDP43",
        "260K",
        "Tl",
        "TDP43/tdp43_protein_single_simulation/condensate/tica_100dimension_10frames_260K",
        "intermediate_tdp43",
        "TDP43/tdp43_protein_single_simulation/condensate/tica_100dimension_10frames_260K/extracted_500ns_bulk_chains/segments_meta.json",
        "ppt_5_11_figures/results/interchain_contact_relaxation_TDP43_260K_parts081-150.npz",
    ),
    condition(
        "TDP43",
        "300K",
        "Th",
        "TDP43/tdp43_protein_single_simulation/condensate/tica_100dimension_10frames_300K",
        "intermediate_tdp43",
        "TDP43/tdp43_condensate_simulation/condensate_300K_200chains/extracted_500ns_bulk_chains/segments_meta.json",
        "ppt_5_11_figures/results/interchain_contact_relaxation_TDP43_300K_parts021-090.npz",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--cutoff-nm", type=float, default=0.66)
    parser.add_argument("--threshold", type=float, default=0.2)
    parser.add_argument("--full-dt-ns", type=float, default=0.01)
    parser.add_argument("--full-to-phase-stride", type=int, default=10)
    parser.add_argument(
        "--phase-stride",
        type=int,
        default=1,
        help="Sample every phase-DCD interval. Default 1 gives 0.1 ns resolution.",
    )
    parser.add_argument("--max-lag-ns", type=float, default=30.0)
    parser.add_argument("--n-bootstrap", type=int, default=400)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def chain_id_from_path(path: Path) -> int:
    match = re.search(r"chain(\d+)_bulk_500ns\.dcd$", path.name)
    if not match:
        raise ValueError(f"Cannot parse chain ID from {path}")
    return int(match.group(1))


def autocorr_fft(values: np.ndarray) -> np.ndarray | None:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size < 4:
        return None
    arr = arr - np.mean(arr)
    var = float(np.dot(arr, arr) / arr.size)
    if not np.isfinite(var) or var <= 1e-12:
        return None
    n = arr.size
    transformed = np.fft.fft(arr, n=2 * n)
    acf = np.fft.ifft(transformed * np.conjugate(transformed))[:n].real
    return acf / (var * np.arange(n, 0, -1, dtype=np.float64))


def crossing_time(time_ns: np.ndarray, acf: np.ndarray, threshold: float) -> float:
    indices = np.where(acf <= threshold)[0]
    if indices.size == 0:
        return float("nan")
    index = int(indices[0])
    if index == 0:
        return float(time_ns[0])
    x0, x1 = float(time_ns[index - 1]), float(time_ns[index])
    y0, y1 = float(acf[index - 1]), float(acf[index])
    if abs(y1 - y0) < 1e-12:
        return x1
    return x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)


def bootstrap_crossing_ci(
    acfs: np.ndarray,
    time_ns: np.ndarray,
    threshold: float,
    n_bootstrap: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_bootstrap):
        indices = rng.integers(0, acfs.shape[0], size=acfs.shape[0])
        value = crossing_time(time_ns, np.mean(acfs[indices], axis=0), threshold)
        if np.isfinite(value):
            values.append(value)
    if not values:
        return float("nan"), float("nan")
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def sampled_local_indices(
    segment: dict,
    trajectory_length: int,
    first_part: int,
    last_part: int,
    frames_per_part: int,
    full_to_phase_stride: int,
    phase_stride: int,
) -> np.ndarray:
    phase_start = int(math.ceil(int(segment["start"]) / full_to_phase_stride))
    phase_end = int(math.floor(int(segment["end"]) / full_to_phase_stride))
    selected_start = (first_part - 1) * frames_per_part
    selected_end = last_part * frames_per_part - 1
    start = max(phase_start, selected_start)
    end = min(phase_end, selected_end)
    if end < start:
        return np.empty(0, dtype=np.int64)
    first = start + (phase_stride - start % phase_stride) % phase_stride
    phase_frames = np.arange(first, end + 1, phase_stride, dtype=np.int64)
    local = phase_frames * full_to_phase_stride - int(segment["start"])
    return local[(local >= 0) & (local < trajectory_length)]


def compute_condition(cond: Condition, args: argparse.Namespace) -> dict:
    out = args.outdir / f"intrachain_contact_relaxation_{cond.system}_{cond.temp}.npz"
    if out.exists() and not args.force:
        loaded = np.load(out, allow_pickle=False)
        expected_dt = (
            float(args.full_dt_ns)
            * int(args.full_to_phase_stride)
            * int(args.phase_stride)
        )
        cache_matches = (
            "dt_sample_ns" in loaded.files
            and np.isclose(float(loaded["dt_sample_ns"]), expected_dt)
            and np.isclose(float(loaded["cutoff_nm"]), float(args.cutoff_nm))
            and np.isclose(float(loaded["threshold"]), float(args.threshold))
        )
        if cache_matches:
            return {key: loaded[key] for key in loaded.files}
        loaded.close()

    for path in (
        cond.feature_path,
        cond.lengths_path,
        cond.metadata_path,
        cond.interchain_cache,
    ):
        if not path.exists():
            raise FileNotFoundError(path)

    feature = np.load(cond.feature_path, mmap_mode="r")
    lengths = np.load(cond.lengths_path).astype(np.int64)
    trajectories = sorted(cond.trajectory_dir.glob("chain*_bulk_500ns.dcd"))
    if len(trajectories) != len(lengths):
        raise RuntimeError(
            f"{cond.system} {cond.temp}: {len(trajectories)} trajectories but "
            f"{len(lengths)} feature lengths"
        )

    with cond.metadata_path.open("r", encoding="utf-8") as handle:
        segments = json.load(handle)
    segment_by_chain = {int(item["chain_id"]): item for item in segments}

    inter = np.load(cond.interchain_cache, allow_pickle=False)
    first_part = int(inter["first_part"])
    last_part = int(inter["last_part"])
    frames_per_part = int(inter["frames_per_part"])
    dt_sample_ns = float(args.full_dt_ns) * int(args.full_to_phase_stride) * int(args.phase_stride)
    max_lag_frames = int(round(float(args.max_lag_ns) / dt_sample_ns))

    acfs = []
    segment_relax = []
    contact_means = []
    used_chain_ids = []
    cursor = 0
    for trajectory, length in zip(trajectories, lengths):
        chain_id = chain_id_from_path(trajectory)
        segment = segment_by_chain.get(chain_id)
        if segment is None:
            cursor += int(length)
            continue

        local_indices = sampled_local_indices(
            segment,
            int(length),
            first_part,
            last_part,
            frames_per_part,
            int(args.full_to_phase_stride),
            int(args.phase_stride),
        )
        if local_indices.size < max(8, max_lag_frames + 1):
            cursor += int(length)
            continue

        global_indices = cursor + local_indices
        distances = np.asarray(feature[global_indices], dtype=np.float32)
        counts = np.count_nonzero(distances < float(args.cutoff_nm), axis=1).astype(np.float64)
        del distances

        acf = autocorr_fft(counts)
        if acf is not None:
            acf = acf[: max_lag_frames + 1]
            if acf.size == max_lag_frames + 1:
                time_ns = np.arange(acf.size, dtype=np.float64) * dt_sample_ns
                acfs.append(acf)
                segment_relax.append(crossing_time(time_ns, acf, float(args.threshold)))
                contact_means.append(float(np.mean(counts)))
                used_chain_ids.append(chain_id)
        cursor += int(length)

    if not acfs:
        raise RuntimeError(f"No usable intra-chain ACFs for {cond.system} {cond.temp}")

    acf_stack = np.stack(acfs, axis=0)
    time_ns = np.arange(acf_stack.shape[1], dtype=np.float64) * dt_sample_ns
    acf_mean = np.mean(acf_stack, axis=0)
    acf_sem = np.std(acf_stack, axis=0, ddof=1) / math.sqrt(acf_stack.shape[0])
    relax_time_ns = crossing_time(time_ns, acf_mean, float(args.threshold))
    ci_low, ci_high = bootstrap_crossing_ci(
        acf_stack,
        time_ns,
        float(args.threshold),
        int(args.n_bootstrap),
        seed=20260626 + len(cond.system) * 10 + int(cond.temp.rstrip("K")),
    )

    payload = {
        "system": np.asarray(cond.system),
        "temp": np.asarray(cond.temp),
        "role": np.asarray(cond.role),
        "time_ns": time_ns.astype(np.float32),
        "acf_mean": acf_mean.astype(np.float32),
        "acf_sem": acf_sem.astype(np.float32),
        "segment_relax_ns": np.asarray(segment_relax, dtype=np.float32),
        "segment_contact_mean": np.asarray(contact_means, dtype=np.float32),
        "used_chain_ids": np.asarray(used_chain_ids, dtype=np.int32),
        "relax_time_ns": np.asarray(relax_time_ns, dtype=np.float32),
        "relax_ci_low_ns": np.asarray(ci_low, dtype=np.float32),
        "relax_ci_high_ns": np.asarray(ci_high, dtype=np.float32),
        "mean_contact_number": np.asarray(np.mean(contact_means), dtype=np.float32),
        "n_segments": np.asarray(acf_stack.shape[0], dtype=np.int32),
        "dt_sample_ns": np.asarray(dt_sample_ns, dtype=np.float32),
        "cutoff_nm": np.asarray(args.cutoff_nm, dtype=np.float32),
        "threshold": np.asarray(args.threshold, dtype=np.float32),
        "sequence_gap_excluded": np.asarray(2, dtype=np.int32),
        "first_part": np.asarray(first_part, dtype=np.int32),
        "last_part": np.asarray(last_part, dtype=np.int32),
        "feature_path": np.asarray(str(cond.feature_path.relative_to(ROOT))),
    }
    np.savez_compressed(out, **payload)
    print(
        f"[save] {cond.system} {cond.temp}: tau={relax_time_ns:.3f} ns "
        f"CI=[{ci_low:.3f}, {ci_high:.3f}], n={acf_stack.shape[0]} -> {out}",
        flush=True,
    )
    return payload


def scalar(payload: dict, key: str) -> float:
    return float(np.asarray(payload[key]).item())


def save_intrachain_tables(results: dict[tuple[str, str], dict], outdir: Path) -> None:
    summary_fields = [
        "system",
        "temp",
        "role",
        "relax_time_ns",
        "relax_ci_low_ns",
        "relax_ci_high_ns",
        "mean_contact_number",
        "n_segments",
        "first_part",
        "last_part",
        "threshold",
        "cutoff_nm",
        "sequence_gap_excluded",
    ]
    with (outdir / "intrachain_contact_relaxation_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        for system in SYSTEMS:
            for role in ROLES:
                item = results[(system, role)]
                writer.writerow(
                    {
                        key: (
                            str(np.asarray(item[key]).item())
                            if key in item
                            else ""
                        )
                        for key in summary_fields
                    }
                )

    ratio_fields = [
        "system",
        "tl_temp",
        "th_temp",
        "tau_tl_ns",
        "tau_th_ns",
        "ratio_th_over_tl",
        "ratio_ci_low",
        "ratio_ci_high",
    ]
    with (outdir / "intrachain_contact_relaxation_ratio.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=ratio_fields)
        writer.writeheader()
        for system in SYSTEMS:
            low = results[(system, "Tl")]
            high = results[(system, "Th")]
            tau_low = scalar(low, "relax_time_ns")
            tau_high = scalar(high, "relax_time_ns")
            writer.writerow(
                {
                    "system": system,
                    "tl_temp": str(np.asarray(low["temp"]).item()),
                    "th_temp": str(np.asarray(high["temp"]).item()),
                    "tau_tl_ns": tau_low,
                    "tau_th_ns": tau_high,
                    "ratio_th_over_tl": tau_high / tau_low,
                    "ratio_ci_low": scalar(high, "relax_ci_low_ns")
                    / scalar(low, "relax_ci_high_ns"),
                    "ratio_ci_high": scalar(high, "relax_ci_high_ns")
                    / scalar(low, "relax_ci_low_ns"),
                }
            )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save_curves(results: dict[tuple[str, str], dict], outdir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.25), sharex=True, sharey=True)
    for ax, system in zip(axes, SYSTEMS):
        for role in ROLES:
            item = results[(system, role)]
            time = np.asarray(item["time_ns"], dtype=float)
            mean = np.asarray(item["acf_mean"], dtype=float)
            sem = np.asarray(item["acf_sem"], dtype=float)
            label = f"{np.asarray(item['temp']).item()} ({scalar(item, 'relax_time_ns'):.2f} ns)"
            ax.plot(time, mean, color=TEMP_COLORS[role], lw=2, label=label)
            ax.fill_between(time, mean - sem, mean + sem, color=TEMP_COLORS[role], alpha=0.18)
            ax.axvline(
                scalar(item, "relax_time_ns"),
                color=TEMP_COLORS[role],
                lw=1,
                ls=":",
                alpha=0.65,
            )
        ax.axhline(0.2, color="#555555", lw=1, ls="--")
        ax.set_title(system, fontweight="bold")
        ax.set_xlabel("Lag time (ns)")
        ax.legend(frameon=False, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Intrachain contact-number ACF")
    fig.suptitle("Condensate intrachain contact-number relaxation", fontweight="bold")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"intrachain_contact_relaxation_curves.{ext}", dpi=300)
    plt.close(fig)


def save_two_bar_figure(outdir: Path) -> None:
    intra_rows = {row["system"]: row for row in read_csv(outdir / "intrachain_contact_relaxation_ratio.csv")}
    inter_rows = {
        row["system"]: row
        for row in read_csv(INTERCHAIN_DIR / "interchain_contact_relaxation_ratio_corrected.csv")
    }

    x = np.arange(len(SYSTEMS), dtype=float)
    width = 0.31
    intra_values = np.array([float(intra_rows[s]["ratio_th_over_tl"]) for s in SYSTEMS])
    inter_values = np.array([float(inter_rows[s]["ratio_th_over_tl"]) for s in SYSTEMS])
    intra_low = np.array([float(intra_rows[s]["ratio_ci_low"]) for s in SYSTEMS])
    intra_high = np.array([float(intra_rows[s]["ratio_ci_high"]) for s in SYSTEMS])
    inter_low = np.array([float(inter_rows[s]["ratio_ci_low"]) for s in SYSTEMS])
    inter_high = np.array([float(inter_rows[s]["ratio_ci_high"]) for s in SYSTEMS])

    fig, ax = plt.subplots(figsize=(8.1, 4.35))
    bars_intra = ax.bar(
        x - width / 2,
        intra_values,
        width,
        color=BAR_COLORS["intrachain"],
        label="Intrachain contact relaxation",
        yerr=np.vstack((intra_values - intra_low, intra_high - intra_values)),
        capsize=3,
        error_kw={"lw": 1.1, "ecolor": "#333333"},
    )
    bars_inter = ax.bar(
        x + width / 2,
        inter_values,
        width,
        color=BAR_COLORS["interchain"],
        label="Interchain contact relaxation",
        yerr=np.vstack((inter_values - inter_low, inter_high - inter_values)),
        capsize=3,
        error_kw={"lw": 1.1, "ecolor": "#333333"},
    )
    ax.axhline(1.0, color="#333333", lw=1.1, ls="--")
    ax.set_xticks(x, SYSTEMS, fontweight="bold")
    ax.set_ylabel(r"Relaxation ratio $\tau(T_h)/\tau(T_l)$")
    ax.set_title(
        "Heating decouples intrachain and interchain relaxation",
        fontweight="bold",
        fontsize=13,
    )
    ax.legend(frameon=False, ncol=2, loc="upper center")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.18)

    for bars in (bars_intra, bars_inter):
        for bar in bars:
            value = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + 0.04,
                f"{value:.2f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ymax = max(float(np.max(intra_high)), float(np.max(inter_high)), 2.0)
    ax.set_ylim(0, ymax * 1.10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"intra_inter_contact_relaxation_ratio.{ext}", dpi=300)
    plt.close(fig)

    source_fields = [
        "system",
        "observable",
        "ratio_th_over_tl",
        "ratio_ci_low",
        "ratio_ci_high",
    ]
    with (outdir / "intra_inter_contact_relaxation_ratio_source.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=source_fields)
        writer.writeheader()
        for system in SYSTEMS:
            for observable, rows in (("intrachain", intra_rows), ("interchain", inter_rows)):
                row = rows[system]
                writer.writerow(
                    {
                        "system": system,
                        "observable": observable,
                        "ratio_th_over_tl": row["ratio_th_over_tl"],
                        "ratio_ci_low": row["ratio_ci_low"],
                        "ratio_ci_high": row["ratio_ci_high"],
                    }
                )


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    results = {}
    for cond in CONDITIONS:
        results[(cond.system, cond.role)] = compute_condition(cond, args)
    save_intrachain_tables(results, args.outdir)
    save_curves(results, args.outdir)
    save_two_bar_figure(args.outdir)
    print(f"[done] Results written to {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
