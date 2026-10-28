#!/usr/bin/env python3
"""Compute inter-chain molecular partner survival in condensate bulk chains.

The observable is the survival of molecular partners around each selected
bulk-chain segment:

    S_partner(tau) = < |P_i(t) & P_i(t + tau)| / |P_i(t)| >

where P_i(t) is the set of other chain IDs contacting chain i at time t.  The
sampling windows, part-DCD selection, cutoff, and bulk-chain metadata are kept
parallel to compute_interchain_contact_relaxation_parts.py so the result can be
compared directly with the inter-chain contact-number ACF.
"""

from __future__ import annotations

import argparse
import copy
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from scipy.spatial import cKDTree

from compute_interchain_contact_relaxation_parts import (
    CONDITIONS,
    ROLE_LABEL,
    ROOT,
    SYSTEM_ORDER,
    TEMP_COLORS,
    bootstrap_crossing_ci,
    build_schedule,
    choose_part_selection,
    crossing_time,
    load_segments,
    text_value,
)


OUTDIR_DEFAULT = ROOT / "ppt_5_11_figures" / "results"
INTERCHAIN_SUMMARY_DEFAULT = ROOT / "storyline_figures_for_ppt_corrected_v3" / "interchain_contact_relaxation_summary_corrected.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=OUTDIR_DEFAULT)
    parser.add_argument("--systems", default="FUS,A1,TDP43")
    parser.add_argument("--roles", default="Tl,Th")
    parser.add_argument("--phase-dt-ns", type=float, default=0.1)
    parser.add_argument("--phase-stride", type=int, default=5)
    parser.add_argument("--full-to-phase-stride", type=int, default=10)
    parser.add_argument("--part-duration-ns", type=float, default=100.0)
    parser.add_argument("--frames-per-part", type=int, default=0)
    parser.add_argument("--last-parts", type=int, default=0)
    parser.add_argument("--start-part", type=int, default=0)
    parser.add_argument("--end-part", type=int, default=0)
    parser.add_argument("--auto-fraction", type=float, default=0.25)
    parser.add_argument("--auto-min-parts", type=int, default=10)
    parser.add_argument("--auto-max-parts", type=int, default=30)
    parser.add_argument("--cutoff-nm", type=float, default=0.66)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--max-lag-ns", type=float, default=30.0)
    parser.add_argument("--min-segment-ns", type=float, default=0.0)
    parser.add_argument("--chunk", type=int, default=20)
    parser.add_argument("--n-chains", type=int, default=200)
    parser.add_argument("--chain-id-base", choices=["auto", "0", "1"], default="auto")
    parser.add_argument("--max-segments", type=int, default=0)
    parser.add_argument("--n-bootstrap", type=int, default=400)
    parser.add_argument(
        "--interchain-summary",
        type=Path,
        default=INTERCHAIN_SUMMARY_DEFAULT,
        help="Use first_part/last_part from the existing Figure 5 interchain summary when available.",
    )
    parser.add_argument(
        "--no-match-interchain-window",
        action="store_true",
        help="Do not override part windows from --interchain-summary.",
    )
    parser.add_argument(
        "--compute-only",
        action="store_true",
        help="Only compute condition caches; skip shared CSV/figure writing.",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def select_conditions(args: argparse.Namespace):
    systems = {x.strip() for x in args.systems.split(",") if x.strip()}
    roles = {x.strip() for x in args.roles.split(",") if x.strip()}
    return [cond for cond in CONDITIONS if cond.system in systems and cond.role in roles]


def interchain_window_for(cond, args: argparse.Namespace) -> tuple[int, int] | None:
    if args.no_match_interchain_window or not args.interchain_summary.exists():
        return None
    with args.interchain_summary.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("system") == cond.system and row.get("temp") == cond.temp and row.get("role") == cond.role:
                try:
                    return int(row["first_part"]), int(row["last_part"])
                except (KeyError, TypeError, ValueError):
                    return None
    return None


def partner_sets_for_chains(
    positions_nm: np.ndarray,
    box_nm: np.ndarray,
    chain_len: int,
    chain_ids: list[int],
    cutoff_nm: float,
    atom_chain_index: np.ndarray,
) -> dict[int, set[int]]:
    if box_nm is None or len(box_nm) != 3 or np.any(np.asarray(box_nm) <= 0):
        raise ValueError(f"Invalid unitcell_lengths: {box_nm}")
    pos = np.asarray(positions_nm, dtype=np.float64) % box_nm
    tree = cKDTree(pos, boxsize=box_nm)
    out: dict[int, set[int]] = {}
    for chain_id in chain_ids:
        a0 = int(chain_id) * chain_len
        a1 = (int(chain_id) + 1) * chain_len
        neighbors = tree.query_ball_point(pos[a0:a1], r=cutoff_nm)
        partners: set[int] = set()
        for nb in neighbors:
            if not nb:
                continue
            nb_arr = np.asarray(nb, dtype=np.int64)
            partner_ids = atom_chain_index[nb_arr]
            for partner in partner_ids[partner_ids != chain_id]:
                partners.add(int(partner))
        out[int(chain_id)] = partners
    return out


def segment_partner_survival(
    partner_series: list[set[int] | None],
    max_lag_frames: int,
) -> tuple[np.ndarray, np.ndarray]:
    values = np.full(max_lag_frames + 1, np.nan, dtype=np.float64)
    counts = np.zeros(max_lag_frames + 1, dtype=np.int32)
    n = len(partner_series)
    for lag in range(max_lag_frames + 1):
        vals = []
        stop = n - lag
        for i in range(stop):
            current = partner_series[i]
            future = partner_series[i + lag]
            if current is None or future is None or len(current) == 0:
                continue
            vals.append(len(current.intersection(future)) / float(len(current)))
        if vals:
            values[lag] = float(np.mean(vals))
            counts[lag] = len(vals)
    return values, counts


def bootstrap_survival_t50(
    segment_survival: np.ndarray,
    time_ns: np.ndarray,
    threshold: float,
    n_boot: int,
    seed: int,
) -> tuple[float, float]:
    finite_segments = segment_survival[np.any(np.isfinite(segment_survival), axis=1)]
    if finite_segments.shape[0] < 2:
        value = crossing_time(time_ns, np.nanmean(finite_segments, axis=0), threshold)
        return value, value
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, finite_segments.shape[0], size=finite_segments.shape[0])
        vals.append(crossing_time(time_ns, np.nanmean(finite_segments[idx], axis=0), threshold))
    vals = np.asarray(vals, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def cache_name(cond, selection, args: argparse.Namespace) -> Path:
    max_lag_tag = str(float(args.max_lag_ns)).replace(".", "p")
    return args.outdir / (
        f"interchain_partner_survival_{cond.system}_{cond.temp}_"
        f"parts{selection.first_part:03d}-{selection.last_part:03d}_"
        f"pstride{int(args.phase_stride):03d}_lag{max_lag_tag}ns.npz"
    )


def compute_condition(cond, args: argparse.Namespace) -> dict[str, np.ndarray]:
    for path in (cond.topology, cond.metadata):
        if not path.exists():
            raise FileNotFoundError(path)

    local_args = copy.copy(args)
    window = interchain_window_for(cond, args)
    if window is not None:
        local_args.start_part, local_args.end_part = window
    selection = choose_part_selection(cond, local_args)
    for part in selection.selected_parts:
        if not part.path.exists():
            raise FileNotFoundError(part.path)

    out = cache_name(cond, selection, args)
    if out.exists() and not args.force:
        data = np.load(out, allow_pickle=False)
        return {key: data[key] for key in data.files}

    segments = load_segments(cond.metadata, args.max_segments, args.n_chains, args.chain_id_base)
    schedule, segment_frames, _clipped_bounds, min_frame, max_frame = build_schedule(
        segments,
        full_to_phase_stride=int(args.full_to_phase_stride),
        phase_stride=int(args.phase_stride),
        selected_global_start=selection.global_start,
        selected_global_end=selection.global_end,
    )
    if not schedule:
        raise RuntimeError(f"No sampled frames for {cond.system} {cond.temp}")

    first_frame = md.load_frame(str(selection.selected_parts[0].path), 0, top=str(cond.topology))
    n_atoms = first_frame.n_atoms
    if n_atoms % int(args.n_chains) != 0:
        raise RuntimeError(f"Atom count {n_atoms} is not divisible by n_chains={args.n_chains}")
    chain_len = n_atoms // int(args.n_chains)
    atom_chain_index = np.repeat(np.arange(int(args.n_chains), dtype=np.int32), chain_len)

    partner_series: list[list[set[int] | None]] = [[None] * len(frames) for frames in segment_frames]
    frame_to_pos = [{frame: i for i, frame in enumerate(frames)} for frames in segment_frames]
    dt_sample_ns = float(args.phase_dt_ns) * float(args.phase_stride)
    max_lag_frames = int(round(float(args.max_lag_ns) / dt_sample_ns))

    print(
        f"[compute] {cond.system} {cond.temp} parts={selection.first_part:03d}-{selection.last_part:03d} "
        f"segments={len(segments)} sampled_frames={len(schedule)} dt={dt_sample_ns:.3f} ns "
        f"chain_len={chain_len}",
        flush=True,
    )

    n_filled = 0
    for part in selection.selected_parts:
        part_global_start = (part.index - 1) * selection.frames_per_part
        part_global_end = part.index * selection.frames_per_part - 1
        if part_global_end < min_frame or part_global_start > max_frame:
            continue

        local_start = (int(args.phase_stride) - (part_global_start % int(args.phase_stride))) % int(args.phase_stride)
        local_counter = 0
        iterator = md.iterload(
            str(part.path),
            top=str(cond.topology),
            chunk=int(args.chunk),
            skip=int(local_start),
            stride=int(args.phase_stride),
        )
        for chunk in iterator:
            if chunk.n_frames == 0:
                continue
            for local_idx in range(chunk.n_frames):
                local_frame = local_start + local_counter * int(args.phase_stride)
                local_counter += 1
                global_frame = part_global_start + local_frame

                if global_frame < min_frame:
                    continue
                if global_frame > max_frame or global_frame > part_global_end:
                    break

                active = schedule.get(global_frame)
                if not active:
                    continue

                chains = sorted({int(segments[seg_idx]["chain_id"]) for seg_idx in active})
                partners = partner_sets_for_chains(
                    chunk.xyz[local_idx],
                    np.asarray(chunk.unitcell_lengths[local_idx], dtype=np.float64),
                    chain_len,
                    chains,
                    float(args.cutoff_nm),
                    atom_chain_index,
                )
                for seg_idx in active:
                    chain_id = int(segments[seg_idx]["chain_id"])
                    pos_idx = frame_to_pos[seg_idx][global_frame]
                    partner_series[seg_idx][pos_idx] = partners[chain_id]
                    n_filled += 1

    min_segment_points = int(round(float(args.min_segment_ns) / dt_sample_ns))
    segment_curves = []
    segment_counts = []
    initial_partner_means = []
    for series in partner_series:
        present = [item for item in series if item is not None]
        if len(present) < max(8, min_segment_points):
            continue
        curve, counts = segment_partner_survival(series, max_lag_frames)
        if np.any(np.isfinite(curve)):
            segment_curves.append(curve)
            segment_counts.append(counts)
            initial_partner_means.append(float(np.mean([len(item) for item in present])))

    if not segment_curves:
        raise RuntimeError(f"No usable partner-survival segments for {cond.system} {cond.temp}")

    segment_survival = np.stack(segment_curves, axis=0)
    segment_counts_arr = np.stack(segment_counts, axis=0)
    time_ns = np.arange(max_lag_frames + 1, dtype=np.float64) * dt_sample_ns
    survival_mean = np.nanmean(segment_survival, axis=0)
    survival_sem = np.nanstd(segment_survival, axis=0, ddof=1) / math.sqrt(max(segment_survival.shape[0], 1))
    t50 = crossing_time(time_ns, survival_mean, float(args.threshold))
    ci_low, ci_high = bootstrap_survival_t50(
        segment_survival,
        time_ns,
        float(args.threshold),
        int(args.n_bootstrap),
        seed=20260713 + len(cond.system) * 10 + int(cond.temp.rstrip("K")),
    )
    auc_0_max = float(np.trapz(np.nan_to_num(survival_mean, nan=0.0), time_ns))

    payload = {
        "system": np.asarray(cond.system),
        "temp": np.asarray(cond.temp),
        "role": np.asarray(cond.role),
        "time_ns": time_ns.astype(np.float32),
        "survival_mean": survival_mean.astype(np.float32),
        "survival_sem": survival_sem.astype(np.float32),
        "segment_survival": segment_survival.astype(np.float32),
        "segment_counts": segment_counts_arr.astype(np.int32),
        "t50_ns": np.asarray(t50, dtype=np.float32),
        "t50_ci_low_ns": np.asarray(ci_low, dtype=np.float32),
        "t50_ci_high_ns": np.asarray(ci_high, dtype=np.float32),
        "auc_0_max_ns": np.asarray(auc_0_max, dtype=np.float32),
        "mean_initial_partner_number": np.asarray(np.mean(initial_partner_means), dtype=np.float32),
        "n_segments": np.asarray(segment_survival.shape[0], dtype=np.int32),
        "n_filled": np.asarray(n_filled, dtype=np.int32),
        "first_part": np.asarray(selection.first_part, dtype=np.int32),
        "last_part": np.asarray(selection.last_part, dtype=np.int32),
        "phase_stride": np.asarray(args.phase_stride, dtype=np.int32),
        "phase_dt_ns": np.asarray(args.phase_dt_ns, dtype=np.float32),
        "cutoff_nm": np.asarray(args.cutoff_nm, dtype=np.float32),
        "threshold": np.asarray(args.threshold, dtype=np.float32),
        "chain_len": np.asarray(chain_len, dtype=np.int32),
    }
    np.savez_compressed(out, **payload)
    print(
        f"[save] {cond.system} {cond.temp}: t50={t50:.3f} ns "
        f"CI=[{ci_low:.3f}, {ci_high:.3f}], AUC={auc_0_max:.3f}, "
        f"n={segment_survival.shape[0]} -> {out}",
        flush=True,
    )
    return payload


def scalar(value) -> float:
    arr = np.asarray(value)
    if arr.shape == ():
        return float(arr)
    return float(arr.ravel()[0])


def write_summary(outdir: Path, results: list[dict[str, np.ndarray]]) -> Path:
    path = outdir / "interchain_partner_survival_summary.csv"
    fieldnames = [
        "system",
        "temp",
        "role",
        "t50_ns",
        "t50_ci_low_ns",
        "t50_ci_high_ns",
        "auc_0_max_ns",
        "mean_initial_partner_number",
        "n_segments",
        "n_filled",
        "first_part",
        "last_part",
        "phase_stride",
        "phase_dt_ns",
        "cutoff_nm",
        "threshold",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in results:
            writer.writerow(
                {
                    "system": text_value(item["system"]),
                    "temp": text_value(item["temp"]),
                    "role": text_value(item["role"]),
                    "t50_ns": f"{scalar(item['t50_ns']):.8g}",
                    "t50_ci_low_ns": f"{scalar(item['t50_ci_low_ns']):.8g}",
                    "t50_ci_high_ns": f"{scalar(item['t50_ci_high_ns']):.8g}",
                    "auc_0_max_ns": f"{scalar(item['auc_0_max_ns']):.8g}",
                    "mean_initial_partner_number": f"{scalar(item['mean_initial_partner_number']):.8g}",
                    "n_segments": str(int(scalar(item["n_segments"]))),
                    "n_filled": str(int(scalar(item["n_filled"]))),
                    "first_part": str(int(scalar(item["first_part"]))),
                    "last_part": str(int(scalar(item["last_part"]))),
                    "phase_stride": str(int(scalar(item["phase_stride"]))),
                    "phase_dt_ns": f"{scalar(item['phase_dt_ns']):.8g}",
                    "cutoff_nm": f"{scalar(item['cutoff_nm']):.8g}",
                    "threshold": f"{scalar(item['threshold']):.8g}",
                }
            )
    return path


def write_curves(outdir: Path, results: list[dict[str, np.ndarray]]) -> Path:
    path = outdir / "interchain_partner_survival_curves.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["system", "temp", "role", "time_ns", "survival_mean", "survival_sem"])
        for item in results:
            system = text_value(item["system"])
            temp = text_value(item["temp"])
            role = text_value(item["role"])
            for time, mean, sem in zip(item["time_ns"], item["survival_mean"], item["survival_sem"]):
                writer.writerow([system, temp, role, f"{float(time):.6f}", f"{float(mean):.8f}", f"{float(sem):.8f}"])
    return path


def write_ratio(outdir: Path, results: list[dict[str, np.ndarray]]) -> Path:
    by_key = {(text_value(item["system"]), text_value(item["role"])): item for item in results}
    path = outdir / "interchain_partner_survival_ratio.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "system",
            "tl_temp",
            "th_temp",
            "t50_tl_ns",
            "t50_th_ns",
            "ratio_th_over_tl",
            "ratio_ci_low",
            "ratio_ci_high",
            "auc_ratio_th_over_tl",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for system in SYSTEM_ORDER:
            low = by_key.get((system, "Tl"))
            high = by_key.get((system, "Th"))
            if low is None or high is None:
                continue
            tau_l = scalar(low["t50_ns"])
            tau_h = scalar(high["t50_ns"])
            low_l = scalar(low["t50_ci_low_ns"])
            high_l = scalar(low["t50_ci_high_ns"])
            low_h = scalar(high["t50_ci_low_ns"])
            high_h = scalar(high["t50_ci_high_ns"])
            ratio = tau_h / tau_l if tau_l > 0 and np.isfinite(tau_l) and np.isfinite(tau_h) else float("nan")
            ci_low = low_h / high_l if high_l > 0 else float("nan")
            ci_high = high_h / low_l if low_l > 0 else float("nan")
            auc_ratio = scalar(high["auc_0_max_ns"]) / scalar(low["auc_0_max_ns"])
            writer.writerow(
                {
                    "system": system,
                    "tl_temp": text_value(low["temp"]),
                    "th_temp": text_value(high["temp"]),
                    "t50_tl_ns": f"{tau_l:.8g}",
                    "t50_th_ns": f"{tau_h:.8g}",
                    "ratio_th_over_tl": f"{ratio:.8g}",
                    "ratio_ci_low": f"{ci_low:.8g}",
                    "ratio_ci_high": f"{ci_high:.8g}",
                    "auc_ratio_th_over_tl": f"{auc_ratio:.8g}",
                }
            )
    return path


def plot_results(outdir: Path, results: list[dict[str, np.ndarray]]) -> None:
    by_key = {(text_value(item["system"]), text_value(item["role"])): item for item in results}
    ratio_rows = []
    ratio_path = outdir / "interchain_partner_survival_ratio.csv"
    if ratio_path.exists():
        with ratio_path.open("r", newline="", encoding="utf-8") as handle:
            ratio_rows = list(csv.DictReader(handle))

    fig = plt.figure(figsize=(13.8, 8.8), dpi=300)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.92], hspace=0.42, wspace=0.34)

    for idx, system in enumerate(SYSTEM_ORDER):
        ax = fig.add_subplot(gs[0, idx])
        for role in ("Tl", "Th"):
            item = by_key.get((system, role))
            if item is None:
                continue
            t = np.asarray(item["time_ns"], dtype=float)
            y = np.asarray(item["survival_mean"], dtype=float)
            e = np.asarray(item["survival_sem"], dtype=float)
            color = TEMP_COLORS[role]
            ax.plot(t, y, color=color, lw=2.0, label=f"{ROLE_LABEL[role]} ({text_value(item['temp'])})")
            ax.fill_between(t, y - e, y + e, color=color, alpha=0.14, linewidth=0)
            t50 = scalar(item["t50_ns"])
            if np.isfinite(t50):
                ax.axvline(t50, color=color, lw=0.9, ls="--", alpha=0.7)
        threshold = scalar(next(iter(by_key.values()))["threshold"])
        ax.axhline(threshold, color="#777777", lw=0.9, ls=":")
        ax.set_title(system)
        ax.set_xlabel("Lag time (ns)")
        if idx == 0:
            ax.set_ylabel("Partner identity survival")
        ax.set_ylim(-0.03, 1.03)
        ax.grid(alpha=0.22, linewidth=0.6)
        ax.legend(frameon=False, fontsize=8)

    ax_ratio = fig.add_subplot(gs[1, :2])
    x = np.arange(len(SYSTEM_ORDER), dtype=float)
    values = []
    lows = []
    highs = []
    for system in SYSTEM_ORDER:
        row = next((r for r in ratio_rows if r["system"] == system), None)
        if row is None:
            values.append(np.nan)
            lows.append(np.nan)
            highs.append(np.nan)
            continue
        value = float(row["ratio_th_over_tl"])
        values.append(value)
        lows.append(float(row["ratio_ci_low"]))
        highs.append(float(row["ratio_ci_high"]))
    values_arr = np.asarray(values, dtype=float)
    lows_arr = np.asarray(lows, dtype=float)
    highs_arr = np.asarray(highs, dtype=float)
    ax_ratio.bar(
        x,
        values_arr,
        width=0.55,
        color="#c96d28",
        edgecolor="#333333",
        yerr=np.vstack((values_arr - lows_arr, highs_arr - values_arr)),
        error_kw={"lw": 0.8, "capsize": 3, "ecolor": "#444444"},
    )
    ax_ratio.axhline(1, color="#444444", lw=1.0, ls="--")
    ax_ratio.set_xticks(x, SYSTEM_ORDER)
    ax_ratio.set_ylabel(r"Partner-survival ratio $t_{50}(T_H)/t_{50}(T_L)$")
    ax_ratio.set_title("Molecular partner persistence under heating")
    ax_ratio.grid(axis="y", alpha=0.22, linewidth=0.6)

    ax_note = fig.add_subplot(gs[1, 2])
    ax_note.axis("off")
    ax_note.text(
        0.0,
        0.88,
        "Definition",
        fontsize=12,
        fontweight="bold",
        transform=ax_note.transAxes,
    )
    ax_note.text(
        0.0,
        0.69,
        r"$S_{\mathrm{partner}}(\tau)=\langle |P_i(t)\cap P_i(t+\tau)|/|P_i(t)|\rangle$",
        fontsize=10,
        transform=ax_note.transAxes,
    )
    ax_note.text(
        0.0,
        0.43,
        "P_i(t): chains contacting the selected bulk chain i.\n"
        "t50: first crossing below 0.5.\n"
        "Same bulk-chain windows and cutoff as Figure 5.",
        fontsize=9,
        linespacing=1.45,
        transform=ax_note.transAxes,
    )

    fig.suptitle(
        "Inter-chain molecular partner survival in condensate bulk",
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.08, top=0.91)
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"interchain_partner_survival.{ext}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    conditions = select_conditions(args)
    if not conditions:
        raise RuntimeError("No selected conditions")
    results = [compute_condition(cond, args) for cond in conditions]
    if args.compute_only:
        print("[done] compute-only caches", flush=True)
        return
    write_summary(args.outdir, results)
    write_curves(args.outdir, results)
    write_ratio(args.outdir, results)
    plot_results(args.outdir, results)
    print(f"[save] {args.outdir / 'interchain_partner_survival.png'}", flush=True)


if __name__ == "__main__":
    main()
