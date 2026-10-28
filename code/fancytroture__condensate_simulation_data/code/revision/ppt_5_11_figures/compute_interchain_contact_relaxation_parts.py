#!/usr/bin/env python3
"""Compute inter-chain contact-number relaxation from segmented slab DCD files.

This version does NOT use the merged *_slab_phase_diagram.dcd, because that file
may be corrupted or unreadable after some frames. Instead, it auto-detects
segmented files like:

    fus_slab_vampnet_part001.dcd ... fus_slab_vampnet_part150.dcd

For each condition, the script automatically selects the last part of the
trajectory. By default it uses roughly the last 25% of available DCD parts,
with a maximum of 30 parts and a minimum of 10 parts. This means:

    150 parts -> last 30 parts, about 3 us if each part is 100 ns
     80 parts -> last 20 parts, about 2 us
     40 parts -> last 10 parts, about 1 us

You can override this with --last-parts, --start-part, and --end-part.

The script still uses segments_meta.json to choose pre-selected bulk-chain
windows. The metadata start/end are mapped to phase-DCD frame indices using
--full-to-phase-stride, then intersected with the selected part range.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
OUTDIR_DEFAULT = ROOT / "ppt_5_11_figures" / "results"

TEMP_COLORS = {"Tl": "#2f7f73", "Th": "#b54b5f"}
ROLE_LABEL = {"Tl": r"$T_l$", "Th": r"$T_h$"}
SYSTEM_ORDER = ["FUS", "A1", "TDP43"]


@dataclass(frozen=True)
class Condition:
    system: str
    temp: str
    role: str
    directory: Path
    topology: Path
    metadata: Path
    part_glob: str


CONDITIONS = [
    Condition(
        "FUS", "280K", "Tl",
        ROOT / "fus/fus_condensate_simulation/condensate_280K_200chains",
        ROOT / "fus/fus_condensate_simulation/condensate_280K_200chains/fus_slab_vampnet_start.pdb",
        ROOT / "fus/fus_condensate_simulation/condensate_280K_200chains/extracted_500ns_bulk_chains/segments_meta.json",
        "fus_slab_vampnet_part*.dcd",
    ),
    Condition(
        "FUS", "325K", "Th",
        ROOT / "fus/fus_condensate_simulation/condensate_325K_200chains",
        ROOT / "fus/fus_condensate_simulation/condensate_325K_200chains/fus_slab_vampnet_start.pdb",
        ROOT / "fus/fus_condensate_simulation/condensate_325K_200chains/extracted_500ns_bulk_chains/segments_meta.json",
        "fus_slab_vampnet_part*.dcd",
    ),
    Condition(
        "A1", "230K", "Tl",
        ROOT / "A1/A1_condensate_simulation/condensate_230K_200chains",
        ROOT / "A1/A1_condensate_simulation/condensate_230K_200chains/a1_slab_vampnet_start.pdb",
        ROOT / "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_230K/extracted_500ns_bulk_chains/segments_meta.json",
        "a1_slab_vampnet_part*.dcd",
    ),
    Condition(
        "A1", "260K", "Th",
        ROOT / "A1/A1_condensate_simulation/condensate_260K_200chains",
        ROOT / "A1/A1_condensate_simulation/condensate_260K_200chains/a1_slab_vampnet_start.pdb",
        ROOT / "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_260K/extracted_500ns_bulk_chains/segments_meta.json",
        "a1_slab_vampnet_part*.dcd",
    ),
    Condition(
        "TDP43", "260K", "Tl",
        ROOT / "TDP43/tdp43_condensate_simulation/condensate_260K_200chains",
        ROOT / "TDP43/tdp43_condensate_simulation/condensate_260K_200chains/tdp43_slab_vampnet_start.pdb",
        ROOT / "TDP43/tdp43_condensate_simulation/condensate_260K_200chains/extracted_500ns_bulk_chains/segments_meta.json",
        "tdp43_slab_vampnet_part*.dcd",
    ),
    Condition(
        "TDP43", "300K", "Th",
        ROOT / "TDP43/tdp43_condensate_simulation/condensate_300K_200chains",
        ROOT / "TDP43/tdp43_condensate_simulation/condensate_300K_200chains/tdp43_slab_vampnet_start.pdb",
        ROOT / "TDP43/tdp43_condensate_simulation/condensate_300K_200chains/extracted_500ns_bulk_chains/segments_meta.json",
        "tdp43_slab_vampnet_part*.dcd",
    ),
]


@dataclass(frozen=True)
class PartFile:
    index: int
    path: Path


@dataclass(frozen=True)
class PartSelection:
    all_parts: list[PartFile]
    selected_parts: list[PartFile]
    first_part: int
    last_part: int
    frames_per_part: int
    global_start: int
    global_end: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=OUTDIR_DEFAULT)
    parser.add_argument("--systems", default="FUS,A1,TDP43")
    parser.add_argument("--roles", default="Tl,Th")

    # Time mapping. Default assumes the full simulation was written every 0.01 ns,
    # while phase/part DCDs are every 10 full frames = 0.1 ns.
    parser.add_argument("--phase-dt-ns", type=float, default=0.1)
    parser.add_argument("--phase-stride", type=int, default=5,
                        help="Analyze every N phase-DCD frames. Default 5 -> 0.5 ns spacing.")
    parser.add_argument("--full-to-phase-stride", type=int, default=10,
                        help="metadata start/end divided by this to get phase-DCD frames.")

    # Part-DCD mapping.
    parser.add_argument("--part-duration-ns", type=float, default=100.0,
                        help="Expected duration of one part DCD in ns. Default: 100 ns.")
    parser.add_argument("--frames-per-part", type=int, default=0,
                        help="Override frames per part. 0 means infer from part-duration-ns / phase-dt-ns.")
    parser.add_argument("--last-parts", type=int, default=0,
                        help="Use the last N part DCDs for every condition. 0 means auto.")
    parser.add_argument("--start-part", type=int, default=0,
                        help="Manual first part index to use. Requires --end-part or uses available end.")
    parser.add_argument("--end-part", type=int, default=0,
                        help="Manual last part index to use. Requires --start-part or uses available start.")
    parser.add_argument("--auto-fraction", type=float, default=0.25,
                        help="Auto policy uses this fraction of available parts, capped by --auto-max-parts.")
    parser.add_argument("--auto-min-parts", type=int, default=10)
    parser.add_argument("--auto-max-parts", type=int, default=30)

    # Contact and ACF settings.
    parser.add_argument("--cutoff-nm", type=float, default=0.66)
    parser.add_argument("--threshold", type=float, default=0.1,
                        help="Relaxation time is first ACF crossing below this value.")
    parser.add_argument("--max-lag-ns", type=float, default=30.0)
    parser.add_argument("--min-segment-ns", type=float, default=0.0, help="Discard segments shorter than this length after part-window overlap. 0 disables filtering.")
    parser.add_argument("--chunk", type=int, default=20)
    parser.add_argument("--n-chains", type=int, default=200)
    parser.add_argument("--chain-id-base", choices=["auto", "0", "1"], default="auto")
    parser.add_argument("--max-segments", type=int, default=0,
                        help="0 means use all selected bulk segments.")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def select_conditions(args: argparse.Namespace) -> list[Condition]:
    systems = {x.strip() for x in args.systems.split(",") if x.strip()}
    roles = {x.strip() for x in args.roles.split(",") if x.strip()}
    return [c for c in CONDITIONS if c.system in systems and c.role in roles]


def part_index_from_name(path: Path) -> int | None:
    match = re.search(r"part(\d+)\.dcd$", path.name)
    if not match:
        return None
    return int(match.group(1))


def discover_parts(cond: Condition) -> list[PartFile]:
    paths = sorted(cond.directory.glob(cond.part_glob))
    parts = []
    for path in paths:
        idx = part_index_from_name(path)
        if idx is not None:
            parts.append(PartFile(idx, path))
    parts.sort(key=lambda x: x.index)
    if not parts:
        raise FileNotFoundError("No part DCDs found: {} / {}".format(cond.directory, cond.part_glob))
    return parts


def choose_part_selection(cond: Condition, args: argparse.Namespace) -> PartSelection:
    all_parts = discover_parts(cond)
    first_available = all_parts[0].index
    last_available = all_parts[-1].index
    part_by_index = {p.index: p for p in all_parts}

    if args.frames_per_part > 0:
        frames_per_part = int(args.frames_per_part)
    else:
        frames_per_part = int(round(float(args.part_duration_ns) / float(args.phase_dt_ns)))
        if frames_per_part <= 0:
            raise ValueError("frames_per_part must be positive")

    if args.start_part or args.end_part:
        start = int(args.start_part) if args.start_part else first_available
        end = int(args.end_part) if args.end_part else last_available
    elif args.last_parts > 0:
        end = last_available
        start = max(first_available, end - int(args.last_parts) + 1)
    else:
        n_total = len(all_parts)
        n_use = int(math.ceil(float(n_total) * float(args.auto_fraction)))
        n_use = max(int(args.auto_min_parts), n_use)
        n_use = min(int(args.auto_max_parts), n_use)
        n_use = min(n_use, n_total)
        end = last_available
        start = max(first_available, end - n_use + 1)

    selected = [part_by_index[i] for i in range(start, end + 1) if i in part_by_index]
    if not selected:
        raise RuntimeError("No selected part files for {} {} with requested range {}..{}".format(
            cond.system, cond.temp, start, end
        ))

    # Global phase-frame convention: part001 covers 0..frames_per_part-1.
    global_start = (selected[0].index - 1) * frames_per_part
    global_end = selected[-1].index * frames_per_part - 1
    return PartSelection(
        all_parts=all_parts,
        selected_parts=selected,
        first_part=selected[0].index,
        last_part=selected[-1].index,
        frames_per_part=frames_per_part,
        global_start=global_start,
        global_end=global_end,
    )


def load_segments(path: Path, max_segments: int, n_chains: int, chain_id_base: str) -> list[dict]:
    with path.open("r", encoding="utf-8") as handle:
        segments = json.load(handle)
    if not segments:
        raise RuntimeError("Empty segments metadata: {}".format(path))

    # Normalize chain_id to 0-based indexing.
    ids = [int(seg["chain_id"]) for seg in segments]
    if chain_id_base == "1":
        offset = 1
    elif chain_id_base == "0":
        offset = 0
    else:
        # Auto: if IDs look like 1..n_chains, convert to 0-based.
        offset = 1 if min(ids) >= 1 and max(ids) <= n_chains and 0 not in ids else 0

    normalized = []
    for seg in segments:
        copied = dict(seg)
        copied["chain_id"] = int(copied["chain_id"]) - offset
        if copied["chain_id"] < 0 or copied["chain_id"] >= n_chains:
            raise ValueError("chain_id out of range after normalization: {} from {}".format(
                copied["chain_id"], path
            ))
        normalized.append(copied)

    if max_segments and len(normalized) > max_segments:
        idx = np.linspace(0, len(normalized) - 1, int(max_segments), dtype=int)
        normalized = [normalized[int(i)] for i in idx]
    return normalized


def segment_phase_bounds(segment: dict, full_to_phase_stride: int) -> tuple[int, int]:
    start = int(math.ceil(float(segment["start"]) / float(full_to_phase_stride)))
    end = int(math.floor(float(segment["end"]) / float(full_to_phase_stride)))
    return start, end


def build_schedule(
    segments: list[dict],
    full_to_phase_stride: int,
    phase_stride: int,
    selected_global_start: int,
    selected_global_end: int,
) -> tuple[dict[int, list[int]], list[list[int]], list[tuple[int, int]], int, int]:
    raw_bounds = [segment_phase_bounds(seg, full_to_phase_stride) for seg in segments]
    clipped_bounds = []
    schedule: dict[int, list[int]] = defaultdict(list)
    segment_frames: list[list[int]] = []

    for seg_idx, (raw_start, raw_end) in enumerate(raw_bounds):
        start = max(int(raw_start), int(selected_global_start))
        end = min(int(raw_end), int(selected_global_end))
        clipped_bounds.append((start, end))
        if end < start:
            segment_frames.append([])
            continue

        # Preserve a common phase-stride grid based on absolute global frame index.
        offset = (phase_stride - (start % phase_stride)) % phase_stride
        first = start + offset
        frames = list(range(first, end + 1, phase_stride))
        segment_frames.append(frames)
        for frame_idx in frames:
            schedule[frame_idx].append(seg_idx)

    if schedule:
        min_frame = min(schedule.keys())
        max_frame = max(schedule.keys())
    else:
        min_frame = selected_global_start
        max_frame = selected_global_end
    return schedule, segment_frames, clipped_bounds, min_frame, max_frame


def autocorr_fft(values: np.ndarray) -> np.ndarray | None:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size < 4:
        return None
    arr = arr - np.mean(arr)
    var = float(np.dot(arr, arr) / arr.size)
    if not np.isfinite(var) or var <= 1e-12:
        return None
    n = arr.size
    f = np.fft.fft(arr, n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[:n].real
    denom = var * np.arange(n, 0, -1, dtype=np.float64)
    return acf / denom


def crossing_time(time_ns: np.ndarray, acf: np.ndarray, threshold: float) -> float:
    idx = np.where(acf <= threshold)[0]
    if idx.size == 0:
        return float("nan")
    i = int(idx[0])
    if i == 0:
        return float(time_ns[0])
    x0, x1 = float(time_ns[i - 1]), float(time_ns[i])
    y0, y1 = float(acf[i - 1]), float(acf[i])
    if abs(y1 - y0) < 1e-12:
        return x1
    return x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)


def bootstrap_crossing_ci(
    acfs: np.ndarray,
    time_ns: np.ndarray,
    threshold: float,
    n_boot: int = 400,
    seed: int = 20240511,
) -> tuple[float, float]:
    if acfs.shape[0] < 2:
        t = crossing_time(time_ns, acfs[0], threshold)
        return t, t
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, acfs.shape[0], size=acfs.shape[0])
        vals.append(crossing_time(time_ns, np.nanmean(acfs[idx], axis=0), threshold))
    vals = np.asarray(vals, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def count_contacts_for_chains(
    positions_nm: np.ndarray,
    box_nm: np.ndarray,
    chain_len: int,
    chain_ids: list[int],
    cutoff_nm: float,
    atom_chain_index: np.ndarray,
) -> dict[int, int]:
    if box_nm is None or len(box_nm) != 3 or np.any(np.asarray(box_nm) <= 0):
        raise ValueError("Invalid unitcell_lengths: {}".format(box_nm))
    pos = np.asarray(positions_nm, dtype=np.float64) % box_nm
    tree = cKDTree(pos, boxsize=box_nm)
    out: dict[int, int] = {}
    for chain_id in chain_ids:
        a0 = int(chain_id) * chain_len
        a1 = (int(chain_id) + 1) * chain_len
        neighbors = tree.query_ball_point(pos[a0:a1], r=cutoff_nm)
        count = 0
        for nb in neighbors:
            if nb:
                nb_arr = np.asarray(nb, dtype=np.int64)
                count += int(np.count_nonzero(atom_chain_index[nb_arr] != chain_id))
        out[int(chain_id)] = count
    return out


def cache_name(cond: Condition, selection: PartSelection, args: argparse.Namespace) -> Path:
    return args.outdir / (
        "interchain_contact_relaxation_{}_{}_parts{:03d}-{:03d}.npz".format(
            cond.system, cond.temp, selection.first_part, selection.last_part
        )
    )


def compute_condition(cond: Condition, args: argparse.Namespace) -> dict:
    for path in (cond.topology, cond.metadata):
        if not path.exists():
            raise FileNotFoundError(path)

    selection = choose_part_selection(cond, args)
    for part in selection.selected_parts:
        if not part.path.exists():
            raise FileNotFoundError(part.path)

    cache_path = cache_name(cond, selection, args)
    if cache_path.exists() and not args.force:
        data = np.load(cache_path, allow_pickle=True)
        return {key: data[key] for key in data.files}

    segments = load_segments(cond.metadata, args.max_segments, args.n_chains, args.chain_id_base)
    schedule, segment_frames, clipped_bounds, min_frame, max_frame = build_schedule(
        segments,
        full_to_phase_stride=int(args.full_to_phase_stride),
        phase_stride=int(args.phase_stride),
        selected_global_start=selection.global_start,
        selected_global_end=selection.global_end,
    )
    if not schedule:
        raw_bounds = [segment_phase_bounds(seg, int(args.full_to_phase_stride)) for seg in segments]
        raise RuntimeError(
            "No sampled frames after intersecting metadata windows with selected parts.\n"
            "Condition: {} {}\n"
            "Selected parts: {}..{} global frames {}..{}\n"
            "Metadata phase-frame range: {}..{}".format(
                cond.system, cond.temp,
                selection.first_part, selection.last_part,
                selection.global_start, selection.global_end,
                min(s for s, _ in raw_bounds), max(e for _, e in raw_bounds),
            )
        )

    # Use the first selected part to determine topology size and chain length.
    first_frame = md.load_frame(str(selection.selected_parts[0].path), 0, top=str(cond.topology))
    if first_frame.n_frames != 1:
        raise RuntimeError("Could not read first frame from {}".format(selection.selected_parts[0].path))
    n_atoms = first_frame.n_atoms
    if n_atoms % int(args.n_chains) != 0:
        raise RuntimeError("Atom count {} is not divisible by n_chains={}".format(n_atoms, args.n_chains))
    chain_len = n_atoms // int(args.n_chains)
    atom_chain_index = np.repeat(np.arange(int(args.n_chains), dtype=np.int32), chain_len)

    series = [np.full(len(frames), np.nan, dtype=np.float64) for frames in segment_frames]
    frame_to_pos = [{frame: i for i, frame in enumerate(frames)} for frames in segment_frames]
    dt_sample_ns = float(args.phase_dt_ns) * float(args.phase_stride)

    print(
        "[compute] {} {} parts={:03d}-{:03d}/{} segments={} "
        "global_frames={}..{} sampled_frames={} dt={:.3f} ns chain_len={}".format(
            cond.system, cond.temp,
            selection.first_part, selection.last_part, len(selection.all_parts),
            len(segments), min_frame, max_frame, len(schedule), dt_sample_ns, chain_len,
        ),
        flush=True,
    )

    n_filled = 0
    n_failed_chunks = 0
    n_empty_chunks = 0

    for part in selection.selected_parts:
        part_global_start = (part.index - 1) * selection.frames_per_part
        part_global_end = part.index * selection.frames_per_part - 1
        if part_global_end < min_frame or part_global_start > max_frame:
            continue

        local_start = (int(args.phase_stride) - (part_global_start % int(args.phase_stride))) % int(args.phase_stride)
        local_counter = 0
        try:
            iterator = md.iterload(
                str(part.path),
                top=str(cond.topology),
                chunk=int(args.chunk),
                skip=int(local_start),
                stride=int(args.phase_stride),
            )
            for chunk in iterator:
                if chunk.n_frames == 0:
                    n_empty_chunks += 1
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
                    box = np.asarray(chunk.unitcell_lengths[local_idx], dtype=np.float64)
                    counts = count_contacts_for_chains(
                        chunk.xyz[local_idx],
                        box,
                        chain_len,
                        chains,
                        float(args.cutoff_nm),
                        atom_chain_index,
                    )
                    for seg_idx in active:
                        chain_id = int(segments[seg_idx]["chain_id"])
                        pos_idx = frame_to_pos[seg_idx][global_frame]
                        series[seg_idx][pos_idx] = counts[chain_id]
                        n_filled += 1
        except Exception as exc:
            n_failed_chunks += 1
            print("[warn] failed reading {}: {}".format(part.path, repr(exc)), flush=True)
            continue

    valid_counts = []
    mean_values = []
    std_values = []
    for values in series:
        finite = values[np.isfinite(values)]
        valid_counts.append(int(finite.size))
        if finite.size:
            mean_values.append(float(np.mean(finite)))
            std_values.append(float(np.std(finite)))

    print(
        "[debug] {} {} filled_values={} segments_with_data={}/{} "
        "valid_points_min={} valid_points_max={} mean_contact_mean={:.3f} mean_contact_std={:.3f}".format(
            cond.system, cond.temp, n_filled,
            sum(v >= 8 for v in valid_counts), len(valid_counts),
            min(valid_counts) if valid_counts else 0,
            max(valid_counts) if valid_counts else 0,
            np.mean(mean_values) if mean_values else np.nan,
            np.mean(std_values) if std_values else np.nan,
        ),
        flush=True,
    )

    max_lag_frames = int(round(float(args.max_lag_ns) / dt_sample_ns))
    acfs = []
    contact_means = []
    contact_means_per_bead = []
    segment_relax = []

    min_segment_points = int(round(float(args.min_segment_ns) / dt_sample_ns))

    for values in series:
        values = values[np.isfinite(values)]

        if values.size < max(8, min_segment_points):
            continue

        contact_means.append(float(np.mean(values)))
        acf = autocorr_fft(values)
        if acf is None:
            continue
        n_keep = min(acf.size, max_lag_frames + 1)
        acf = acf[:n_keep]
        acfs.append(acf)
        time = np.arange(n_keep, dtype=np.float64) * dt_sample_ns
        segment_relax.append(crossing_time(time, acf, float(args.threshold)))

    if not acfs:
        raise RuntimeError("No usable contact-number ACFs for {} {}".format(cond.system, cond.temp))

    n_common = min(acf.size for acf in acfs)
    acf_stack = np.stack([acf[:n_common] for acf in acfs], axis=0)
    time_ns = np.arange(n_common, dtype=np.float64) * dt_sample_ns
    acf_mean = np.nanmean(acf_stack, axis=0)
    acf_sem = np.nanstd(acf_stack, axis=0, ddof=1) / math.sqrt(max(acf_stack.shape[0], 1))
    relax_time_ns = crossing_time(time_ns, acf_mean, float(args.threshold))
    ci_low, ci_high = bootstrap_crossing_ci(acf_stack, time_ns, float(args.threshold))

    payload = {
        "system": np.asarray(cond.system),
        "temp": np.asarray(cond.temp),
        "role": np.asarray(cond.role),
        "time_ns": time_ns.astype(np.float32),
        "segment_acfs": acf_stack.astype(np.float32),
        "acf_mean": acf_mean.astype(np.float32),
        "acf_sem": acf_sem.astype(np.float32),
        "segment_relax_ns": np.asarray(segment_relax, dtype=np.float32),
        "segment_contact_mean": np.asarray(contact_means, dtype=np.float32),
        "segment_contact_mean_per_bead": np.asarray(contact_means_per_bead, dtype=np.float32),
        "relax_time_ns": np.asarray(relax_time_ns, dtype=np.float32),
        "relax_ci_low_ns": np.asarray(ci_low, dtype=np.float32),
        "relax_ci_high_ns": np.asarray(ci_high, dtype=np.float32),
        "mean_contact_number": np.asarray(float(np.mean(contact_means)), dtype=np.float32),
        "mean_contact_number_per_bead": np.asarray(float(np.mean(contact_means_per_bead)), dtype=np.float32),
        "n_segments": np.asarray(acf_stack.shape[0], dtype=np.int32),
        "n_filled": np.asarray(n_filled, dtype=np.int32),
        "phase_stride": np.asarray(args.phase_stride, dtype=np.int32),
        "phase_dt_ns": np.asarray(args.phase_dt_ns, dtype=np.float32),
        "cutoff_nm": np.asarray(args.cutoff_nm, dtype=np.float32),
        "threshold": np.asarray(args.threshold, dtype=np.float32),
        "chain_len": np.asarray(chain_len, dtype=np.int32),
        "first_part": np.asarray(selection.first_part, dtype=np.int32),
        "last_part": np.asarray(selection.last_part, dtype=np.int32),
        "n_available_parts": np.asarray(len(selection.all_parts), dtype=np.int32),
        "frames_per_part": np.asarray(selection.frames_per_part, dtype=np.int32),
        "selected_global_start": np.asarray(selection.global_start, dtype=np.int32),
        "selected_global_end": np.asarray(selection.global_end, dtype=np.int32),
        "actual_min_frame": np.asarray(min_frame, dtype=np.int32),
        "actual_max_frame": np.asarray(max_frame, dtype=np.int32),
        "n_failed_part_reads": np.asarray(n_failed_chunks, dtype=np.int32),
        "n_empty_chunks": np.asarray(n_empty_chunks, dtype=np.int32),
    }
    np.savez(cache_path, **payload)
    return payload


def scalar(value) -> float:
    arr = np.asarray(value)
    return float(arr.reshape(-1)[0])


def text_value(value) -> str:
    arr = np.asarray(value)
    return str(arr.reshape(-1)[0])


def write_summary(outdir: Path, rows: list[dict]) -> Path:
    path = outdir / "interchain_contact_relaxation_summary.csv"
    fieldnames = [
        "system", "temp", "role",
        "first_part", "last_part", "n_available_parts",
        "actual_min_frame", "actual_max_frame",
        "chain_len", "mean_contact_number", "mean_contact_number_per_bead",
        "relax_time_ns", "relax_ci_low_ns", "relax_ci_high_ns",
        "n_segments", "n_filled", "phase_stride", "phase_dt_ns", "cutoff_nm", "threshold",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return path


def write_curves(outdir: Path, results: list[dict]) -> Path:
    path = outdir / "interchain_contact_relaxation_curves.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["system", "temp", "role", "time_ns", "acf_mean", "acf_sem"])
        for item in results:
            system = text_value(item["system"])
            temp = text_value(item["temp"])
            role = text_value(item["role"])
            for t, y, e in zip(item["time_ns"], item["acf_mean"], item["acf_sem"]):
                writer.writerow([system, temp, role, "{:.6f}".format(float(t)), "{:.8f}".format(float(y)), "{:.8f}".format(float(e))])
    return path


def plot_results(outdir: Path, results: list[dict], threshold: float) -> None:
    systems = ["FUS", "A1", "TDP43"]

    for system in systems:
        system_results = [
            item for item in results
            if text_value(item["system"]) == system
        ]

        if not system_results:
            continue

        # 1. ACF curve figure for one system
        fig, ax = plt.subplots(figsize=(5.2, 4.0), dpi=300)

        for item in system_results:
            role = text_value(item["role"])
            temp = text_value(item["temp"])
            color = TEMP_COLORS[role]

            t = np.asarray(item["time_ns"], dtype=float)
            y = np.asarray(item["acf_mean"], dtype=float)
            e = np.asarray(item["acf_sem"], dtype=float)

            ax.plot(
                t,
                y,
                color=color,
                lw=2.0,
                label=f"{ROLE_LABEL[role]} ({temp})"
            )
            ax.fill_between(
                t,
                y - e,
                y + e,
                color=color,
                alpha=0.16,
                linewidth=0
            )

            rt = scalar(item["relax_time_ns"])
            if np.isfinite(rt):
                ax.axvline(
                    rt,
                    color=color,
                    lw=1.0,
                    ls="--",
                    alpha=0.75
                )

        ax.axhline(
            threshold,
            color="#777777",
            lw=0.9,
            ls=":",
            zorder=0
        )

        ax.set_title(f"{system}: inter-chain contact relaxation")
        ax.set_xlabel("Lag time (ns)")
        ax.set_ylabel("Contact-number autocorrelation")
        ax.grid(alpha=0.20, linestyle="--", linewidth=0.6)
        ax.legend(frameon=False, fontsize=8)

        fig.tight_layout()

        for ext in ("png", "pdf"):
            fig.savefig(
                outdir / f"interchain_contact_relaxation_curves_{system.lower()}.{ext}",
                bbox_inches="tight"
            )
        plt.close(fig)

        # 2. Relaxation time summary figure for one system
        fig, ax = plt.subplots(figsize=(4.4, 4.0), dpi=300)

        x = np.arange(len(system_results), dtype=float)
        labels = []
        vals = []
        lows = []
        highs = []
        colors = []

        # Sort Tl before Th
        system_results_sorted = sorted(
            system_results,
            key=lambda item: 0 if text_value(item["role"]) == "Tl" else 1
        )

        for item in system_results_sorted:
            role = text_value(item["role"])
            temp = text_value(item["temp"])

            val = scalar(item["relax_time_ns"])
            ci_low = scalar(item["relax_ci_low_ns"])
            ci_high = scalar(item["relax_ci_high_ns"])

            labels.append(f"{ROLE_LABEL[role]}\n{temp}")
            vals.append(val)
            lows.append(max(0.0, val - ci_low))
            highs.append(max(0.0, ci_high - val))
            colors.append(TEMP_COLORS[role])

        ax.bar(
            x,
            vals,
            color=colors,
            alpha=0.82,
            edgecolor="#303030",
            linewidth=0.8
        )

        ax.errorbar(
            x,
            vals,
            yerr=np.vstack([lows, highs]),
            fmt="none",
            ecolor="#303030",
            elinewidth=1.0,
            capsize=3
        )

        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Inter-chain contact relaxation time (ns)")
        ax.set_title(f"{system}: relaxation time")
        ax.grid(axis="y", alpha=0.22, linestyle="--", linewidth=0.6)

        fig.tight_layout()

        for ext in ("png", "pdf"):
            fig.savefig(
                outdir / f"interchain_contact_relaxation_summary_{system.lower()}.{ext}",
                bbox_inches="tight"
            )
        plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    conditions = select_conditions(args)
    results = []
    summary_rows = []

    for cond in conditions:
        payload = compute_condition(cond, args)
        results.append(payload)
        summary_rows.append(
            {
                "system": text_value(payload["system"]),
                "temp": text_value(payload["temp"]),
                "role": text_value(payload["role"]),
                "first_part": str(int(scalar(payload["first_part"]))),
                "last_part": str(int(scalar(payload["last_part"]))),
                "n_available_parts": str(int(scalar(payload["n_available_parts"]))),
                "actual_min_frame": str(int(scalar(payload["actual_min_frame"]))),
                "actual_max_frame": str(int(scalar(payload["actual_max_frame"]))),
                "chain_len": str(int(scalar(payload["chain_len"]))),
                "mean_contact_number": "{:.6f}".format(scalar(payload["mean_contact_number"])),
                "mean_contact_number_per_bead": "{:.6f}".format(scalar(payload["mean_contact_number_per_bead"])),
                "relax_time_ns": "{:.6f}".format(scalar(payload["relax_time_ns"])),
                "relax_ci_low_ns": "{:.6f}".format(scalar(payload["relax_ci_low_ns"])),
                "relax_ci_high_ns": "{:.6f}".format(scalar(payload["relax_ci_high_ns"])),
                "n_segments": str(int(scalar(payload["n_segments"]))),
                "n_filled": str(int(scalar(payload["n_filled"]))),
                "phase_stride": str(int(scalar(payload["phase_stride"]))),
                "phase_dt_ns": "{:.6f}".format(scalar(payload["phase_dt_ns"])),
                "cutoff_nm": "{:.6f}".format(scalar(payload["cutoff_nm"])),
                "threshold": "{:.6f}".format(scalar(payload["threshold"])),
            }
        )

    summary_csv = write_summary(args.outdir, summary_rows)
    curves_csv = write_curves(args.outdir, results)
    plot_results(args.outdir, results, float(args.threshold))
    print("[save] {}".format(summary_csv))
    print("[save] {}".format(curves_csv))
    print("[save] {}".format(args.outdir / "interchain_contact_relaxation_curves.png"))
    print("[save] {}".format(args.outdir / "interchain_contact_relaxation_summary.png"))


if __name__ == "__main__":
    main()
