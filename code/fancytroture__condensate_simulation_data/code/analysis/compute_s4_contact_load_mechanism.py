#!/usr/bin/env python3
"""Cache per-partner bead-contact loads for the Figure S4 mechanism analysis."""

from __future__ import annotations

import argparse
import concurrent.futures
import copy
import sys
from pathlib import Path

import mdtraj as md
import numpy as np
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
PPT_DIR = ROOT / "ppt_5_11_figures"
sys.path.insert(0, str(PPT_DIR))

from compute_interchain_contact_relaxation_parts import (  # noqa: E402
    CONDITIONS,
    build_schedule,
    choose_part_selection,
    load_segments,
)
from compute_interchain_partner_survival import (  # noqa: E402
    INTERCHAIN_SUMMARY_DEFAULT,
    interchain_window_for,
)


DEFAULT_OUTDIR = ROOT / "source_data/supplementary/s4_contact_load_cache"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
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
    parser.add_argument("--min-segment-ns", type=float, default=0.0)
    parser.add_argument("--chunk", type=int, default=20)
    parser.add_argument("--n-chains", type=int, default=200)
    parser.add_argument("--chain-id-base", choices=["auto", "0", "1"], default="auto")
    parser.add_argument("--max-segments", type=int, default=0)
    parser.add_argument("--interchain-summary", type=Path, default=INTERCHAIN_SUMMARY_DEFAULT)
    parser.add_argument("--no-match-interchain-window", action="store_true")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def selected_conditions(args: argparse.Namespace):
    systems = {value.strip() for value in args.systems.split(",") if value.strip()}
    roles = {value.strip() for value in args.roles.split(",") if value.strip()}
    return [cond for cond in CONDITIONS if cond.system in systems and cond.role in roles]


def cache_path(cond, selection, args: argparse.Namespace) -> Path:
    return args.outdir / (
        f"partner_contact_load_{cond.system}_{cond.temp}_"
        f"parts{selection.first_part:03d}-{selection.last_part:03d}_"
        f"pstride{int(args.phase_stride):03d}.npz"
    )


def partner_contact_loads(
    positions_nm: np.ndarray,
    box_nm: np.ndarray,
    chain_len: int,
    chain_ids: list[int],
    cutoff_nm: float,
    atom_chain_index: np.ndarray,
    n_chains: int,
) -> dict[int, np.ndarray]:
    """Count bead contacts from each focal chain to every molecular partner."""
    if box_nm is None or len(box_nm) != 3 or np.any(np.asarray(box_nm) <= 0):
        raise ValueError(f"Invalid unitcell_lengths: {box_nm}")
    pos = np.asarray(positions_nm, dtype=np.float64) % box_nm
    tree = cKDTree(pos, boxsize=box_nm)
    output: dict[int, np.ndarray] = {}
    for chain_id in chain_ids:
        start = int(chain_id) * chain_len
        stop = (int(chain_id) + 1) * chain_len
        neighbors = tree.query_ball_point(pos[start:stop], r=cutoff_nm)
        nonempty = [np.asarray(indices, dtype=np.int64) for indices in neighbors if indices]
        if not nonempty:
            output[int(chain_id)] = np.zeros(n_chains, dtype=np.uint16)
            continue
        neighbor_atoms = np.concatenate(nonempty)
        partner_ids = atom_chain_index[neighbor_atoms]
        partner_ids = partner_ids[partner_ids != int(chain_id)]
        counts = np.bincount(partner_ids, minlength=n_chains)
        if counts.max(initial=0) > np.iinfo(np.uint16).max:
            raise OverflowError("Per-partner contact load exceeds uint16 capacity")
        output[int(chain_id)] = counts.astype(np.uint16, copy=False)
    return output


def compute_condition(cond, args: argparse.Namespace) -> Path:
    local_args = copy.copy(args)
    window = interchain_window_for(cond, args)
    if window is not None:
        local_args.start_part, local_args.end_part = window
    selection = choose_part_selection(cond, local_args)
    output = cache_path(cond, selection, args)
    if output.exists() and not args.force:
        print(f"[cached] {output}", flush=True)
        return output

    for path in (cond.topology, cond.metadata):
        if not path.exists():
            raise FileNotFoundError(path)
    segments = load_segments(
        cond.metadata,
        int(args.max_segments),
        int(args.n_chains),
        args.chain_id_base,
    )
    schedule, segment_frames, _bounds, min_frame, max_frame = build_schedule(
        segments,
        full_to_phase_stride=int(args.full_to_phase_stride),
        phase_stride=int(args.phase_stride),
        selected_global_start=selection.global_start,
        selected_global_end=selection.global_end,
    )
    if not schedule:
        raise RuntimeError(f"No sampled frames for {cond.system} {cond.temp}")

    first_frame = md.load_frame(
        str(selection.selected_parts[0].path),
        0,
        top=str(cond.topology),
    )
    n_atoms = first_frame.n_atoms
    if n_atoms % int(args.n_chains) != 0:
        raise RuntimeError(f"{n_atoms} atoms are not divisible by {args.n_chains} chains")
    chain_len = n_atoms // int(args.n_chains)
    atom_chain_index = np.repeat(
        np.arange(int(args.n_chains), dtype=np.int32),
        chain_len,
    )

    max_length = max(len(frames) for frames in segment_frames)
    contact_load = np.zeros(
        (len(segments), max_length, int(args.n_chains)),
        dtype=np.uint16,
    )
    valid = np.zeros((len(segments), max_length), dtype=bool)
    frame_to_pos = [
        {frame: index for index, frame in enumerate(frames)}
        for frames in segment_frames
    ]
    dt_ns = float(args.phase_dt_ns) * int(args.phase_stride)

    print(
        f"[compute] {cond.system} {cond.temp} "
        f"parts={selection.first_part:03d}-{selection.last_part:03d} "
        f"segments={len(segments)} frames={len(schedule)} dt={dt_ns:.3f} ns",
        flush=True,
    )

    n_filled = 0
    for part in selection.selected_parts:
        part_global_start = (part.index - 1) * selection.frames_per_part
        part_global_end = part.index * selection.frames_per_part - 1
        if part_global_end < min_frame or part_global_start > max_frame:
            continue

        local_start = (
            int(args.phase_stride) - part_global_start % int(args.phase_stride)
        ) % int(args.phase_stride)
        local_counter = 0
        iterator = md.iterload(
            str(part.path),
            top=str(cond.topology),
            chunk=int(args.chunk),
            skip=int(local_start),
            stride=int(args.phase_stride),
        )
        for chunk in iterator:
            for local_index in range(chunk.n_frames):
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
                chain_ids = sorted(
                    {int(segments[segment_index]["chain_id"]) for segment_index in active}
                )
                loads = partner_contact_loads(
                    chunk.xyz[local_index],
                    np.asarray(chunk.unitcell_lengths[local_index], dtype=np.float64),
                    chain_len,
                    chain_ids,
                    float(args.cutoff_nm),
                    atom_chain_index,
                    int(args.n_chains),
                )
                for segment_index in active:
                    chain_id = int(segments[segment_index]["chain_id"])
                    position = frame_to_pos[segment_index][global_frame]
                    contact_load[segment_index, position] = loads[chain_id]
                    valid[segment_index, position] = True
                    n_filled += 1

    min_points = max(
        8,
        int(round(float(args.min_segment_ns) / dt_ns)),
    )
    usable = np.flatnonzero(np.sum(valid, axis=1) >= min_points)
    if usable.size == 0:
        raise RuntimeError(f"No usable contact-load segments for {cond.system} {cond.temp}")

    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        system=np.asarray(cond.system),
        temp=np.asarray(cond.temp),
        role=np.asarray(cond.role),
        contact_load=contact_load[usable],
        valid=valid[usable],
        original_segment_index=usable.astype(np.int16),
        chain_id=np.asarray(
            [int(segments[index]["chain_id"]) for index in usable],
            dtype=np.int16,
        ),
        dt_ns=np.asarray(dt_ns, dtype=np.float32),
        cutoff_nm=np.asarray(args.cutoff_nm, dtype=np.float32),
        chain_len=np.asarray(chain_len, dtype=np.int16),
        n_chains=np.asarray(args.n_chains, dtype=np.int16),
        first_part=np.asarray(selection.first_part, dtype=np.int16),
        last_part=np.asarray(selection.last_part, dtype=np.int16),
        n_filled=np.asarray(n_filled, dtype=np.int32),
    )
    means = []
    for index in range(usable.size):
        mask = valid[usable[index]]
        means.append(float(np.mean(np.sum(contact_load[usable[index], mask], axis=1))))
    print(
        f"[save] {output} usable_segments={usable.size} "
        f"mean_contact={np.mean(means):.6f}",
        flush=True,
    )
    return output


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    conditions = selected_conditions(args)
    if int(args.workers) <= 1:
        for condition in conditions:
            compute_condition(condition, args)
        return

    with concurrent.futures.ProcessPoolExecutor(max_workers=int(args.workers)) as executor:
        futures = {
            executor.submit(compute_condition, condition, args): condition
            for condition in conditions
        }
        for future in concurrent.futures.as_completed(futures):
            condition = futures[future]
            try:
                future.result()
            except Exception as error:
                raise RuntimeError(
                    f"Failed {condition.system} {condition.temp}"
                ) from error


if __name__ == "__main__":
    main()
