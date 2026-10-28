#!/usr/bin/env python3
"""Shared post-processing tools for TDP43 CTD condensate trajectories."""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import time
from pathlib import Path

import numpy as np

# Keep post-processing independent from OpenMM/OpenABC.  The simulation helper
# module imports OpenMM at module import time, but extraction/density analysis
# only needs these constants.
CHAIN_LEN = 148
N_MOL_DEFAULT = 200
RES_MASS = {
    "ALA": 71.0788,
    "ARG": 156.1875,
    "ASN": 114.1038,
    "ASP": 115.0886,
    "CYS": 103.1388,
    "GLN": 128.1307,
    "GLU": 129.1155,
    "GLY": 57.0519,
    "HIS": 137.1411,
    "ILE": 113.1594,
    "LEU": 113.1594,
    "LYS": 128.1741,
    "MET": 131.1926,
    "PHE": 147.1766,
    "PRO": 97.1167,
    "SER": 87.0782,
    "THR": 101.1051,
    "TRP": 186.2132,
    "TYR": 163.1760,
    "VAL": 99.1326,
}

NA = 6.022e23
TWOPI = 2.0 * np.pi


def natural_key(path: str):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", path)]


def unwrap_chain_inplace(pos, box):
    for i in range(1, pos.shape[0]):
        delta = pos[i] - pos[i - 1]
        delta -= np.round(delta / box) * box
        pos[i] = pos[i - 1] + delta


def wrap_by_com_inplace(pos, box):
    com = pos.mean(axis=0)
    shift = -np.floor(com / box) * box
    pos += shift
    pos[:] = pos - np.floor(pos / box) * box


def main_deal_dcd(argv: list[str] | None = None) -> None:
    import MDAnalysis as mda

    parser = argparse.ArgumentParser(description="Fix periodic imaging in a TDP43 slab trajectory.")
    parser.add_argument("-s", "--top", default="tdp43_slab_vampnet_start.pdb")
    parser.add_argument("-f", "--traj", default="tdp43_slab_phase_diagram.dcd")
    parser.add_argument("-o", "--out", default="fix_stride.dcd")
    parser.add_argument("--out-pdb", default=None)
    parser.add_argument("--n-mol", type=int, default=N_MOL_DEFAULT)
    parser.add_argument("--chain-len", type=int, default=CHAIN_LEN)
    parser.add_argument("--mode", choices=("nojump", "center_wrap"), default="center_wrap")
    parser.add_argument("--center-sel", default="all")
    parser.add_argument("--force-box", nargs=3, type=float, default=None, help="Force box lengths in nm: Lx Ly Lz")
    parser.add_argument("--stride", type=int, default=20)
    parser.add_argument("--log-every", type=int, default=50)
    args = parser.parse_args(argv)

    if not os.path.exists(args.top):
        raise FileNotFoundError(args.top)
    if not os.path.exists(args.traj):
        raise FileNotFoundError(args.traj)

    u = mda.Universe(args.top, args.traj)
    expected_atoms = args.n_mol * args.chain_len
    if u.atoms.n_atoms != expected_atoms:
        raise RuntimeError(
            f"Atom count mismatch: got {u.atoms.n_atoms}, expected {expected_atoms} (= {args.n_mol}*{args.chain_len})."
        )

    center_ag = u.select_atoms(args.center_sel)
    written = 0
    with mda.Writer(args.out, n_atoms=u.atoms.n_atoms) as writer:
        for ts in u.trajectory[:: args.stride]:
            if args.force_box is not None:
                box_nm = np.array(args.force_box, dtype=float)
                ts.dimensions = np.array([box_nm[0] * 10.0, box_nm[1] * 10.0, box_nm[2] * 10.0, 90.0, 90.0, 90.0])
            else:
                if ts.dimensions is None or ts.dimensions.shape[0] < 3 or ts.dimensions[0] == 0:
                    raise RuntimeError("No unit cell in trajectory. Use --force-box Lx Ly Lz in nm.")
                box_nm = ts.dimensions[:3].astype(float) * 0.1

            pos_nm = u.atoms.positions * 0.1
            for mol_idx in range(args.n_mol):
                a0 = mol_idx * args.chain_len
                a1 = (mol_idx + 1) * args.chain_len
                unwrap_chain_inplace(pos_nm[a0:a1], box_nm)

            if args.mode == "center_wrap":
                sel_center = (center_ag.positions * 0.1).mean(axis=0)
                pos_nm[:] = pos_nm + (0.5 * box_nm - sel_center)
                for mol_idx in range(args.n_mol):
                    a0 = mol_idx * args.chain_len
                    a1 = (mol_idx + 1) * args.chain_len
                    wrap_by_com_inplace(pos_nm[a0:a1], box_nm)

            u.atoms.positions = pos_nm * 10.0
            writer.write(u.atoms)
            written += 1
            if written % args.log_every == 0:
                print(f"[progress] written={written} last_frame={ts.frame}")

    if args.out_pdb is not None:
        u.atoms.write(args.out_pdb)
        print(f"[save] {args.out_pdb}")
    print(f"[save] {args.out}")


def main_rebuild_phase_dcd(argv: list[str] | None = None) -> None:
    import MDAnalysis as mda

    parser = argparse.ArgumentParser(description="Rebuild the low-resolution phase DCD from chunked kinetic DCDs.")
    parser.add_argument("--top", default="tdp43_slab_vampnet_start.pdb")
    parser.add_argument("--part-pattern", default="tdp43_slab_vampnet_part*.dcd")
    parser.add_argument("--out", default="tdp43_slab_phase_diagram_rebuilt.dcd")
    parser.add_argument("--stride", type=int, default=10, help="Take one frame every N high-resolution frames.")
    args = parser.parse_args(argv)

    if not os.path.exists(args.top):
        raise FileNotFoundError(args.top)
    parts = sorted(glob.glob(args.part_pattern), key=natural_key)
    if not parts:
        raise FileNotFoundError(args.part_pattern)

    u = mda.Universe(args.top, parts)
    written = 0
    with mda.Writer(args.out, u.atoms.n_atoms) as writer:
        try:
            for ts in u.trajectory[:: args.stride]:
                writer.write(u.atoms)
                written += 1
                if written % 1000 == 0:
                    print(f"[progress] written={written}")
        except OSError as exc:
            print(f"[warn] stop early because of trajectory read error: {exc!r}")
    print(f"[save] {args.out} frames={written}")


def main_check_bond(argv: list[str] | None = None) -> None:
    import mdtraj as md

    parser = argparse.ArgumentParser(description="Check nearest-neighbor CA bond distances in a TDP43 PDB.")
    parser.add_argument("--pdb", default="compressed_relaxed_10ns.pdb")
    parser.add_argument("--warn-raw-nm", type=float, default=0.5)
    args = parser.parse_args(argv)

    traj = md.load(args.pdb)
    top = traj.topology
    bonds = []
    for chain in top.chains:
        ca = [atom.index for atom in chain.atoms if atom.name == "CA"]
        for i in range(len(ca) - 1):
            bonds.append((ca[i], ca[i + 1]))

    d_pbc = md.compute_distances(traj, bonds, periodic=True)
    d_raw = md.compute_distances(traj, bonds, periodic=False)

    def report(name, dist_nm):
        print(
            f"[{name}] mean={dist_nm.mean() * 10:.3f} A std={dist_nm.std() * 10:.3f} A max={dist_nm.max() * 10:.3f} A"
        )

    print(f"[info] unitcell_lengths_nm={traj.unitcell_lengths}")
    print(f"[info] pairs={len(bonds)} chains={top.n_chains} atoms={top.n_atoms}")
    report("PBC", d_pbc)
    report("RAW", d_raw)
    if d_raw.max() > args.warn_raw_nm:
        print("[warn] RAW(no-PBC) max distance exceeds threshold. Coordinates may have been wrapped incorrectly.")
    else:
        print("[ok] RAW(no-PBC) distances look continuous.")


class DSU:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]

    def largest_component(self) -> list[int]:
        groups: dict[int, list[int]] = {}
        for i in range(len(self.parent)):
            root = self.find(i)
            groups.setdefault(root, []).append(i)
        return max(groups.values(), key=len) if groups else []


def get_com_nm(pos_angstrom: np.ndarray, n_chains: int) -> np.ndarray:
    n_atoms = pos_angstrom.shape[0]
    beads_per_chain = n_atoms // n_chains
    pos = pos_angstrom.reshape(n_chains, beads_per_chain, 3).astype(np.float32)
    box_angstrom = np.array([200.0, 200.0, 2800.0], dtype=np.float32)
    ref = pos[:, :1, :]
    delta = pos - ref
    delta -= np.round(delta / box_angstrom) * box_angstrom
    com_angstrom = (ref + delta).mean(axis=1)
    return (com_angstrom / 10.0) % np.array([20.0, 20.0, 280.0], dtype=np.float32)


def solve_pbc_z_range(z_coords_nm: np.ndarray, Lz_nm: float, low_pct: float, high_pct: float):
    z_sorted = np.sort(z_coords_nm)
    gaps = np.diff(z_sorted, append=z_sorted[0] + Lz_nm)
    max_gap_idx = np.argmax(gaps)
    split_point = z_sorted[(max_gap_idx + 1) % len(z_sorted)]
    z_shifted = (z_coords_nm - split_point) % Lz_nm
    return split_point, np.percentile(z_shifted, low_pct), np.percentile(z_shifted, high_pct)


def pick_diverse_segments(
    in_cluster,
    n_chains: int,
    n_target: int,
    time_bins: int,
    window_coarse: int,
    min_bulk_fraction: float,
):
    n_start = in_cluster.shape[0] - window_coarse + 1
    if n_start <= 0:
        return []
    threshold = int(np.ceil(float(min_bulk_fraction) * float(window_coarse)))
    threshold = min(max(threshold, 1), window_coarse)
    candidates = [[] for _ in range(time_bins)]
    bin_size = max(n_start // time_bins, 1)
    kernel = np.ones(window_coarse, dtype=np.int32)
    for chain_idx in range(n_chains):
        series = in_cluster[:, chain_idx].astype(np.int32)
        counts = np.convolve(series, kernel, mode="valid")
        starts = np.where(counts >= threshold)[0]
        for start in starts:
            bin_idx = min(int(start // bin_size), time_bins - 1)
            candidates[bin_idx].append((chain_idx, int(start)))

    chosen = []
    used = set()
    for bin_candidates in candidates:
        for chain_idx, start in bin_candidates:
            if chain_idx not in used and len(chosen) < n_target:
                chosen.append((chain_idx, start))
                used.add(chain_idx)
                break
    return chosen


def main_extract_bulk_segments(argv: list[str] | None = None) -> None:
    import MDAnalysis as mda
    from MDAnalysis.lib.nsgrid import FastNS

    parser = argparse.ArgumentParser(description="Extract 0.5 us chain trajectories from the slab dense phase.")
    parser.add_argument("--pdb", default="tdp43_slab_vampnet_start.pdb")
    parser.add_argument("--dcd-pattern", default="tdp43_slab_vampnet_part*.dcd")
    parser.add_argument("--n-chains", type=int, default=N_MOL_DEFAULT)
    parser.add_argument("--dt-ns", type=float, default=0.01)
    parser.add_argument("--stride", type=int, default=20)
    parser.add_argument("--skip-us", type=float, default=5.0)
    parser.add_argument("--window-us", type=float, default=0.5)
    parser.add_argument("--n-target", type=int, default=40)
    parser.add_argument(
        "--min-target-chains",
        type=int,
        default=20,
        help="Abort extraction if fewer than this many unique dense-phase chains are found.",
    )
    parser.add_argument(
        "--min-bulk-fraction",
        type=float,
        default=1.0,
        help="Minimum fraction of the extraction window that a chain must stay in dense bulk.",
    )
    parser.add_argument("--time-bins", type=int, default=40)
    parser.add_argument("--cutoff-nm", type=float, default=4.0)
    parser.add_argument("--cluster-min-size", type=int, default=30)
    parser.add_argument("--low-percentile", type=float, default=20.0)
    parser.add_argument("--high-percentile", type=float, default=80.0)
    parser.add_argument("--out-dir", default="extracted_500ns_bulk_chains")
    parser.add_argument(
        "--scan-log-every",
        type=int,
        default=200,
        help="Print dense-phase scan progress every N sampled frames.",
    )
    parser.add_argument(
        "--write-log-every",
        type=int,
        default=5000,
        help="Print DCD writing progress every N full-resolution frames.",
    )
    args = parser.parse_args(argv)

    parts = sorted(glob.glob(args.dcd_pattern), key=natural_key)
    if not parts:
        raise FileNotFoundError(args.dcd_pattern)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    u = mda.Universe(args.pdb, parts)

    skip_frames = int(round(args.skip_us * 1000.0 / args.dt_ns))
    dt_coarse_ns = args.dt_ns * args.stride
    window_coarse = int(round(args.window_us * 1000.0 / dt_coarse_ns))
    window_full = int(round(args.window_us * 1000.0 / args.dt_ns))
    min_bulk_fraction = min(max(float(args.min_bulk_fraction), 0.0), 1.0)
    min_bulk_count = int(np.ceil(min_bulk_fraction * float(window_coarse)))
    frame_indices = list(range(skip_frames, len(u.trajectory), args.stride))
    pseudo_box = np.array([20.0, 20.0, 280.0, 90.0, 90.0, 90.0], dtype=np.float32)
    in_cluster = np.zeros((len(frame_indices), args.n_chains), dtype=np.bool_)

    print(
        "[extract] "
        f"parts={len(parts)} total_frames={len(u.trajectory)} skip_frames={skip_frames} "
        f"stride={args.stride} sampled_frames={len(frame_indices)} "
        f"window_full_frames={window_full} min_bulk_fraction={min_bulk_fraction:.3f} "
        f"min_bulk_count={min_bulk_count}/{window_coarse}",
        flush=True,
    )
    scan_t0 = time.monotonic()
    for ti, frame_idx in enumerate(frame_indices):
        u.trajectory[frame_idx]
        com_nm = get_com_nm(u.atoms.positions, args.n_chains)
        ns = FastNS(args.cutoff_nm, com_nm, pseudo_box, True)
        pairs = ns.self_search().get_pairs()
        dsu = DSU(args.n_chains)
        for a, b in pairs:
            dsu.union(int(a), int(b))
        largest = dsu.largest_component()
        if len(largest) < args.cluster_min_size:
            continue

        z_largest = com_nm[largest, 2]
        split_point, z_low, z_high = solve_pbc_z_range(
            z_largest, 280.0, args.low_percentile, args.high_percentile
        )
        z_shifted_all = (com_nm[:, 2] - split_point) % 280.0
        in_bulk = (z_shifted_all >= z_low) & (z_shifted_all <= z_high)
        largest_set = set(largest)
        for chain_idx in range(args.n_chains):
            if (chain_idx in largest_set) and in_bulk[chain_idx]:
                in_cluster[ti, chain_idx] = True

        if args.scan_log_every > 0 and (
            ti == 0 or (ti + 1) % args.scan_log_every == 0 or (ti + 1) == len(frame_indices)
        ):
            elapsed = max(time.monotonic() - scan_t0, 1e-9)
            rate = (ti + 1) / elapsed
            pct = 100.0 * (ti + 1) / max(len(frame_indices), 1)
            print(
                "[scan] "
                f"{ti + 1}/{len(frame_indices)} ({pct:.1f}%) "
                f"frame={frame_idx} rate={rate:.2f} sampled_frames/s "
                f"elapsed={elapsed / 60.0:.1f} min",
                flush=True,
            )

    chosen = pick_diverse_segments(
        in_cluster,
        n_chains=args.n_chains,
        n_target=args.n_target,
        time_bins=args.time_bins,
        window_coarse=window_coarse,
        min_bulk_fraction=min_bulk_fraction,
    )
    if not chosen:
        raise RuntimeError("No dense-phase segments found. Inspect the slab trajectory first.")
    unique_chains = sorted({int(chain_idx) for chain_idx, _ in chosen})
    if len(unique_chains) < args.min_target_chains:
        raise RuntimeError(
            "Not enough dense-phase chains found: "
            f"selected_segments={len(chosen)}, unique_chains={len(unique_chains)}, "
            f"required_min={args.min_target_chains}. "
            "Inspect the slab trajectory or relax extraction criteria."
        )
    print(
        f"[extract] selected_segments={len(chosen)} unique_chains={len(unique_chains)}",
        flush=True,
    )

    metadata = []
    for seg_id, (chain_idx, start_coarse) in enumerate(chosen):
        start_frame = frame_indices[start_coarse]
        end_frame = start_frame + window_full - 1
        metadata.append({"id": seg_id, "chain_id": int(chain_idx), "start": int(start_frame), "end": int(end_frame)})

    with (out_dir / "segments_meta.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    atoms_per_chain = u.atoms.n_atoms // args.n_chains
    for item in metadata:
        seg_id = int(item["id"])
        chain_idx = item["chain_id"]
        atom_group = u.atoms[chain_idx * atoms_per_chain : (chain_idx + 1) * atoms_per_chain]
        out_dcd = out_dir / f"chain{chain_idx:03d}_bulk_500ns.dcd"
        n_write = int(item["end"] - item["start"] + 1)
        print(
            "[write] "
            f"segment={seg_id + 1}/{len(metadata)} chain={chain_idx:03d} "
            f"frames={item['start']}..{item['end']} ({n_write}) -> {out_dcd}",
            flush=True,
        )
        write_t0 = time.monotonic()
        with mda.coordinates.DCD.DCDWriter(str(out_dcd), atom_group.n_atoms) as writer:
            for local_i, frame_idx in enumerate(range(item["start"], item["end"] + 1), start=1):
                u.trajectory[frame_idx]
                writer.write(atom_group)
                if args.write_log_every > 0 and (
                    local_i == 1 or local_i % args.write_log_every == 0 or local_i == n_write
                ):
                    elapsed = max(time.monotonic() - write_t0, 1e-9)
                    rate = local_i / elapsed
                    pct = 100.0 * local_i / max(n_write, 1)
                    print(
                        "[write] "
                        f"segment={seg_id + 1}/{len(metadata)} "
                        f"{local_i}/{n_write} ({pct:.1f}%) "
                        f"rate={rate:.2f} frames/s elapsed={elapsed / 60.0:.1f} min",
                        flush=True,
                    )
        print(f"[save] {out_dcd}", flush=True)


def circular_center_z(z_in_0L, Lz, weights):
    theta = TWOPI * (z_in_0L / Lz)
    C = np.sum(weights * np.cos(theta))
    S = np.sum(weights * np.sin(theta))
    theta0 = np.arctan2(S, C)
    return (theta0 % TWOPI) * (Lz / TWOPI)


def make_edges(z_window_angstrom: float, n_bins: int):
    edges = np.linspace(-float(z_window_angstrom), float(z_window_angstrom), n_bins + 1)
    centers_angstrom = 0.5 * (edges[:-1] + edges[1:])
    return edges, centers_angstrom


def hist_to_density_factor(Lx_angstrom, Ly_angstrom, bin_width_angstrom):
    bin_vol_angstrom3 = float(Lx_angstrom) * float(Ly_angstrom) * float(bin_width_angstrom)
    return (1e27 / NA) / bin_vol_angstrom3


def compute_profile(u, protein, bead_mass, edges, start_frame, stride, frames_limit=None):
    Lx0, Ly0, _, *_ = u.dimensions
    bin_width_angstrom = edges[1] - edges[0]
    total_hist = np.zeros(len(edges) - 1, dtype=np.float64)
    n_used = 0

    for ts in u.trajectory[start_frame::stride]:
        if (frames_limit is not None) and (n_used >= frames_limit):
            break
        Lx, Ly, Lz, *_ = ts.dimensions
        z = protein.positions[:, 2].astype(np.float64) % Lz
        z0 = circular_center_z(z, Lz, bead_mass)
        zc = ((z - z0 + 0.5 * Lz) % Lz) - 0.5 * Lz
        hist, _ = np.histogram(zc, bins=edges, weights=bead_mass)
        total_hist += hist
        n_used += 1

    if n_used == 0:
        raise RuntimeError("No frames selected for density profile.")

    mean_hist = total_hist / float(n_used)
    factor_mgml = hist_to_density_factor(Lx0, Ly0, bin_width_angstrom)
    return 0.5 * (edges[:-1] + edges[1:]), mean_hist * factor_mgml, n_used


def compute_blocks(u, protein, bead_mass, edges, frames_per_block, stride, n_blocks):
    u.trajectory[0]
    Lx, Ly, _, *_ = u.dimensions
    bin_width_angstrom = edges[1] - edges[0]
    factor_mgml = hist_to_density_factor(Lx, Ly, bin_width_angstrom)
    blocks = []
    for block_idx in range(n_blocks):
        start = int(block_idx * frames_per_block * stride)
        if start >= len(u.trajectory):
            break
        total_hist = np.zeros(len(edges) - 1, dtype=np.float64)
        used = 0
        for k in range(frames_per_block):
            frame_idx = start + k * stride
            if frame_idx >= len(u.trajectory):
                break
            u.trajectory[frame_idx]
            _, _, Lz, *_ = u.dimensions
            z = protein.positions[:, 2].astype(np.float64) % Lz
            z0 = circular_center_z(z, Lz, bead_mass)
            zc = ((z - z0 + 0.5 * Lz) % Lz) - 0.5 * Lz
            hist, _ = np.histogram(zc, bins=edges, weights=bead_mass)
            total_hist += hist
            used += 1
        if used == 0:
            break
        blocks.append((total_hist / float(used)) * factor_mgml)
    return 0.5 * (edges[:-1] + edges[1:]), np.asarray(blocks, dtype=np.float64)


def plot_two_col(z_nm, y, unit_label, title, out_png, log_floor):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), dpi=300)
    axes[0].plot(z_nm, y, lw=1.8)
    axes[0].set_xlabel("Z (nm)")
    axes[0].set_ylabel(f"Density ({unit_label})")
    axes[0].set_title(title + " (linear)")
    axes[0].grid(True, alpha=0.3)

    yy = np.clip(y, log_floor, None)
    axes[1].plot(z_nm, yy, lw=1.8)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Z (nm)")
    axes[1].set_ylabel(f"Density ({unit_label})")
    axes[1].set_title(title + " (log)")
    axes[1].grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(out_png)
    plt.close(fig)


def format_us(value: float) -> str:
    return f"{value:g}"


def block_label(block_idx: int, block_us: float) -> str:
    start_us = block_idx * float(block_us)
    end_us = (block_idx + 1) * float(block_us)
    return f"{format_us(start_us)}-{format_us(end_us)} us"


def plot_blocks_two_col(z_nm, blocks_mgml, out_png, log_floor, block_us):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), dpi=300)
    for i in range(blocks_mgml.shape[0]):
        label = block_label(i, block_us)
        axes[0].plot(z_nm, blocks_mgml[i], lw=1.2, label=label)
    axes[0].set_xlabel("Z (nm)")
    axes[0].set_ylabel("Density (mg/mL)")
    axes[0].set_title(f"rho(z) per {format_us(block_us)} us block (linear)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8, ncol=2, frameon=False)

    for i in range(blocks_mgml.shape[0]):
        axes[1].plot(z_nm, np.clip(blocks_mgml[i], log_floor, None), lw=1.2)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Z (nm)")
    axes[1].set_ylabel("Density (mg/mL)")
    axes[1].set_title(f"rho(z) per {format_us(block_us)} us block (log)")
    axes[1].grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(out_png)
    plt.close(fig)


def main_plot_density_profile(argv: list[str] | None = None) -> None:
    import MDAnalysis as mda

    parser = argparse.ArgumentParser(description="Compute TDP43 slab density profiles.")
    parser.add_argument("--topology", default="tdp43_slab_vampnet_start.pdb")
    parser.add_argument("--trajectory", default="tdp43_slab_phase_diagram.dcd")
    parser.add_argument("--frame-dt-ns", type=float, default=0.1)
    parser.add_argument("--stride", type=int, default=50)
    parser.add_argument("--skip-ns-overall", type=float, default=3000.0)
    parser.add_argument("--block-us", type=float, default=1.0)
    parser.add_argument("--z-window-ang", type=float, default=400.0)
    parser.add_argument("--n-bins", type=int, default=500)
    parser.add_argument("--overall-dat", default="overall_rhoz_gml.dat")
    parser.add_argument("--block-dat", default="block_rhoz_mgml.dat")
    parser.add_argument("--overall-fig", default="overall_rhoz_gml_2col.png")
    parser.add_argument("--block-fig", default="block_rhoz_mgml_2col.png")
    args = parser.parse_args(argv)

    if not (os.path.exists(args.topology) and os.path.exists(args.trajectory)):
        raise FileNotFoundError("Missing topology or trajectory.")

    u = mda.Universe(args.topology, args.trajectory)
    u.trajectory[0]
    protein = u.select_atoms("protein")
    bead_mass = np.array([RES_MASS.get(atom.resname.strip().upper(), 100.0) for atom in protein.atoms], dtype=float)
    edges, centers_angstrom = make_edges(args.z_window_ang, args.n_bins)
    z_nm = centers_angstrom / 10.0

    start_frame_overall = int(round(args.skip_ns_overall / args.frame_dt_ns))
    _, density_mgml_overall, n_used = compute_profile(u, protein, bead_mass, edges, start_frame_overall, args.stride)
    density_gml_overall = density_mgml_overall / 1000.0
    np.savetxt(
        args.overall_dat,
        np.column_stack([z_nm, density_gml_overall]),
        header="Z(nm) rho(z) (g/mL)",
        fmt="%.6f",
    )
    plot_two_col(z_nm, density_gml_overall, "g/mL", "Overall rho(z)", args.overall_fig, log_floor=1e-6)
    print(f"[save] {args.overall_dat} frames_used={n_used}")
    print(f"[save] {args.overall_fig}")

    eff_dt_ns = args.frame_dt_ns * args.stride
    frames_per_block = max(1, int(round((args.block_us * 1000.0) / eff_dt_ns)))
    total_ns = len(u.trajectory) * args.frame_dt_ns
    n_blocks = max(1, int(total_ns // (args.block_us * 1000.0)))
    print(
        f"[info] trajectory_frames={len(u.trajectory)} total_us={total_ns / 1000.0:.4f} "
        f"stride={args.stride} eff_dt_ns={eff_dt_ns:.4f} block_us={args.block_us:g} "
        f"frames_per_block={frames_per_block} full_blocks={n_blocks}",
        flush=True,
    )
    _, blocks_mgml = compute_blocks(u, protein, bead_mass, edges, frames_per_block, args.stride, n_blocks)
    output_matrix = [z_nm]
    headers = ["Z(nm)"]
    for i in range(blocks_mgml.shape[0]):
        output_matrix.append(blocks_mgml[i])
        headers.append(f"block_{format_us(i * args.block_us)}_{format_us((i + 1) * args.block_us)}us_mgml")
    np.savetxt(args.block_dat, np.column_stack(output_matrix), header=" ".join(headers), fmt="%.6f")
    plot_blocks_two_col(z_nm, blocks_mgml, args.block_fig, log_floor=1e-6, block_us=args.block_us)
    print(f"[save] {args.block_dat} blocks={blocks_mgml.shape[0]}")
    print(f"[save] {args.block_fig}")


__all__ = [
    "main_check_bond",
    "main_deal_dcd",
    "main_extract_bulk_segments",
    "main_plot_density_profile",
    "main_rebuild_phase_dcd",
]
