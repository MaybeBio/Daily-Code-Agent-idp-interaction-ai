#!/usr/bin/env python3
"""Extract a real local condensate environment around one tagged FUS-LCD chain."""

from pathlib import Path
import sys
import numpy as np


def main(source: str, target: int, output: str, n_neighbors: int = 20) -> None:
    lines = Path(source).read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("ATOM  ")]
    xyz = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    for line in atom_lines], dtype=float).reshape(-1, 163, 3)
    box = np.array([200.0, 200.0, 2800.0])

    unwrapped = xyz.copy()
    for chain in range(len(unwrapped)):
        for residue in range(1, 163):
            delta = xyz[chain, residue] - xyz[chain, residue - 1]
            delta -= box * np.round(delta / box)
            unwrapped[chain, residue] = unwrapped[chain, residue - 1] + delta

    centers = unwrapped.mean(axis=1)
    delta = centers - centers[target]
    delta -= box * np.round(delta / box)
    distances = np.linalg.norm(delta, axis=1)
    neighbors = [i for i in np.argsort(distances) if i != target][:n_neighbors]
    selected = [target] + neighbors

    local = []
    target_center = centers[target]
    for chain in selected:
        shift = delta[chain] - (centers[chain] - target_center)
        local.append(unwrapped[chain] + shift - target_center)
    local = np.asarray(local)

    # Use the tagged-chain principal axes for a reproducible, readable view;
    # apply the same rotation to its actual neighboring chains.
    _, _, axes = np.linalg.svd(local[0], full_matrices=False)
    local = local @ axes.T

    output_lines = ["REMARK REAL LOCAL ENVIRONMENT FROM EQUILIBRATED 200-CHAIN SLAB SNAPSHOT"]
    serial = 1
    # Keep the tagged chain intact.  For surrounding chains, retain only the
    # chain segments close to the tagged chain.  This produces a dense local
    # environment without filling the panel with irrelevant remote tails.
    displayed = [(selected[0], local[0], np.arange(163))]
    for source_chain, coords in zip(selected[1:], local[1:]):
        nearest = np.linalg.norm(coords[:, None, :] - local[0][None, :, :], axis=2).min(axis=1)
        near = nearest < 100.0
        keep = near.copy()
        for offset in range(1, 6):
            keep[offset:] |= near[:-offset]
            keep[:-offset] |= near[offset:]
        residue_ids = np.where(keep)[0]
        if residue_ids.size:
            displayed.append((source_chain, coords[residue_ids], residue_ids))

    for displayed_chain, (source_chain, coords, residue_ids) in enumerate(displayed):
        chain_id = chr(65 + displayed_chain % 26)
        segid = "TARG" if displayed_chain == 0 else f"N{displayed_chain:03d}"
        first_serial = serial
        atom_serials = []
        for residue0, coord in zip(residue_ids, coords):
            residue = int(residue0) + 1
            template = atom_lines[(source_chain * 163) + residue - 1]
            resname = template[17:20]
            output_lines.append(
                f"ATOM  {serial:5d}  CA  {resname:>3s} {chain_id}{residue:4d}    "
                f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00      {segid:>4s} C"
            )
            atom_serials.append(serial)
            serial += 1
        for k in range(len(atom_serials) - 1):
            if residue_ids[k + 1] == residue_ids[k] + 1:
                output_lines.append(f"CONECT{atom_serials[k]:5d}{atom_serials[k + 1]:5d}")
        output_lines.append("TER")
    output_lines.append("END")
    Path(output).write_text("\n".join(output_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 20)
