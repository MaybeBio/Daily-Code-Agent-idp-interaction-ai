#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
OPENABC_ROOT = REPO_ROOT / "fus/fus_condensate_simulation/OpenABC"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(OPENABC_ROOT) not in sys.path:
    sys.path.insert(0, str(OPENABC_ROOT))

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as unit

from openabc.forcefields import HPSModel
from openabc.forcefields.parsers import HPSParser

from TDP43.common.tdp43_config import CTD_HELIX_RANGE, CTD_RANGE, CTD_SEQUENCE

# ================= 1) 配置区 =================
n_mol = 200
temp = 260 * unit.kelvin

dt = 5 * unit.femtosecond
friction = 1.0 / unit.picosecond

initial_L = 300.0
target_L = 20.0
compress_factor = 0.99
relax_steps_each_stage = 5000
minimize_max_iter = 5000

final_relax_time_ns = 10.0
final_relax_steps = int(final_relax_time_ns * 1000 / dt.value_in_unit(unit.picosecond))

reference_pdb = "tdp43_ctd_monomer.pdb"
monomer_pdb = "tdp43_ctd_monomer_collapsed.pdb"
start_pdb = f"start_L{int(initial_L)}.pdb"
out_pdb = "compressed_relaxed_10ns.pdb"
log_csv = "compress_log.csv"

CUDA_DEVICE_INDEX = os.environ.get("CUDA_DEVICE", "3")
CUDA_PRECISION = os.environ.get("CUDA_PRECISION", "mixed")

urry_mu = 1.0
urry_delta = 0.08
dh_ldby_nm = 1.0

helix_k = 50000.0
helix_cutoff_nm = 1.2
helix_min_seq_sep = 2

chain_len = len(CTD_SEQUENCE)
ctd_start, _ = CTD_RANGE
helix_start, helix_end = CTD_HELIX_RANGE


# ================= 2) TDP43 helix elastic =================
def parse_ca_pdb(pdb_path):
    residues = []
    coords_angstrom = []
    with open(pdb_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            residues.append(int(line[22:26]))
            coords_angstrom.append(
                [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ]
            )
    if not residues:
        raise RuntimeError(f"No CA atoms found in {pdb_path}")
    return residues, np.asarray(coords_angstrom, dtype=float) * 0.1


def build_helix_elastic_pairs(reference_ca_pdb):
    residues, coords_nm = parse_ca_pdb(reference_ca_pdb)
    helix_indices = [i for i, resseq in enumerate(residues) if helix_start <= resseq <= helix_end]
    if len(helix_indices) != (helix_end - helix_start + 1):
        raise RuntimeError(
            f"Expected {helix_end - helix_start + 1} helix residues in {reference_ca_pdb}, "
            f"found {len(helix_indices)}"
        )

    pairs = []
    for pos_i, atom_i in enumerate(helix_indices):
        for atom_j in helix_indices[pos_i + 1 :]:
            seq_sep = abs(residues[atom_j] - residues[atom_i])
            if seq_sep < helix_min_seq_sep:
                continue
            r0_nm = float(np.linalg.norm(coords_nm[atom_j] - coords_nm[atom_i]))
            if r0_nm <= helix_cutoff_nm:
                pairs.append((atom_i, atom_j, r0_nm))
    return pairs


def add_helix_elastic_network(system, reference_ca_pdb, n_chain):
    pairs = build_helix_elastic_pairs(reference_ca_pdb)
    force = mm.HarmonicBondForce()
    force.setUsesPeriodicBoundaryConditions(True)
    force.setForceGroup(4)
    k = helix_k * unit.kilojoule_per_mole / unit.nanometer**2
    for mol_idx in range(n_chain):
        offset = mol_idx * chain_len
        for atom_i, atom_j, r0_nm in pairs:
            force.addBond(offset + atom_i, offset + atom_j, r0_nm * unit.nanometer, k)
    system.addForce(force)
    return pairs


def write_helix_pair_table(pairs, out_csv="helix_elastic_network.csv"):
    with open(out_csv, "w", encoding="utf-8") as handle:
        handle.write("local_a1,local_a2,resseq_a1,resseq_a2,r0_nm,k_kj_mol_nm2\n")
        for atom_i, atom_j, r0_nm in pairs:
            handle.write(
                f"{atom_i},{atom_j},{ctd_start + atom_i},{ctd_start + atom_j},{r0_nm:.8f},{helix_k:.3f}\n"
            )


# ================= 3) 初始结构生成 =================
def random_rotation_matrix(rng):
    try:
        from scipy.spatial.transform import Rotation

        return Rotation.random(random_state=rng).as_matrix()
    except Exception:
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        w, x, y, z = q
        return np.array(
            [
                [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
            ],
            dtype=float,
        )


def safe_place_molecules(monomer_ca_pdb, output_pdb, n_chain, box_size, seed=42, buffer_nm=1.5):
    print(f"[init] Reading monomer: {monomer_ca_pdb}", flush=True)
    pdb_mono = app.PDBFile(monomer_ca_pdb)
    pos_mono = np.asarray(pdb_mono.positions.value_in_unit(unit.nanometer), dtype=float)
    pos_mono -= pos_mono.mean(axis=0)

    max_r = np.max(np.linalg.norm(pos_mono, axis=1))
    spacing = 2.0 * max_r + buffer_nm
    print(f"[init] Monomer Radius: {max_r:.2f} nm | Spacing: {spacing:.2f} nm", flush=True)

    n_side = int(np.ceil(n_chain ** (1.0 / 3.0)))
    needed = (n_side - 1) * spacing + 2.0 * (max_r + buffer_nm)
    if needed > box_size:
        print(f"[warn] Suggested box > {needed:.1f} nm, current {box_size:.1f} nm", flush=True)

    grid = []
    start = max_r + buffer_nm
    for ix in range(n_side):
        for iy in range(n_side):
            for iz in range(n_side):
                grid.append([start + ix * spacing, start + iy * spacing, start + iz * spacing])

    rng = np.random.RandomState(seed)
    indices = np.arange(len(grid))
    rng.shuffle(indices)
    centers = [grid[i] for i in indices[:n_chain]]

    modeller = app.Modeller(pdb_mono.topology, pdb_mono.positions)
    for _ in range(n_chain - 1):
        modeller.add(pdb_mono.topology, pdb_mono.positions)

    all_positions = []
    for center in centers:
        rot_matrix = random_rotation_matrix(rng)
        new_pos = np.dot(pos_mono, rot_matrix.T) + center
        all_positions.extend([mm.Vec3(x, y, z) for x, y, z in new_pos])

    final_positions = all_positions * unit.nanometer
    top = modeller.topology
    top.setUnitCellDimensions(mm.Vec3(box_size, box_size, box_size))

    with open(output_pdb, "w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(top, final_positions, handle)
    print(f"[save] {output_pdb}", flush=True)


# ================= 4) COM 缩放 =================
def scale_positions_by_chain_com_and_wrap(pos_nm, topology, old_L, new_L):
    pos = np.asarray(pos_nm, dtype=float)
    old_center = np.array([old_L / 2.0] * 3, dtype=float)
    new_center = np.array([new_L / 2.0] * 3, dtype=float)
    scale = float(new_L) / float(old_L)
    L_new = np.array([new_L, new_L, new_L], dtype=float)

    new_pos = pos.copy()
    for chain in topology.chains():
        idx = np.array([atom.index for atom in chain.atoms()], dtype=int)
        coords = pos[idx].copy()

        com_old = coords.mean(axis=0)
        rel_vec = com_old - old_center
        new_com = new_center + rel_vec * scale
        scaled_coords = new_com + (coords - com_old)

        c = scaled_coords.copy()
        for i in range(1, len(c)):
            d = c[i] - c[i - 1]
            d -= np.round(d / L_new) * L_new
            c[i] = c[i - 1] + d

        chain_com = c.mean(axis=0)
        shift = -np.floor(chain_com / L_new) * L_new
        c = c + shift
        new_pos[idx] = c

    return new_pos


def make_whole_and_wrap(positions, topology, box_vectors):
    pos = positions.value_in_unit(unit.nanometer)
    pos = np.asarray(pos, dtype=float)
    box_v = [v.value_in_unit(unit.nanometer) for v in box_vectors]
    L = np.array([box_v[0][0], box_v[1][1], box_v[2][2]], dtype=float)

    new_pos = pos.copy()
    for chain in topology.chains():
        idx = np.array([atom.index for atom in chain.atoms()], dtype=int)
        c = new_pos[idx].copy()
        for i in range(1, len(c)):
            d = c[i] - c[i - 1]
            d -= np.round(d / L) * L
            c[i] = c[i - 1] + d
        center = c.mean(axis=0)
        shift = -np.floor(center / L) * L
        c = c + shift
        new_pos[idx] = c
    return new_pos * unit.nanometer


def build_L_schedule(L0, L1, factor=0.99):
    values = [float(L0)]
    current = float(L0)
    while current * factor > L1:
        current *= factor
        values.append(float(current))
    if abs(values[-1] - float(L1)) > 1e-12:
        values.append(float(L1))
    return values


# ================= 5) 系统搭建 =================
def build_tdp43_system(start_ca_pdb):
    parser = HPSParser(reference_pdb)
    model = HPSModel()
    for _ in range(n_mol):
        model.append_mol(parser)

    pdb_ref = app.PDBFile(start_ca_pdb)
    model.create_system(pdb_ref.getTopology(), box_a=initial_L, box_b=initial_L, box_c=initial_L)
    model.add_protein_bonds(force_group=1)
    model.add_contacts("Urry", mu=urry_mu, delta=urry_delta, force_group=2)
    model.add_dh_elec(ldby=dh_ldby_nm * unit.nanometer, force_group=3)
    helix_pairs = add_helix_elastic_network(model.system, reference_pdb, n_mol)
    return model, pdb_ref, helix_pairs


# ================= 6) 主流程 =================
def main():
    if not os.path.exists(reference_pdb):
        raise FileNotFoundError(f"Missing reference PDB: {reference_pdb}")
    if not os.path.exists(monomer_pdb):
        raise RuntimeError(
            f"Missing {monomer_pdb}.\n"
            "Suggested order:\n"
            "1. python make_monomer_from_sequence.py\n"
            "2. prepare tdp43_ctd_monomer_collapsed.pdb\n"
            "   or use python copy_monomer_to_collapsed.py as a temporary placeholder\n"
        )

    print(f"[config] temp={temp} | cuda_device={CUDA_DEVICE_INDEX}", flush=True)
    safe_place_molecules(monomer_pdb, start_pdb, n_mol, initial_L)

    print("[step] Building HPS-Urry + helix elastic system...", flush=True)
    model, pdb_ref, helix_pairs = build_tdp43_system(start_pdb)
    write_helix_pair_table(helix_pairs, "helix_elastic_network.csv")
    print(f"[save] helix_elastic_network.csv | pairs={len(helix_pairs)}", flush=True)

    integrator = mm.LangevinMiddleIntegrator(temp, friction, dt)

    try:
        platform = mm.Platform.getPlatformByName("CUDA")
        props = {"DeviceIndex": CUDA_DEVICE_INDEX, "Precision": CUDA_PRECISION}
    except Exception:
        platform = mm.Platform.getPlatformByName("CPU")
        props = {}
    print(f"[platform] {platform.getName()} | props={props}", flush=True)

    print("[step] Creating Simulation...", flush=True)
    simulation = app.Simulation(model.top, model.system, integrator, platform, props)
    simulation.context.setPositions(pdb_ref.getPositions())
    simulation.context.setVelocitiesToTemperature(temp)

    print("[run] Initial Minimize...", flush=True)
    simulation.minimizeEnergy(maxIterations=minimize_max_iter)

    print("[run] Pre-equilibration (10000 steps)...", flush=True)
    simulation.step(10000)
    simulation.minimizeEnergy(maxIterations=max(1000, minimize_max_iter // 5))

    L_schedule = build_L_schedule(initial_L, target_L, factor=compress_factor)
    print(f"[plan] Compression schedule: {len(L_schedule)} stages", flush=True)

    cur_L = float(initial_L)
    with open(log_csv, "w", encoding="utf-8") as handle:
        handle.write("stage,old_L,new_L,step,pot_E\n")

        for stage, new_L in enumerate(L_schedule[1:], start=1):
            state = simulation.context.getState(getPositions=True)
            pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            pos_nm2 = scale_positions_by_chain_com_and_wrap(pos_nm, simulation.topology, cur_L, new_L)

            a = mm.Vec3(new_L, 0, 0) * unit.nanometer
            b = mm.Vec3(0, new_L, 0) * unit.nanometer
            c = mm.Vec3(0, 0, new_L) * unit.nanometer
            simulation.context.setPeriodicBoxVectors(a, b, c)
            simulation.context.setPositions(pos_nm2 * unit.nanometer)
            simulation.context.setVelocitiesToTemperature(temp)

            simulation.minimizeEnergy(maxIterations=0)
            st = simulation.context.getState(getEnergy=True)
            e_check = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            if e_check > 1e6:
                integrator.setStepSize(0.1 * unit.femtosecond)
                simulation.step(1000)
                integrator.setStepSize(1.0 * unit.femtosecond)
                simulation.step(1000)
                integrator.setStepSize(dt)

            simulation.step(relax_steps_each_stage)

            state = simulation.context.getState(getEnergy=True)
            pot = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            if not np.isfinite(pot) or pot > 1e9:
                raise RuntimeError(f"Simulation exploded at stage {stage}, L={new_L} nm. E={pot:.2e}")

            print(f"[stage {stage:03d}] {cur_L:.2f} -> {new_L:.2f} nm | E={pot:.2e} kJ/mol", flush=True)

            cur_L = float(new_L)
            handle.write(f"{stage},{L_schedule[stage-1]:.3f},{cur_L:.3f},{simulation.currentStep},{pot}\n")
            handle.flush()

    print("[run] Compression Done. Starting Long Relaxation...", flush=True)
    simulation.context.setVelocitiesToTemperature(temp)
    simulation.minimizeEnergy(maxIterations=0)

    print("[run] Warming up...", flush=True)
    integrator.setStepSize(1.0 * unit.femtosecond)
    simulation.step(5000)
    integrator.setStepSize(dt)

    print(f"[run] Long Relax ({final_relax_time_ns} ns)...", flush=True)
    n_chunks = 10
    steps_chunk = max(1, final_relax_steps // n_chunks)
    for i in range(n_chunks):
        simulation.step(steps_chunk)
        st = simulation.context.getState(getEnergy=True)
        e = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"[relax {i+1}/{n_chunks}] E = {e:.2e}", flush=True)

    state = simulation.context.getState(getPositions=True)
    clean_pos = make_whole_and_wrap(state.getPositions(), simulation.topology, state.getPeriodicBoxVectors())

    simulation.topology.setUnitCellDimensions(mm.Vec3(cur_L, cur_L, cur_L))
    with open(out_pdb, "w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(simulation.topology, clean_pos, handle)

    print(f"[save] Saved to {out_pdb}", flush=True)


if __name__ == "__main__":
    main()
