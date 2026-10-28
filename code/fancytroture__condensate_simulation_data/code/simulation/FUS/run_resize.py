#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import openmm as mm
import openmm.app as app
from openmm import unit
from openabc.forcefields.parsers import HPSParser
from openabc.forcefields import HPSModel

# ================= 1) 配置区 =================
n_mol = 200
temp = 280 * unit.kelvin

# 压缩阶段用小步长
dt = 5 * unit.femtosecond
friction = 1.0 / unit.picosecond

initial_L = 300.0  # 初始大盒子
# ✅ 修正：必须是 20.0 nm (200条链的安全密度)
target_L  = 20.0   

# ✅ 改进 1: 压缩因子改为 0.99 (非常温和，防止硬冲突)
compress_factor = 0.99 
relax_steps_each_stage = 5000
minimize_max_iter      = 5000 # 增加最小化能力

final_relax_time_ns    = 10.0
final_relax_steps      = int(final_relax_time_ns * 1000 / (dt.value_in_unit(unit.picosecond)))

monomer_pdb = "fus_monomer_collapsed.pdb"
start_pdb   = f"start_L{int(initial_L)}.pdb"
out_pdb     = f"compressed_L{int(target_L)}_relaxed_10ns.pdb"
log_csv     = "compress_log.csv"

CUDA_DEVICE_INDEX = "0"
CUDA_PRECISION    = "mixed"

urry_mu    = 1.0
urry_delta = 0.08
dh_ldby_nm = 1.0

# ================= 2) 初始结构生成 (晶格 + 随机旋转) =================
def safe_place_molecules(monomer_pdb, output_pdb, n_mol, box_size, seed=42, buffer_nm=1.5):
    print(f"[init] Reading monomer: {monomer_pdb} ...")
    pdb_mono = app.PDBFile(monomer_pdb)
    pos_mono = np.asarray(pdb_mono.positions.value_in_unit(unit.nanometer), dtype=float)
    # 居中单体
    pos_mono -= pos_mono.mean(axis=0)

    # 计算半径
    max_r = np.max(np.linalg.norm(pos_mono, axis=1))
    spacing = 2.0 * max_r + buffer_nm
    print(f"[init] Monomer Radius: {max_r:.2f} nm | Spacing: {spacing:.2f} nm")

    # 晶格计算
    n_side = int(np.ceil(n_mol ** (1.0 / 3.0)))
    needed = (n_side - 1) * spacing + 2.0 * (max_r + buffer_nm)
    
    if needed > box_size:
        print(f"⚠️ 警告: 建议盒子 > {needed:.1f} nm，当前 {box_size} nm。")

    print(f"[init] Generating grid with random rotations...")
    grid = []
    start = max_r + buffer_nm
    for ix in range(n_side):
        x = start + ix * spacing
        for iy in range(n_side):
            y = start + iy * spacing
            for iz in range(n_side):
                z = start + iz * spacing
                grid.append([x, y, z])

    rng = np.random.RandomState(seed)
    indices = np.arange(len(grid))
    rng.shuffle(indices)
    centers = [grid[i] for i in indices[:n_mol]]

    all_positions = []
    modeller = app.Modeller(pdb_mono.topology, pdb_mono.positions)
    
    # 扩展拓扑
    for i in range(n_mol - 1):
        modeller.add(pdb_mono.topology, pdb_mono.positions)

    # ✅ 改进 2: 添加随机旋转
    for i, center in enumerate(centers):
        # 随机旋转矩阵
        from scipy.spatial.transform import Rotation
        # 随机生成一个旋转
        rot_matrix = Rotation.random(random_state=rng).as_matrix()
        
        # 旋转单体 + 平移到晶格点
        new_pos = np.dot(pos_mono, rot_matrix.T) + center
        all_positions.extend([mm.Vec3(x,y,z) for x,y,z in new_pos])

    final_positions = all_positions * unit.nanometer
    top = modeller.topology
    top.setUnitCellDimensions(mm.Vec3(box_size, box_size, box_size))
    
    with open(output_pdb, "w") as f:
        app.PDBFile.writeFile(top, final_positions, f)
    print(f"[init] Success! Saved to {output_pdb}")

# ================= 3) COM 缩放 =================
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

        # Unwrap & Wrap
        c = scaled_coords.copy()
        for i in range(1, len(c)):
            d = c[i] - c[i-1]
            d -= np.round(d / L_new) * L_new
            c[i] = c[i-1] + d

        chain_com = c.mean(axis=0)
        shift = -np.floor(chain_com / L_new) * L_new
        c = c + shift
        new_pos[idx] = c

    return new_pos

# ================= 4) 最终 Wrap =================
def make_whole_and_wrap(positions, topology, box_vectors):
    pos = positions.value_in_unit(unit.nanometer)
    pos = np.asarray(pos, dtype=float)
    if hasattr(box_vectors, 'value_in_unit'):
        box_v = box_vectors.value_in_unit(unit.nanometer)
    else:
        box_v = [v.value_in_unit(unit.nanometer) for v in box_vectors]
    Lx, Ly, Lz = box_v[0][0], box_v[1][1], box_v[2][2]
    L = np.array([Lx, Ly, Lz], dtype=float)

    new_pos = pos.copy()
    for chain in topology.chains():
        idx = np.array([atom.index for atom in chain.atoms()], dtype=int)
        c = new_pos[idx].copy()
        for i in range(1, len(c)):
            d = c[i] - c[i-1]
            d -= np.round(d / L) * L
            c[i] = c[i-1] + d
        center = c.mean(axis=0)
        shift = -np.floor(center / L) * L
        c = c + shift
        new_pos[idx] = c
    return new_pos * unit.nanometer

# ================= 5) 缩盒计划 =================
def build_L_schedule(initial_L, target_L, factor=0.99):
    Ls = [float(initial_L)]
    L = float(initial_L)
    while L * factor > target_L:
        L = L * factor
        Ls.append(float(L))
    if abs(Ls[-1] - float(target_L)) > 1e-12:
        Ls.append(float(target_L))
    return Ls

# ================= 6) 主流程 =================
def main():
    if not os.path.exists(monomer_pdb):
        raise RuntimeError(f"❌ 找不到 {monomer_pdb}！请先运行 make_collapsed_monomer.py")

    # 自动重新生成 start_pdb 以应用随机旋转
    if os.path.exists(start_pdb):
        print(f"Removing old {start_pdb} to apply random rotation...")
        os.remove(start_pdb)
    
    safe_place_molecules(monomer_pdb, start_pdb, n_mol, initial_L)

    parser = HPSParser(monomer_pdb)
    model = HPSModel()
    for _ in range(n_mol):
        model.append_mol(parser)

    pdb_ref = app.PDBFile(start_pdb)
    model.create_system(pdb_ref.getTopology(), box_a=initial_L, box_b=initial_L, box_c=initial_L)
    model.add_protein_bonds(force_group=1)
    model.add_contacts("Urry", mu=urry_mu, delta=urry_delta, force_group=2)
    model.add_dh_elec(ldby=dh_ldby_nm * unit.nanometer, force_group=3)

    integrator = mm.LangevinMiddleIntegrator(temp, friction, dt)

    try:
        platform = mm.Platform.getPlatformByName("CUDA")
        props = {"DeviceIndex": CUDA_DEVICE_INDEX, "Precision": CUDA_PRECISION}
    except:
        platform = mm.Platform.getPlatformByName("CPU")
        props = {}
    print(f"[run] Platform: {platform.getName()}")

    simulation = app.Simulation(model.top, model.system, integrator, platform, props)
    simulation.context.setPositions(pdb_ref.getPositions())
    simulation.context.setVelocitiesToTemperature(temp)

    print("[run] Initial Minimize...")
    simulation.minimizeEnergy(maxIterations=5000)

    print("[run] Pre-equilibration (10000 steps)...")
    simulation.step(10000)
    simulation.minimizeEnergy(maxIterations=1000)

    L_schedule = build_L_schedule(initial_L, target_L, factor=compress_factor)
    print(f"[plan] Compression schedule: {len(L_schedule)} stages (target={target_L} nm, factor={compress_factor})")

    cur_L = float(initial_L)

    with open(log_csv, "w") as f:
        f.write("stage,old_L,new_L,step,pot_E\n")
        
        for stage, new_L in enumerate(L_schedule[1:], start=1):
            print(f"[stage {stage}/{len(L_schedule)-1}] {cur_L:.2f} -> {new_L:.2f} nm", end=' | ')

            state = simulation.context.getState(getPositions=True)
            pos_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            pos_nm2 = scale_positions_by_chain_com_and_wrap(pos_nm, simulation.topology, cur_L, new_L)

            a = mm.Vec3(new_L, 0, 0) * unit.nanometer
            b = mm.Vec3(0, new_L, 0) * unit.nanometer
            c = mm.Vec3(0, 0, new_L) * unit.nanometer
            simulation.context.setPeriodicBoxVectors(a, b, c)
            simulation.context.setPositions(pos_nm2 * unit.nanometer)
            simulation.context.setVelocitiesToTemperature(temp)

            # Safety Valve Logic
            simulation.minimizeEnergy(maxIterations=0) # 强力最小化
            
            # 检查能量，如果还高，用极小步长缓冲
            st = simulation.context.getState(getEnergy=True)
            e_check = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
            if e_check > 1e6:
                integrator.setStepSize(0.1 * unit.femtosecond)
                simulation.step(1000) # 0.1 ps
                integrator.setStepSize(1.0 * unit.femtosecond)
                simulation.step(1000) # 1 ps
                integrator.setStepSize(dt)
            
            simulation.step(relax_steps_each_stage)

            state = simulation.context.getState(getEnergy=True)
            pot = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            
            if not np.isfinite(pot) or pot > 1e9:
                raise RuntimeError(f"💥 Simulation exploded at stage {stage}, L={new_L} nm. E={pot:.2e}")

            print(f"E={pot:.2e} kJ/mol")

            cur_L = float(new_L)
            step = simulation.currentStep
            f.write(f"{stage},{L_schedule[stage-1]:.3f},{cur_L:.3f},{step},{pot}\n")
            f.flush()

    print("\n[run] Compression Done. Starting Long Relaxation...")
    
    simulation.context.setVelocitiesToTemperature(temp)
    simulation.minimizeEnergy(maxIterations=0)
    
    print("[run] Warming up...")
    integrator.setStepSize(1.0 * unit.femtosecond)
    simulation.step(5000)
    integrator.setStepSize(dt)

    print(f"[run] Long Relax ({final_relax_time_ns} ns)...")
    n_chunks = 10
    steps_chunk = final_relax_steps // n_chunks

    for i in range(n_chunks):
        simulation.step(steps_chunk)
        st = simulation.context.getState(getEnergy=True)
        e = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  {i+1}/{n_chunks}: E = {e:.2e}")

    state = simulation.context.getState(getPositions=True)
    clean_pos = make_whole_and_wrap(state.getPositions(), simulation.topology, state.getPeriodicBoxVectors())

    simulation.topology.setUnitCellDimensions(mm.Vec3(cur_L, cur_L, cur_L))
    with open(out_pdb, "w") as f:
        app.PDBFile.writeFile(simulation.topology, clean_pos, f)

    print(f"[save] Saved to {out_pdb}")

if __name__ == "__main__":
    main()