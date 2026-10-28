import os
import sys
import numpy as np
import openmm as mm
import openmm.app as app
from openmm import unit
from openabc.forcefields.parsers import HPSParser
from openabc.forcefields import HPSModel
from openabc.utils.helper_functions import build_straight_CA_chain, write_pdb

# ================= 1. 配置区域 =================
raw_sequence = """MASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQSSNFGPMKGGNFGGRSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF"""
sequence_str = raw_sequence.replace(" ", "").replace("\n", "").strip().upper()
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

n_mol = 200
temp = 230 * unit.kelvin
dt = 10 * unit.femtosecond
friction = 1.0 / unit.picosecond

total_time_ns = 15000 
total_steps = int(total_time_ns * 1000 / (dt.value_in_unit(unit.picosecond)))

# ✅ 修正 1: 步数计算
# 1 ns = 100,000 steps
steps_per_ns = int(1.0 / dt.value_in_unit(unit.nanosecond)) 
chunk_size_ns = 100
steps_per_chunk = chunk_size_ns * steps_per_ns # 100 ns = 10,000,000 步
n_chunks = total_time_ns // chunk_size_ns

interval_kinetic = 1000 
interval_thermo = 10000
file_thermo = "a1_slab_phase_diagram.dcd"

slab_Lx = 20.0
slab_Ly = 20.0
slab_Lz = 280.0

job_name = "a1_slab_vampnet"
input_compressed_pdb = 'compressed_relaxed_10ns.pdb'
start_slab_pdb = f'{job_name}_start.pdb'
monomer_pdb = 'a1_monomer.pdb'
chk_file = f'{job_name}_checkpoint.chk'

CUDA_DEVICE = os.environ.get("CUDA_DEVICE_INDEX", os.environ.get("OPENMM_CUDA_DEVICE", "0"))


def validate_sequence_or_raise(seq):
    if (not seq) or ("REPLACE_WITH_" in seq):
        raise RuntimeError(
            "A1 sequence placeholder is still present.\n"
            "Replace raw_sequence in run_a1_condensate.py, or run:\n"
            "  python make_monomer_from_sequence.py --sequence YOUR_A1_LCD_SEQUENCE"
        )

    bad = sorted(set(seq) - VALID_AA)
    if bad:
        raise RuntimeError(
            "Sequence contains unsupported residue codes: {}\n"
            "Please keep only standard one-letter amino-acid codes.".format(
                ",".join(bad)
            )
        )

# ================= 2. 构建 Slab (带 Unwrap 功能) =================
def build_centered_slab(input_pdb, output_pdb, Lx, Ly, Lz):
    print(f"Building Slab with unwrap: {input_pdb} -> {output_pdb}")
    pdb = app.PDBFile(input_pdb)
    pos = pdb.positions.value_in_unit(unit.nanometer)
    pos = np.array(pos)
    old_vecs = pdb.topology.getUnitCellDimensions().value_in_unit(unit.nanometer)
    L_old = np.array([old_vecs.x, old_vecs.y, old_vecs.z])
    
    # ✅ 修正 3: 先 make whole (缝合跨界链)
    new_pos = pos.copy()
    for chain in pdb.topology.chains():
        idx = [a.index for a in chain.atoms()]
        c = new_pos[idx]
        # 沿所有轴 unwrap
        for i in range(1, len(c)):
            delta = c[i] - c[i-1]
            delta -= np.round(delta / L_old) * L_old
            c[i] = c[i-1] + delta
        new_pos[idx] = c
    
    # 再整体平移
    # 计算 unwrap 后的几何中心
    center = np.mean(new_pos, axis=0)
    # 目标中心
    target_center = np.array([Lx/2.0, Ly/2.0, Lz/2.0])
    shift = target_center - center
    
    new_pos += shift
    
    top = pdb.topology
    top.setUnitCellDimensions(mm.Vec3(Lx, Ly, Lz))
    with open(output_pdb, 'w') as f:
        app.PDBFile.writeFile(top, new_pos * unit.nanometer, f)
    print(f"✅ Slab generated. Centered.")

if not os.path.exists(start_slab_pdb):
    if not os.path.exists(input_compressed_pdb):
        raise FileNotFoundError(
            "Missing input: {}\n"
            "Prepare inputs in this order:\n"
            "1. python make_monomer_from_sequence.py\n"
            "2. prepare a1_monomer_collapsed.pdb\n"
            "   or use python copy_monomer_to_collapsed.py as a temporary placeholder\n"
            "3. python run_resize.py\n".format(input_compressed_pdb)
        )
    build_centered_slab(input_compressed_pdb, start_slab_pdb, slab_Lx, slab_Ly, slab_Lz)

# ================= 3. 系统设置 =================
if not os.path.exists(monomer_pdb):
    validate_sequence_or_raise(sequence_str)
    print(f"Generating straight monomer from sequence -> {monomer_pdb}")
    atoms = build_straight_CA_chain(sequence_str, r0=0.38)
    write_pdb(atoms, monomer_pdb)

try:
    protein_parser = HPSParser(monomer_pdb)
except:
    content = open(monomer_pdb, 'r').read().replace('\r\n', '\n')
    with open(monomer_pdb, 'w') as f:
        f.write(content)
    protein_parser = HPSParser(monomer_pdb)

model = HPSModel()
for _ in range(n_mol):
    model.append_mol(protein_parser)

pdb_ref = app.PDBFile(start_slab_pdb)
model.create_system(pdb_ref.topology, box_a=slab_Lx, box_b=slab_Ly, box_c=slab_Lz)
model.add_protein_bonds(force_group=1)
model.add_contacts('Urry', mu=1.0, delta=0.08, force_group=2)
model.add_dh_elec(ldby=1.0 * unit.nanometer, force_group=3)

# ================= 4. 模拟环境 =================
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': CUDA_DEVICE, 'Precision': 'mixed'}
except:
    platform = mm.Platform.getPlatformByName('CPU')
    properties = {}
print(f"Platform: {platform.getName()}")

integrator = mm.LangevinMiddleIntegrator(temp, friction, dt)
simulation = app.Simulation(model.top, model.system, integrator, platform, properties)

# 断点续传逻辑
start_chunk_idx = 0
append_mode = False

# 正确的盒子向量定义
box_vectors = (
    mm.Vec3(slab_Lx, 0, 0) * unit.nanometer,
    mm.Vec3(0, slab_Ly, 0) * unit.nanometer,
    mm.Vec3(0, 0, slab_Lz) * unit.nanometer
)

if os.path.exists(chk_file):
    print(f"Loading checkpoint: {chk_file}")
    try:
        with open(chk_file, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())
        simulation.context.setPeriodicBoxVectors(*box_vectors) # 强制重设
        current_step = simulation.context.getState().getStepCount()
        start_chunk_idx = int(current_step / steps_per_chunk)
        print(f"Resuming from Chunk {start_chunk_idx + 1} / {n_chunks}")
        append_mode = True
    except Exception as e:
        print(f"Checkpoint failed ({e}). Restarting.")
        # ✅ 修正 2: `except` 块中的 `setPeriodicBoxVectors`
        simulation.context.setPositions(pdb_ref.getPositions())
        simulation.context.setPeriodicBoxVectors(*box_vectors)
        simulation.context.setVelocitiesToTemperature(temp)
        print("Minimizing energy...")
        simulation.minimizeEnergy()
else:
    print("Starting new simulation...")
    simulation.context.setPositions(pdb_ref.getPositions())
    simulation.context.setPeriodicBoxVectors(*box_vectors)
    simulation.context.setVelocitiesToTemperature(temp)
    print("Minimizing energy...")
    simulation.minimizeEnergy()

# ================= 5. 配置 Reporters =================
simulation.reporters.append(app.StateDataReporter(f"{job_name}_global.csv", interval_thermo, step=True, potentialEnergy=True, temperature=True, speed=True, remainingTime=True, totalSteps=total_steps, append=append_mode))
simulation.reporters.append(app.DCDReporter(file_thermo, interval_thermo, append=append_mode))
simulation.reporters.append(app.CheckpointReporter(chk_file, 1000000))
# ================= 6. 运行分段循环 =================
print(f"\n🚀 Simulation Started")
print(f"  Stream 1 (VAMPnet): {chunk_size_ns} ns files @ 10 ps interval")
print(f"  Stream 2 (PhaseDiagram): Single file @ 100 ps interval")
print(f"  Steps per chunk: {steps_per_chunk}")
print("-" * 60)
for chunk_idx in range(start_chunk_idx, n_chunks):
    kinetic_traj_name = f"{job_name}_part{chunk_idx+1:03d}.dcd"
    print(f"▶️ Running Chunk {chunk_idx+1}/{n_chunks} -> {kinetic_traj_name}")
    hi_res_reporter = app.DCDReporter(kinetic_traj_name, interval_kinetic)
    simulation.reporters.append(hi_res_reporter)
    simulation.step(steps_per_chunk)
    simulation.reporters.remove(hi_res_reporter)
    with open(chk_file, 'wb') as f:
        f.write(simulation.context.createCheckpoint())
    print(f"✅ Chunk {chunk_idx+1} Done.")
print("-" * 60)
print("🎉 All Done!")
