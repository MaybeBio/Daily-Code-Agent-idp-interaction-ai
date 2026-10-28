#!/usr/bin/env python3
"""Paired slab temperature jump. All times are elapsed since velocity initialization."""
import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
import numpy as np
import openmm as mm
from openmm import app, unit
from MDAnalysis.coordinates.DCD import DCDReader
from mdtraj.formats import DCDTrajectoryFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREVIOUS = ROOT / 'revision_2026_09_30'
BOX = np.array([20., 20., 280.])
EDGES = np.arange(-140., 140.001, .5)
GPU = 'GPU-29d16388-8ec1-2739-6b28-26cbc21152d3'
DT_FS = 10.
CHUNK_NS = .5
CHUNK_STEPS = 50000
INTERVAL_STEPS = 1000
THERMO_STEPS = 100


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def signatures(system):
    return [hashlib.sha256(mm.XmlSerializer.serialize(f).encode()).hexdigest() for f in system.getForces()]


def whole(x):
    x = np.asarray(x, float).reshape(200, 163, 3).copy()
    dx = np.diff(x, axis=1)
    dx -= np.rint(dx / BOX) * BOX
    x[:, 1:] = x[:, :1] + np.cumsum(dx, axis=1)
    return x


def measure(x, mass):
    z = np.remainder(x[:, 2], BOX[2])
    phase = 2 * np.pi * z / BOX[2]
    phasor = np.sum(mass * np.exp(1j * phase))
    center = np.angle(phasor) * BOX[2] / (2 * np.pi)
    centered = np.remainder(z - center + BOX[2] / 2, BOX[2]) - BOX[2] / 2
    rho = np.histogram(centered, bins=EDGES, weights=mass)[0] * (1e21 / 6.02214076e23) / 200
    np.testing.assert_allclose(rho.sum() * 200 / (1e21 / 6.02214076e23), mass.sum(), rtol=1e-12)
    zc = (EDGES[:-1] + EDGES[1:]) / 2
    core = [float(rho[abs(zc) < radius].mean()) for radius in [3, 5, 7]]
    return rho, core, float(center), float(abs(phasor) / mass.sum())


def prepare():
    if (HERE / 'manifest.json').exists():
        raise RuntimeError('Already prepared; verify existing inputs instead of replacing them')
    inputs = HERE / 'inputs'
    inputs.mkdir(exist_ok=False)
    reference = PREVIOUS / 'matched_density_FUS'
    original_manifest = json.loads((reference / 'prepared.json').read_text())
    assert sha(reference / 'inputs/system.xml') == original_manifest['inputs']['system.xml']
    system = mm.XmlSerializer.deserialize((reference / 'inputs/system.xml').read_text())
    original_signature = signatures(system)
    assert original_signature == [f['sha256'] for f in original_manifest['force_signature']]
    system.setDefaultPeriodicBoxVectors(*[mm.Vec3(*r) * unit.nanometer for r in np.diag(BOX)])
    assert signatures(system) == original_signature
    assert [type(f).__name__ for f in system.getForces()] == [
        'CMMotionRemover', 'HarmonicBondForce', 'CustomNonbondedForce', 'CustomNonbondedForce']
    (inputs / 'system.xml').write_text(mm.XmlSerializer.serialize(system))
    source = ROOT / 'fus/fus_condensate_simulation/condensate_325K_200chains'
    shutil.copy2(source / 'fus_slab_vampnet_start.pdb', inputs / 'topology.pdb')
    topology = app.PDBFile(str(inputs / 'topology.pdb')).topology
    assert len(list(topology.chains())) == 200 and topology.getNumAtoms() == 32600
    assert all(len(list(c.atoms())) == 163 for c in topology.chains())
    mass = np.array([system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(32600)])
    inventory = json.loads((PREVIOUS / 'fus_critical_temperature/profiles/FUS_325_sources.json').read_text())
    snapshots = []
    # Prespecified times, not chosen by their later response or by density closeness.
    for number, desired in enumerate(range(10000, 14000, 500), 1):
        part = next(p for p in inventory['parts'] if p['first_header_time_ns'] <= desired <= p['last_header_time_ns'])
        path = ROOT / part['path']
        assert path.stat().st_size == part['bytes']
        local = round((desired - part['first_header_time_ns']) / part['raw_dt_ns'])
        with DCDReader(str(path)) as reader:
            frame = reader[local]
            np.testing.assert_allclose(frame.dimensions[:3] / 10, BOX)
            np.testing.assert_allclose(frame.dimensions[3:], 90)
            x = frame.positions.astype(float) / 10
            actual_time = float(frame.time / 1000)
            baseline_profiles = [measure(reader[i].positions.astype(float) / 10, mass)[0]
                                 for i in range(local - 50, local)]
        with DCDTrajectoryFile(str(path)) as reader:
            reader.seek(local)
            xyz, lengths, angles = reader.read(n_frames=1)
        np.testing.assert_allclose(x, xyz[0].astype(float) / 10, rtol=0, atol=0)
        chain = whole(x)
        bonds = np.linalg.norm(np.diff(chain, axis=1), axis=2)
        assert bonds.min() > .2 and bonds.max() < .7
        rng = np.random.default_rng(20261008 + number)
        # R in kJ/(mol K), mass in daltons: velocity unit is nm/ps.
        velocity = rng.normal(size=x.shape) * np.sqrt(.00831446261815324 * 325 / mass[:, None])
        velocity -= np.sum(velocity * mass[:, None], axis=0) / mass.sum()
        rho, core, center, order = measure(x, mass)
        com = np.sum(chain * mass.reshape(200, 163, 1), axis=1) / mass.reshape(200, 163).sum(1)[:, None]
        center_z = np.remainder(com[:, 2] - center + 140, 280) - 140
        cohort = np.flatnonzero(abs(center_z) < 5)
        baseline = np.mean(baseline_profiles, axis=0)
        np.savez_compressed(inputs / f'snapshot{number:02d}.npz', positions_nm=x,
            velocities_nm_ps=velocity, masses_da=mass, initial_profile=rho,
            baseline_profile=baseline, baseline_profiles=baseline_profiles,
            initial_center_chain_ids=cohort, edges_nm=EDGES)
        snapshots.append(dict(snapshot=number, desired_source_time_ns=desired, source_time_ns=actual_time,
            source_path=str(path), source_local_frame=local, source_part=part['part'],
            source_file_bytes=path.stat().st_size, source_mtime_ns=path.stat().st_mtime_ns,
            selected_coordinates_sha256=hashlib.sha256(x.tobytes()).hexdigest(),
            independent_raw_reader_exact_match=True, core_density_g_ml=core,
            baseline_500ps_core5_g_ml=float(baseline[abs((EDGES[:-1] + EDGES[1:]) / 2) < 5].mean()),
            mass_weighted_initial_temperature_K=float(np.sum(mass[:, None] * velocity**2) / ((3 * len(mass) - 3) * .00831446261815324)),
            initial_center_chains=len(cohort), circular_order=order))
    manifest = dict(openmm_version=mm.__version__, source_temperature_K=325,
        branch_temperatures_K=[325, 280], snapshots=snapshots, box_nm=BOX.tolist(),
        dt_fs=DT_FS, friction_per_ps=1., initial_velocity_temperature_K=325,
        initial_velocity_policy='One Maxwell draw per snapshot, COM removed, shared by both branches; not original DCD velocities',
        minimization=False, position_scaling=False, barostat=False, position_restraints=False,
        max_elapsed_ns=20., chunk_ns=CHUNK_NS, trajectory_interval_ps=10., thermo_interval_ps=1.,
        analysis_start_ns=.1, density_block_ns=.5, density_tolerances=[.02, .05],
        reference_force_signatures=original_signature,
        inputs={p.name: sha(p) for p in inputs.iterdir()},
        limitations=['Single source slab; separated snapshots are not automatically independent',
                    'Transient response, not low-temperature equilibrium kinetics'])
    write_json(HERE / 'manifest.json', manifest)
    print(json.dumps(snapshots, indent=2))


def validate():
    manifest = json.loads((HERE / 'manifest.json').read_text())
    assert manifest['openmm_version'] == mm.__version__
    for name, expected in manifest['inputs'].items():
        assert sha(HERE / 'inputs' / name) == expected, f'Input changed: {name}'
    system = mm.XmlSerializer.deserialize((HERE / 'inputs/system.xml').read_text())
    assert signatures(system) == manifest['reference_force_signatures']
    box = np.array([v.value_in_unit(unit.nanometer) for v in system.getDefaultPeriodicBoxVectors()])
    np.testing.assert_allclose(box, np.diag(BOX))
    assert system.getNumParticles() == 32600 and system.getNumConstraints() == 0
    return manifest, system


def check_state(sim, temperature):
    s = sim.context.getState(getPositions=True, getEnergy=True)
    np.testing.assert_allclose(s.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer), np.diag(BOX))
    assert np.isfinite(s.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).all()
    energy = s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    assert np.isfinite(energy)
    np.testing.assert_allclose(sim.integrator.getTemperature().value_in_unit(unit.kelvin), temperature)
    np.testing.assert_allclose(sim.integrator.getFriction().value_in_unit(unit.picosecond**-1), 1.)
    np.testing.assert_allclose(sim.integrator.getStepSize().value_in_unit(unit.femtosecond), DT_FS)
    return dict(step=sim.currentStep, time_ns=s.getTime().value_in_unit(unit.nanosecond), potential_kj_mol=energy)


def completed(path, identity):
    m = json.loads((path / 'complete.json').read_text())
    assert m['identity'] == identity
    for name, expected in m['files'].items():
        assert sha(path / name) == expected, f'Completed file changed: {path/name}'
    return m


def run(snapshot, temperature, until_ns):
    manifest, system = validate()
    if os.environ.get('CUDA_VISIBLE_DEVICES') != GPU:
        raise RuntimeError('Set CUDA_VISIBLE_DEVICES to the verified physical GPU 1 UUID')
    if not 0 < until_ns <= manifest['max_elapsed_ns'] or until_ns / CHUNK_NS != round(until_ns / CHUNK_NS):
        raise ValueError('until-ns must be a positive multiple of 0.5, at most 20')
    seed = 202610080 + temperature + snapshot
    identity = dict(manifest_sha256=sha(HERE / 'manifest.json'), code_sha256=sha(__file__),
                    snapshot=snapshot, temperature_K=temperature, integrator_seed=seed)
    root = HERE / 'runs' / f'snapshot{snapshot:02d}_{temperature}K'
    root.mkdir(parents=True, exist_ok=True)
    with (root / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if list(root.glob('*.incomplete')):
            raise RuntimeError('Interrupted output retained; archive it before retrying from last complete segment')
        marker = root / 'identity.json'
        if marker.exists():
            assert json.loads(marker.read_text()) == identity
        else:
            write_json(marker, identity)
        pdb = app.PDBFile(str(HERE / 'inputs/topology.pdb'))
        integrator = mm.LangevinMiddleIntegrator(temperature * unit.kelvin, 1. / unit.picosecond, DT_FS * unit.femtosecond)
        integrator.setRandomNumberSeed(seed)
        sim = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatformByName('CUDA'),
                             {'DeviceIndex': '0', 'Precision': 'mixed'})
        data = np.load(HERE / 'inputs' / f'snapshot{snapshot:02d}.npz')
        mass = data['masses_da']
        sim.context.setPositions(data['positions_nm'] * unit.nanometer)
        sim.context.setVelocities(data['velocities_nm_ps'] * unit.nanometer / unit.picosecond)
        sim.context.setTime(0 * unit.picosecond)
        sim.context.setStepCount(0)
        check_state(sim, temperature)
        write_json(root / f'launch_{time.time_ns()}.json', dict(pid=os.getpid(), gpu_uuid=GPU, until_ns=until_ns, identity=identity))
        for n in range(1, round(until_ns / CHUNK_NS) + 1):
            final = root / f'part{n:04d}'
            if final.exists():
                meta = completed(final, identity)
                assert meta['begin']['step'] == (n - 1) * CHUNK_STEPS
                assert meta['end']['step'] == n * CHUNK_STEPS
                sim.context.loadCheckpoint((final / 'checkpoint.chk').read_bytes())
                state = check_state(sim, temperature)
                assert state['step'] == meta['end']['step']
                np.testing.assert_allclose(state['time_ns'], n * CHUNK_NS, atol=1e-6)
                continue
            work = root / f'part{n:04d}.incomplete'
            work.mkdir()
            begin = check_state(sim, temperature)
            assert begin['step'] == (n - 1) * CHUNK_STEPS
            dcd = app.DCDReporter(str(work / 'trajectory.dcd'), INTERVAL_STEPS, enforcePeriodicBox=True)
            start = time.monotonic()
            with (work / 'thermo.csv').open('w') as f:
                thermo = app.StateDataReporter(f, THERMO_STEPS, step=True, time=True, temperature=True,
                    potentialEnergy=True, kineticEnergy=True, volume=True, density=True, speed=True)
                sim.reporters.extend([dcd, thermo])
                try:
                    sim.step(CHUNK_STEPS)
                finally:
                    sim.reporters.clear(); dcd._out.close()
            end = check_state(sim, temperature)
            np.testing.assert_allclose(end['time_ns'], n * CHUNK_NS, atol=1e-6)
            records, profiles = [], []
            with DCDReader(str(work / 'trajectory.dcd')) as reader:
                assert len(reader) == 50 and reader.n_atoms == 32600
                np.testing.assert_allclose(reader.dt, 10., rtol=1e-5)
                np.testing.assert_allclose(reader[0].time / 1000, begin['time_ns'] + .01, atol=1e-5)
                np.testing.assert_allclose(reader[-1].time / 1000, end['time_ns'], atol=1e-5)
                for frame in reader:
                    np.testing.assert_allclose(frame.dimensions[:3] / 10, BOX)
                    rho, core, center, order = measure(frame.positions.astype(float) / 10, mass)
                    profiles.append(rho)
                    records.append([frame.time / 1000, *core, order])
            np.savez_compressed(work / 'density.npz', records=records, profiles=profiles, edges_nm=EDGES)
            (work / 'checkpoint.chk').write_bytes(sim.context.createCheckpoint())
            state = sim.context.getState(getPositions=True, getVelocities=True, getParameters=True)
            (work / 'state.xml').write_text(mm.XmlSerializer.serialize(state))
            meta = dict(identity=identity, begin=begin, end=end, elapsed_s=time.monotonic() - start,
                        files={p.name: sha(p) for p in work.iterdir() if p.is_file()})
            write_json(work / 'complete.json', meta)
            work.rename(final)
            print(json.dumps(dict(snapshot=snapshot, temperature_K=temperature, elapsed_ns=end['time_ns'],
                                   core5_mean_g_ml=float(np.mean(np.asarray(records)[:, 2])), wall_seconds=meta['elapsed_s'])), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    sub.add_parser('prepare'); sub.add_parser('validate')
    r = sub.add_parser('run')
    r.add_argument('--snapshot', type=int, choices=range(1, 9), required=True)
    r.add_argument('--temperature', type=int, choices=[280, 325], required=True)
    r.add_argument('--until-ns', type=float, default=2.)
    args = p.parse_args()
    if args.command == 'prepare':
        prepare()
    elif args.command == 'validate':
        validate(); print('Static validation passed; no Context created')
    else:
        run(args.snapshot, args.temperature, args.until_ns)


if __name__ == '__main__':
    main()
