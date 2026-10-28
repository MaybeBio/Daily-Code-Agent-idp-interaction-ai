#!/usr/bin/env python3
"""A1/TDP43 paired transient slabs, adapted from the frozen FUS pilot protocol."""
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

CODE = Path(__file__).resolve().parent
REVISION = CODE.parent
ROOT = REVISION.parent
CASES = {
    'A1': dict(chain_length=135, high=260, low=230, seed_base=202610180,
               reference='A1_267K', source='A1/A1_condensate_simulation/condensate_260K_200chains', prefix='a1'),
    'TDP43': dict(chain_length=148, high=300, low=260, seed_base=202610280,
                  reference='TDP43_315K', source='TDP43/tdp43_condensate_simulation/condensate_300K_200chains', prefix='tdp43'),
}
HERE = None
CASE = None
SETTINGS = None


def configure(case):
    global HERE, CASE, SETTINGS
    SETTINGS = CASES[case]
    CASE = case
    HERE = REVISION / f'{case}_slab_temperature_jump'
    return HERE

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
    x = np.asarray(x, float).reshape(200, SETTINGS["chain_length"], 3).copy()
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


def validate():
    manifest = json.loads((HERE / 'manifest.json').read_text())
    assert manifest['openmm_version'] == mm.__version__
    assert manifest['case'] == CASE
    assert manifest['branch_temperatures_K'] == [SETTINGS['high'], SETTINGS['low']]
    assert manifest['protocol_sha256'] == sha(__file__)
    assert manifest['dt_fs'] == DT_FS and manifest['friction_per_ps'] == 1.
    assert manifest['max_elapsed_ns'] == 20. and manifest['chunk_ns'] == CHUNK_NS
    assert not any(manifest[k] for k in ('minimization', 'position_scaling', 'barostat', 'position_restraints'))
    for name, expected in manifest['inputs'].items():
        assert sha(HERE / 'inputs' / name) == expected, f'Input changed: {name}'
    system = mm.XmlSerializer.deserialize((HERE / 'inputs/system.xml').read_text())
    assert signatures(system) == manifest['reference_force_signatures']
    box = np.array([v.value_in_unit(unit.nanometer) for v in system.getDefaultPeriodicBoxVectors()])
    np.testing.assert_allclose(box, np.diag(BOX))
    assert system.getNumParticles() == 200 * SETTINGS["chain_length"] and system.getNumConstraints() == 0
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
    assert temperature in (SETTINGS['high'], SETTINGS['low'])
    assert 1 <= snapshot <= 8
    seed = SETTINGS['seed_base'] + temperature + snapshot
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
                assert len(reader) == 50 and reader.n_atoms == 200 * SETTINGS["chain_length"]
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
    p.add_argument('--case', choices=CASES, required=True)
    sub = p.add_subparsers(dest='command', required=True)
    sub.add_parser('validate')
    r = sub.add_parser('run')
    r.add_argument('--snapshot', type=int, choices=range(1, 9), required=True)
    r.add_argument('--temperature', type=int, required=True)
    r.add_argument('--until-ns', type=float, default=2.)
    args = p.parse_args()
    configure(args.case)
    if args.command == 'validate':
        validate(); print('Static validation passed; no Context created')
    else:
        run(args.snapshot, args.temperature, args.until_ns)


if __name__ == '__main__':
    main()
