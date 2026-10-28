"""Frozen input validation for the three additional coexistence temperatures."""
import hashlib
import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BOX = np.array([20., 20., 280.])
CASES = ("FUS_335K", "A1_267K", "TDP43_315K")
SETTINGS = dict(n_chains=200, box_nm=BOX.tolist(), dt_fs=10., friction_per_ps=1.,
                integrator="LangevinMiddleIntegrator", total_time_ns=15000,
                chunk_ns=100, kinetic_interval_steps=1000, thermo_interval_steps=10000,
                checkpoint_interval_steps=1000000, cuda_precision="mixed",
                integrator_seed=0, initial_velocity_seed=None,
                initial_velocity_policy="OpenMM setVelocitiesToTemperature at target T; omit seed, original osrngseed default",
                minimize_tolerance_kj_mol_nm=10., minimize_max_iterations=0,
                hps_scale="Urry", mu=1., delta=.08, epsilon_kj_mol=.8368,
                debye_length_nm=1., dielectric_water=80., dh_cutoff_nm=3.5,
                barostat=False, positional_restraints=False)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def signatures(system):
    return [dict(type=type(force).__name__, group=force.getForceGroup(),
                 sha256=hashlib.sha256(mm.XmlSerializer.serialize(force).encode()).hexdigest())
            for force in system.getForces()]


def integrator(config):
    result = mm.LangevinMiddleIntegrator(config["temperature_K"] * unit.kelvin,
        config["friction_per_ps"] / unit.picosecond, config["dt_fs"] * unit.femtosecond)
    result.setRandomNumberSeed(config["integrator_seed"])
    return result


def load_case(case):
    if case not in CASES:
        raise ValueError(case)
    folder = HERE / case
    manifest = json.loads((folder / "manifest.json").read_text())
    for name, expected in manifest["files"].items():
        if sha(folder / name) != expected:
            raise RuntimeError(f"Frozen input changed: {folder / name}")
    config = json.loads((folder / "config.json").read_text())
    for name, value in SETTINGS.items():
        if config[name] != value:
            raise RuntimeError(f"Physical setting mismatch: {name}")
    if config["temperature_K"] != dict(FUS_335K=335, A1_267K=267, TDP43_315K=315)[case]:
        raise RuntimeError("Unexpected target temperature")
    if manifest["openmm_version"] != mm.__version__:
        raise RuntimeError("OpenMM version changed; revalidate before running")
    system = mm.XmlSerializer.deserialize((folder / "inputs/system.xml").read_text())
    pdb = app.PDBFile(str(folder / "inputs/start_slab.pdb"))
    assert signatures(system) == manifest["forces"]
    assert system.getNumParticles() == pdb.topology.getNumAtoms() == 200 * config["chain_length"]
    assert system.getNumConstraints() == 0
    np.testing.assert_allclose(np.array([v.value_in_unit(unit.nanometer)
        for v in system.getDefaultPeriodicBoxVectors()]), np.diag(BOX), atol=1e-10)
    np.testing.assert_allclose(pdb.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer), np.diag(BOX))
    return folder, manifest, config, system, pdb


def check_state(simulation, config, expected_step):
    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    assert simulation.currentStep == expected_step == state.getStepCount()
    np.testing.assert_allclose(state.getTime().value_in_unit(unit.picosecond),
                               expected_step * config["dt_fs"] / 1000, rtol=1e-7, atol=1e-5)
    np.testing.assert_allclose(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer), np.diag(BOX))
    for data in (state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
                 state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)):
        assert np.isfinite(data).all()
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    assert np.isfinite(energy)
    it = simulation.integrator
    assert isinstance(it, mm.LangevinMiddleIntegrator)
    np.testing.assert_allclose(it.getTemperature().value_in_unit(unit.kelvin), config["temperature_K"])
    np.testing.assert_allclose(it.getStepSize().value_in_unit(unit.femtosecond), config["dt_fs"])
    np.testing.assert_allclose(it.getFriction().value_in_unit(unit.picosecond**-1), config["friction_per_ps"])
    return dict(step=expected_step, time_ns=state.getTime().value_in_unit(unit.nanosecond),
                potential_kj_mol=energy)
