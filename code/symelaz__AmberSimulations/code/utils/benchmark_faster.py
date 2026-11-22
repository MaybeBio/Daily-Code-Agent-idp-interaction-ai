"""
benchmark_faster.py <parm7> <coordinates> [--steps N] [--stability-ps PS] [--label TEXT]

Classic (stock OpenMM) vs --faster (NVIDIA cuPotential) on a built system,
with the settings run.sh uses for production (inputs/*/step5_production.inp
with hmr = yes): PME, 0.9 nm cutoff, ewaldErrorTolerance 5e-4, LJ dispersion
correction, HBonds + rigid water, HMR 3.024 Da (water excluded), LangevinMiddle
4 fs, 303.15 K, 1/ps, isotropic MonteCarloBarostat 1 bar every 25 steps
(OpenMM's default, which omm_barostat.py uses for isotropic coupling).

<coordinates>: an equilibrated OpenMM state (step4.*.rst / step5_*.rst, XML;
recommended, it has velocities and the equilibrated box) or step3_input.rst7.

Reports
  1. accuracy: energy and force error of each engine vs a converged reference
     (stock OpenMM, double precision, ewaldErrorTolerance 1e-6, same LJ model)
  2. speed: ns/day, classic (single, = production) and faster (mixed, = run.sh
     --faster; also single for information), each with its own integrator/barostat
  3. stability: NPT from the same state, classic vs faster: temperature,
     potential energy, box volume

Run on a GPU node in the openmm_faster env. Takes ~5-15 min depending on size.
"""
import argparse
import time

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u

import nv_openmm_plugin as nv

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("parm7")
ap.add_argument("coordinates", help="OpenMM XML state (.rst) or Amber .rst7")
ap.add_argument("--steps", type=int, default=5000, help="timed steps per engine (default 5000)")
ap.add_argument("--stability-ps", type=float, default=200, help="NPT length per engine, ps (default 200, 0 = skip)")
ap.add_argument("--label", default=None)
a = ap.parse_args()

T, FRIC, DT, P, PFREQ, WARMUP = 303.15 * u.kelvin, 1.0 / u.picosecond, 0.004 * u.picoseconds, 1.0 * u.bar, 25, 500
WATER = {'HOH', 'WAT', 'TIP3', 'T4E', 'TP4', 'TIP4', 'OPC', 'SOL', 'SPC'}

prmtop = app.AmberPrmtopFile(a.parm7)
N = prmtop.topology.getNumAtoms()
if a.coordinates.endswith((".rst7", ".inpcrd")):
    crd = app.AmberInpcrdFile(a.coordinates)
    pos0, vel0, box0 = crd.positions, None, crd.boxVectors
    print("NOTE: .rst7 input (not equilibrated): accuracy is fine to compare, but use an "
          "equilibrated .rst state for speed/stability numbers.")
else:
    st = mm.XmlSerializer.deserialize(open(a.coordinates).read())
    pos0, vel0, box0 = st.getPositions(), st.getVelocities(), st.getPeriodicBoxVectors()
plat = mm.Platform.getPlatformByName("CUDA")
print(f"### {a.label or a.parm7}: {N} atoms, OpenMM {mm.__version__}", flush=True)


def hmr(system, hmass=3.024):  # = utils/omm_hmr.py
    for a1, a2 in prmtop.topology.bonds():
        if a1.residue.name.upper() in WATER:
            continue
        if a1.element == app.element.hydrogen:
            a1, a2 = a2, a1
        if a2.element != app.element.hydrogen or a1.element in (app.element.hydrogen, None):
            continue
        extra = hmass * u.dalton - system.getParticleMass(a2.index)
        system.setParticleMass(a2.index, hmass * u.dalton)
        system.setParticleMass(a1.index, system.getParticleMass(a1.index) - extra)
    return system


def make_system(ewald_tol=5e-4):
    s = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.9 * u.nanometer,
                            constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=ewald_tol)
    for f in s.getForces():
        if isinstance(f, mm.NonbondedForce):
            f.setUseDispersionCorrection(True)
    return hmr(s)


def to_nv(s):
    i = [k for k, f in enumerate(s.getForces()) if isinstance(f, mm.NonbondedForce)][0]
    s.addForce(nv.NonbondedForce(s.getForce(i)))
    s.removeForce(i)
    for k, f in enumerate(s.getForces()):
        if isinstance(f, mm.CMMotionRemover):
            s.removeForce(k)
            break
    return s


def context(s, integ, precision):
    ctx = mm.Context(s, integ, plat, {"Precision": precision})
    ctx.setPeriodicBoxVectors(*box0)
    ctx.setPositions(pos0)
    if vel0 is not None:
        ctx.setVelocities(vel0)
    else:
        ctx.setVelocitiesToTemperature(T)
    return ctx


def energy_forces(s, precision):
    ctx = context(s, mm.VerletIntegrator(0.001), precision)
    st = ctx.getState(getEnergy=True, getForces=True)
    return (st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole),
            np.asarray(st.getForces(asNumpy=True).value_in_unit(u.kilojoule_per_mole / u.nanometer)))


# ---------------------------------------------------------------- 1. accuracy
print("\n## 1. Accuracy vs converged reference (stock, double, ewaldErrorTolerance 1e-6, same LJ)", flush=True)
e_ref, f_ref = energy_forces(make_system(1e-6), "double")
real = np.linalg.norm(f_ref, axis=1) > 0  # massless virtual sites carry no force
frms = np.sqrt(np.mean(np.sum(f_ref[real] ** 2, axis=1)))
print(f"{'engine':32s} {'dE (kJ/mol)':>12s} {'rel F error':>12s} {'max |dF|':>9s}")
for name, mk, prec in [("classic single (production)", make_system, "single"),
                       ("faster mixed (run.sh --faster)", lambda: to_nv(make_system()), "mixed")]:
    e, f = energy_forces(mk(), prec)
    d = f[real] - f_ref[real]
    print(f"{name:32s} {e - e_ref:12.2f} {np.sqrt(np.mean(np.sum(d ** 2, axis=1))) / frms:12.2e} "
          f"{np.abs(d).max():9.1f}", flush=True)


# ---------------------------------------------------------------- 2. speed
def sim(engine, precision):
    s = make_system()
    if engine == "faster":
        s = to_nv(s)
        s.addForce(nv.MonteCarloBarostat(P, T, frequency=PFREQ))
        integ = nv.LangevinMiddleIntegrator(T, FRIC, DT)
    else:
        s.addForce(mm.MonteCarloBarostat(P, T, PFREQ))
        integ = mm.LangevinMiddleIntegrator(T, FRIC, DT)
    return context(s, integ, precision), integ, s


print(f"\n## 2. Speed: {WARMUP} warm-up + {a.steps} timed steps, 4 fs, NPT", flush=True)
speed = {}
for engine, prec in [("classic", "single"), ("faster", "mixed"), ("faster", "single")]:
    ctx, integ, _ = sim(engine, prec)
    integ.step(WARMUP)
    ctx.getState(getEnergy=True)
    t0 = time.perf_counter()
    integ.step(a.steps)
    ctx.getState(getEnergy=True)  # waits for the GPU
    speed[(engine, prec)] = a.steps * 0.004e-3 / (time.perf_counter() - t0) * 86400
    note = " (run.sh --faster)" if (engine, prec) == ("faster", "mixed") else \
           " (run.sh, production)" if engine == "classic" else " (for information)"
    print(f"{engine:8s} {prec:7s} {speed[(engine, prec)]:9.1f} ns/day{note}", flush=True)
    del ctx, integ
base = speed[("classic", "single")]
print(f"speedup run.sh --faster vs classic: {speed[('faster', 'mixed')] / base:.2f}x")

# ---------------------------------------------------------------- 3. stability
if a.stability_ps > 0:
    nstep = int(a.stability_ps / 0.004)
    every = max(1, nstep // 10)
    print(f"\n## 3. Stability: {a.stability_ps:g} ps NPT from the same state", flush=True)
    for engine, prec in [("classic", "single"), ("faster", "mixed")]:
        ctx, integ, s = sim(engine, prec)
        nd = 3 * sum(1 for i in range(s.getNumParticles()) if s.getParticleMass(i).value_in_unit(u.dalton) > 0) \
            - s.getNumConstraints() - 3
        rows = []
        for _ in range(nstep // every):
            integ.step(every)
            st = ctx.getState(getEnergy=True)
            rows.append((st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole),
                         2 * st.getKineticEnergy().value_in_unit(u.kilojoule_per_mole) / (nd * 0.00831446261815324),
                         st.getPeriodicBoxVolume().value_in_unit(u.nanometer ** 3)))
        r = np.array(rows)
        finite = np.isfinite(ctx.getState(getPositions=True).getPositions(asNumpy=True)
                             .value_in_unit(u.nanometer)).all()
        print(f"{engine:8s} T = {r[:, 1].mean():6.1f} +- {r[:, 1].std():.1f} K | E_pot mean {r[:, 0].mean():.0f} kJ/mol"
              f" | V mean {r[:, 2].mean():.2f} nm^3 | finite = {finite}", flush=True)
        del ctx, integ
print("DONE")
