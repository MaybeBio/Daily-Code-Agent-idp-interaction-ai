"""
check_faster.py [parm7 rst7|state.xml]

Is this machine + environment ready for `run.sh --faster` (NVIDIA cuPotential)?
Run it on a GPU node, inside the environment made by utils/install_faster.sh.

  1. versions: OpenMM, cuPotential plugin, the runtime's CUDA, the driver's CUDA,
     the NVRTC OpenMM will use (must not be newer than the driver)
  2. GPU: compute capability covered by the cuPotential runtime?
  3. accuracy: classic and faster vs a converged reference (stock OpenMM, double
     precision, ewaldErrorTolerance 1e-6, same LJ model). Both approximate the
     same electrostatics, so they are NOT compared with each other: the faster
     engine must be at least as accurate as the classic one (within 1.5x, plus
     a small floor), which is what we measured on our systems (~2x better).
  4. 1000 steps of NPT at 4 fs with the cuPotential integrator and barostat:
     finite, temperature near target

Without arguments it uses a small generated water box; give a built system
(step3_input.parm7 + step3_input.rst7, or an equilibrated .rst XML state) to
test the real thing. Exit code 0 = all PASS.
"""
import ctypes
import glob
import json
import os
import sys

import numpy as np

results = []


def report(name, ok, detail, warn=False):
    tag = "PASS" if ok else ("WARN" if warn else "FAIL")
    print(f"{tag:4s}  {name:44s} {detail}", flush=True)
    if not ok and not warn:
        results.append(name)


def cuda_ver(v):
    return f"{v // 1000}.{(v % 1000) // 10}"


# ---------------------------------------------------------------- versions
import openmm as mm
import openmm.app as app
import openmm.unit as u

report("OpenMM >= 8.4", tuple(map(int, mm.__version__.split('.')[:2])) >= (8, 4), mm.__version__)
try:
    import nv_openmm_plugin as nv
    pkg = os.path.dirname(nv.__file__)
    report("cuPotential plugin importable", True, pkg)
except ImportError as e:
    report("cuPotential plugin importable", False, f"{e} -> run utils/install_faster.sh")
    sys.exit(1)

runtime = {}
for f in glob.glob(os.path.join(pkg, "**", "*.json"), recursive=True):
    try:
        d = json.load(open(f))
        if isinstance(d, dict) and "cuda_min" in d:
            runtime = d
            break
    except Exception:
        pass

try:
    drv = ctypes.c_int()
    ctypes.CDLL("libcuda.so.1").cuDriverGetVersion(ctypes.byref(drv))
    driver = drv.value
    report("GPU driver present", driver > 0, f"supports CUDA {cuda_ver(driver)}")
except OSError as e:
    report("GPU driver present", False, f"{e} (not a GPU node?)")
    sys.exit(1)

if runtime:
    need = tuple(map(int, runtime["cuda_min"].split('.')))
    have = (driver // 1000, (driver % 1000) // 10)
    report("driver vs cuPotential runtime CUDA", have[0] == need[0],
           f"runtime built for CUDA {runtime['cuda_min']}-{runtime.get('cuda_max', '?')}, driver {cuda_ver(driver)}"
           + (" (older minor: OK via minor-version compatibility)" if have < need and have[0] == need[0] else ""))

try:
    nvrtc = ctypes.CDLL("libnvrtc.so.13")
    a, b = ctypes.c_int(), ctypes.c_int()
    nvrtc.nvrtcVersion(ctypes.byref(a), ctypes.byref(b))
    ok = (a.value, b.value) <= (driver // 1000, (driver % 1000) // 10)
    report("NVRTC (OpenMM kernel compiler) <= driver", ok,
           f"NVRTC {a.value}.{b.value}" + ("" if ok else
           " is newer than the driver -> CUDA_ERROR_UNSUPPORTED_PTX_VERSION; re-run install_faster.sh "
           "(it installs an NVRTC matching the driver) and re-activate the env"))
except OSError:
    report("NVRTC (OpenMM kernel compiler) <= driver", True, "not preloaded; checked by the Context test below", warn=True)

# ---------------------------------------------------------------- GPU
plat = mm.Platform.getPlatformByName("CUDA")
dev = os.environ.get("CUDA_VISIBLE_DEVICES", "0").split(",")[0]
try:
    import subprocess
    out = subprocess.run(["nvidia-smi", "--query-gpu=name,compute_cap", "--format=csv,noheader", "-i", dev],
                         capture_output=True, text=True).stdout.strip()
    name, cc = [x.strip() for x in out.split(",")]
    archs = runtime.get("compact_architectures", [])
    report("GPU architecture supported", not archs or int(cc.replace('.', '')) in archs,
           f"{name}, compute capability {cc}" + (f" (runtime: {archs})" if archs else ""))
except Exception as e:
    report("GPU architecture supported", True, f"could not query nvidia-smi ({e})", warn=True)

# ---------------------------------------------------------------- system
if len(sys.argv) >= 3:
    prmtop = app.AmberPrmtopFile(sys.argv[1])
    topology = prmtop.topology
    if sys.argv[2].endswith(".rst7") or sys.argv[2].endswith(".inpcrd"):
        crd = app.AmberInpcrdFile(sys.argv[2])
        positions, box = crd.positions, crd.boxVectors
    else:
        st = mm.XmlSerializer.deserialize(open(sys.argv[2]).read())
        positions, box = st.getPositions(), st.getPeriodicBoxVectors()

    def create():
        return prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.9 * u.nanometer,
                                   constraints=app.HBonds, ewaldErrorTolerance=5e-4)
    label = os.path.basename(sys.argv[1])
else:
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    mod = app.Modeller(app.Topology(), [])
    mod.addSolvent(ff, boxSize=mm.Vec3(4.0, 4.0, 4.0) * u.nanometer)
    topology, positions, box = mod.topology, mod.positions, mod.topology.getPeriodicBoxVectors()

    def create():
        return ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=0.9 * u.nanometer,
                               constraints=app.HBonds, ewaldErrorTolerance=5e-4)
    label = "generated 4 nm TIP3P water box"
print(f"\nsystem: {label}, {topology.getNumAtoms()} atoms")


def to_nv(s):
    i = [k for k, f in enumerate(s.getForces()) if isinstance(f, mm.NonbondedForce)][0]
    s.addForce(nv.NonbondedForce(s.getForce(i)))
    s.removeForce(i)
    return s


def create_ref():
    s = create()
    for f in s.getForces():
        if isinstance(f, mm.NonbondedForce):
            f.setEwaldErrorTolerance(1e-6)
    return s


def energy_forces(s, precision):
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(s, integ, plat, {"Precision": precision})
    ctx.setPeriodicBoxVectors(*box)
    ctx.setPositions(positions)
    st = ctx.getState(getEnergy=True, getForces=True)
    return (st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole),
            np.asarray(st.getForces(asNumpy=True).value_in_unit(u.kilojoule_per_mole / u.nanometer)))


def hint(ex):
    s = str(ex)
    if "PTX" in s:
        return s.splitlines()[0] + " -> NVRTC newer than the driver; re-run install_faster.sh and re-activate"
    return s.splitlines()[0][:150]


try:
    e0, f0 = energy_forces(create(), "single")
    report("classic engine Context", True, f"E_pot {e0:.1f} kJ/mol")
except Exception as ex:
    report("classic engine Context", False, hint(ex))
    sys.exit(1)
try:
    e1, f1 = energy_forces(to_nv(create()), "mixed")
    report("faster engine Context", True, f"E_pot {e1:.1f} kJ/mol")
except Exception as ex:
    report("faster engine Context", False, hint(ex))
    sys.exit(1)

try:
    er, fr = energy_forces(create_ref(), "double")
except Exception as ex:
    report("reference Context (double precision)", False, hint(ex))
    sys.exit(1)
real = np.linalg.norm(fr, axis=1) > 0
frms = np.sqrt(np.mean(np.sum(fr[real] ** 2, axis=1)))
def ferr(f):
    return np.sqrt(np.mean(np.sum((f - fr)[real] ** 2, axis=1))) / frms
de0, de1 = abs(e0 - er) / abs(er), abs(e1 - er) / abs(er)
df0, df1 = ferr(f0), ferr(f1)
report("energy error vs converged reference", de1 <= max(1.5 * de0, 1e-5),
       f"faster {de1:.1e}, classic {de0:.1e} (faster must be <= 1.5x classic, or <= 1e-5)")
report("force error vs converged reference", df1 <= max(1.5 * df0, 2e-4),
       f"faster {df1:.1e}, classic {df0:.1e} (faster must be <= 1.5x classic, or <= 2e-4)")

# ---------------------------------------------------------------- dynamics
T = 303.15
s = to_nv(create())
for k, f in enumerate(s.getForces()):
    if isinstance(f, mm.CMMotionRemover):
        s.removeForce(k)
        break
s.addForce(nv.MonteCarloBarostat(1.0 * u.bar, T * u.kelvin, frequency=25))
for a, b in topology.bonds():  # HMR 3.024 Da on non-water hydrogens, as omm_hmr.py
    if a.residue.name in ('HOH', 'WAT'):
        continue
    if a.element == app.element.hydrogen:
        a, b = b, a
    if b.element == app.element.hydrogen and a.element not in (app.element.hydrogen, None):
        extra = 3.024 * u.dalton - s.getParticleMass(b.index)
        s.setParticleMass(b.index, 3.024 * u.dalton)
        s.setParticleMass(a.index, s.getParticleMass(a.index) - extra)
integ = nv.LangevinMiddleIntegrator(T * u.kelvin, 1.0 / u.picosecond, 0.004 * u.picoseconds)
try:
    ctx = mm.Context(s, integ, plat, {"Precision": "mixed"})
    ctx.setPeriodicBoxVectors(*box)
    ctx.setPositions(positions)
    mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=200)
    ctx.setVelocitiesToTemperature(T * u.kelvin)
    integ.step(1000)
    st = ctx.getState(getEnergy=True, getPositions=True)
    nd = 3 * sum(1 for i in range(s.getNumParticles()) if s.getParticleMass(i).value_in_unit(u.dalton) > 0) \
        - s.getNumConstraints() - 3
    temp = 2 * st.getKineticEnergy().value_in_unit(u.kilojoule_per_mole) / (nd * 0.00831446261815324)
    finite = np.isfinite(st.getPositions(asNumpy=True).value_in_unit(u.nanometer)).all()
    report("faster engine 1000 steps NPT, 4 fs", finite and abs(temp - T) < 15,
           f"T = {temp:.1f} K (target {T}), finite = {finite}")
except Exception as ex:
    report("faster engine 1000 steps NPT, 4 fs", False, hint(ex))

print("\nRESULT:", "READY for run.sh --faster" if not results else "NOT READY: " + ", ".join(results))
sys.exit(1 if results else 0)
