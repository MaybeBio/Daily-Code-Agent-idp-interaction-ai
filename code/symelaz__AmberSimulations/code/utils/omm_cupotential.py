"""
omm_cupotential.py

The --faster engine: NVIDIA cuPotential (https://github.com/NVIDIA/cuPotential)
in place of OpenMM's own nonbonded force, LangevinMiddle integrator and
isotropic Monte Carlo barostat. Everything else (bonded terms, CMAP,
restraints, constraints, HMR, reporters) stays stock OpenMM.

openmm_run.py builds the System exactly as for the classic engine and then
hands it to faster(), which swaps the pieces and returns the integrator. The
physics is the same model: the same cutoff, Ewald tolerance target,
dispersion correction, 1-4 exceptions, temperature, friction, pressure and
barostat frequency. cuPotential evaluates the long-range electrostatics with
its own splitting/spreading functions instead of standard PME; on our systems
its force error vs a converged reference was ~half that of stock PME at the
same settings (see README, "Faster engine").

Anything the plugin does not support, or that has not been tested in this
pipeline, is REFUSED with an explanation instead of being run differently
from what the .inp asks for.

Setup and compatibility check: Installation_readme.md, "Faster engine
(cuPotential)", utils/install_faster.sh and utils/check_faster.py.
"""
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

PRECISION = 'mixed'   # cuPotential's validated configuration (NVIDIA's benchmarks, our tests)


def _import_plugin():
    try:
        import nv_openmm_plugin as nv
        return nv
    except ImportError as e:
        raise SystemExit(
            "--faster: the cuPotential OpenMM plugin (nv_openmm_plugin) is not installed in this "
            "Python environment (%s).\nInstall it with utils/install_faster.sh and run with "
            "`conda activate openmm_faster`; see Installation_readme.md, 'Faster engine'." % e)


def _refuse(msg):
    raise SystemExit("--faster: %s\nRun without --faster for this input (classic engine)." % msg)


def validate(inputs, fftype):
    """Refuse, before anything is built, what the plugin cannot do or what was
    never tested in this pipeline. Called by openmm_run.py right after the .inp is read."""
    _import_plugin()
    if fftype.upper() != 'AMBER':
        _refuse("only tested with -ff amber (CHARMM force-switch uses a CustomNonbondedForce).")
    if inputs.coulomb != PME:
        _refuse("coulomb must be PME (the plugin replaces PME electrostatics).")
    if inputs.vdw not in ('CutoffPeriodic',):
        _refuse("vdw = %s is not supported or not tested; use CutoffPeriodic." % inputs.vdw)
    if inputs.implicitSolvent:
        _refuse("implicit solvent is not supported (periodic explicit solvent only).")
    if inputs.annealing == 'yes':
        _refuse("simulated annealing is not tested with the cuPotential integrator.")
    if inputs.pcouple == 'yes' and inputs.p_type != 'isotropic':
        _refuse("p_type = %s; only the isotropic MonteCarloBarostat is supported." % inputs.p_type)
    if inputs.e14scale != 1.0:
        _refuse("e14scale != 1 is not tested with the cuPotential force.")


def faster(system, inputs, fftype):
    """Swap cuPotential into an already built System. Returns the integrator."""
    nv = _import_plugin()
    nonbonded = [i for i, f in enumerate(system.getForces()) if isinstance(f, NonbondedForce)]
    if len(nonbonded) != 1:
        _refuse("expected exactly one NonbondedForce, found %d." % len(nonbonded))
    if any(isinstance(f, CustomNonbondedForce) for f in system.getForces()):
        _refuse("the system has a CustomNonbondedForce (not tested).")

    # ---- nonbonded force: always (works with stock and cuPotential integrators) ----
    i = nonbonded[0]
    src = system.getForce(i)
    system.addForce(nv.NonbondedForce(src))   # copies charges, LJ, exceptions, cutoff,
    system.removeForce(i)                     # Ewald tolerance and dispersion correction
    print("--faster: NonbondedForce -> cuPotential NonbondedForce "
          "(cutoff %s, ewaldErrorTolerance %s, dispersion correction %s)"
          % (src.getCutoffDistance(), src.getEwaldErrorTolerance(), src.getUseDispersionCorrection()))

    # ---- integrator + barostat: cuPotential provides LangevinMiddle only ----
    if inputs.integrator == 'langevin_middle':
        for k, f in enumerate(system.getForces()):
            if isinstance(f, CMMotionRemover):    # the cuPotential integrator removes
                system.removeForce(k)             # centre-of-mass motion itself
                break
        for k, f in enumerate(system.getForces()):
            if isinstance(f, MonteCarloBarostat):  # same parameters the classic path used
                p, t, freq = f.getDefaultPressure(), f.getDefaultTemperature(), f.getFrequency()
                system.removeForce(k)
                system.addForce(nv.MonteCarloBarostat(p, t, frequency=freq))
                print("--faster: MonteCarloBarostat -> cuPotential MonteCarloBarostat "
                      "(%s, %s, every %d steps)" % (p, t, freq))
                break
        integrator = nv.LangevinMiddleIntegrator(inputs.temp*kelvin, inputs.fric_coeff/picosecond,
                                                 inputs.dt*picoseconds)
        print("--faster: integrator -> cuPotential LangevinMiddleIntegrator")
    else:
        # e.g. the 1 fs restrained equilibration (integrator = langevin): cuPotential has no
        # plain Langevin integrator, so only the nonbonded force is accelerated here.
        integrator = LangevinIntegrator(inputs.temp*kelvin, inputs.fric_coeff/picosecond,
                                        inputs.dt*picoseconds)
        print("--faster: integrator = langevin -> stock LangevinIntegrator "
              "(cuPotential nonbonded force only)")
    return integrator
