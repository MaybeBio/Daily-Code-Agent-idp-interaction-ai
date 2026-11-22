#!/usr/bin/env python
"""
run_calvados_cal.py — CALVADOS-2 run with the capped hydration forcing (saturating_forcing.py) built into the energy.
Written 2026-10-24 for the capped ER-NTD runs (BENDER-MT/code/make_era_calibrated_inputs.py). The run uses the plain
residues file; the whole forcing (TEA and hydrophilic terms, both capped) comes from one added force:
    E = -eps/2 sum_i w_i c_sat (1 - exp(-c_i / c_sat)),   c_i = -sum_{j, |i-j| >= 2} [U_AH(r_ij; 1) - U_AH(r_ij; 0)] / eps
c_i is residue i's contact number in full contacts with the Ashbaugh-Hatch shape of the run itself (sigma_ij = mean of
the run's per-particle s, shifted at the run's cutoff), and w_i = gamma_TEA w_TEA + gamma_h w_hyd at the run temperature
(dimensionless, lambda units), read per residue letter from --cap (header w_lambda and a '# c_sat X' line). A
CustomGBForce: the per-residue sum is a computed value, the cap a single-particle energy term; same force group,
cutoff, periodicity and exclusions as the run's Ashbaugh-Hatch force.
Everything else is the run's own run.py (written by GAIVI_PAGE4/calvados/prepare_page4.py), unchanged: this wrapper adds
the force when CALVADOS creates its Simulation, then runs run.py.
    python run_calvados_cal.py --cap w_S1cap_277K.tsv --path RUNDIR
    python run_calvados_cal.py --cap w_S1cap_277K.tsv --path RUNDIR --inspect   (list the forces, add and check the capped
                                                                                 force, stop before simulating)
The Ashbaugh-Hatch force is found by its form, not its name (find_ah, 2026-10-24: the GAIVI CALVADOS version does not
name it AH), and the reading of its parameters is checked against the force's own energy before anything runs.
"""
import glob
import os
import re
import sys

FORCE_NAME = 'BENDER capped hydration forcing'


def read_cap(path):
    lines = open(path).read().splitlines()
    head = next(x for x in lines if x and not x.startswith('#'))
    if head.split('\t')[1] != 'w_lambda':
        raise ValueError(f'{path} is not a CALVADOS capped weight file')
    csat = float(next(x for x in lines if x.startswith('# c_sat')).split()[2])
    w = {}
    for x in lines:
        if x and not x.startswith('#') and x != head:
            a, v = x.split('\t')[:2]
            w[a] = float(v)
    return w, csat


def run_sequence(path):
    fa = [f for f in glob.glob(os.path.join(path, '*.fasta'))]
    if len(fa) != 1:
        raise RuntimeError(f'expected one FASTA in {path}')
    return ''.join(x.strip() for x in open(fa[0]) if not x.startswith('>'))


def describe_forces(system):
    """One line per force: class, name, group, per-particle and global parameters, cutoff, start of the expression."""
    import openmm as mm
    out = []
    for k, f in enumerate(system.getForces()):
        line = f'  [{k}] {type(f).__name__} name={f.getName()!r} group={f.getForceGroup()}'
        if isinstance(f, mm.CustomNonbondedForce):
            line += (f' params={[f.getPerParticleParameterName(i) for i in range(f.getNumPerParticleParameters())]}'
                     f' globals={[(f.getGlobalParameterName(i), f.getGlobalParameterDefaultValue(i)) for i in range(f.getNumGlobalParameters())]}'
                     f' cutoff={f.getCutoffDistance()} exclusions={f.getNumExclusions()} expr={f.getEnergyFunction()[:200]!r}')
        out.append(line)
    return '\n'.join(out)


SIGMA_NAMES, LAMBDA_NAMES = ('s', 'sigma', 'sig'), ('l', 'lambda', 'lam', 'lambdas')


def find_ah(system):
    """The Ashbaugh-Hatch force, by its form rather than its name (CALVADOS versions differ: the GAIVI version does not
    name it AH, 2026-10-24): a CustomNonbondedForce with a sigma-like and a lambda-like per-particle parameter whose
    expression switches at 2^(1/6) sigma. Returns the force, the parameter indices of sigma and lambda, and epsilon."""
    import openmm as mm
    cands = []
    for f in system.getForces():
        if not isinstance(f, mm.CustomNonbondedForce):
            continue
        names = [f.getPerParticleParameterName(i) for i in range(f.getNumPerParticleParameters())]
        si = next((names.index(p) for p in SIGMA_NAMES if p in names), None)
        li = next((names.index(p) for p in LAMBDA_NAMES if p in names), None)
        if si is not None and li is not None and '2^(1/6)' in f.getEnergyFunction().replace(' ', ''):
            cands.append((f, si, li))
    if len(cands) != 1:
        raise RuntimeError(f'expected one Ashbaugh-Hatch force, found {len(cands)}; forces:\n' + describe_forces(system))
    f, si, li = cands[0]
    expr = f.getEnergyFunction().replace(' ', '')
    glob_ = {f.getGlobalParameterName(i): f.getGlobalParameterDefaultValue(i) for i in range(f.getNumGlobalParameters())}
    eps = glob_.get('eps')
    if eps is None:
        m = re.match(r'([0-9.eE+-]+)\*', expr) or re.search(r'(?:^|;)eps=([0-9.eE+-]+)', expr)
        eps = float(m.group(1)) if m else None
    if eps is None or abs(eps - 0.2 * 4.184) > 1e-6:
        raise RuntimeError(f'could not read the Ashbaugh-Hatch epsilon (got {eps}); forces:\n' + describe_forces(system))
    return f, si, li, float(eps)


def check_ah_reading(system, ah, si, li, eps, rc, excl):
    """Recompute the run's own Ashbaugh-Hatch energy on a test chain from the parameters as read (sigma and lambda
    mixed by arithmetic means, shift at the cutoff, bonded pairs excluded) and require agreement with OpenMM's value of
    that force alone; this confirms the reading that the capped force relies on, whatever the CALVADOS version."""
    import numpy as np
    import openmm as mm
    n = system.getNumParticles()
    test = mm.System()
    for i in range(n):
        test.addParticle(system.getParticleMass(i))
    test.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    test.addForce(mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(ah)))
    rng = np.random.default_rng(7)                 # self-avoiding chain: bonds 0.38 nm, no two beads closer than 0.5 nm
    x = np.zeros((n, 3))
    for k in range(1, n):
        for _ in range(1000):
            v = rng.normal(size=3)
            t = x[k - 1] + 0.38 * v / np.linalg.norm(v)
            if k < 2 or np.min(np.linalg.norm(x[:k - 1] - t, axis=1)) > 0.5:
                break
        x[k] = t
    x = x - x.mean(0) + 0.5 * np.array([system.getDefaultPeriodicBoxVectors()[k][k].value_in_unit(mm.unit.nanometer)
                                        for k in range(3)])
    ctx = mm.Context(test, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    ctx.setPositions(x)
    e_omm = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
    p = np.array([[float(v) for v in ah.getParticleParameters(i)] for i in range(n)])
    s_i, l_i = p[:, si], p[:, li]
    i, j = np.triu_indices(n, k=1)
    keep = np.array([(a, b) not in excl for a, b in zip(i, j)])
    i, j = i[keep], j[keep]
    r = np.linalg.norm(x[i] - x[j], axis=1)
    s, l = 0.5 * (s_i[i] + s_i[j]), 0.5 * (l_i[i] + l_i[j])
    lj = (s / r) ** 12 - (s / r) ** 6
    shift = (s / rc) ** 12 - (s / rc) ** 6
    v = np.where(r < 2 ** (1 / 6) * s, 4 * eps * (lj - l * shift) + eps * (1 - l), 4 * eps * l * (lj - shift))
    e_np = float(np.where(r < rc, v, 0.0).sum())
    if abs(e_omm - e_np) > 1e-4 * max(1.0, abs(e_omm)):
        raise RuntimeError(f'Ashbaugh-Hatch energy as read ({e_np:.6f}) differs from the force ({e_omm:.6f}); forces:\n'
                           + describe_forces(system))
    return e_omm, e_np


def add_capped_forcing_calvados(system, w, csat, seq):
    """Add E to a CALVADOS system (idempotent). w: {one-letter: lambda units}; seq: the run sequence."""
    import openmm as mm
    from openmm import unit
    for f in system.getForces():
        if f.getName() == FORCE_NAME:
            return f
    ah, si, li, eps = find_ah(system)
    n = system.getNumParticles()
    if n != len(seq):
        raise RuntimeError(f'{n} particles but the run sequence has {len(seq)} residues')
    s_i = [float(ah.getParticleParameters(i)[si]) for i in range(n)]
    rc = ah.getCutoffDistance()
    rc = rc.value_in_unit(unit.nanometer) if hasattr(rc, 'value_in_unit') else float(rc)
    excl = {tuple(sorted(ah.getExclusionParticles(k))) for k in range(ah.getNumExclusions())}
    if excl != {(i, i + 1) for i in range(n - 1)}:
        raise RuntimeError('Ashbaugh-Hatch exclusions are not exactly the bonded (i, i+1) pairs; forces:\n'
                           + describe_forces(system))
    e_omm, e_np = check_ah_reading(system, ah, si, li, eps, rc, excl)
    f = mm.CustomGBForce()
    f.setName(FORCE_NAME)
    f.addPerParticleParameter('sg')
    f.addPerParticleParameter('w')
    f.addComputedValue('c', f'step({rc:.10g}-r)*select(step(r-2^(1/6)*s), -4*(lj-shift), 1+4*shift);'
                            f'lj=(s/r)^12-(s/r)^6; shift=(s/{rc:.10g})^12-(s/{rc:.10g})^6; s=0.5*(sg1+sg2)',
                       mm.CustomGBForce.ParticlePair)
    f.addEnergyTerm(f'-0.5*{eps:.10g}*w*{csat:.10g}*(1-exp(-c/{csat:.10g}))', mm.CustomGBForce.SingleParticle)
    for i in range(n):
        f.addParticle([s_i[i], w[seq[i]]])
    for a, b in sorted(excl):
        f.addExclusion(a, b)
    f.setNonbondedMethod(mm.CustomGBForce.CutoffPeriodic if ah.usesPeriodicBoundaryConditions()
                         else mm.CustomGBForce.CutoffNonPeriodic)
    f.setCutoffDistance(rc)
    f.setForceGroup(ah.getForceGroup())                    # same group and cutoff as AH: one cutoff per group
    system.addForce(f)
    print(f'  {FORCE_NAME} (c_sat {csat:g}): {n} beads, cutoff {rc} nm (group {f.getForceGroup()}), w range '
          f'{min(w[a] for a in seq):+.4f} to {max(w[a] for a in seq):+.4f}; AH force {ah.getName()!r} read and checked '
          f'({e_omm:.4f} vs {e_np:.4f} kJ/mol on a test chain)', flush=True)
    return f


def main():
    from openmm import app
    argv = sys.argv[1:]
    if '--cap' not in argv or '--path' not in argv:
        raise SystemExit('usage: run_calvados_cal.py --cap FILE --path RUNDIR [--inspect]')
    w, csat = read_cap(argv[argv.index('--cap') + 1])
    path = argv[argv.index('--path') + 1]
    inspect = '--inspect' in argv
    seq = run_sequence(path)
    orig = app.Simulation.__init__
    calls = []

    def init(self, topology, system, *a, **kw):
        if inspect:
            print(describe_forces(system), flush=True)
        add_capped_forcing_calvados(system, w, csat, seq)
        calls.append(1)
        if inspect:
            print('inspect: the capped force can be added to this run; no simulation started', flush=True)
            raise SystemExit(0)
        orig(self, topology, system, *a, **kw)
    app.Simulation.__init__ = init
    import runpy
    sys.argv = [os.path.join(path, 'run.py'), '--path', path]
    runpy.run_path(os.path.join(path, 'run.py'), run_name='__main__')
    if not calls:
        raise RuntimeError('run.py created no Simulation; the forcing was not applied')


if __name__ == '__main__':
    main()
