#!/usr/bin/env python
"""
run_openmpipi_cal.py — Mpipi-Recharged run with the calibrated hydration forcing built into the energy, so the ensemble
at T is the calibrated model's own ensemble rather than a reweighted plain one. Written 2026-10-23 for the direct
calibrated PNt runs (BENDER-MT/code/make_pnt_calibrated_inputs.py).

The forcing is the energy the reweighting adds (tea_calibrate.recharged_class_energies, heldout_analysis.bias):
    E_cal(x; T) = sum_{i<j, |i-j|>=2} (delta_i(T) + delta_j(T)) / 2 * s(r_ij; sigma_ij, mu_ij)
with s the attractive part of the Wang-Frenkel shape (-1 inside the minimum r_min, the Wang-Frenkel curve from r_min to
rc = 3 sigma, 0 beyond; nu = 1, sigma and mu from the run's own Wang-Frenkel table) and delta_i(T) in kJ/mol per
residue type, read from --delta (written by make_pnt_calibrated_inputs.py with the code of the reweighting).
Minimum image in the run's periodic box; the bonded (i, i+1) pairs are excluded, as in the Wang-Frenkel force.

Everything else is GAIVI_PAGE4/mpipi/run_openmpipi_page4.py, unchanged: this wrapper adds the force to the System when
that script creates its Simulation, then runs it with the remaining arguments.
    python run_openmpipi_cal.py --delta delta_A0_283K.tsv <run_openmpipi_page4.py arguments>
--delta may also name a capped weight file (header w_kJmol, '# c_sat X'), which adds the capped forcing of
saturating_forcing.py instead (add_capped_forcing, 2026-10-24).
"""
import os
import sys

import numpy as np

FORCE_NAME = 'BENDER calibrated hydration forcing'
N0 = 21                                  # OpenMpipi type table size
ORDER = 'MGKTRADEYVLQWFSHNPCI'           # OpenMpipi mapping_dict indices 0..19 (20 = RNA U)


def read_delta(path):
    d = {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        a, v = line.split('\t')[:2]
        if a == 'residue':
            continue
        d[a] = float(v)
    if set(d) != set(ORDER):
        raise ValueError(f'{path} must list all 20 residues')
    return d


def wf_tables():
    """sigma and mu tables from OpenMpipi's recharged_params.txt, [y, x] = [type_j, type_i]."""
    import OpenMpipi
    P = np.loadtxt(os.path.join(os.path.dirname(OpenMpipi.__file__), 'data_files', 'recharged_params.txt'))
    w = P[:-N0 * N0].reshape(3, N0, N0)
    return w, w[1], np.floor(w[2])


def _wf_setup(system):
    """The run's Wang-Frenkel force and what both forcings need from it: type table size, per-particle types, the amino
    acid each type copies, sigma and mu tables, the residue letters, the bonded exclusions and the largest 3 sigma."""
    import openmm as mm
    wf = [f for f in system.getForces() if isinstance(f, mm.CustomNonbondedForce) and f.getNumTabulatedFunctions()
          and f.getTabulatedFunctionName(0) == 'wf_table']
    if len(wf) != 1:
        raise RuntimeError('expected one Wang-Frenkel force with a wf_table')
    wf = wf[0]
    w, _, _ = wf_tables()
    xs, ys, zs, vals = wf.getTabulatedFunction(0).getFunctionParameters()
    nt = xs
    if (ys, zs) != (nt, 3) or nt < N0:
        raise RuntimeError('unexpected Wang-Frenkel table shape')
    wt = np.asarray(vals).reshape(3, nt, nt)
    if not np.allclose(wt[:, :N0, :N0], w):
        raise RuntimeError('the Wang-Frenkel table differs from recharged_params.txt')
    # phosphate-like types added by apply_phospho (2026-10-24): each copies one amino acid's Wang-Frenkel row (Glu for
    # pSer/pThr, Tyr for pTyr) and takes that residue's hydration weight, as the reweighting reads the run sequence
    base = list(range(N0))
    for t in range(N0, nt):
        match = [b for b in range(len(ORDER)) if np.allclose(wt[:, t, :N0], wt[:, b, :N0])
                 and np.allclose(wt[:, :N0, t], wt[:, :N0, b])]
        if len(match) != 1:
            raise RuntimeError(f'type {t} does not copy exactly one amino-acid Wang-Frenkel row')
        base.append(match[0])
    sig, mu = wt[1], np.floor(wt[2])
    n = system.getNumParticles()
    types = [int(wf.getParticleParameters(i)[0]) for i in range(n)]
    if max(base[t] for t in types) >= len(ORDER):
        raise RuntimeError('a particle is not one of the 20 amino-acid types')
    seq = ''.join(ORDER[base[t]] for t in types)
    excl = {tuple(sorted(wf.getExclusionParticles(k))) for k in range(wf.getNumExclusions())}
    if excl != {(i, i + 1) for i in range(n - 1)}:
        raise RuntimeError('Wang-Frenkel exclusions are not exactly the bonded (i, i+1) pairs')
    present = sorted(set(types))
    rc_max = 3.0 * max(sig[a, b] for a in present for b in present)
    return {'wf': wf, 'nt': nt, 'sig': sig, 'mu': mu, 'n': n, 'types': types, 'seq': seq, 'excl': excl, 'rc_max': rc_max}


def _place(system, f, info):
    """Periodicity, exclusions, cutoff and force group of a new force. OpenMM needs one cutoff per force group. The
    expressions are already 0 beyond each pair's 3 sigma, so any cutoff >= rc_max gives the same energy: share the
    Wang-Frenkel group and cutoff when that is long enough, otherwise use the shortest exact cutoff in a force group of
    its own (fixed 2026-10-24: the first runs failed with "All Forces in a single force group must use the same cutoff
    distance")."""
    import openmm as mm
    from openmm import unit
    wf = info['wf']
    for a, b in sorted(info['excl']):
        f.addExclusion(a, b)
    periodic = wf.usesPeriodicBoundaryConditions()
    if isinstance(f, mm.CustomGBForce):
        f.setNonbondedMethod(mm.CustomGBForce.CutoffPeriodic if periodic else mm.CustomGBForce.CutoffNonPeriodic)
    else:
        f.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic if periodic else mm.CustomNonbondedForce.CutoffNonPeriodic)
        f.setUseSwitchingFunction(False)
        f.setUseLongRangeCorrection(False)
    wf_cut = wf.getCutoffDistance()
    wf_cut = wf_cut.value_in_unit(unit.nanometer) if hasattr(wf_cut, 'value_in_unit') else float(wf_cut)
    if wf_cut >= info['rc_max']:
        f.setCutoffDistance(wf_cut)
        f.setForceGroup(wf.getForceGroup())
    else:
        used = {g.getForceGroup() for g in system.getForces()}
        free = [g for g in range(32) if g not in used]
        if not free:
            raise RuntimeError('no free OpenMM force group for the hydration forcing')
        f.setCutoffDistance(info['rc_max'] * 1.0001)
        f.setForceGroup(free[0])
    system.addForce(f)


WF_SHAPE = ('select(step(r-rmin), select(step(r-rc), 0, alpha*((s/r)^(2*m)-1)*((rc/r)^(2*m)-1)^2), -1);'
            'rmin=rc*(3/(1+2*(rc/s)^(2*m)))^(1/(2*m));'
            'alpha=2*3^(2*m)*(3/(2*(3^(2*m)-1)))^3;'
            'rc=3*s; s=sig_t(t1, t2); m=mu_t(t1, t2)')


def add_hydration_forcing(system, delta):
    """Add E_cal to system (idempotent). delta: {one-letter: kJ/mol}. Returns the new force."""
    import openmm as mm
    for f in system.getForces():
        if f.getName() == FORCE_NAME:
            return f
    info = _wf_setup(system)
    f = mm.CustomNonbondedForce('0.5*(d1+d2)*' + WF_SHAPE)
    f.setName(FORCE_NAME)
    f.addPerParticleParameter('t')
    f.addPerParticleParameter('d')
    f.addTabulatedFunction('sig_t', mm.Discrete2DFunction(info['nt'], info['nt'], info['sig'].ravel()))
    f.addTabulatedFunction('mu_t', mm.Discrete2DFunction(info['nt'], info['nt'], info['mu'].ravel()))
    for i in range(info['n']):
        f.addParticle([info['types'][i], delta[info['seq'][i]]])
    _place(system, f, info)
    seq = info['seq']
    print(f'  {FORCE_NAME}: {info["n"]} beads, cutoff {f.getCutoffDistance()} (group {f.getForceGroup()}), delta range '
          f'{min(delta[a] for a in seq):+.3f} to {max(delta[a] for a in seq):+.3f} kJ/mol', flush=True)
    return f


def add_capped_forcing(system, w, csat):
    """Add the capped forcing (saturating_forcing.py, 2026-10-24), idempotent:
        E = -1/2 sum_i w_i c_sat (1 - exp(-c_i / c_sat)),   c_i = -sum_{j, |i-j| >= 2} s(r_ij)
    with s the attractive Wang-Frenkel shape of the linear forcing, so c_i is residue i's contact number in full
    contacts; w: {one-letter: kJ/mol} (gamma_TEA w_TEA + gamma_h w_hyd at the run temperature). A CustomGBForce: the
    per-residue sum is a computed value, the cap a single-particle energy term."""
    import openmm as mm
    for f in system.getForces():
        if f.getName() == FORCE_NAME:
            return f
    info = _wf_setup(system)
    f = mm.CustomGBForce()
    f.setName(FORCE_NAME)
    f.addPerParticleParameter('t')
    f.addPerParticleParameter('w')
    f.addTabulatedFunction('sig_t', mm.Discrete2DFunction(info['nt'], info['nt'], info['sig'].ravel()))
    f.addTabulatedFunction('mu_t', mm.Discrete2DFunction(info['nt'], info['nt'], info['mu'].ravel()))
    f.addComputedValue('c', '-(' + WF_SHAPE.split(';', 1)[0] + ');' + WF_SHAPE.split(';', 1)[1],
                       mm.CustomGBForce.ParticlePair)
    f.addEnergyTerm(f'-0.5*w*{csat:.10g}*(1-exp(-c/{csat:.10g}))', mm.CustomGBForce.SingleParticle)
    for i in range(info['n']):
        f.addParticle([info['types'][i], w[info['seq'][i]]])
    _place(system, f, info)
    seq = info['seq']
    print(f'  {FORCE_NAME} (capped, c_sat {csat:g}): {info["n"]} beads, cutoff {f.getCutoffDistance()} '
          f'(group {f.getForceGroup()}), w range {min(w[a] for a in seq):+.3f} to {max(w[a] for a in seq):+.3f} kJ/mol',
          flush=True)
    return f


def read_forcing(path):
    """('linear', delta) from a delta file, or ('capped', (w, c_sat)) from a capped weight file (header w_kJmol and a
    '# c_sat X' line)."""
    lines = open(path).read().splitlines()
    head = next(x for x in lines if x and not x.startswith('#'))
    if head.split('\t')[1] == 'w_kJmol':
        csat = float(next(x for x in lines if x.startswith('# c_sat')).split()[2])
        d = {}
        for x in lines:
            if x and not x.startswith('#') and x != head:
                a, v = x.split('\t')[:2]
                d[a] = float(v)
        if set(d) != set(ORDER):
            raise ValueError(f'{path} must list all 20 residues')
        return 'capped', (d, csat)
    return 'linear', read_delta(path)


def main():
    from openmm import app
    argv = sys.argv[1:]
    if '--delta' not in argv:
        raise SystemExit('usage: run_openmpipi_cal.py --delta FILE <run_openmpipi_page4.py arguments>')
    k = argv.index('--delta')
    kind, forcing = read_forcing(argv[k + 1])
    rest = argv[:k] + argv[k + 2:]
    page4 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_openmpipi_page4.py')
    orig = app.Simulation.__init__
    calls = []

    def init(self, topology, system, *a, **kw):
        add_hydration_forcing(system, forcing) if kind == 'linear' else add_capped_forcing(system, *forcing)
        calls.append(1)
        orig(self, topology, system, *a, **kw)
    app.Simulation.__init__ = init
    import runpy
    sys.argv = [page4] + rest
    runpy.run_path(page4, run_name='__main__')
    if not calls:
        raise RuntimeError('run_openmpipi_page4.py created no Simulation; the forcing was not applied')


if __name__ == '__main__':
    main()
