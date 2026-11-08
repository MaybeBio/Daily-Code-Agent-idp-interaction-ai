"""
phospho_recharged_ess.py — phospho_direct_plan.md Amendment 3 (hash 5ee23c3b): a fluctuation measure native to
Mpipi-Recharged, and the secondary analyses M1 to M3. Written and hashed before it was run.

  --verify   the systems built here must serialize identically to those of run_openmpipi_idp.py (snapshot 683290f1) for
             one unmodified (PHOS-O00193-WT-r1) and one phosphorylated entry (PHOS-O00193-LS2-r1); the run script is executed
             up to its Simulation constructor with the same arguments and its System is captured
  (default)  for every site, dU_f = U_P(x_f) - U_W(x_f) on the 1,000-frame rule of bender_common.load_calvados (second half,
             stride n // 1000), evaluated by OpenMM (CPU) on systems built with the run script's code path (OpenMpipi IDP,
             get_mpipi_system at 300 K and 150 mM, then apply_phospho with Yukawa x2 for the phosphorylated site); weights
             w_f proportional to exp(-dU_f / kT); ESS = 1 / sum(w^2) / n; Recharged log ESS of a site = mean over the three
             new unmodified replicates; the stored BENDER Recharged run as a check; reweighted Rg change with the same weights;
             direct Rg change from the 3 + 3 Recharged runs on the same frame rule. Then M1 to M3 of the amendment.

  ~/miniforge3/envs/openmpipi/bin/python phospho_recharged_ess.py --verify
  ~/miniforge3/envs/openmpipi/bin/python phospho_recharged_ess.py   ->  phospho_direct_recharged/{ess_recharged.csv, summary_amend3.txt}
"""
import argparse
import os
import runpy
import sys

import mdtraj as md
import numpy as np
import openmm as mm
import pandas as pd
from openmm import app, unit
from scipy.stats import binomtest, spearmanr

from OpenMpipi import IDP, get_mpipi_system

EV = os.path.dirname(os.path.abspath(__file__))
RUN_SCRIPT = os.path.expanduser((__import__('os').environ['BSS_DATA'] + '/Scripts/OpenMPIPI-runs/run_openmpipi_idp.py'))
FASTA = os.path.join(EV, 'gaivi', 'phospho_direct', 'OpenMPIPI_snapshot_amend2', 'input_phospho', 'phospho_direct_300K.fasta')
RUNS = (__import__('os').environ['BSS_DATA'] + '/BENDER-DYN/phospho_direct_recharged/300K')
STORED = (__import__('os').environ['BSS_DATA'] + '/MPIPI Recharged/300K')
SEL = os.path.join(EV, 'gaivi', 'phospho_direct', 'selection_amend1.csv')
SITES_R = os.path.join(EV, 'phospho_direct_recharged', 'sites_recharged_final.csv')
MIMIC = (__import__('os').environ['BSS_DATA'] + '/BENDER-MT/Traj Level Analysis/Analysis_FDT_RECHARGED/traj_pass/phosphomimetic.csv')
OUT = os.path.join(EV, 'phospho_direct_recharged')
T, CSX, SCALE, K = 300.0, 150.0, 2.0, 1000
KT = 0.0083144626 * T
FORM = {'A': 'annotated', 'U': 'partner', 'L': 'low'}


def read_fasta():
    seqs, name = {}, None
    for line in open(FASTA):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split('|')[0]
            seqs[name] = ''
        elif name:
            seqs[name] += line
    return seqs


def saw_coords(n, bond=0.381, min_dist=0.6, rng=None, max_tries=200):
    """As run_openmpipi_idp.py (the positions do not enter the System; checked by --verify)."""
    rng = rng or np.random.default_rng()
    while True:
        xyz = np.zeros((n, 3))
        ok = True
        for i in range(1, n):
            for _ in range(max_tries):
                v = rng.normal(size=3)
                trial = xyz[i - 1] + bond * v / np.linalg.norm(v)
                if i < 2 or np.min(np.linalg.norm(xyz[:i - 1] - trial, axis=1)) > min_dist:
                    xyz[i] = trial
                    break
            else:
                ok = False
                break
        if ok:
            return xyz - xyz.mean(axis=0)


def apply_phospho(system, sites, seq):
    """Copied unchanged from run_openmpipi_idp.py (683290f1)."""
    for i, _ in sites:
        if seq[i] != 'E':
            raise ValueError(f'residue {i + 1} is {seq[i]}, expected a Glu bead')
    wf = next(f for f in system.getForces() if isinstance(f, mm.CustomNonbondedForce) and f.getNumTabulatedFunctions()
              and f.getTabulatedFunctionName(0) == 'wf_table')
    yu = next(f for f in system.getForces() if isinstance(f, mm.CustomNonbondedForce) and f.getNumTabulatedFunctions()
              and f.getTabulatedFunctionName(0) == 'A_table')
    n0, glu = 21, 7
    nb = np.loadtxt(os.path.join(os.path.dirname(__import__('OpenMpipi').__file__), 'data_files', 'recharged_params.txt'))
    W0, A0 = nb[:-n0 * n0].reshape(3, n0, n0), nb[-n0 * n0:].reshape(n0, n0)
    xs, ys, zs, wv = wf.getTabulatedFunction(0).getFunctionParameters()
    ax, ay, av = yu.getTabulatedFunction(0).getFunctionParameters()
    if (xs, ys, zs, ax, ay) != (n0, n0, 3, n0, n0) or not (np.allclose(wv, W0.ravel()) and np.allclose(av, A0.ravel())):
        raise RuntimeError('OpenMpipi tables differ from recharged_params.txt (21 types expected); check the version')
    if any(wf.getParticleParameters(i)[0] != glu or yu.getParticleParameters(i)[0] != glu for i, _ in sites):
        raise RuntimeError('a phosphosite is not on a Glu-type particle')
    n = n0 + len(sites)
    m = np.r_[np.arange(n0), np.full(len(sites), glu)]
    sc = np.r_[np.ones(n0), [s for _, s in sites]]
    W = W0[:, m][:, :, m]
    A = A0[np.ix_(m, m)] * np.outer(sc, sc)
    wf.getTabulatedFunction(0).setFunctionParameters(n, n, 3, W.ravel())
    yu.getTabulatedFunction(0).setFunctionParameters(n, n, A.ravel())
    for k, (i, s) in enumerate(sites):
        _, g = wf.getParticleParameters(i)
        wf.setParticleParameters(i, [n0 + k, g])
        yu.setParticleParameters(i, [n0 + k])
    return sites


def build(seq, site=None, name='X'):
    """System of the run script for seq (site: 0-based index of a phosphate-like bead, seq[site] must be 'E')."""
    n = len(seq)
    L = int(n / 5) + 20 if n > 300 else 50
    chain = IDP(name, seq)
    top = chain.topology
    pos = saw_coords(n, rng=np.random.default_rng(0)) + L / 2
    box = [mm.Vec3(L, 0, 0), mm.Vec3(0, L, 0), mm.Vec3(0, 0, L)] * unit.nanometer
    top.setPeriodicBoxVectors(box)
    system = get_mpipi_system(pos, top, {name: []}, T, CSX, CM_remover=False, periodic=True)
    system.setDefaultPeriodicBoxVectors(*box)
    if site is not None:
        apply_phospho(system, [(site, SCALE)], seq)
    return system


def run_script_system(name):
    """Execute run_openmpipi_idp.py with --phospho up to the Simulation constructor and capture its System."""
    captured = {}

    class Stop(Exception):
        pass

    class Fake:
        def __init__(self, topology, system, *a, **k):
            captured['system'] = system
            raise Stop()
    real = app.Simulation
    app.Simulation = Fake
    argv = sys.argv
    sys.argv = [RUN_SCRIPT, '--name', name, '--fasta', FASTA, '--temp', '300', '--seed', '7', '--platform', 'CPU',
                '--outdir', os.path.join('/private/tmp', 'amend3_verify'), '--phospho']
    try:
        runpy.run_path(RUN_SCRIPT, run_name='__main__')
    except Stop:
        pass
    finally:
        app.Simulation, sys.argv = real, argv
    return captured['system']


def verify():
    seqs = read_fasta()
    ok = True
    for name, site in (('PHOS-O00193-WT-r1', None), ('PHOS-O00193-LS2-r1', 1)):
        seq = seqs[name].replace('B', 'E').replace('J', 'E')
        a = mm.XmlSerializer.serialize(run_script_system(name))
        b = mm.XmlSerializer.serialize(build(seq, site))
        print(f'{name}: identical System XML: {a == b}')
        ok &= a == b
    if not ok:
        sys.exit('VERIFY FAILED')
    print('VERIFY OK')


def even_indices(n, k):
    if k is None or k >= n:
        return np.arange(n)
    return np.arange(0, n, max(1, n // k))


def frames(dcd, top):
    traj = md.load(dcd, top=top)
    start = traj.n_frames // 2
    return traj.xyz[start:][even_indices(traj.n_frames - start, K)].astype(np.float64)


def rg(X):
    return np.sqrt(((X - X.mean(1, keepdims=True)) ** 2).sum(2).mean(1))


def energies(ctx, X):
    out = np.empty(len(X))
    for f, x in enumerate(X):
        ctx.setPositions(x)
        out[f] = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return out


def context(system):
    return mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'), {'Threads': '2'})


def compute():
    seqs = read_fasta()
    sel = pd.read_csv(SEL)
    rows = []
    for uid in sorted(sel.UniProt_ID.unique()):
        s = sel[sel.UniProt_ID == uid].set_index('role')
        wt = seqs[f'PHOS-{uid}-WT-r1']
        ctxW = context(build(wt))
        names = {f: f'{f}{s.loc[r, "residue"]}{int(s.loc[r, "pos"])}' for f, r in FORM.items()}
        ctxP = {}
        for f, r in FORM.items():
            i = int(s.loc[r, 'pos']) - 1
            assert wt[i] == s.loc[r, 'residue'], (uid, r)
            ctxP[f] = context(build(wt[:i] + 'E' + wt[i + 1:], i))
        ens = [(f'r{k}', os.path.join(RUNS, f'PHOS-{uid}-WT-r{k}', f'PHOS-{uid}-WT-r{k}.dcd'),
                os.path.join(RUNS, f'PHOS-{uid}-WT-r{k}', 'top.pdb')) for k in (1, 2, 3)]
        ens.append(('stored', os.path.join(STORED, uid, f'{uid}.dcd'), os.path.join(STORED, uid, 'top.pdb')))
        res = {f: {} for f in FORM}
        rgW = []
        for tag, dcd, top in ens:
            X = frames(dcd, top)
            r_ = rg(X)
            if tag != 'stored':
                rgW.append(r_.mean())
            uW = energies(ctxW, X)
            for f in FORM:
                du = energies(ctxP[f], X) - uW
                x = -du / KT
                w = np.exp(x - x.max())
                w /= w.sum()
                res[f][tag] = (float(np.log(1.0 / (w ** 2).sum() / len(w))), float(w @ r_ / r_.mean() - 1) * 100)
        rgWm = np.mean(rgW)
        for f, r in FORM.items():
            rgP = np.mean([rg(frames(os.path.join(RUNS, f'PHOS-{uid}-{names[f]}-r{k}', f'PHOS-{uid}-{names[f]}-r{k}.dcd'),
                                     os.path.join(RUNS, f'PHOS-{uid}-{names[f]}-r{k}', 'top.pdb'))).mean() for k in (1, 2, 3)])
            rows.append({'UniProt_ID': uid, 'role': r, 'pos': int(s.loc[r, 'pos']), 'residue': s.loc[r, 'residue'],
                         **{f'log_ess_{t}': res[f][t][0] for t in ('r1', 'r2', 'r3', 'stored')},
                         'log_ess_rech': float(np.mean([res[f][t][0] for t in ('r1', 'r2', 'r3')])),
                         'dRg_rw_rech_pct': float(np.mean([res[f][t][1] for t in ('r1', 'r2', 'r3')])),
                         'dRg_direct_rech_pct': (rgP / rgWm - 1) * 100})
        print(f'{uid} done', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, 'ess_recharged.csv'), index=False)
    return df


def centre(v, g):
    v = pd.Series(np.asarray(v, float))
    return (v - v.groupby(g).transform('mean')).to_numpy()


def one_sided_neg(x, y):
    r = spearmanr(x, y)
    return r.correlation, (r.pvalue / 2 if r.correlation < 0 else 1 - r.pvalue / 2)


def analyse(df):
    sr = pd.read_csv(SITES_R)
    sel = pd.read_csv(SEL)
    m = df.merge(sr, on=['UniProt_ID', 'role']).merge(sel[['UniProt_ID', 'role', 'ess']], on=['UniProt_ID', 'role'])
    m['log_ess_cal'] = np.log(m.ess)
    g = m.UniProt_ID.to_numpy()
    L = ['Amendment 3 (5ee23c3b): Mpipi-Recharged-native fluctuation and M1 to M3', '',
         f'Recharged log ESS: median {m.log_ess_rech.median():+.3f} (CALVADOS-2 {m.log_ess_cal.median():+.3f}); agreement of the '
         f'three replicate values: median pairwise Spearman over sites '
         f'{np.median([spearmanr(m[a], m[b]).correlation for a, b in (("log_ess_r1", "log_ess_r2"), ("log_ess_r1", "log_ess_r3"), ("log_ess_r2", "log_ess_r3"))]):+.3f}; '
         f'new replicates vs stored BENDER run Spearman {spearmanr(m.log_ess_rech, m.log_ess_stored).correlation:+.3f}', '']
    r, p = one_sided_neg(m.log_ess_rech, m.E_occ)
    L.append(f'M1  Recharged E (hub occupancy) vs Recharged log ESS over 60 sites: rho {r:+.3f}, one-sided P {p:.4g} -> '
             f'{"rewiring follows the Recharged fluctuation" if p < 0.05 else "not shown"}')
    r2, p2 = one_sided_neg(m.log_ess_rech, m.E_deg)
    L.append(f'    degree: rho {r2:+.3f}, one-sided P {p2:.4g}; within protein (centred): occupancy '
             f'{spearmanr(centre(m.log_ess_rech, g), centre(m.E_occ, g)).correlation:+.3f}, degree '
             f'{spearmanr(centre(m.log_ess_rech, g), centre(m.E_deg, g)).correlation:+.3f}')
    L.append('')
    piv = m.pivot(index='UniProt_ID', columns='role', values='log_ess_rech')
    for other in ('partner', 'annotated'):
        d = (piv['low'] - piv[other]).to_numpy()
        L.append(f'M2  Recharged log ESS, low minus {other}: median {np.median(d):+.3f}, positive in {int((d > 0).sum())}/20, '
                 f'one-sided sign P {binomtest(int((d > 0).sum()), 20, 0.5, alternative="greater").pvalue:.4g}')
    L.append(f'    Recharged vs CALVADOS-2 log ESS: Spearman over sites {spearmanr(m.log_ess_rech, m.log_ess_cal).correlation:+.3f}, '
             f'within protein {spearmanr(centre(m.log_ess_rech, g), centre(m.log_ess_cal, g)).correlation:+.3f}')
    L.append('    (extended reading applies only if R1 and R2 were not both positive; they were, so the replication stands)')
    L.append('')
    r3 = spearmanr(m.dRg_rw_rech_pct, m.dRg_direct_rech_pct)
    agree = int((np.sign(m.dRg_rw_rech_pct) == np.sign(m.dRg_direct_rech_pct)).sum())
    L.append(f'M3  Recharged reweighted vs direct Rg change: Spearman {r3.correlation:+.3f} (P {r3.pvalue:.3g}), same sign {agree}/60; '
             f'median |dRg| reweighted {m.dRg_rw_rech_pct.abs().median():.2f}%, direct {m.dRg_direct_rech_pct.abs().median():.2f}%')
    mim = pd.read_csv(MIMIC)
    mim = mim[mim.mimic == 'E'][['UniProt_ID', 'pos', 'ess']].rename(columns={'ess': 'ess_mimic'})
    mm_ = m.merge(mim, on=['UniProt_ID', 'pos'], how='left')
    ok = mm_.ess_mimic.notna()
    if ok.sum() > 5:
        L.append(f'Comparison only: BENDER-MT Recharged Ser/Thr->Glu mimic log ESS vs Recharged log ESS, {int(ok.sum())} sites: Spearman '
                 f'{spearmanr(np.log(mm_.ess_mimic[ok]), mm_.log_ess_rech[ok]).correlation:+.3f}')
    txt = '\n'.join(L) + '\n'
    open(os.path.join(OUT, 'summary_amend3.txt'), 'w').write(txt)
    print(txt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    a = ap.parse_args()
    if a.verify:
        verify()
        return
    analyse(compute())


if __name__ == '__main__':
    main()
