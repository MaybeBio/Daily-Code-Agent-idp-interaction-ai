"""
page4_direct.py — the direct simulations of the PAGE4 phosphoforms (GAIVI_PAGE4; CALVADOS-2 and Mpipi-Recharged, 300 ns
per run, 10 ps frames, first half discarded), against the reweighting predictions of page4_phosphoforms.py, the BENDER
WT ensembles and the measurements (Kulkarni et al. PNAS 2017).

prepare   per run: 3,000 evenly spaced production frames (every 50 ps, the BENDER stored-frame spacing); Rg, distances
          18-63 and 63-102 and FRET efficiencies (page4_phosphoforms.observables); H, dH/dT and G_T at the run's T; the
          solvation forcing D at 300 K (calibrated strengths); the class energies E_c(T) of the Wuttke calibration.
          Mpipi-Recharged: per-force energies at T and T +/- 1 K from the run's own system
            (mpipi_recharged_fdt.build_system plus apply_phospho, taken verbatim from the run script), G_T as
            mpipi_recharged_fdt.per_frame; D from mutation_scan_recharged.pair_terms (gamma_TEA 0.4, gamma_hyd 0.209);
            E_c from tea_calibrate.recharged_class_energies. Check: total energy against the potential energy the
            simulation logged (100 logged steps per run, all on saved frames).
          CALVADOS-2: energies and G_T from TODO_T2_fdt_computation.per_frame_energies with the run's residues file
            (phosphosites as custom letters), pH, ionic strength and T; D as heat_crosses / page4_phosphoforms
            (gamma_TEA 2.0, gamma_hyd 0.358, on the unmodified Ser/Thr letters, whose lambda and sigma the phosphosites
            keep); E_c from tea_reweight.class_energies and tea_calibrate.calvados_hyd_energies. Check (no energy log
            is written): H and dH/dT against CALVADOS' own OpenMM system built from the run's config.yaml and
            components.yaml at T and T +/- 1 K, on 20 frames.
summary   per run: means and 10-block SEs; FDT heating responses, plain (G_T) at every T and calibrated
          (G_T + D / RT) at 300 K, 10-block jackknife SEs. Per form, the heating curve over the simulated temperatures,
          plain and calibrated (each T reweighted with exp(-sum_c gamma_c E_c(T) / RT), as the Wuttke calibration), and
          its slope at 300 K from a quadratic fit (as tea_calibrate.report), jackknife over blocks at all temperatures.
          Comparisons: (1) direct WT against the stored BENDER WT ensemble, (2) every form against its reweighting
          prediction, (3) the FDT response against the finite-temperature slope, (4) the heat x phosphorylation change,
          (5) the measurement-condition runs against SAXS and smFRET. Runs missing from the output root are skipped.
    BENDER_REANALYSIS_DIR=... python page4_direct.py prepare --model recharged|calvados [--workers 8]   (openmpipi env)
    BENDER_REANALYSIS_DIR=... python page4_direct.py summary --model recharged|calvados
Output: Traj Level Analysis/Analysis_FDT/total_response/page4/direct_<model>/{npz/, runs.csv, curves.csv, summary.txt}
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import ast  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import tempfile  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import page4_phosphoforms as P4  # noqa: E402
import tea_reweight as W  # noqa: E402

HC, MS, TF, fdt, L = P4.HC, P4.MS, P4.TF, P4.fdt, P4.L
BM = '/Volumes/User Homes/Nu Project/BENDER-MT'
GAIVI = os.path.expanduser('~/Downloads/Scripts/GAIVI_PAGE4')
K = 3000                      # production frames kept (of 15,000)
T0 = 300.0
RT0 = fdt.R_GAS * T0
CLS = ('aliph', 'pro', 'arom', 'arg', 'rest')
MODELS = {
    'recharged': {'label': 'Mpipi-Recharged', 'root': f'{BM}/page4_MPIPI', 'gaivi': f'{GAIVI}/mpipi',
                  'gam': dict({c: 0.4 for c in CLS}, hyd=0.209),        # A0 strengths, pinned (not the scan env)
                  'pred': {'HIPK1occ': 'HIPK1 (occupancy) (phosphate-like)', 'HIPK1T51': 'HIPK1 T51 only (phosphate-like)',
                           'CLK2occ': 'CLK2 (occupancy) (phosphate-like)', 'CLK2full': 'CLK2 (phosphate-like)',
                           'Nmotif': 'S7+S9 (phosphate-like)', 'HIPK1glu': 'HIPK1 (Glu mimic)', 'CLK2glu': 'CLK2 (Glu mimic)'}},
    'calvados': {'label': 'CALVADOS-2', 'root': f'{BM}/page4_CALVADOS', 'gaivi': f'{GAIVI}/calvados',
                 'gam': dict({c: 2.0 for c in CLS}, hyd=0.358),        # A0 strengths, pinned (not the scan env)
                 'pred': {'HIPK1occ': 'HIPK1 (occupancy)', 'HIPK1T51': 'HIPK1 T51 only', 'CLK2occ': 'CLK2 (occupancy)',
                          'CLK2full': 'CLK2', 'Nmotif': 'S7+S9'}}}
MEAS = {'WT': {'saxs': 36.2, 'E18_63': 0.55, 'E63_102': 0.65},
        'HIPK1occ': {'saxs': 34.7, 'E18_63': 0.52, 'E63_102': 0.65},
        'CLK2occ': {'saxs': 49.8, 'E18_63': 0.35, 'E63_102': 0.58}}
OBS = ['Rg', 'E18_63', 'E63_102']
PHOS = 'BJOUXZ'               # CALVADOS custom phosphoresidue letters (residues_<form>.csv)


def out_dir(model):
    return os.path.join(P4.OUT, f'direct_{model}')


def runs_table(model):
    m = MODELS[model]
    r = pd.read_csv(os.path.join(m['gaivi'], 'runs.tsv'), sep='\t')
    fa = {}
    for block in open(os.path.join(m['gaivi'], 'forms.fasta')).read().split('>')[1:]:
        h, *s = block.split('\n')
        fa[h.strip()] = ''.join(s).strip()
    r['seq'] = r.form.map(fa)
    r['wt_seq'] = fa['WT']
    r['dir'] = [os.path.join(m['root'], f'{int(t)}K', n) for n, t in zip(r['name'], r['temp'])]
    r['present'] = [os.path.exists(os.path.join(d, f'{n}.dcd')) for d, n in zip(r['dir'], r['name'])]
    r['model'] = model
    return r


def load_frames(row):
    """All frames of the run, and K evenly spaced production frames (whole chain). The DCD header's frame count must
    equal the frames in the file: a mismatch means two jobs wrote the same file (found once, CALVADOS-2 Nmotif 300 K)."""
    import struct
    dcd = os.path.join(row['dir'], f"{row['name']}.dcd")
    xyz_all, top_seq, info = HC.bc.load_calvados(dcd, os.path.join(row['dir'], 'top.pdb'), k=None, production=False)
    n = len(xyz_all)
    nset = struct.unpack('<i', open(dcd, 'rb').read(12)[8:12])[0]
    if nset != n:
        raise ValueError(f'DCD header lists {nset} frames but the file holds {n}; rerun this simulation')
    idx = n // 2 + HC.bc.even_indices(n - n // 2, K)
    xyz = xyz_all[idx]
    mb = HC.bc.max_bond(xyz)
    if mb > 1.0:
        xyz = HC.bc.unwrap_chain(xyz, 50.0)
    return xyz_all, xyz, idx, mb, top_seq


# ----------------------------------------------------------------------------- Mpipi-Recharged
def apply_phospho_from_run_script():
    """apply_phospho exactly as run_openmpipi_page4.py defines it (the script runs on import, so the function is
    extracted from its source)."""
    import openmm as mm
    src = open(os.path.join(MODELS['recharged']['gaivi'], 'run_openmpipi_page4.py')).read()
    fn = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == 'apply_phospho')
    ns = {'mm': mm, 'np': np, 'os': os, '__import__': __import__}
    exec(compile(ast.Module([fn], []), 'run_openmpipi_page4.apply_phospho', 'exec'), ns)
    return ns['apply_phospho']


def build(seq, T, ionic, spec, xyz0, name='PAGE4'):
    import mpipi_recharged_fdt as R
    system = R.build_system(name, seq, T, ionic, xyz0)
    if spec not in ('-', '', None):
        apply_phospho_from_run_script()(system, spec, seq)
    return system


def prepare_recharged(row):
    import mpipi_recharged_fdt as R
    import tea_calibrate as TCal
    T, ionic, spec, seq = float(row['temp']), float(row['ionic']), row['phospho'], row['seq']
    xyz_all, xyz, idx, mb, top_seq = load_frames(row)
    if top_seq != seq:
        raise ValueError('top.pdb sequence differs from forms.fasta')
    box = R.box_nm(len(seq))
    ev = {t: R.Evaluator(build(seq, t, ionic, spec, xyz[0]), 'CPU') for t in (T, T - 1.0, T + 1.0)}
    E = {t: e.energies(xyz) for t, e in ev.items()}
    pf = R.per_frame(ev[T].labels, E[T], (E[T + 1.0] - E[T - 1.0]) / 2.0, T)
    # check against the logged potential energy (frame f is step (f + 1) * N_save; the last logged step is the last frame)
    pe = pd.read_csv(os.path.join(row['dir'], f"{row['name']}.log"), sep='\t')
    steps, logged = pe.iloc[:, 0].to_numpy(), pe.iloc[:, 1].to_numpy()
    n_save = int(steps.max() // len(xyz_all))
    f = steps // n_save - 1
    ok = (steps % n_save == 0) & (f >= 0) & (f < len(xyz_all))
    dlog = np.array([ev[T].total(x) for x in xyz_all[f[ok]]]) - logged[ok]
    tb = MS.tables()
    ix = np.array([tb['ix'][a] for a in seq])
    ft = {a: MS.G_TEA * TF.DL['charmm36'][a] + MS.G_HYD * TF.DL_CLASSES['hyd'][a] for a in MS.AA}
    fv = np.array([ft[a] for a in seq])
    i, j = np.triu_indices(len(seq), k=2)
    dv = xyz[:, i] - xyz[:, j]
    dv -= box * np.round(dv / box)
    D = MS.pair_terms(np.linalg.norm(dv, axis=2), ix[i], ix[j], tb, MS.kappa_of(T0), fv[i], fv[j])[3]
    Ec = TCal.recharged_class_energies(xyz, seq, T, box)
    return xyz, idx, mb, pf['G_T'], pf['H'], D, Ec, {
        'check': 'logged potential energy', 'n_check': int(ok.sum()), 'check_H_maxabs': float(np.abs(dlog).max()),
        'check_dHdT_maxabs': float(np.abs(pf['dHdT_non_yukawa']).max())}


# ----------------------------------------------------------------------------- CALVADOS-2
def calvados_own_energies(row, xyz, T):
    """Total potential energy of frames under CALVADOS' own system, built from the run's YAML files at temperature T
    (local residues and FASTA paths; files the build writes go to a temporary folder)."""
    import openmm as mm
    import yaml
    from calvados.sim import Sim
    cfg = yaml.safe_load(open(os.path.join(row['dir'], 'config.yaml')))
    comp = yaml.safe_load(open(os.path.join(row['dir'], 'components.yaml')))
    cfg.update(temp=float(T), platform='CPU', verbose=False)
    cfg.pop('runtime', None)                     # run length only; the local CALVADOS version rejects runtime=0
    comp['defaults'].update(ffasta=os.path.join(row['dir'], f"{row['name']}.fasta"),
                            fresidues=os.path.join(MODELS['calvados']['gaivi'], row['residues']))
    with tempfile.TemporaryDirectory() as tmp:
        sim = Sim(tmp, cfg, comp)
        sim.build_system()
        ctx = mm.Context(sim.system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        out = []
        for x in xyz:
            ctx.setPositions(x)
            out.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole))
    return np.array(out)


def prepare_calvados(row):
    import tea_calibrate as TCal
    G = MODELS['calvados']['gaivi']
    T, ionic, ph, seq, base = float(row['temp']), float(row['ionic']), float(row['ph']), row['seq'], row['wt_seq']
    fseq = ''.join(open(os.path.join(row['dir'], f"{row['name']}.fasta")).read().split('\n')[1:]).strip()
    if fseq != seq:
        raise ValueError("run FASTA differs from forms.fasta")
    res = pd.read_csv(os.path.join(G, row['residues'])).set_index('one')
    res0 = pd.read_csv(os.path.join(G, 'residues_WT.csv')).set_index('one')
    for k, (a, b) in enumerate(zip(seq, base)):       # phosphosites: custom letters on Ser/Thr with their lambda, sigma
        if a != b and not (a in PHOS and b in 'ST' and np.allclose(res.loc[a, ['lambdas', 'sigmas']].astype(float),
                                                                   res0.loc[b, ['lambdas', 'sigmas']].astype(float))):
            raise ValueError(f'position {k + 1}: {a} is not a phosphorylated {b}')
    xyz_all, xyz, idx, mb, _ = load_frames(row)
    q, lam, sig, b0 = fdt.calvados_residue_params(seq, res, ph=ph)
    pf = fdt.per_frame_energies(xyz, q, lam, sig, b0, fdt.get_dielectric_params(T, ionic))
    i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
    dl = {k: np.array([v[a] for a in base]) for k, v in (('tea', TF.DL['charmm36']), ('hyd', TF.DL_CLASSES['hyd']))}
    fpair = HC.G_TEA * 0.5 * (dl['tea'][i] + dl['tea'][j]) + HC.G_HYD * 0.5 * (dl['hyd'][i] + dl['hyd'][j])
    D = np.zeros(len(xyz))
    step = max(1, int(fdt.MAX_ELEMS // len(i)))
    for a in range(0, len(xyz), step):
        r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
        D[a:a + step] = ((fdt.ah_pair_energy(r, s, np.ones_like(l), shift)
                          - fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)) * fpair).sum(1)
    Ec = W.class_energies(xyz, base, T, residues=res0)
    Ec.update(TCal.calvados_hyd_energies(xyz, base, T))
    # check against CALVADOS' own system
    sel = np.linspace(0, len(xyz) - 1, 20).astype(int)
    own = {t: calvados_own_energies(row, xyz[sel], t) for t in (T, T - 1.0, T + 1.0)}
    dH = own[T] - pf['H'][sel]
    dT = (own[T + 1.0] - own[T - 1.0]) / 2.0 - pf['dVdT_DH'][sel]
    return xyz, idx, mb, pf['G_T'], pf['H'], D, Ec, {
        'check': "CALVADOS' own system", 'n_check': len(sel), 'check_H_maxabs': float(np.abs(dH).max()),
        'check_dHdT_maxabs': float(np.abs(dT).max())}


# ----------------------------------------------------------------------------- prepare
def _prepare_one(row):
    out = os.path.join(out_dir(row['model']), 'npz', f"{row['name']}_{int(row['temp'])}K.npz")
    tag = f"{row['name']} {int(row['temp'])} K"
    if os.path.exists(out):
        return tag, 'exists'
    try:
        fn = prepare_recharged if row['model'] == 'recharged' else prepare_calvados
        xyz, idx, mb, G_T, H, D, Ec, chk = fn(row)
        obs = P4.observables(xyz)
        np.savez_compressed(out, T=float(row['temp']), ionic=float(row['ionic']), form=row['form'], frame_index=idx,
                            max_bond_raw=mb, G_T=G_T, H=H, D300=D, **chk, **{f'obs_{k}': v for k, v in obs.items()},
                            **{f'E_{c}': v for c, v in Ec.items()})
        return tag, (f"ok: {chk['check']}, {chk['n_check']} frames, max |d| H {chk['check_H_maxabs']:.4f} kJ/mol, "
                     f"dH/dT {chk['check_dHdT_maxabs']:.2e}")
    except Exception as ex:  # noqa: BLE001
        return tag, f'FAILED {ex!r} | {traceback.format_exc()[-400:]}'.replace('\n', ' ')


def prepare(model, workers):
    os.makedirs(os.path.join(out_dir(model), 'npz'), exist_ok=True)
    r = runs_table(model)
    for _, x in r[~r.present].iterrows():
        print(f"{x['name']} {int(x['temp'])} K: not in {MODELS[model]['root']}, skipped")
    with mp.get_context('spawn').Pool(workers) as pool:
        for tag, msg in pool.imap_unordered(_prepare_one, r[r.present].to_dict('records')):
            print(f'{tag}: {msg}', flush=True)


# ----------------------------------------------------------------------------- summary
def block_mean_se(a, n=10):
    b = [x.mean() for x in np.array_split(a, n)]
    return float(a.mean()), float(np.std(b, ddof=1) / np.sqrt(n))


def fdt_total(A, G):
    f = lambda ix: -HC.wcov(A[ix], G[ix], np.ones(len(ix))) / A[ix].mean() * 100  # noqa: E731
    return f(np.arange(len(A))), HC.jk(f, len(A))


def curve_slope(series, gam, drop=None):
    """Relative curve <Rg>(T) / <Rg>(300) and its slope at 300 K (% per K, quadratic fit); drop = block to leave out."""
    m = {}
    for T, z in series.items():
        A = z['obs_Rg']
        keep = np.ones(len(A), bool)
        if drop is not None:
            keep[np.array_split(np.arange(len(A)), HC.N_BLOCKS)[drop]] = False
        m[T] = W.reweight(A[keep], {c: z[f'E_{c}'][keep] for c in gam}, gam, float(T))[0]
    Ts = np.array(sorted(m))
    rel = np.array([m[T] / m[300] for T in Ts])
    return Ts, rel, 100 * np.polyfit(Ts - T0, rel, 2)[1]


def slope_with_se(series, gam):
    """Curve, slope and their jackknife SEs (blocks left out at all temperatures together)."""
    Ts, rel, s = curve_slope(series, gam)
    jk = [curve_slope(series, gam, b) for b in range(HC.N_BLOCKS)]
    f = (HC.N_BLOCKS - 1) / HC.N_BLOCKS
    v = np.array([x[2] for x in jk])
    r = np.array([x[1] for x in jk])
    return Ts, rel, s, float(np.sqrt(f * ((v - v.mean()) ** 2).sum())), np.sqrt(f * ((r - r.mean(0)) ** 2).sum(0))


def pct(m1, s1, m0, s0):
    return 100 * (m1 / m0 - 1), 100 * m1 / m0 * np.sqrt((s1 / m1) ** 2 + (s0 / m0) ** 2)


def stored_wt(model):
    """The stored BENDER WT ensemble (300 K, 150 mM, 1,000 frames): <Rg>, E means and the calibrated total of each, as
    page4_phosphoforms builds them."""
    uid = P4.UID
    if model == 'recharged':
        import mpipi_recharged_fdt as R
        d = MS.TC.rc_dirs()[uid]
        z = np.load(os.path.join(MS.TC.RC_NPZ, f'{uid}.npz'))
        xyz, seq, info = R.load_run(d, uid)
        if not np.array_equal(info['production_start'] + info['stride'] * np.arange(len(xyz)), z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        tb, box = MS.tables(), R.box_nm(len(seq))
        ix = np.array([tb['ix'][a] for a in seq])
        ft = {a: MS.G_TEA * TF.DL['charmm36'][a] + MS.G_HYD * TF.DL_CLASSES['hyd'][a] for a in MS.AA}
        fv = np.array([ft[a] for a in seq])
        i, j = np.triu_indices(len(seq), k=2)
        dv = xyz[:, i] - xyz[:, j]
        dv -= box * np.round(dv / box)
        D = MS.pair_terms(np.linalg.norm(dv, axis=2), ix[i], ix[j], tb, MS.kappa_of(T0), fv[i], fv[j])[3]
    else:
        d = TF.index_dirs([uid])[uid]
        z = np.load(os.path.join(TF.NPZ, f'{uid}.npz'))
        xyz_all, seq, info = HC.bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=None)
        idx = HC.bc.even_indices(len(xyz_all), L.K_FRAMES)
        if not np.array_equal(info['production_start'] + idx, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        xyz = xyz_all[idx].astype(float)
        res0 = pd.read_csv(os.path.join(MODELS['calvados']['gaivi'], 'residues_WT.csv')).set_index('one')
        q, lam, sig, _ = fdt.calvados_residue_params(seq, res0)
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        dl = {k: np.array([v[a] for a in seq]) for k, v in (('tea', TF.DL['charmm36']), ('hyd', TF.DL_CLASSES['hyd']))}
        fpair = HC.G_TEA * 0.5 * (dl['tea'][i] + dl['tea'][j]) + HC.G_HYD * 0.5 * (dl['hyd'][i] + dl['hyd'][j])
        r = np.linalg.norm(xyz[:, i] - xyz[:, j], axis=2)
        D = ((fdt.ah_pair_energy(r, s, np.ones_like(l), shift) - fdt.ah_pair_energy(r, s, np.zeros_like(l), shift))
             * fpair).sum(1)
    G = z['G_T'].astype(float) + D / RT0
    obs = P4.observables(xyz)
    out = {}
    for k in OBS:
        out[f'{k}_mean'], out[f'{k}_se'] = block_mean_se(obs[k])
        out[f'{k}_total'], out[f'{k}_total_se'] = fdt_total(obs[k], G)
    return out


def summary(model):
    M, od = MODELS[model], out_dir(model)
    runs = runs_table(model)
    runs['present'] = runs['present'].to_numpy() & np.array([os.path.exists(os.path.join(od, 'npz', f"{n}_{int(t)}K.npz")) for n, t in zip(runs['name'], runs.temp)])
    missing = runs[~runs.present]
    runs = runs[runs.present]
    rows, data = [], {}
    for _, r in runs.iterrows():
        z = dict(np.load(os.path.join(od, 'npz', f"{r['name']}_{int(r['temp'])}K.npz")))
        data[(r.form, int(r.temp), float(r.ionic))] = z
        row = {'name': r['name'], 'form': r.form, 'set': r.set, 'T': int(r.temp), 'ionic': float(r.ionic),
               'check_H_maxabs': float(z['check_H_maxabs']), 'check_dHdT_maxabs': float(z['check_dHdT_maxabs']),
               'n_check': int(z['n_check']), 'max_bond_raw': float(z['max_bond_raw'])}
        for k in OBS + ['d18_63_r2', 'd63_102_r2']:
            row[f'{k}_mean'], row[f'{k}_se'] = block_mean_se(z[f'obs_{k}'])
        for k in OBS:
            row[f'{k}_plain'], row[f'{k}_plain_se'] = fdt_total(z[f'obs_{k}'], z['G_T'])
            if int(r.temp) == 300:
                row[f'{k}_total'], row[f'{k}_total_se'] = fdt_total(z[f'obs_{k}'], z['G_T'] + z['D300'] / RT0)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(od, 'runs.csv'), index=False)

    crow = []
    for form in runs[runs.set == 'heat'].form.unique():
        series = {T: data[(form, T, 0.15)] for T in sorted(runs[(runs.form == form) & (runs.set == 'heat')].temp.astype(int))}
        if 300 not in series or len(series) < 4:
            continue
        for lab, gam in (('plain', {c: 0.0 for c in M['gam']}), ('calibrated', M['gam'])):
            Ts, rel, s, se, rse = slope_with_se(series, gam)
            crow.append({'form': form, 'curve': lab, 'temps': ','.join(map(str, Ts)), 'slope300': s, 'slope300_se': se,
                         **{f'rel_{T}K': v for T, v in zip(Ts, rel)}, **{f'rel_se_{T}K': v for T, v in zip(Ts, rse)}})
    cv = pd.DataFrame(crow)
    cv.to_csv(os.path.join(od, 'curves.csv'), index=False)

    rw = pd.read_csv(os.path.join(P4.TAB, 'phosphoforms.csv'))      # the table of the active forcing (A0: page4_A0_linear)
    rw = rw[rw.model == M['label']].set_index('pattern')
    h = df[(df['T'] == 300) & (df.ionic == 0.15)].set_index('form')
    wt = h.loc['WT']
    Lns = [f"PAGE4 direct {M['label']} runs (300 ns each, 3,000 production frames), 10-block SEs; {len(df)} runs",
           f"check ({df.n_check.iloc[0]} frames per run): max |d| H {df.check_H_maxabs.max():.4f} kJ/mol, "
           f"dH/dT {df.check_dHdT_maxabs.max():.2e} kJ/mol/K"
           + (" (recharged: dH/dT outside the Yukawa term)" if model == 'recharged' else ''),
           'not available: ' + (', '.join(f"{n} {int(t)} K" for n, t in zip(missing['name'], missing.temp)) or 'none'), '']
    b = stored_wt(model)
    if abs(b['Rg_total'] - rw.loc['WT', 'Rg_total_wt']) > 1e-6:
        raise ValueError('stored WT total differs from phosphoforms.csv')
    Lns += ['(1) WT at 300 K, 150 mM: direct run (300 ns) vs stored BENDER ensemble (100 ns, 1,000 frames); direct - stored in SEs']
    for k, f in (('Rg', 10), ('E18_63', 1), ('E63_102', 1)):
        for q, e in (('mean', 'se'), ('total', 'total_se')):
            a_, b_ = wt[f'{k}_{q}'], b[f'{k}_{q}']
            se = np.hypot(wt[f'{k}_{e}'], b[f'{k}_{e}'])
            sc = f if q == 'mean' else 1
            Lns.append(f"  {k} {q}: {sc * a_:+.3f} +- {sc * wt[f'{k}_{e}']:.3f} vs {sc * b_:+.3f} +- {sc * b[f'{k}_{e}']:.3f}"
                       f"  ({(a_ - b_) / se:+.1f} SE)")
    Lns.append('')
    Lns += ['(2) phosphoforms at 300 K, 150 mM: direct vs reweighted from the BENDER WT ensemble (ESS = Kish fraction)']
    t2 = []
    for form, pat in M['pred'].items():
        if form not in h.index:
            continue
        z, p = h.loc[form], rw.loc[pat]
        dr, dse = pct(z.Rg_mean, z.Rg_se, wt.Rg_mean, wt.Rg_se)
        t2.append({'form': form, 'ESS': p.ess, 'dRg% direct': dr, 'se': dse, 'dRg% reweighted': p.Rg_pct, 'se_rw': p.Rg_pct_se,
                   'E18 direct': z.E18_63_mean, 'E18 rw': p.E18_63_mean, 'E63 direct': z.E63_102_mean, 'E63 rw': p.E63_102_mean,
                   'dtotal direct': z.Rg_total - wt.Rg_total, 'dtotal_se': np.hypot(z.Rg_total_se, wt.Rg_total_se),
                   'dtotal rw': p.Rg_dtotal, 'dtotal_rw_se': p.Rg_dtotal_se})
    Lns += [pd.DataFrame(t2).round(3).to_string(index=False), '']
    Lns += ['(3) heating response of Rg at 300 K (% per K): FDT at 300 K vs slope of the direct curve (temperatures listed)']
    t3 = []
    for form in cv.form.unique() if len(cv) else []:
        c = cv[cv.form == form].set_index('curve')
        z = h.loc[form]
        t3.append({'form': form, 'temps': c.loc['plain', 'temps'], 'FDT plain': z.Rg_plain, 'se': z.Rg_plain_se,
                   'curve plain': c.loc['plain', 'slope300'], 'se_c': c.loc['plain', 'slope300_se'],
                   'FDT calibrated': z.Rg_total, 'se_t': z.Rg_total_se,
                   'curve calibrated': c.loc['calibrated', 'slope300'], 'se_cc': c.loc['calibrated', 'slope300_se']})
    Lns += [pd.DataFrame(t3).round(3).to_string(index=False), '']
    Lns += ['(4) change of the calibrated heating response from WT (FDT, 300 K), direct; per 10 K for E']
    for form in [f for f in ('HIPK1occ', 'HIPK1T51', 'CLK2occ', 'CLK2full', 'Nmotif', 'CLK2glu') if f in h.index]:
        z = h.loc[form]
        s = []
        for k in OBS:
            if k == 'Rg':
                s.append(f"Rg {z.Rg_total:+.3f} vs {wt.Rg_total:+.3f}: {z.Rg_total - wt.Rg_total:+.3f} +- "
                         f"{np.hypot(z.Rg_total_se, wt.Rg_total_se):.3f} % per K")
            else:
                per10 = lambda q, c='total': q[f'{k}_{c}'] / 100 * q[f'{k}_mean'] * 10  # noqa: E731
                s.append(f'{k} dE/10K {per10(wt):+.4f} -> {per10(z):+.4f} ({per10(z) - per10(wt):+.4f} +- '
                         f"{np.hypot(per10(z, 'total_se'), per10(wt, 'total_se')):.4f})")
        Lns.append(f'  {form}: ' + '; '.join(s))
    Lns.append('')
    m = df[df['T'] == 295].set_index(['form', 'ionic'])
    Lns += ['(5) measurement conditions, 295 K (Rg: SAXS buffer, I 0.19 M; E: smFRET buffer, I 0.128 M); measured in brackets']
    for form, v in MEAS.items():
        part = []
        if (form, 0.19) in m.index:
            s = m.loc[(form, 0.19)]
            part.append(f"Rg {10 * s.Rg_mean:.1f} +- {10 * s.Rg_se:.1f} A [{v['saxs']}]")
        if (form, 0.128) in m.index:
            f_ = m.loc[(form, 0.128)]
            part.append(f"E18-63 {f_.E18_63_mean:.3f} +- {f_.E18_63_se:.3f} [{v['E18_63']}]; "
                        f"E63-102 {f_.E63_102_mean:.3f} +- {f_.E63_102_se:.3f} [{v['E63_102']}]")
        Lns.append(f'  {form}: ' + ('; '.join(part) if part else 'no runs'))
    if ('WT', 0.19) in m.index:
        for form in ['HIPK1occ', 'CLK2occ']:
            s, s0 = m.loc[(form, 0.19)], m.loc[('WT', 0.19)]
            d_, se = pct(s.Rg_mean, s.Rg_se, s0.Rg_mean, s0.Rg_se)
            Lns.append(f"  {form} Rg change from WT: {d_:+.1f} +- {se:.1f} % [SAXS {100 * (MEAS[form]['saxs'] / 36.2 - 1):+.1f} %]")
    open(os.path.join(od, 'summary.txt'), 'w').write('\n'.join(Lns) + '\n')
    print('\n'.join(Lns))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['prepare', 'summary'])
    ap.add_argument('--model', choices=list(MODELS), required=True)
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    prepare(a.model, a.workers) if a.mode == 'prepare' else summary(a.model)
