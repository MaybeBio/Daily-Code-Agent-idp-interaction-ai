"""
era_direct.py — the direct CALVADOS-2 and Mpipi-Recharged simulations of the ER-NTD Ser118 forms (GAIVI_ERA,
make_era_inputs.py), as
page4_direct.py does for PAGE4: size and its temperature response for WT, S118A, S118D, pS118 and E56Q-pS118, against
the SEC-SAXS data of Du et al. 2025 (SASBDB SASDSK2, SASDSF5, SASDVC5).

prepare   per run: 3,000 evenly spaced production frames (page4_direct.load_frames; DCD header checked); sequence checked
          against GAIVI_ERA forms.fasta from the residue names of top.pdb (pS118 = SEP); Rg, the Kirkwood 1/r mean
          (Rh = 1/<invR>) and contacts between the two aromatic clusters of Du et al. (cluster I = UniProt 52-62,
          cluster II = 117-125, bead centres closer than 1.2 nm); energies, dH/dT and G_T at the run's T with the run's
          residues file, pH and ionic strength (TODO_T2_fdt_computation.per_frame_energies, as page4_direct, validated there
          against CALVADOS' own system; the run folders here hold only the DCD and top.pdb, so that check is not repeated);
          the solvation forcing D (gamma_TEA 2.0, gamma_hyd 0.358) and the class energies E_c(T) on the form's sequence
          with the phosphosite read as Ser (it keeps the Ser hydration, as in the PAGE4 runs).
          Mpipi-Recharged: per-force energies at T and T +/- 1 K from the run's own system (mpipi_recharged_fdt.build_system
          plus apply_phospho, the function of GAIVI_PAGE4/mpipi/run_openmpipi_page4.py recovered verbatim from the session
          that wrote it, since the script is not on this machine); D from mutation_scan_recharged.pair_terms (gamma_TEA 0.4,
          gamma_hyd 0.209) and E_c from tea_calibrate.recharged_class_energies, both on the run sequence (Glu bead at the
          phosphosite), as page4_direct.prepare_recharged. The run folders hold no energy log, so the logged-energy check
          of page4_direct is not repeated.
summary   per form: <Rg> at every temperature, plain and calibrated (each T reweighted with exp(-sum_c gamma_c E_c / RT)),
          slope at 300 K from a quadratic fit with jackknife SE (page4_direct.slope_with_se); the fluctuation-response
          prediction at 300 K, plain (G_T) and calibrated (G_T + D/RT), for Rg, Rh and the cluster contacts; and the 277 K
          runs under the two SAXS buffers against the measured Rg.
    BENDER_REANALYSIS_DIR=... python era_direct.py prepare --model calvados|recharged [--workers 6]   (openmpipi env)
    BENDER_REANALYSIS_DIR=... python era_direct.py summary --model calvados|recharged
Output: Traj Level Analysis/Analysis_FDT/total_response/era/direct_<model>/{npz/, runs.csv, curves.csv, summary.txt}
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import page4_direct as PD  # noqa: E402

HC, TF, fdt, W = PD.HC, PD.TF, PD.fdt, PD.W
BM = PD.BM
PKG = os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_ERA'))
OUTROOT = f'{BM}/Traj Level Analysis/Analysis_FDT/total_response/era'
CFG = {'calvados': {'root': f'{BM}/era_CALVADOS2', 'gaivi': f'{PKG}/calvados', 'label': 'CALVADOS-2'},
       'recharged': {'root': f'{BM}/era_MPIPI', 'gaivi': f'{PKG}/mpipi', 'label': 'Mpipi-Recharged'}}
OFFSET = 3                                     # SNA tag: UniProt residue n is construct position n + 3
CL1 = np.arange(52, 63) + OFFSET - 1           # 0-based construct indices of UniProt 52-62
CL2 = np.arange(117, 126) + OFFSET - 1         # UniProt 117-125
THREE = {'SEP': 'B'}
TEMPS = [278, 288, 300, 310, 320]
MEAS = {('WT', 'SAXS-I195'): 33.9, ('pS118', 'SAXS-I095'): 33.9 + 5.2}   # Du et al. form-factor Rg (A)


def out_dir(model):
    return os.path.join(OUTROOT, f'direct_{model}')


def runs(model):
    c = CFG[model]
    r = pd.read_csv(os.path.join(c['gaivi'], 'runs.tsv'), sep='\t')
    fa = {}
    for block in open(os.path.join(c['gaivi'], 'forms.fasta')).read().split('>')[1:]:
        h, *s = block.split('\n')
        fa[h.strip()] = ''.join(s).strip()
    r['seq'] = r.form.map(fa)
    r['dir'] = [os.path.join(c['root'], f'{int(t)}K', n) for n, t in zip(r['name'], r['temp'])]
    r['present'] = [os.path.exists(os.path.join(d, f'{n}.dcd')) for d, n in zip(r['dir'], r['name'])]
    r['model'] = model
    return r


def apply_phospho(system, spec, seq):
    """Give each listed residue (a Glu bead) its own type: Glu's Wang-Frenkel parameters, Yukawa amplitudes x scale.
    Verbatim from GAIVI_PAGE4/mpipi/run_openmpipi_page4.py (the function the ER-NTD Mpipi runs used), except the
    progress print."""
    import openmm as mm
    sites = [(int(p) - 1, float(s)) for p, s in (x.split(':') for x in spec.split(','))]
    for i, _ in sites:
        if seq[i] != 'E':
            raise ValueError(f'residue {i + 1} is {seq[i]}, expected the Glu bead from forms.fasta')
    wf = next(f for f in system.getForces() if isinstance(f, mm.CustomNonbondedForce) and f.getNumTabulatedFunctions()
              and f.getTabulatedFunctionName(0) == 'wf_table')
    yu = next(f for f in system.getForces() if isinstance(f, mm.CustomNonbondedForce) and f.getNumTabulatedFunctions()
              and f.getTabulatedFunctionName(0) == 'A_table')
    n0, glu = 21, 7                                # OpenMpipi type table size and Glu's index (mapping_dict 'pE')
    nb = np.loadtxt(os.path.join(os.path.dirname(__import__('OpenMpipi').__file__), 'data_files', 'recharged_params.txt'))
    W0, A0 = nb[:-n0 * n0].reshape(3, n0, n0), nb[-n0 * n0:].reshape(n0, n0)   # [z, y, x] and [y, x]
    n = n0 + len(sites)
    m = np.r_[np.arange(n0), np.full(len(sites), glu)]           # new type -> Glu
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


def top_sequence(d):
    import mdtraj as md
    t = md.load_topology(os.path.join(d, 'top.pdb'))
    one = {**{k.upper(): v for k, v in md.core.residue_names._AMINO_ACID_CODES.items() if v}, **THREE}
    return ''.join(one[r.name] for r in t.residues)


def observables(xyz):
    n = xyz.shape[1]
    rg = np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    iu = np.triu_indices(n, 1)
    inv = np.zeros(len(xyz))
    for a in range(0, len(xyz), 200):
        dd = np.linalg.norm(xyz[a:a + 200, iu[0]] - xyz[a:a + 200, iu[1]], axis=2)
        inv[a:a + 200] = 2.0 * (1.0 / dd).sum(1) / n ** 2
    dc = np.linalg.norm(xyz[:, CL1, None, :] - xyz[:, None, CL2, :], axis=3)
    return {'obs_Rg': rg, 'obs_invR': inv, 'obs_c12': (dc < 1.2).sum((1, 2)).astype(float)}


def prepare_recharged(row):
    import mpipi_recharged_fdt as R
    import tea_calibrate as TCal
    MS = PD.MS
    T, ionic, spec, seq = float(row['temp']), float(row['ionic']), row['phospho'], row['seq']
    xyz_all, xyz, idx, mb, top_seq = PD.load_frames(row)
    if top_seq != seq:
        raise ValueError('top.pdb sequence differs from forms.fasta')
    box = R.box_nm(len(seq))

    def build(t):
        system = R.build_system('ERa', seq, t, ionic, xyz[0])
        if spec not in ('-', '', None):
            apply_phospho(system, spec, seq)
        return system
    ev = {t: R.Evaluator(build(t), 'CPU') for t in (T, T - 1.0, T + 1.0)}
    E = {t: e.energies(xyz) for t, e in ev.items()}
    pf = R.per_frame(ev[T].labels, E[T], (E[T + 1.0] - E[T - 1.0]) / 2.0, T)
    tb = MS.tables()
    ix = np.array([tb['ix'][a] for a in seq])
    ft = {a: MS.G_TEA * TF.DL['charmm36'][a] + MS.G_HYD * TF.DL_CLASSES['hyd'][a] for a in MS.AA}
    fv = np.array([ft[a] for a in seq])
    i, j = np.triu_indices(len(seq), k=2)
    dv = xyz[:, i] - xyz[:, j]
    dv -= box * np.round(dv / box)
    D = MS.pair_terms(np.linalg.norm(dv, axis=2), ix[i], ix[j], tb, MS.kappa_of(PD.T0), fv[i], fv[j])[3]
    Ec = TCal.recharged_class_energies(xyz, seq, T, box)
    return xyz_all, xyz, idx, mb, pf, D, Ec


def prepare_one(row):
    out = os.path.join(out_dir(row['model']), 'npz', f"{row['name']}_{int(row['temp'])}K.npz")
    tag = f"{row['name']} {int(row['temp'])} K"
    if os.path.exists(out):
        return tag, 'exists'
    try:
        if row['model'] == 'recharged':
            xyz_all, xyz, idx, mb, pf, D, Ec = prepare_recharged(row)
            os.makedirs(os.path.dirname(out), exist_ok=True)
            np.savez_compressed(out, T=float(row['temp']), ionic=float(row['ionic']), form=row['form'], set=row['set'],
                                frame_index=idx, max_bond_raw=mb, G_T=pf['G_T'], H=pf['H'], D300=D, **observables(xyz),
                                **{f'E_{c}': v for c, v in Ec.items()})
            return tag, f'ok, {len(xyz_all)} frames, max bond {mb:.2f} nm'
        import tea_calibrate as TCal
        GAIVI = CFG['calvados']['gaivi']
        T, ionic, ph, seq = float(row['temp']), float(row['ionic']), float(row['ph']), row['seq']
        ts = top_sequence(row['dir'])
        if ts != seq:
            raise ValueError(f'top.pdb sequence differs from forms.fasta ({len(ts)} vs {len(seq)} residues)')
        base = seq.replace('B', 'S')                      # the phosphosite keeps the Ser hydration
        res = pd.read_csv(os.path.join(GAIVI, row['residues'])).set_index('one')
        res0 = pd.read_csv(os.path.join(GAIVI, 'residues_WT.csv')).set_index('one')
        xyz_all, xyz, idx, mb, _ = PD.load_frames(row)
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
        os.makedirs(os.path.dirname(out), exist_ok=True)
        np.savez_compressed(out, T=T, ionic=ionic, form=row['form'], set=row['set'], frame_index=idx, max_bond_raw=mb,
                            G_T=pf['G_T'], H=pf['H'], D300=D, **observables(xyz), **{f'E_{c}': v for c, v in Ec.items()})
        return tag, f'ok, {len(xyz_all)} frames, max bond {mb:.2f} nm'
    except Exception as ex:  # noqa: BLE001
        return tag, f'FAILED {ex!r} | {traceback.format_exc()[-400:]}'.replace('\n', ' ')


def prepare(model, workers):
    r = runs(model)
    for _, x in r[~r.present].iterrows():
        print(f"{x['name']} {int(x['temp'])} K: missing, skipped")
    with mp.get_context('spawn').Pool(workers) as pool:
        for tag, msg in pool.imap_unordered(prepare_one, r[r.present].to_dict('records')):
            print(f'{tag}: {msg}', flush=True)


def load(model, name, T):
    return dict(np.load(os.path.join(out_dir(model), 'npz', f'{name}_{T}K.npz')))


def summary(model):
    r = runs(model)
    GAM = PD.MODELS[model]['gam']
    RT = PD.RT0 if model == 'calvados' else PD.MS.RT
    OUT = out_dir(model)
    L, rows, cv = [f"ER-NTD direct {CFG[model]['label']} runs (Du et al. 2025 construct, SNA + 1-184)", ''], [], []
    L.append('== heating series 278-320 K, 150 mM, pH 7.4: <Rg> slope at 300 K (% per K, quadratic fit, jackknife SE) and')
    L.append('   the fluctuation-response prediction at 300 K (plain G_T; calibrated G_T + D/RT) for Rg, Rh and cluster I-II contacts')
    heat = r[r.set == 'heat']
    for form in heat.form.unique():
        name = heat[heat.form == form]['name'].iloc[0]
        ser = {T: load(model, name, T) for T in TEMPS}
        out = {'form': form}
        for lab, gam in (('plain', {c: 0.0 for c in GAM}), ('calibrated', GAM)):
            Ts, rel, s, se, rse = PD.slope_with_se(ser, gam)
            out[f'slope_{lab}'], out[f'slope_{lab}_se'] = s, se
            for T, v, e in zip(Ts, rel, rse):
                cv.append({'form': form, 'hydration': lab, 'T': int(T), 'Rg_rel': v, 'se': e})
        z = ser[300]
        Gt = z['G_T'] + z['D300'] / RT
        for obs in ('Rg', 'invR', 'c12'):
            for lab, G in (('plain', z['G_T']), ('calibrated', Gt)):
                f, se = PD.fdt_total(z[f'obs_{obs}'], G)
                key = 'Rh' if obs == 'invR' else obs
                out[f'fdt_{key}_{lab}'] = -f if obs == 'invR' else f
                out[f'fdt_{key}_{lab}_se'] = se
        out['Rg300_A'] = 10 * z['obs_Rg'].mean()
        out['c12_300'] = z['obs_c12'].mean()
        rows.append(out)
    df = pd.DataFrame(rows)
    w = df.set_index('form')
    for form, x in w.iterrows():
        L.append(f"  {form:10s} Rg {x.Rg300_A:5.1f} A, contacts {x.c12_300:5.2f} | direct slope plain {x.slope_plain:+.3f}+-{x.slope_plain_se:.3f}, "
                 f"calibrated {x.slope_calibrated:+.3f}+-{x.slope_calibrated_se:.3f} | FDT Rg {x.fdt_Rg_plain:+.3f}/{x.fdt_Rg_calibrated:+.3f}, "
                 f"Rh {x.fdt_Rh_plain:+.3f}/{x.fdt_Rh_calibrated:+.3f}, contacts {x.fdt_c12_plain:+.2f}/{x.fdt_c12_calibrated:+.2f}")
    L.append('  form minus WT (calibrated direct slope; FDT Rg calibrated):')
    for form, x in w.iterrows():
        if form != 'WT':
            d = x.slope_calibrated - w.loc['WT', 'slope_calibrated']
            se = np.hypot(x.slope_calibrated_se, w.loc['WT', 'slope_calibrated_se'])
            L.append(f"    {form:10s} {d:+.3f} +- {se:.3f}   FDT {x.fdt_Rg_calibrated - w.loc['WT', 'fdt_Rg_calibrated']:+.3f}")
    L += ['', '== SAXS conditions, 277 K: <Rg> (A) +- 10-block SE; change from WT in the same buffer']
    sx = r[r.set.str.startswith('SAXS')]
    vals = {}
    for _, x in sx.iterrows():
        z = load(model, x['name'], 277)
        m, se = PD.block_mean_se(10 * z['obs_Rg'])
        vals[(x.form, x.set)] = (m, se)
    for (form, st), (m, se) in sorted(vals.items(), key=lambda k: (k[0][1], k[0][0])):
        wt = vals.get(('WT', st))
        d = f', {m - wt[0]:+.2f} +- {np.hypot(se, wt[1]):.2f} from WT' if wt and form != 'WT' else ''
        meas = f' [measured {MEAS[(form, st)]:.1f}]' if (form, st) in MEAS else ''
        L.append(f'  {st:9s} {form:10s} {m:5.2f} +- {se:.2f}{d}{meas}')
    L.append('  measured (Du et al.): WT 33.9 +- 0.2 A; pS118 +5.2 +- 0.4 A (pS118 in 0.095 M, WT in 0.195 M); S118D about +3 A')
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, 'runs.csv'), index=False)
    pd.DataFrame(cv).to_csv(os.path.join(OUT, 'curves.csv'), index=False)
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['prepare', 'summary'])
    ap.add_argument('--model', choices=list(CFG), required=True)
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    prepare(a.model, a.workers) if a.mode == 'prepare' else summary(a.model)
