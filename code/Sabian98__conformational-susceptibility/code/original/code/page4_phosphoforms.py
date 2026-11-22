"""
page4_phosphoforms.py — PAGE4 (O60829) phosphoforms by exact multisite reweighting of the stored 300 K ensembles, in
CALVADOS-2 and Mpipi-Recharged, for comparison with the HIPK1-PAGE4 and CLK2-PAGE4 measurements (SAXS, smFRET; Kulkarni
et al. PNAS 2017 114:E2644; occupancies from mass spectrometry/NMR there, mean-field charges as Lin et al. J Mol Biol 2018):
  HIPK1  T51 (>95%), S9 (30-40%)
  CLK2   S7, S9, T51, T85 (>95%), S73 (~75%), T94 (~60%), T71 (~50%), S79 (~50%)
Measured (Kulkarni et al. PNAS 2017; also Lin et al. J Mol Biol 2018 Table 1), WT / HIPK1 / CLK2:
  SAXS Rg 36.2+-1.1 / 34.7+-1.2 / 49.8+-1.9 A (100 mM KPi, 10 mM NaN3, 0.5 mM DTT, pH 7.0; temperature not stated)
  smFRET <E> A18C:C63 0.55 / 0.52 / 0.35, C63:P102 0.65 / 0.65 / 0.58 (Alexa 555/647, R0 5.24 nm, gamma 1; 20 mM Tris
  pH 7.8, 100 mM NaOAc, 5 mM MgCl2); RMS distances 56 / 59 / 75 A and 50 / 50 / 55 A; smFRET-derived Rg 34 / 35 / 43 A.

Patterns, each with full occupancy and with occupancy-weighted charge (mean field, as Lin et al.):
  CALVADOS-2       phosphate = charge -2 on the Ser/Thr bead (as heat_crosses.py), dH = E_DH(q') - E_DH(q) per frame
                   (mutual terms between phosphosites included), dG from the dielectric derivative of the same sum
  Mpipi-Recharged  (a) Glu mimic at every site (full occupancy only); (b) phosphate-like: Glu bead with its Yukawa
                   amplitudes scaled by 2p (pairs between two sites by 2p_i x 2p_j); full pair energies recomputed
                   (Wang-Frenkel + Yukawa + hydration forcing of the new residues), as mutation_scan_recharged.py
For each pattern: Kish effective sample fraction; reweighted <Rg>, distances 18-63 and 63-102, their squares (for RMS
distances) and FRET efficiencies E = 1/(1 + (r/R0)^6) with R0 5.24 nm (bead positions, no dye linkers); % change of each
mean; the total heating response of each (% per K, calibrated hydration) and its change; 10-block jackknife SEs.
Checks (2026-10-10): Rg from the frames equals the npz Rg (2e-7 nm), chains unbroken by the periodic box; T51-only rows
equal heat_crosses.py (CALVADOS-2) and recharged_traj_pass.py (Glu mimic) single-site values; CLK2 (occupancy) Debye-Hueckel
energy equals CALVADOS' own system built from GAIVI_PAGE4 inputs (4 decimals) and its dV/dT equals a finite difference;
Mpipi phosphate-like Wang-Frenkel/Yukawa energies and the Yukawa dT term equal OpenMpipi with run_openmpipi_page4's patch.
    BENDER_REANALYSIS_DIR=... python page4_phosphoforms.py        (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/page4/phosphoforms.csv, summary.txt (capped; A0: page4_A0_linear)
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import heat_crosses as HC  # noqa: E402
import mutation_scan_recharged as MS  # noqa: E402

TF, TC, fdt, L = HC.TF, MS.TC, HC.fdt, HC.L
UID = 'O60829'
OUT = os.path.join(HC.TLA, 'Analysis_FDT', 'total_response', 'page4')
# BENDER_SCAN_FORCING=capped (2026-10-25): the capped main model (heat_crosses.FORCING). CALVADOS-2 phosphates keep the
# Ser/Thr hydration, so only the baseline forcing changes (capped_frames.term); in Mpipi-Recharged the Glu beads take the
# charged-group coefficient, so the capped forcing of each pattern's sequence is recomputed in full. Tables (TAB): page4
# for the capped main model, page4_A0_linear for the linear A0 forcing (folder layout since 2026-10-25); OUT stays the
# shared page4 folder, which also holds the direct runs (page4_direct.out_dir).
TAB = OUT if HC.FORCING == 'capped' else OUT + '_A0_linear'
PAIRS = {'d18_63': (17, 62), 'd63_102': (62, 101)}
R0 = 5.24                                   # nm, Forster radius used by Kulkarni et al. / Lin et al.
FULL = {'HIPK1': {51: 1.0, 9: 1.0}, 'HIPK1 T51 only': {51: 1.0},
        'CLK2': {7: 1.0, 9: 1.0, 51: 1.0, 85: 1.0, 73: 1.0, 94: 1.0, 71: 1.0, 79: 1.0}, 'S7+S9': {7: 1.0, 9: 1.0}}
OCC = {'HIPK1 (occupancy)': {51: 1.0, 9: 0.35},
       'CLK2 (occupancy)': {7: 1.0, 9: 1.0, 51: 1.0, 85: 1.0, 73: 0.75, 94: 0.6, 71: 0.5, 79: 0.5}}


def observables(xyz):
    rg = np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    obs = {'Rg': rg}
    for k, (a, b) in PAIRS.items():
        r = np.linalg.norm(xyz[:, a] - xyz[:, b], axis=1)
        obs[k] = r
        obs[f'{k}_r2'] = r ** 2
        obs[f'E{k[1:]}'] = 1.0 / (1.0 + (r / R0) ** 6)
    return obs


def evaluate(model, pattern, obs, Gtot, dH, dG, RT):
    nf = len(dH)
    x = -dH / RT
    w = np.exp(x - x.max())
    row = {'model': model, 'pattern': pattern, 'ess': float(w.sum() ** 2 / (w ** 2).sum() / nf)}
    for k, A in obs.items():
        m0 = A.mean()
        f_mean = lambda ix: ((w[ix] * A[ix]).sum() / w[ix].sum() / A[ix].mean() - 1) * 100  # noqa: E731
        f_tot = lambda ix: -HC.wcov(A[ix], Gtot[ix] + dG[ix], w[ix]) / ((w[ix] * A[ix]).sum() / w[ix].sum()) * 100  # noqa: E731
        f_base = lambda ix: -HC.wcov(A[ix], Gtot[ix], np.ones(len(ix))) / A[ix].mean() * 100  # noqa: E731
        f_d = lambda ix: f_tot(ix) - f_base(ix)  # noqa: E731
        allx = np.arange(nf)
        row.update({f'{k}_wt_nm': m0, f'{k}_mean': (w * A).sum() / w.sum(), f'{k}_pct': f_mean(allx), f'{k}_pct_se': HC.jk(f_mean, nf),
                    f'{k}_total_wt': f_base(allx), f'{k}_total': f_tot(allx), f'{k}_dtotal': f_d(allx), f'{k}_dtotal_se': HC.jk(f_d, nf)})
    return row


def calvados():
    d = TF.index_dirs([UID])[UID]
    z = np.load(os.path.join(TF.NPZ, f'{UID}.npz'))
    xyz_all, seq, info = HC.bc.load_calvados(os.path.join(d, f'{UID}.dcd'), os.path.join(d, 'top.pdb'), k=None)
    idx = HC.bc.even_indices(len(xyz_all), L.K_FRAMES)
    if not np.array_equal(info['production_start'] + idx, z['frame_index']):
        raise ValueError('frames differ from the stored FDT frames')
    xyz = xyz_all[idx].astype(float)
    T0, RT, CAL = HC.T0, HC.RT, HC.CAL
    rp = L.RESIDUES_CSV_PATH if os.path.exists(L.RESIDUES_CSV_PATH) else \
        os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv')
    q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params(rp))
    i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
    dl = {k: np.array([v[a] for a in seq]) for k, v in (('tea', TF.DL['charmm36']), ('hyd', TF.DL_CLASSES['hyd']))}
    fpair = HC.G_TEA * 0.5 * (dl['tea'][i] + dl['tea'][j]) + HC.G_HYD * 0.5 * (dl['hyd'][i] + dl['hyd'][j])
    eps, deps, rc, k1, dk1 = CAL['eps_yu'], CAL['deps_yu_dT'], fdt.RC_YU, CAL['kappa'], CAL['dkappa_dT']
    r = np.linalg.norm(xyz[:, i] - xyz[:, j], axis=2)
    u1, u0 = fdt.ah_pair_energy(r, s, np.ones_like(l), shift), fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
    inside = r < rc
    er, ec = np.exp(-k1 * r), np.exp(-k1 * rc)
    Y = np.where(inside, er / r - ec / rc, 0.0)
    dY = deps * Y + eps * dk1 * np.where(inside, ec - er, 0.0)
    V = lambda qq: (qq * eps * Y).sum(1)  # noqa: E731
    dV = lambda qq: (qq * dY).sum(1)  # noqa: E731
    if np.abs(V(q[i] * q[j]) - z['V_DH'].astype(float)).max() > 1e-3 * max(1.0, np.abs(z['V_DH']).max()):
        raise ValueError('V_DH differs from the npz')
    Gtot = z['G_T'].astype(float) + ((u1 - u0) * fpair).sum(1) / RT
    if HC.FORCING == 'capped':
        import capped_frames as CF
        Gtot = z['G_T'].astype(float) + CF.term('C', UID, z['rg'].astype(float))
    obs = observables(xyz)
    rows = []
    for name, pat in {'WT': {}, **FULL, **OCC}.items():
        q2 = q.copy()
        for p, occ in pat.items():
            assert seq[p - 1] in 'ST', (p, seq[p - 1])
            q2[p - 1] += -2.0 * occ
        dH = V(q2[i] * q2[j]) - V(q[i] * q[j])
        dG = -dH / (fdt.R_GAS * T0 ** 2) + (dV(q2[i] * q2[j]) - dV(q[i] * q[j])) / RT
        rows.append(evaluate('CALVADOS-2', name, obs, Gtot, dH, dG, RT))
    return rows


def recharged():
    import mpipi_recharged_fdt as R
    d = TC.rc_dirs()[UID]
    z = np.load(os.path.join(TC.RC_NPZ, f'{UID}.npz'))
    xyz, seq, info = R.load_run(d, UID)
    if not np.array_equal(info['production_start'] + info['stride'] * np.arange(len(xyz)), z['frame_index']):
        raise ValueError('frames differ from the stored FDT frames')
    xyz = xyz.astype(float)
    tb, T0, RT = MS.tables(), MS.T0, MS.RT
    k = MS.kappa_of(T0)
    dk = MS.kappa_of(T0 + 0.5) - MS.kappa_of(T0 - 0.5)
    N, Lb = len(seq), R.box_nm(len(seq))
    ft = {a: MS.G_TEA * TF.DL['charmm36'][a] + MS.G_HYD * TF.DL_CLASSES['hyd'][a] for a in MS.AA}
    i, j = np.triu_indices(N, k=2)
    dv = xyz[:, i] - xyz[:, j]
    dv -= Lb * np.round(dv / Lb)
    r = np.linalg.norm(dv, axis=2)

    def energies(sq, scale=None):
        ix = np.array([tb['ix'][a] for a in sq])
        fv = np.array([ft[a] for a in sq])
        tbs = tb
        W, Yd, Ye, D = MS.pair_terms(r, ix[i], ix[j], tbs, k, fv[i], fv[j])
        if HC.FORCING == 'capped':                  # capped: -1/2 sum_i f_i h(c_i), c_i = -sum_j WF_att_ij (kJ: x 4.184)
            Sres = (HC.CP.incidence(N, i, j).T @ MS.wf_att(r, ix[i], ix[j], tb).T).T
            fc = np.array([4.184 * (MS.G_TEA * HC.CP.F_TEA[a] + MS.G_HYD * HC.CP.F_CHG[a]) for a in sq])
            D = -0.5 * (HC.CP.h(-Sres) @ fc)
        if scale is not None:                       # rescale the Yukawa part pair by pair
            sc = scale[i] * scale[j]
            Aij = tb['A'][ix[j], ix[i]] * sc
            ex = np.where(r < 3.5, np.exp(-k * r), 0.0)
            Yd, Ye = (Aij * ex / r).sum(1), (Aij * ex).sum(1)
        return W, Yd, Ye, D

    W0, Y0, E0, D0 = energies(seq)
    Gtot = z['G_T'].astype(float) + D0 / RT
    obs = observables(xyz)
    m = obs['Rg'].mean()
    rows = []
    variants = [('WT', {}, False)] + [(f'{n} (Glu mimic)', p, False) for n, p in FULL.items()] + \
               [(f'{n} (phosphate-like)', p, True) for n, p in {**FULL, **OCC}.items()]
    for name, pat, phos in variants:
        sq = list(seq)
        for p in pat:
            assert seq[p - 1] in 'ST'
            sq[p - 1] = 'E'
        scale = None
        if phos:
            scale = np.ones(N)
            for p, occ in pat.items():
                scale[p - 1] = 2.0 * occ
        W1, Y1, E1, D1 = energies(''.join(sq), scale)
        dH = (W1 - W0) + (Y1 - Y0)
        dG = -dH / (MS.R_KJ * T0 ** 2) + (-dk * (E1 - E0)) / RT + (D1 - D0) / RT
        rows.append(evaluate('Mpipi-Recharged', name, obs, Gtot, dH, dG, RT))
    base = -HC.fdt._cov(obs['Rg'], Gtot) / m * 100
    return rows, base


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = calvados()
    rr, base_r = recharged()
    df = pd.DataFrame(rows + rr)
    os.makedirs(TAB, exist_ok=True)
    df.to_csv(os.path.join(TAB, 'phosphoforms.csv'), index=False)
    cols = ['model', 'pattern', 'ess', 'Rg_pct', 'Rg_pct_se', 'd18_63_pct', 'd18_63_pct_se', 'd63_102_pct', 'd63_102_pct_se',
            'Rg_total_wt', 'Rg_dtotal', 'Rg_dtotal_se', 'd18_63_dtotal', 'd18_63_dtotal_se']
    lines = ['PAGE4 (O60829) phosphoforms at 300 K, exact reweighting; % change from the unphosphorylated ensemble',
             'measured (Kulkarni et al. 2017), WT / HIPK1 / CLK2: SAXS Rg 36.2 / 34.7 / 49.8 A; E18-63 0.55 / 0.52 / 0.35; '
             'E63-102 0.65 / 0.65 / 0.58; RMS 18-63 56 / 59 / 75 A, 63-102 50 / 50 / 55 A; smFRET-derived Rg 34 / 35 / 43 A',
             f"WT ensembles: d18-63 {df.d18_63_wt_nm.iloc[0]:.2f} nm (C2), {df.d18_63_wt_nm.iloc[-1]:.2f} nm (Recharged); "
             f"d63-102 {df.d63_102_wt_nm.iloc[0]:.2f} / {df.d63_102_wt_nm.iloc[-1]:.2f} nm; Recharged total(Rg) {base_r:.3f} % per K", '',
             df[cols].round(3).to_string(index=False)]
    absc = df[['model', 'pattern', 'ess', 'Rg_mean', 'd18_63_r2_mean', 'd63_102_r2_mean', 'E18_63_mean', 'E63_102_mean',
               'E18_63_total_wt', 'E18_63_total']].copy()
    absc['Rg_A'] = absc.pop('Rg_mean') * 10
    absc['rms18_63_A'] = np.sqrt(absc.pop('d18_63_r2_mean')) * 10
    absc['rms63_102_A'] = np.sqrt(absc.pop('d63_102_r2_mean')) * 10
    absc['dE18_63_per10K'] = absc.E18_63_total / 100 * absc.E18_63_mean * 10
    lines += ['', 'absolute values (WT rows: unreweighted ensemble; E with R0 5.24 nm, bead positions, no dye linkers)',
              absc.round(3).to_string(index=False)]
    open(os.path.join(TAB, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
