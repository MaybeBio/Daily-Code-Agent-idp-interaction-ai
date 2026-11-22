"""
phospho_heat_direct_capped.py — the direct check of phospho_heat_direct.py under the capped hydration forcing (the main
model from 2026-10-24), on the same 240 CALVADOS-2 runs and 60 sites.

Per run (same frames as phospho_heat_direct.prepare_one: BENDER frame rule, 1,000 frames of the second half, unwrapped):
per-residue contact numbers c_i = -sum_j S_ij / eps_AH with the run's own pair shapes (phosphoresidues B, J have the
Ser/Thr lambda and sigma), and the capped forcing per unit strength
    D_chg = -eps/2 sum_i f_chg(a_i) h(c_i),  D_tea = -eps/2 sum_i f_tea(a_i) h(c_i),  h(c) = 2 (1 - exp(-c / 2))
with the coefficients of the unmodified residue (phosphosites keep the Ser/Thr hydration, as in the reweighting).
Check per run: Rg equals the stored phospho_heat/runs/<name>.npz to 1e-6 and the linear forcing recomputed on the same
frames reproduces its D300 to 1e-6 relative. Totals: G_T + (g_h D_chg + g_tea D_tea) / RT with the capped strengths
(capped_frames.strengths). Predictions: the capped reweighting of the same sites (total_response/heat_crosses/
sites.csv; the A0 run is in heat_crosses_A0_linear). Same criteria P1-P3 as phospho_heat_direct.py; WT replicates against the capped proteome total.
    BENDER_REANALYSIS_DIR=... python phospho_heat_direct_capped.py prepare [--workers 4]
    BENDER_REANALYSIS_DIR=... python phospho_heat_direct_capped.py evaluate
Output: heldout_temperature/phospho_heat/capped/{runs/, runs.csv, sites_direct.csv, summary.txt}
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import phospho_heat_direct as PHD  # noqa: E402

OUT = os.path.join(PHD.OUT, 'capped')
HCX = os.path.join(PHD.BM, 'Traj Level Analysis', 'Analysis_FDT', 'total_response', 'heat_crosses')


def prepare_one(row):
    import capped_proteome as CP
    import page4_direct as PD
    HC, fdt, TF = PD.HC, PD.fdt, PD.TF
    out = os.path.join(OUT, 'runs', f"{row['name']}.npz")
    if os.path.exists(out):
        return row['name'], 'exists'
    try:
        seq = row['seq']
        base = ''.join(PHD.PHOS.get(a, a) for a in seq)
        res = pd.read_csv(os.path.join(PHD.INP, 'residues_phospho.csv')).set_index('one')
        dcd = os.path.join(row['dir'], f"{row['name']}.dcd")
        xyz, _, info = HC.bc.load_calvados(dcd, os.path.join(row['dir'], 'top.pdb'), k=PHD.K, production=True)
        if HC.bc.max_bond(xyz) > 1.0:
            xyz = HC.bc.unwrap_chain(xyz, 50.0)
        st = np.load(os.path.join(PHD.OUT, 'runs', f"{row['name']}.npz"))
        rg = np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        if rg.shape != st['rg'].shape or np.abs(rg - st['rg']).max() > 1e-6 * rg.mean():
            raise ValueError('frames differ from the stored run')
        q, lam, sig, _ = fdt.calvados_residue_params(seq, res, ph=PHD.PH)
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        N, nf = len(seq), len(xyz)
        Mi = CP.incidence(N, i, j)
        dl = {k: np.array([v[a] for a in base]) for k, v in (('tea', TF.DL['charmm36']), ('hyd', TF.DL_CLASSES['hyd']))}
        fpair = HC.G_TEA * 0.5 * (dl['tea'][i] + dl['tea'][j]) + HC.G_HYD * 0.5 * (dl['hyd'][i] + dl['hyd'][j])
        Sres, Dlin = np.zeros((nf, N)), np.zeros(nf)
        step = max(1, int(fdt.MAX_ELEMS // len(i)))
        for a in range(0, nf, step):
            r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            S = fdt.ah_pair_energy(r, s, np.ones_like(l), shift) - fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
            Sres[a:a + step] = (Mi.T @ S.T).T
            Dlin[a:a + step] = (S * fpair).sum(1)
        dev = np.abs(Dlin - st['D300']).max() / max(np.abs(st['D300']).max(), 1e-12)
        if dev > 1e-6:
            raise ValueError(f'linear forcing differs from the stored D300 ({dev:.1e})')
        hc = CP.h(-Sres / fdt.EPS_AH)
        e0 = fdt.EPS_AH
        D_chg = -0.5 * e0 * (hc @ np.array([CP.F_CHG[a] for a in base]))
        D_tea = -0.5 * e0 * (hc @ np.array([CP.F_TEA[a] for a in base]))
        os.makedirs(os.path.dirname(out), exist_ok=True)
        np.savez_compressed(out, rg=rg, G_T=st['G_T'], D_chg_cap=D_chg, D_tea_cap=D_tea, lin_check=dev)
        return row['name'], f'ok {nf} frames, linear check {dev:.1e}'
    except Exception as ex:  # noqa: BLE001
        return row['name'], f'FAILED {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def prepare(workers):
    r = PHD.runs()
    print(len(r), 'runs', flush=True)
    with mp.get_context('spawn').Pool(workers) as pool:
        for name, msg in pool.imap_unordered(prepare_one, r.to_dict('records')):
            print(name, msg, flush=True)


def evaluate():
    import capped_frames as CF
    gh, gt = CF.strengths()['C']
    RT0 = 8.314462618e-3 * PHD.T0
    r = PHD.runs()
    rows = []
    for _, x in r.iterrows():
        z = np.load(os.path.join(OUT, 'runs', f"{x['name']}.npz"))
        tot, tot_se = PHD.response(z['rg'], z['G_T'] + (gh * z['D_chg_cap'] + gt * z['D_tea_cap']) / RT0)
        pl, pl_se = PHD.response(z['rg'], z['G_T'])
        rows.append({**x.drop(['seq', 'dir']).to_dict(), 'total': tot, 'total_se': tot_se, 'plain': pl, 'plain_se': pl_se,
                     'rg_mean': float(z['rg'].mean())})
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, 'runs.csv'), index=False)

    def form_mean(g, col):
        m = g[col].mean()
        se = max(g[col].std(ddof=1) / np.sqrt(len(g)), np.sqrt((g[f'{col}_se'] ** 2).mean() / len(g)))
        return m, se

    frozen = pd.read_csv(os.path.join(PHD.OUT, 'predictions.csv'))
    cap = pd.read_csv(os.path.join(HCX, 'sites.csv')).drop_duplicates(['UniProt_ID', 'pos']).set_index(['UniProt_ID', 'pos'])
    out = []
    for _, p in frozen.iterrows():
        c = cap.loc[(p.UniProt_ID, p.pos)]
        code = {'annotated': 'A', 'partner': 'U', 'low': 'L'}[p.role] + p.residue + str(int(p.pos))
        ph = R[(R.UniProt_ID == p.UniProt_ID) & (R.form == code)]
        wt = R[(R.UniProt_ID == p.UniProt_ID) & (R.form == 'WT')]
        if len(ph) != 3 or len(wt) != 3:
            raise ValueError(f'{p.UniProt_ID} {code}: {len(ph)} phosphorylated and {len(wt)} WT runs')
        row = {'UniProt_ID': p.UniProt_ID, 'role': p.role, 'pos': p.pos, 'residue': p.residue,
               'pred_dtotal': c.phospho_dtotal, 'pred_se': c.phospho_dtotal_se, 'pred_dRg_pct': c.phospho_dRg_pct}
        for col in ('total', 'plain'):
            m1, s1 = form_mean(ph, col)
            m0, s0 = form_mean(wt, col)
            row[f'direct_d{col}'], row[f'direct_d{col}_se'] = m1 - m0, float(np.hypot(s1, s0))
        row['direct_dRg_pct'] = 100 * (ph.rg_mean.mean() / wt.rg_mean.mean() - 1)
        out.append(row)
    S = pd.DataFrame(out)
    S.to_csv(os.path.join(OUT, 'sites_direct.csv'), index=False)
    L = ['Direct check of the reweighted change of the heat response on single-site phosphorylation, capped forcing',
         f'(CALVADOS-2, 300 K, g_h {gh:.4f}, g_tea {gt:.1f}; % per K)', '']
    rho1, p1 = PHD.perm_p(S.pred_dtotal.values, S.direct_dtotal.values)
    L.append(f'P1  Spearman rho(predicted, direct change of total response) = {rho1:+.3f}, one-sided permutation P = {p1:.4f} '
             f'(n = {len(S)})')
    both = S[(S.pred_dtotal.abs() > 2 * S.pred_se) & (S.direct_dtotal.abs() > 2 * S.direct_dtotal_se)]
    agree = int((np.sign(both.pred_dtotal) == np.sign(both.direct_dtotal)).sum())
    pb = float(stats.binomtest(agree, len(both), 0.5, alternative='greater').pvalue) if len(both) else float('nan')
    L.append(f'P2  sign agreement where both exceed 2 SE: {agree}/{len(both)}, one-sided binomial P = {pb:.4f}')
    rho3, p3 = PHD.perm_p(S.direct_dRg_pct.values, S.direct_dtotal.values)
    L.append(f'P3  compaction rule in the direct data, rho(direct dRg, direct d total) = {rho3:+.3f}, one-sided P = {p3:.4f}')
    L.append(f'  rho(predicted dRg, direct dRg) = {stats.spearmanr(S.pred_dRg_pct, S.direct_dRg_pct).statistic:+.3f}')
    L.append(f'  median |predicted| {S.pred_dtotal.abs().median():.3f}, median |direct| {S.direct_dtotal.abs().median():.3f}, '
             f'median direct SE {S.direct_dtotal_se.median():.3f}; |direct| > 2 SE at '
             f'{int((S.direct_dtotal.abs() > 2 * S.direct_dtotal_se).sum())} sites')
    import total_response_evidence as TRE
    assert TRE.FORCING == 'capped'
    D, _ = TRE.load()
    wt = R[R.form == 'WT'].groupby('UniProt_ID').total.mean()
    ref = D['C'].total.reindex(wt.index)
    L.append(f'  check: WT replicate mean against the capped proteome total, Spearman {stats.spearmanr(wt, ref).statistic:.3f}, '
             f'median difference {float((wt - ref).median()):+.3f} (n = {len(wt)})')
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('step', choices=['prepare', 'evaluate'])
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    prepare(a.workers) if a.step == 'prepare' else evaluate()
