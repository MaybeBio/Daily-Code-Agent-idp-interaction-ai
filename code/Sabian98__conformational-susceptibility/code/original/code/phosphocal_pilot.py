"""
phosphocal_pilot.py — the direct CALVADOS-2 PAGE4 runs at 300 K reweighted to the Phospho-CALVADOS phosphoresidues of
Rauh et al. (Biophys. J. 2026; external/rauh2025_phospho_calvados/residues.csv), a pilot for the reviewer question why the
paper's charge-only phosphosites were used (2026-10-27).

Our runs: phosphosites keep the stickiness and size of Ser/Thr and carry charge -2 x occupancy (page4_direct.py).
Phospho-CALVADOS: pSer/pThr stickiness of Ser/Thr - 0.37, larger sigma, charge -(1 + 1/(1 + 10^(pKa - pH))). Partial
occupancy o is applied as a mean field, lambda, sigma and charge each (1 - o) x Ser/Thr + o x phosphoresidue.
Variants: 'charge' (only the charge, their charge-only baseline) and 'phospho' (their fitted model).
Per form: frame weights w = exp(-(H_new - H_old) / RT) from TODO_T2_fdt_computation.per_frame_energies; <Rg> and its
change from WT; the calibrated total response at 300 K, -100 Cov_w(Rg, G_T,new + D / RT) / <Rg>_w, with the saturating
hydration forcing D of page4_capped_fdt (contacts on the unmodified chain, as there); 10-block jackknife SEs; Kish
effective sample fraction. Check: the stored H is reproduced by the old parameters, and the unchanged variant gives
weights of one.
    BENDER_REANALYSIS_DIR=... python phosphocal_pilot.py      (openmpipi env)
Output: heldout_temperature/capped_forcing/page4/phosphocal/{pilot.csv, pilot.txt}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import page4_capped_fdt as CF  # noqa: E402
import page4_capped_predictions as PC  # noqa: E402

PD, SF, fdt = PC.PD, PC.SF, PC.PD.fdt
RES_PC = os.path.join(PD.BM, 'external/rauh2025_phospho_calvados/residues.csv')
OUT = os.path.join(PC.OUT, 'phosphocal')
FORMS = ['WT', 'HIPK1occ', 'CLK2occ', 'CLK2full']
T, NB = 300, 10
RT = fdt.R_GAS * T
PCODE = {'S': 'B', 'T': 'O'}          # Phospho-CALVADOS letters for pSer and pThr


def run_row(form):
    r = PD.runs_table('calvados')
    r = r[(r.form == form) & (r.temp == T) & (r.set == 'heat')]
    assert len(r) == 1, form
    return r.iloc[0]


def frames(z, form):
    name = f'PAGE4-{form}{PC.FS.TAG["calvados"]}'
    d = os.path.join(PD.MODELS['calvados']['root'], f'{T}K', name)
    xyz_all, _, _ = PD.HC.bc.load_calvados(os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'), k=None,
                                           production=False)
    xyz = xyz_all[z['frame_index']].astype(float)
    return PD.HC.bc.unwrap_chain(xyz, 50.0) if PD.HC.bc.max_bond(xyz) > 1.0 else xyz


def new_residues(res, seq, base, ph, variant):
    pc = pd.read_csv(RES_PC).set_index('one')
    new = res.copy()
    seen = {}
    for a, b in zip(seq, base):
        if a not in PD.PHOS:
            continue
        assert seen.setdefault(a, b) == b, f'letter {a} stands for both Ser and Thr'
        occ = -float(res.loc[a, 'q']) / 2.0
        p = pc.loc[PCODE[b]]
        q = -(1.0 + 1.0 / (1.0 + 10 ** (float(p['pKa']) - ph)))
        new.loc[a, 'q'] = occ * q
        if variant == 'phospho':
            for col in ('lambdas', 'sigmas'):
                new.loc[a, col] = (1 - occ) * float(res.loc[b, col]) + occ * float(p[col])
    return new


def wstats(rg, force, u, keep=None):
    k = np.ones(len(rg), bool) if keep is None else keep
    w = np.exp(-(u[k] - u[k].min()))
    w /= w.sum()
    m = (w * rg[k]).sum()
    cov = (w * (rg[k] - m) * (force[k] - (w * force[k]).sum())).sum()
    return m, -100.0 * cov / m, 1.0 / (w ** 2).sum() / k.sum()


def with_se(rg, force, u):
    m, tot, ess = wstats(rg, force, u)
    n = len(rg)
    jk = np.array([wstats(rg, force, u, PC.FS.keep_mask(n, b))[:2] for b in range(NB)])
    f = (NB - 1) / NB
    se = np.sqrt(f * ((jk - jk.mean(0)) ** 2).sum(0))
    return m, se[0], tot, se[1], ess


def main():
    g_tea, g_h = PC.capped_gamma('calvados')
    rows = []
    for form in FORMS:
        row = run_row(form)
        z, _ = PC.run_frames('calvados', form, T)
        xyz = frames(z, form)
        rg = z['obs_Rg'] * 10.0
        ct, ch = CF.dEdT('calvados', z, SF.HYDSET['S1'], PC.CSAT)
        res = pd.read_csv(os.path.join(PD.MODELS['calvados']['gaivi'], row['residues'])).set_index('one')
        ph, seq, base = float(row['ph']), row['seq'], row['wt_seq']
        prm = fdt.get_dielectric_params(float(T), float(row['ionic']))
        old = fdt.per_frame_energies(xyz, *fdt.calvados_residue_params(seq, res, ph=ph), prm)
        dev = np.abs(old['H'] - z['H']).max()
        if dev > 1e-3:
            raise ValueError(f'{form}: stored H not reproduced ({dev:.1e} kJ/mol)')
        for variant in ('ours', 'charge', 'phospho'):
            if variant == 'ours':
                new = old
            else:
                new = fdt.per_frame_energies(xyz, *fdt.calvados_residue_params(seq, new_residues(res, seq, base, ph, variant),
                                                                                ph=ph), prm)
            u = (new['H'] - old['H']) / RT
            if form == 'WT' and np.abs(u).max() > 1e-9:
                raise ValueError('WT has no phosphosites but the energies changed')
            m, mse, tot, tse, ess = with_se(rg, new['G_T'] + (g_tea * ct + g_h * ch) / RT, u)
            rows.append({'form': form, 'variant': variant, 'Rg_A': m, 'Rg_se': mse, 'total': tot, 'total_se': tse,
                         'ess': ess, 'H_check': dev})
        print(form, 'done', flush=True)
    df = pd.DataFrame(rows)
    wt = df[df.form == 'WT'].iloc[0]
    df['dRg_pct'] = 100 * (df.Rg_A / wt.Rg_A - 1)
    df['dtotal'] = df.total - wt.total
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, 'pilot.csv'), index=False)
    L = ['PAGE4 CALVADOS-2 direct runs, 300 K, 150 mM, pH 7.5, reweighted to Phospho-CALVADOS phosphoresidues',
         'ours = charge -2 x occupancy, Ser/Thr lambda and sigma; charge = their pKa charge only; phospho = their fitted model',
         'calibrated total response of Rg (% per K, saturating forcing), 10-block SEs; SAXS (0.19 M): HIPK1 -4%, CLK2 +38%', '']
    for r in df.itertuples():
        L.append(f'{r.form:9s} {r.variant:8s} Rg {r.Rg_A:6.2f} +- {r.Rg_se:.2f} A ({r.dRg_pct:+5.1f}% vs WT)  total {r.total:+.3f}'
                 f' +- {r.total_se:.3f} (change {r.dtotal:+.3f})  ESS {r.ess:.2f}')
    L.append(f'\ncheck: stored H reproduced to {df.H_check.max():.1e} kJ/mol; WT weights all one')
    open(os.path.join(OUT, 'pilot.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
