"""
phospho_heat_direct.py — direct check of the reweighted change of the heat response on single-site phosphorylation
(registered in heldout_temperature/phospho_heat/REGISTRATION.md before any direct value was computed).

Data: the 240 CALVADOS-2 runs of BENDER-DYN (BENDER-DYN/phospho_direct/300K/PHOS-<UniProt>-<form>-r<k>/), 20 BENDER
proteins, each unmodified (WT) and with one Ser/Thr phosphorylated at an annotated (A), a matched unannotated (U) and a
low-fluctuation (L) site, three independent replicates each, BENDER protocol (300 K, 150 mM, pH 7.5, BENDER run length).
Phosphoresidues B (pSer) and J (pThr): charge -2, Ser/Thr lambda and sigma, phosphate mass
(Scripts/CALVADOS IDP/input_phospho/residues_phospho.csv), as in the PAGE4 phosphoform runs.

Per run: the BENDER frame rule (second half of the trajectory, 1,000 evenly spaced frames, bender_common.load_calvados),
per-frame CALVADOS-2 energies and thermal force G_T with the run's residues and sequence
(TODO_T2_fdt_computation.per_frame_energies, pH 7.5, 150 mM, 300 K), and the calibrated hydration forcing D on the
unmodified sequence (phosphosites keep the Ser/Thr hydration, as in the reweighting), with the frozen strengths
(gamma_TEA 2.0, gamma_hyd 0.358), exactly as page4_direct.prepare_calvados. Responses at 300 K by the
fluctuation-response relation, plain (G_T) and total (G_T + D/RT), % per K, with 10-block jackknife SEs, and mean Rg.

Per site: direct change = mean over the three phosphorylated replicates minus mean over the three WT replicates of the
same protein; SE = sqrt(SE_phos^2 + SE_wt^2) with SE of a replicate mean = max(SD across replicates / sqrt(3),
sqrt(mean jackknife SE^2 / 3)). Compared with the stored reweighted predictions (predictions.csv, frozen copy of
total_response/heat_crosses/sites.csv for these 60 sites).
Criteria (REGISTRATION.md): P1 Spearman rho(predicted, direct change of the total response) > 0 with one-sided
permutation P < 0.05 (10,000 permutations); P2 sign agreement above 50% among sites where both the prediction and the
direct change exceed two SE (one-sided binomial P reported); P3 compaction rule in the direct data, Spearman
rho(direct change of the total response, direct change of Rg) > 0 with one-sided permutation P < 0.05.
Check: the WT replicates against the stored BENDER-MT total response of the same protein (CALVADOS-2, 300 K).
    BENDER_REANALYSIS_DIR=... python phospho_heat_direct.py prepare [--workers 6]    (openmpipi env)
    BENDER_REANALYSIS_DIR=... python phospho_heat_direct.py evaluate
Output: heldout_temperature/phospho_heat/{runs/, runs.csv, sites_direct.csv, summary.txt}
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

BM = '/Volumes/User Homes/Nu Project/BENDER-MT'
DYN = '/Volumes/User Homes/Nu Project/BENDER-DYN/phospho_direct/300K'
INP = os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP/input_phospho')
OUT = os.path.join(BM, 'heldout_temperature', 'phospho_heat')
T0, IONIC, PH, K, NB = 300.0, 0.15, 7.5, 1000, 10
PHOS = {'B': 'S', 'J': 'T'}


def fasta():
    seqs, name = {}, None
    for line in open(os.path.join(INP, 'phospho_direct_300K.fasta')):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split('|')[0]
        elif line:
            seqs[name] = seqs.get(name, '') + line
    return seqs


def runs():
    rows = []
    for name, seq in fasta().items():
        _, uid, form, rep = name.split('-')
        rows.append({'name': name, 'UniProt_ID': uid, 'form': form, 'rep': rep, 'seq': seq,
                     'dir': os.path.join(DYN, name)})
    return pd.DataFrame(rows)


def prepare_one(row):
    import page4_direct as PD           # frozen code, imported read-only for its CALVADOS-2 helpers
    HC, fdt, TF = PD.HC, PD.fdt, PD.TF
    out = os.path.join(OUT, 'runs', f"{row['name']}.npz")
    if os.path.exists(out):
        return row['name'], 'exists'
    try:
        seq = row['seq']
        base = ''.join(PHOS.get(a, a) for a in seq)
        res = pd.read_csv(os.path.join(INP, 'residues_phospho.csv')).set_index('one')
        dcd = os.path.join(row['dir'], f"{row['name']}.dcd")
        xyz, top_seq, info = HC.bc.load_calvados(dcd, os.path.join(row['dir'], 'top.pdb'), k=K, production=True)
        if xyz.shape[1] != len(seq):
            raise ValueError(f'trajectory has {xyz.shape[1]} beads, sequence {len(seq)}')
        if HC.bc.max_bond(xyz) > 1.0:
            xyz = HC.bc.unwrap_chain(xyz, 50.0)
        q, lam, sig, b0 = fdt.calvados_residue_params(seq, res, ph=PH)
        pf = fdt.per_frame_energies(xyz, q, lam, sig, b0, fdt.get_dielectric_params(T0, IONIC))
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        dl = {k: np.array([v[a] for a in base]) for k, v in (('tea', TF.DL['charmm36']), ('hyd', TF.DL_CLASSES['hyd']))}
        fpair = HC.G_TEA * 0.5 * (dl['tea'][i] + dl['tea'][j]) + HC.G_HYD * 0.5 * (dl['hyd'][i] + dl['hyd'][j])
        D = np.zeros(len(xyz))
        step = max(1, int(fdt.MAX_ELEMS // len(i)))
        for a in range(0, len(xyz), step):
            r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            D[a:a + step] = ((fdt.ah_pair_energy(r, s, np.ones_like(l), shift)
                              - fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)) * fpair).sum(1)
        rg = np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        os.makedirs(os.path.dirname(out), exist_ok=True)
        np.savez_compressed(out, rg=rg, G_T=pf['G_T'], H=pf['H'], D300=D, n_frames_total=info['n_frames_total'],
                            stride=info['stride'])
        return row['name'], f"ok {len(xyz)} frames of {info['n_frames_total']}"
    except Exception as ex:  # noqa: BLE001
        return row['name'], f'FAILED {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def prepare(workers):
    r = runs()
    print(len(r), 'runs')
    with mp.get_context('spawn').Pool(workers) as pool:
        for name, msg in pool.imap_unordered(prepare_one, r.to_dict('records')):
            print(name, msg, flush=True)


def jk(f, n):
    blocks = np.array_split(np.arange(n), NB)
    vals = np.array([f(np.concatenate([b for k, b in enumerate(blocks) if k != m])) for m in range(NB)])
    return float(np.sqrt((NB - 1) / NB * ((vals - vals.mean()) ** 2).sum()))


def response(A, G):
    f = lambda ix: -np.cov(A[ix], G[ix], ddof=0)[0, 1] / A[ix].mean() * 100  # noqa: E731
    return f(np.arange(len(A))), jk(f, len(A))


def perm_p(x, y, n=10000, seed=0):
    rho = stats.spearmanr(x, y).statistic
    rng = np.random.default_rng(seed)
    null = np.array([stats.spearmanr(rng.permutation(x), y).statistic for _ in range(n)])
    return float(rho), float((1 + (null >= rho).sum()) / (n + 1))


def evaluate():
    r = runs()
    RT0 = 8.314462618e-3 * T0
    rows = []
    for _, x in r.iterrows():
        z = np.load(os.path.join(OUT, 'runs', f"{x['name']}.npz"))
        tot, tot_se = response(z['rg'], z['G_T'] + z['D300'] / RT0)
        pl, pl_se = response(z['rg'], z['G_T'])
        rows.append({**x.drop(['seq', 'dir']).to_dict(), 'total': tot, 'total_se': tot_se, 'plain': pl, 'plain_se': pl_se,
                     'rg_mean': float(z['rg'].mean())})
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, 'runs.csv'), index=False)

    def form_mean(g, col):
        m = g[col].mean()
        se = max(g[col].std(ddof=1) / np.sqrt(len(g)), np.sqrt((g[f'{col}_se'] ** 2).mean() / len(g))) \
            if f'{col}_se' in g else g[col].std(ddof=1) / np.sqrt(len(g))
        return m, se

    pred = pd.read_csv(os.path.join(OUT, 'predictions.csv'))
    out = []
    for _, p in pred.iterrows():
        code = {'annotated': 'A', 'partner': 'U', 'low': 'L'}[p.role] + p.residue + str(int(p.pos))
        ph = R[(R.UniProt_ID == p.UniProt_ID) & (R.form == code)]
        wt = R[(R.UniProt_ID == p.UniProt_ID) & (R.form == 'WT')]
        if len(ph) != 3 or len(wt) != 3:
            raise ValueError(f'{p.UniProt_ID} {code}: {len(ph)} phosphorylated and {len(wt)} WT runs')
        row = {'UniProt_ID': p.UniProt_ID, 'role': p.role, 'pos': p.pos, 'residue': p.residue,
               'pred_dtotal': p.phospho_dtotal, 'pred_se': p.phospho_dtotal_se, 'pred_dRg_pct': p.phospho_dRg_pct}
        for col in ('total', 'plain'):
            m1, s1 = form_mean(ph, col)
            m0, s0 = form_mean(wt, col)
            row[f'direct_d{col}'], row[f'direct_d{col}_se'] = m1 - m0, float(np.hypot(s1, s0))
        g1, g0 = ph.rg_mean.mean(), wt.rg_mean.mean()
        row['direct_dRg_pct'] = 100 * (g1 / g0 - 1)
        out.append(row)
    S = pd.DataFrame(out)
    S.to_csv(os.path.join(OUT, 'sites_direct.csv'), index=False)

    L = ['Direct check of the reweighted change of the heat response on single-site phosphorylation (CALVADOS-2, 300 K,',
         'total response with the frozen calibrated hydration, % per K)', '']
    rho1, p1 = perm_p(S.pred_dtotal.values, S.direct_dtotal.values)
    L.append(f'P1  Spearman rho(predicted, direct change of total response) = {rho1:+.3f}, one-sided permutation P = '
             f'{p1:.4f} (n = {len(S)}) -> {"PASS" if (rho1 > 0 and p1 < 0.05) else "FAIL"}')
    both = S[(S.pred_dtotal.abs() > 2 * S.pred_se) & (S.direct_dtotal.abs() > 2 * S.direct_dtotal_se)]
    agree = int((np.sign(both.pred_dtotal) == np.sign(both.direct_dtotal)).sum())
    pb = float(stats.binomtest(agree, len(both), 0.5, alternative='greater').pvalue) if len(both) else float('nan')
    L.append(f'P2  sign agreement where both exceed 2 SE: {agree}/{len(both)}, one-sided binomial P = {pb:.4f} -> '
             f'{"PASS" if (len(both) and agree / len(both) > 0.5) else "FAIL"}')
    rho3, p3 = perm_p(S.direct_dRg_pct.values, S.direct_dtotal.values)
    L.append(f'P3  compaction rule in the direct data, Spearman rho(direct dRg, direct d total) = {rho3:+.3f}, one-sided '
             f'permutation P = {p3:.4f} -> {"PASS" if (rho3 > 0 and p3 < 0.05) else "FAIL"}')
    L += ['', 'secondary (not judged)']
    L.append(f'  rho(predicted, direct change of the plain response) = '
             f'{stats.spearmanr(S.pred_dtotal, S.direct_dplain).statistic:+.3f}')
    L.append(f'  rho(predicted dRg, direct dRg) = {stats.spearmanr(S.pred_dRg_pct, S.direct_dRg_pct).statistic:+.3f}')
    L.append(f'  median |predicted| {S.pred_dtotal.abs().median():.3f}, median |direct| {S.direct_dtotal.abs().median():.3f},'
             f' median direct SE {S.direct_dtotal_se.median():.3f}; sites with |direct| > 2 SE: '
             f'{int((S.direct_dtotal.abs() > 2 * S.direct_dtotal_se).sum())}')
    for role, g in S.groupby('role'):
        L.append(f'  {role:9s} n={len(g)} rho(pred, direct) = {stats.spearmanr(g.pred_dtotal, g.direct_dtotal).statistic:+.3f}')
    # check: WT replicates against the stored BENDER-MT total response
    import total_response_evidence as TRE
    D, _ = TRE.load()
    st = D['C'].total
    w = R[R.form == 'WT'].groupby('UniProt_ID').total.mean()
    common = w.index.intersection(st.index)
    L += ['', f'check: WT replicate mean against stored BENDER-MT total ({len(common)} proteins): Spearman '
          f'{stats.spearmanr(w[common], st[common]).statistic:+.3f}, median difference '
          f'{(w[common] - st[common]).median():+.4f} % per K']
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['prepare', 'evaluate'])
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    prepare(a.workers) if a.cmd == 'prepare' else evaluate()
