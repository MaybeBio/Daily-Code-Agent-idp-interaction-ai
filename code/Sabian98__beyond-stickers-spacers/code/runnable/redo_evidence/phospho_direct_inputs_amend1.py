"""
phospho_direct_inputs_amend1.py — phospho_direct_plan.md (hash 5a69b002) with Amendment 1 (hash 807887ae): reselect 20
proteins that have an annotated site, an ESS-matched unannotated partner and a low-fluctuation unannotated site, write
the frozen reweighting predictions, and regenerate the CALVADOS IDP inputs (4 forms x 3 replicates = 240 runs).
Uses annotation, residue type, position, protein length and the effective sample fraction only.

  ~/.venvs/finches/bin/python -I phospho_direct_inputs_amend1.py
      ->  gaivi/phospho_direct/{selection_amend1.csv, frozen_predictions_amend1.npz},
          $BSS_DATA/Scripts/CALVADOS IDP/input_phospho/{phospho_direct_300K.fasta, residues_phospho.csv, README.md}
"""
import os
import random

import numpy as np
import pandas as pd

DYN = (__import__('os').environ['BSS_REPO'] + '/redo_evidence')
MTS = (__import__('os').environ['BSS_DATA'] + '/BENDER-MT/Traj Level Analysis/Analysis_FDT/phosphosite_test/sites.csv')
CIDP = os.path.expanduser((__import__('os').environ['BSS_DATA'] + '/Scripts/CALVADOS IDP'))
PAGE4_RES = os.path.expanduser((__import__('os').environ['BSS_DATA'] + '/Scripts/GAIVI_PAGE4/calvados/residues_CLK2full.csv'))
OUT = os.path.join(DYN, 'gaivi', 'phospho_direct')
SEED, N_PROT, N_REP, MATCH, LOW = 20261012, 20, 3, 0.10, 0.20


def far(df, ref):
    if df.empty:
        return df
    return df.loc[np.array([min(abs(p - q) for q in ref) >= 4 for p in df.pos], dtype=bool)]


def select():
    S = pd.read_csv(os.path.join(DYN, 'phospho_network', 'sites_network.csv'))
    M = pd.read_csv(MTS)
    M = M[M.model == 'CALVADOS-2'][['UniProt_ID', 'pos', 'annotated', 'evidence']]
    S = S.merge(M, on=['UniProt_ID', 'pos'], how='left', validate='one_to_one')
    S['annotated'] = S['annotated'].astype(str) == 'True'
    S['logess'] = np.log(S['ess_frac_net'])
    trip = {}
    for uid, g in S.groupby('UniProt_ID'):
        if g.N.iloc[0] > 200 or not g.annotated.any() or g.annotated.all():
            continue
        ann = set(g.loc[g.annotated, 'pos'])
        best = None
        for _, a in g[g.annotated & (g.ess_frac_net >= 0.2)].sort_values('pos').iterrows():
            cand = far(g[(~g.annotated) & (g.ess_frac_net >= 0.2) & (g.residue == a.residue)], ann)
            if cand.empty:
                continue
            cand = cand.assign(gap=(cand.logess - a.logess).abs()).sort_values(['gap', 'pos'])
            c = cand.iloc[0]
            if c.gap > MATCH:
                continue
            low = far(g[(~g.annotated) & (g.pos != c.pos) & (g.residue == a.residue)], ann | {c.pos})
            low = low[low.logess - a.logess >= LOW] if not low.empty else low
            if low.empty:
                continue
            lo = low.sort_values(['logess', 'pos'], ascending=[False, True]).iloc[0]
            if best is None or c.gap < best[3]:
                best = (a, c, lo, c.gap)
        if best is not None:
            trip[uid] = best
    qual = sorted(trip)
    chosen = sorted(random.Random(SEED).sample(qual, N_PROT)) if len(qual) > N_PROT else qual
    rows = []
    for uid in chosen:
        a, c, lo, gap = trip[uid]
        for role, r in (('annotated', a), ('partner', c), ('low', lo)):
            rows.append(dict(UniProt_ID=uid, N=int(r.N), role=role, pos=int(r.pos), residue=r.residue,
                             evidence=r.evidence if role == 'annotated' else '', ess=float(r.ess_frac_net),
                             logess_vs_annotated=float(r.logess - a.logess), R_hub_rw=float(r.R_hub),
                             R_mag_rw=float(r.R_mag), R_shape_rw=float(r.R_shape), dRg_rw=float(r.dRg_pct_rw)))
    return pd.DataFrame(rows), len(qual)


def main():
    sel, n_qual = select()
    sel.to_csv(os.path.join(OUT, 'selection_amend1.csv'), index=False)
    prof, seqs = {}, {}
    for uid, g in sel.groupby('UniProt_ID'):
        z = np.load(os.path.join(DYN, 'phospho_network', 'profiles', f'{uid}.npz'))
        seqs[uid] = str(z['seq'])
        pos = list(z['positions'])
        for _, r in g.iterrows():
            k = pos.index(r.pos)
            prof[f'{uid}_{r.role}_ddeg'] = z['ddeg'][k]
            prof[f'{uid}_{r.role}_docc'] = z['docc'][k]
    np.savez_compressed(os.path.join(OUT, 'frozen_predictions_amend1.npz'), **prof)

    idir = os.path.join(CIDP, 'input_phospho')
    os.makedirs(idir, exist_ok=True)
    base = pd.read_csv(os.path.join(CIDP, 'input', 'residues_CALVADOS2.csv'))
    p4 = pd.read_csv(PAGE4_RES)
    res = pd.concat([base, p4[p4.one.isin(['B', 'J'])]], ignore_index=True)
    assert res.one.is_unique
    res.to_csv(os.path.join(idir, 'residues_phospho.csv'), index=False)
    n = 0
    tagmap = {'annotated': 'A', 'partner': 'U', 'low': 'L'}
    with open(os.path.join(idir, 'phospho_direct_300K.fasta'), 'w') as fh:
        for uid, g in sel.groupby('UniProt_ID'):
            s = seqs[uid]
            forms = [('WT', s)]
            for _, r in g.sort_values('role').iterrows():
                i = r.pos - 1
                assert s[i] == r.residue and r.residue in 'ST', (uid, r.pos, s[i])
                forms.append((f'{tagmap[r.role]}{r.residue}{r.pos}', s[:i] + {'S': 'B', 'T': 'J'}[r.residue] + s[i + 1:]))
            for tag, sq in forms:
                for k in range(1, N_REP + 1):
                    fh.write(f'>PHOS-{uid}-{tag}-r{k}|{uid}\n{sq}\n')
                    n += 1
    readme = open(os.path.join(idir, 'README.md')).read()
    readme = readme.replace(
        "Inputs for `BENDER-DYN/redo_evidence/phospho_direct_plan.md` (hash 5a69b002): 20 BENDER proteins, each as the\n"
        "unmodified chain (`-WT-`), with one UniProt-annotated phosphosite phosphorylated (`-A<res><pos>-`), and with a matched\n"
        "unannotated Ser/Thr phosphorylated (`-U<res><pos>-`); three independent replicates each (`-r1`, `-r2`, `-r3`), 180 runs.\n"
        "Generated by `BENDER-DYN/redo_evidence/phospho_direct_inputs.py`; regenerate there rather than editing by hand.",
        "Inputs for `BENDER-DYN/redo_evidence/phospho_direct_plan.md` (hash 5a69b002) with Amendment 1 (hash 807887ae): 20 BENDER\n"
        "proteins, each as the unmodified chain (`-WT-`), with one UniProt-annotated phosphosite phosphorylated (`-A<res><pos>-`),\n"
        "with an unannotated Ser/Thr of matched potential fluctuation phosphorylated (`-U<res><pos>-`), and with an unannotated\n"
        "Ser/Thr of low potential fluctuation phosphorylated (`-L<res><pos>-`); three independent replicates each (`-r1`, `-r2`,\n"
        f"`-r3`), {n} runs. Generated by `BENDER-DYN/redo_evidence/phospho_direct_inputs_amend1.py`; regenerate there rather\n"
        "than editing by hand.")
    assert f'{n} runs' in readme
    open(os.path.join(idir, 'README.md'), 'w').write(readme)
    print(f'qualifying proteins {n_qual}; selected {sel.UniProt_ID.nunique()}; runs {n}')
    print(sel[['UniProt_ID', 'N', 'role', 'pos', 'residue', 'ess', 'logess_vs_annotated', 'R_hub_rw']].to_string(index=False))


if __name__ == '__main__':
    main()
