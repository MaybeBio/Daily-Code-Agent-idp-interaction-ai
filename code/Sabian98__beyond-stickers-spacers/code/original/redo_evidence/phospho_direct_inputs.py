"""
phospho_direct_inputs.py — phospho_direct_plan.md (hash 5a69b002): select annotated / ESS-matched Ser/Thr pairs, write the
frozen reweighting predictions, and generate the GAIVI inputs (CALVADOS-2, BENDER protocol, 3 replicates per form) from
the PAGE4 machinery in ~/Downloads/Scripts/GAIVI_PAGE4. Selection uses annotation, residue type, protein length,
distance to annotated sites and the effective sample fraction only, never a rewiring quantity.

  ~/.venvs/finches/bin/python -I phospho_direct_inputs.py
      ->  gaivi/phospho_direct/{selection.csv, frozen_predictions.npz},  ~/Downloads/Scripts/GAIVI_PHOSPHO/
"""
import os
import random
import re
import shutil

import numpy as np
import pandas as pd

DYN = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
MTS = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT/phosphosite_test/sites.csv'
PAGE4 = os.path.expanduser('~/Downloads/Scripts/GAIVI_PAGE4')
GOUT = os.path.expanduser('~/Downloads/Scripts/GAIVI_PHOSPHO')
OUT = os.path.join(DYN, 'gaivi', 'phospho_direct')
BASE = os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv')
SEED, N_PROT, N_REP = 20261012, 20, 3


def select():
    S = pd.read_csv(os.path.join(DYN, 'phospho_network', 'sites_network.csv'))
    M = pd.read_csv(MTS)
    M = M[M.model == 'CALVADOS-2'][['UniProt_ID', 'pos', 'annotated', 'evidence']]
    S = S.merge(M, on=['UniProt_ID', 'pos'], how='left', validate='one_to_one')
    S['annotated'] = S['annotated'].astype(str) == 'True'
    S['logess'] = np.log(S['ess_frac_net'])
    pairs = []
    for uid, g in S.groupby('UniProt_ID'):
        if g.N.iloc[0] > 200 or not g.annotated.any() or g.annotated.all():
            continue
        ann_pos = set(g.loc[g.annotated, 'pos'])
        best = None
        for _, a in g[g.annotated & (g.ess_frac_net >= 0.2)].sort_values('pos').iterrows():
            cand = g[(~g.annotated) & (g.ess_frac_net >= 0.2) & (g.residue == a.residue)]
            cand = cand[[min(abs(p - q) for q in ann_pos) >= 4 for p in cand.pos]]
            if cand.empty:
                continue
            cand = cand.assign(gap=(cand.logess - a.logess).abs()).sort_values(['gap', 'pos'])
            c = cand.iloc[0]
            if c.gap <= 0.10 and (best is None or c.gap < best[2]):
                best = (a, c, c.gap)
        if best is not None:
            pairs.append(best)
    qual = sorted(p[0].UniProt_ID for p in pairs)
    chosen = sorted(random.Random(SEED).sample(qual, N_PROT)) if len(qual) > N_PROT else qual
    rows = []
    for a, c, gap in pairs:
        if a.UniProt_ID not in chosen:
            continue
        for role, r in (('annotated', a), ('partner', c)):
            rows.append(dict(UniProt_ID=r.UniProt_ID, N=int(r.N), role=role, pos=int(r.pos), residue=r.residue,
                             evidence=r.evidence if role == 'annotated' else '', ess=float(r.ess_frac_net),
                             logess_gap=float(gap), R_hub_rw=float(r.R_hub), R_mag_rw=float(r.R_mag),
                             R_shape_rw=float(r.R_shape), dRg_rw=float(r.dRg_pct_rw)))
    return pd.DataFrame(rows).sort_values(['UniProt_ID', 'role']), len(qual)


def main():
    os.makedirs(OUT, exist_ok=True)
    sel, n_qual = select()
    sel.to_csv(os.path.join(OUT, 'selection.csv'), index=False)
    prof = {}
    seqs = {}
    for uid, g in sel.groupby('UniProt_ID'):
        z = np.load(os.path.join(DYN, 'phospho_network', 'profiles', f'{uid}.npz'))
        seqs[uid] = str(z['seq'])
        pos = list(z['positions'])
        for _, r in g.iterrows():
            k = pos.index(r.pos)
            prof[f'{uid}_{r.role}_ddeg'] = z['ddeg'][k]
            prof[f'{uid}_{r.role}_docc'] = z['docc'][k]
    np.savez_compressed(os.path.join(OUT, 'frozen_predictions.npz'), **prof)

    cdir = os.path.join(GOUT, 'calvados')
    os.makedirs(cdir, exist_ok=True)
    base = pd.read_csv(BASE)
    page4 = pd.read_csv(os.path.join(PAGE4, 'calvados', 'residues_CLK2full.csv'))
    res = pd.concat([base, page4[page4.one.isin(['B', 'J'])]], ignore_index=True)
    assert res.one.is_unique
    res.to_csv(os.path.join(cdir, 'residues_phospho.csv'), index=False)
    runs, forms = [], []
    with open(os.path.join(cdir, 'forms.fasta'), 'w') as fh:
        for uid, g in sel.groupby('UniProt_ID'):
            s = seqs[uid]
            fl = [(f'{uid}_WT', 'WT', s)]
            for _, r in g.iterrows():
                i = r.pos - 1
                assert s[i] == r.residue in 'ST', (uid, r.pos, s[i])
                tag = ('A' if r.role == 'annotated' else 'U') + f'{r.residue}{r.pos}'
                fl.append((f'{uid}_{tag}', tag, s[:i] + {'S': 'B', 'T': 'J'}[r.residue] + s[i + 1:]))
            for form, tag, sq in fl:
                fh.write(f'>{form}\n{sq}\n')
                forms.append(dict(form=form, UniProt_ID=uid, tag=tag, N=len(sq)))
                for k in range(1, N_REP + 1):
                    name = f'PHOS-{uid}-{tag}-r{k}'
                    seed = int.from_bytes(name.encode(), 'little') % 2_000_000_000 + 1
                    runs.append(dict(name=name, form=form, set='direct', temp=300, ionic=0.15, ph=7.5,
                                     residues='residues_phospho.csv', seed=seed))
    pd.DataFrame(runs).to_csv(os.path.join(cdir, 'runs.tsv'), sep='\t', index=False)
    fseq = {f['form']: None for f in forms}
    for line in open(os.path.join(cdir, 'forms.fasta')).read().split('>')[1:]:
        h, sq = line.split('\n', 1)
        fseq[h.strip()] = sq.replace('\n', '')
    with open(os.path.join(GOUT, 'labels_input_phospho.fasta'), 'w') as fh:      # one entry per run, for hub_labels_gaivi.py
        for r in runs:
            fh.write(f">{r['name']}|{r['form'].split('_')[0]}\n{fseq[r['form']]}\n")
    pd.DataFrame(forms).to_csv(os.path.join(GOUT, 'forms.tsv'), sep='\t', index=False)

    # prepare script: the PAGE4 one plus --seed and time factor 1 by default
    src = open(os.path.join(PAGE4, 'calvados', 'prepare_page4.py')).read()
    src = src.replace('prepare_page4.py — CALVADOS-2 setup for one PAGE4 phosphoform run.',
                      'prepare_phospho.py — CALVADOS-2 setup for one single-site phosphorylation run (from prepare_page4.py;\n'
                      'phospho_direct_plan.md). Adds --seed (used when the CALVADOS version accepts random_number_seed).')
    src = src.replace("parser.add_argument('--time-factor', type=int, default=3)",
                      "parser.add_argument('--time-factor', type=int, default=1)\nparser.add_argument('--seed', type=int, default=None)")
    src = src.replace("restart='checkpoint', frestart='restart.chk', verbose=True)",
                      "restart='checkpoint', frestart='restart.chk', verbose=True)\nif args.seed is not None:\n    cfg['random_number_seed'] = args.seed")
    src = src.replace("except Exception:                                # newer CALVADOS versions reject runtime=0; steps then set the length\n    cfg.pop('runtime')\n    config = Config(**cfg)",
                      "except Exception:                                # newer CALVADOS versions reject runtime=0; older ones may reject the seed\n    cfg.pop('runtime')\n    try:\n        config = Config(**cfg)\n    except Exception:\n        cfg.pop('random_number_seed', None)          # OpenMM then picks a random seed per context\n        cfg['runtime'] = 0\n        try:\n            config = Config(**cfg)\n        except Exception:\n            cfg.pop('runtime')\n            config = Config(**cfg)\nprint(f\"  seed: {cfg.get('random_number_seed', 'OpenMM random per context')}\")")
    assert '--seed' in src and 'random_number_seed' in src
    open(os.path.join(cdir, 'prepare_phospho.py'), 'w').write(src)

    # submit script: the PAGE4 one, CALVADOS-2 only, new job names, outroot, time factor 1, sets 'direct', seed column
    sh = open(os.path.join(PAGE4, 'submit_page4.sh')).read()
    reps = [('submit_page4.sh', 'submit_phospho.sh'), ('PAGE4 phosphoform runs on GAIVI (CALVADOS-2 and Mpipi-Recharged)',
             'single-site phosphorylation runs on GAIVI (CALVADOS-2; phospho_direct_plan.md)'),
            ('PAGE4_cal', 'PHOS_cal'), ('PAGE4_mp', 'PHOS_mp'), ('PAGE4_resub_', 'PHOS_resub_'),
            ('OUT_CAL="/data/taseefr/page4"', 'OUT_CAL="/data/taseefr/phospho_direct"'),
            ('OUT_MP="/general/taseefr/page4"', 'OUT_MP="/general/taseefr/phospho_direct_unused"'),
            ('SETS="heat,smFRET,SAXS"', 'SETS="direct"'), ('TF=3', 'TF=1'),
            ('RES=\\$(echo "\\$LINE" | cut -f7)', 'RES=\\$(echo "\\$LINE" | cut -f7); SEED=\\$(echo "\\$LINE" | cut -f8)'),
            ('python "$SCRIPT_DIR/calvados/prepare_page4.py"', 'python "$SCRIPT_DIR/calvados/prepare_phospho.py"'),
            ('--outdir "$OUT" --time-factor $TF\npython "\\$RUNDIR/run.py"', '--outdir "$OUT" --time-factor $TF --seed "\\$SEED"\npython "\\$RUNDIR/run.py"')]
    for a, b in reps:
        assert a in sh, a
        sh = sh.replace(a, b)
    sh = sh.replace('[[ "$MODEL" =~ ^(calvados|mpipi|both)$ ]] || { echo "usage: bash submit_page4.sh calvados|mpipi|both [options]"; exit 1; }',
                    '[[ "$MODEL" == calvados ]] || { echo "usage: bash submit_phospho.sh calvados [options] (CALVADOS-2 only)"; exit 1; }')
    sh = sh.replace('[[ "$MODEL" =~ ^(calvados|mpipi|both)$ ]] || { echo "usage: bash submit_phospho.sh calvados|mpipi|both [options]"; exit 1; }',
                    '[[ "$MODEL" == calvados ]] || { echo "usage: bash submit_phospho.sh calvados [options] (CALVADOS-2 only)"; exit 1; }')
    open(os.path.join(GOUT, 'submit_phospho.sh'), 'w').write(sh)
    n_runs = len(runs)
    readme = f"""# Single-site phosphorylation runs on GAIVI

Direct CALVADOS-2 simulations for `BENDER-DYN/redo_evidence/phospho_direct_plan.md`. They test whether annotated
phosphosites rewire the contact network more than unannotated Ser/Thr of the same protein matched on the fluctuation
of the electrostatic potential at the site (the confound found by the charge-model reweighting). Generated by
`BENDER-DYN/redo_evidence/phospho_direct_inputs.py`; regenerate there rather than editing by hand.

- {sel.UniProt_ID.nunique()} proteins (of {n_qual} that qualified), each with three forms: unmodified (`_WT`), the annotated site
  phosphorylated (`_A<res><pos>`) and the matched partner phosphorylated (`_U<res><pos>`). `forms.tsv` lists them.
- Three independent replicates per form: {n_runs} runs in `calvados/runs.tsv`.
- Phosphosites as in the PAGE4 runs: B = pSer, J = pThr (charge -2, Ser/Thr lambda and sigma, phosphate mass), defined
  in `calvados/residues_phospho.csv`.
- BENDER protocol: 300 K, 150 mM, pH 7.5, time factor 1 (BENDER run length, 10 fs step, first half discarded).
- Seeds: column 8 of `runs.tsv`, passed as `random_number_seed` when the installed CALVADOS accepts it; otherwise OpenMM
  draws a random seed per run, which also gives independent replicates. The prepare log prints which one was used.

## Commands on GAIVI

```bash
rsync -av ~/Downloads/Scripts/GAIVI_PHOSPHO/ <gaivi>:~/GAIVI_PHOSPHO/   # from the Mac
source ~/calvados_env/bin/activate                                       # on GAIVI
cd ~/GAIVI_PHOSPHO
bash submit_phospho.sh calvados                                          # {n_runs} GPU array tasks
bash submit_phospho.sh calvados                                          # later: resubmits only unfinished runs
```

Output: `/data/taseefr/phospho_direct/300K/<name>/` (`<name>.dcd`, `top.pdb`), analysis in `.../300K/data/<name>/`.

## Hub labels on GAIVI (after the runs finish)

The validated label script of the hub study, on this output layout (the analysis script is hashed before this is run):

```bash
python3 <path to BENDER-DYN gaivi copy>/hub_labels_gaivi.py --layout openmm --root /data/taseefr/phospho_direct/300K \
    --fasta ~/GAIVI_PHOSPHO/labels_input_phospho.fasta --out ~/labels_phospho_direct
```

Copy `labels_phospho_direct/` back to `BENDER-DYN/redo_evidence/` (per_residue.csv and per_protein.csv are small).
Nodes GPU4, GPU6 and GPU8 are excluded by default, as in the PAGE4 script.
"""
    open(os.path.join(GOUT, 'README.md'), 'w').write(readme)
    print(f'qualifying proteins {n_qual}; selected {sel.UniProt_ID.nunique()}; runs {n_runs}')
    print(sel[['UniProt_ID', 'N', 'role', 'pos', 'residue', 'evidence', 'ess', 'logess_gap', 'R_hub_rw', 'dRg_rw']].to_string(index=False))


if __name__ == '__main__':
    main()
