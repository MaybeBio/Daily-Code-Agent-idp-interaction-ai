"""
make_page4_inputs.py — GAIVI inputs for direct simulations of the PAGE4 phosphoforms (the test of
page4_phosphoforms.py's reweighting predictions and of the heat x phosphorylation prediction).

Forms (sites and mass-spectrometry occupancies from Kulkarni et al. PNAS 2017; mean-field charge -2p as Lin et al. 2018):
  WT        unphosphorylated
  HIPK1occ  T51 (1.0), S9 (0.35)
  HIPK1T51  T51 only
  CLK2occ   S7, S9, T51, T85 (1.0), S73 (0.75), T94 (0.6), T71, S79 (0.5)
  CLK2full  all eight sites at 1.0
  Nmotif    S7, S9 (the N-motif sites; tests the N-terminal loop mechanism)
CALVADOS-2: each phosphosite is a custom residue letter with the Ser/Thr lambda and sigma, charge -2p and the phosphate
            mass added (as page4_phosphoforms.py's charge model); one residues file per form.
Mpipi-Recharged: phosphate-like = Glu bead with its Yukawa amplitudes scaled by 2p (run_openmpipi_page4.py --phospho);
            Glu mimics (scale 1) for HIPK1 and CLK2 as a second representation.
Conditions: heating series 278-320 K at 150 mM (CALVADOS-2 pH 7.5), as the proteome runs; and the measured conditions at
            295 K (temperature not stated in the papers; room temperature assumed): smFRET buffer (20 mM Tris pH 7.8,
            100 mM NaOAc, 5 mM MgCl2; ionic strength ~0.128 M) and SAXS buffer (100 mM KPi pH 7.0, 10 mM NaN3; ~0.19 M),
            for WT, HIPK1occ and CLK2occ.
    python make_page4_inputs.py [--out ~/Downloads/Scripts/GAIVI_PAGE4]
Writes: calvados/{forms.fasta, residues_<form>.csv, runs.tsv}, mpipi/{forms.fasta, runs.tsv}, forms.tsv
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402

import pandas as pd  # noqa: E402

SEQ = 'MSARVRSRSRGRGDGQEAPDVVAFVAPGESQQEEPPTDNQDIEPGQEREGTPPIEERKVEGDCQEMDLEKTRSERGDGSDVKEKTPPNPKHAKTKEAGDGQP'
CLK2 = {7: 1.0, 9: 1.0, 51: 1.0, 85: 1.0, 73: 0.75, 94: 0.6, 71: 0.5, 79: 0.5}
FORMS = {'WT': {}, 'HIPK1occ': {51: 1.0, 9: 0.35}, 'HIPK1T51': {51: 1.0}, 'CLK2occ': CLK2,
         'CLK2full': {p: 1.0 for p in CLK2}, 'Nmotif': {7: 1.0, 9: 1.0}}
GLU_MIMIC = {'HIPK1glu': {51: 1.0, 9: 1.0}, 'CLK2glu': {p: 1.0 for p in CLK2}}
HPO3 = 79.98
TEMPS = [278, 288, 300, 310, 320]
HEAT = {'ionic': 0.15, 'ph': 7.5}
EXP = {'smFRET': {'ionic': 0.128, 'ph': 7.8}, 'SAXS': {'ionic': 0.19, 'ph': 7.0}}
EXP_T, EXP_FORMS = 295, ['WT', 'HIPK1occ', 'CLK2occ']
RESIDUES = os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv')


def calvados_form(name, sites, base):
    seq, rows, letter = list(SEQ), [], {}
    for p, occ in sorted(sites.items()):
        aa = SEQ[p - 1]
        assert aa in 'ST', (name, p, aa)
        key = (aa, round(occ, 3))
        if key not in letter:
            full = {('S', 1.0): 'B', ('T', 1.0): 'J'}           # same letter for a full phosphate in every form
            letter[key] = full.get(key) or [c for c in 'OUXZ' if c not in letter.values()][0]
            r = base.loc[aa].copy()
            r['q'] = -2.0 * occ
            r['MW'] = r['MW'] + HPO3 * occ
            r['three'] = {'S': 'SEP', 'T': 'TPO'}[aa] if occ == 1.0 else f'{aa}{round(20 * occ):02d}'
            rows.append(r.rename(letter[key]))
        seq[p - 1] = letter[key]
    res = pd.concat([base, pd.DataFrame(rows)]) if rows else base.copy()
    return ''.join(seq), res, {f'{SEQ[p - 1]}{p}': f'{letter[(SEQ[p - 1], round(o, 3))]} (q {-2 * o:+.2f})' for p, o in sorted(sites.items())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.expanduser('~/Downloads/Scripts/GAIVI_PAGE4'))
    a = ap.parse_args()
    cdir, mdir = os.path.join(a.out, 'calvados'), os.path.join(a.out, 'mpipi')
    os.makedirs(cdir, exist_ok=True)
    os.makedirs(mdir, exist_ok=True)
    base = pd.read_csv(RESIDUES).set_index('one')
    manifest, crun, mrun = [], [], []
    with open(os.path.join(cdir, 'forms.fasta'), 'w') as fc:
        for f, sites in FORMS.items():
            seq, res, how = calvados_form(f, sites, base)
            res.reset_index().rename(columns={'index': 'one'})[['one', 'three', 'MW', 'lambdas', 'sigmas', 'q', 'bondlength']] \
                .to_csv(os.path.join(cdir, f'residues_{f}.csv'), index=False)
            fc.write(f'>{f}\n{seq}\n')
            manifest.append({'model': 'CALVADOS-2', 'form': f, 'sites': '; '.join(f'{k} {v}' for k, v in how.items()) or 'none'})
    for f in FORMS:
        conds = [('heat', T, HEAT) for T in TEMPS] + ([(c, EXP_T, v) for c, v in EXP.items()] if f in EXP_FORMS else [])
        for tag, T, c in conds:
            crun.append({'name': f'PAGE4-{f}-I{round(c["ionic"] * 1000)}-pH{round(c["ph"] * 10)}', 'form': f, 'set': tag,
                         'temp': T, 'ionic': c['ionic'], 'ph': c['ph'], 'residues': f'residues_{f}.csv'})
    with open(os.path.join(mdir, 'forms.fasta'), 'w') as fm:
        for f, sites in {**FORMS, **GLU_MIMIC}.items():
            seq = ''.join('E' if (i + 1) in sites else c for i, c in enumerate(SEQ))
            fm.write(f'>{f}\n{seq}\n')
            spec = '' if f in GLU_MIMIC else ','.join(f'{p}:{2 * o:g}' for p, o in sorted(sites.items()))
            manifest.append({'model': 'Mpipi-Recharged', 'form': f,
                             'sites': (('Glu mimic at ' if f in GLU_MIMIC else 'Glu bead, Yukawa x2p at ') + ', '.join(f'{SEQ[p - 1]}{p}' for p in sorted(sites))) if sites else 'none',
                             'phospho': spec})
            conds = [('heat', T, HEAT) for T in TEMPS] + ([(c, EXP_T, v) for c, v in EXP.items()] if f in EXP_FORMS else [])
            for tag, T, c in conds:
                mrun.append({'name': f'PAGE4-{f}-I{round(c["ionic"] * 1000)}', 'form': f, 'set': tag, 'temp': T,
                             'ionic': c['ionic'], 'phospho': spec or '-'})
    pd.DataFrame(crun).to_csv(os.path.join(cdir, 'runs.tsv'), sep='\t', index=False)
    pd.DataFrame(mrun).to_csv(os.path.join(mdir, 'runs.tsv'), sep='\t', index=False)
    pd.DataFrame(manifest).to_csv(os.path.join(a.out, 'forms.tsv'), sep='\t', index=False)
    print(f'CALVADOS-2: {len(crun)} runs; Mpipi-Recharged: {len(mrun)} runs; written to {a.out}')
    print(pd.DataFrame(manifest).to_string(index=False))


if __name__ == '__main__':
    main()
