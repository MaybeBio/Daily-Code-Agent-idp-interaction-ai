"""
make_era_calibrated_inputs.py — GAIVI inputs for direct simulations of the ER-NTD Ser118 forms at 277 K under the
calibrated models, replacing the reweighting of the plain 277 K runs (era_direct.py), whose effective sample fraction is
about 0.1 (CALVADOS-2) and 0.01 (Mpipi-Recharged) at 23 K below the reference temperature. Written 2026-10-24;
exploratory, nothing frozen is changed.

Construct, forms and phosphate model as make_era_inputs.py (SNA + P03372 1-184, 187 residues; sequences and residues
files read from GAIVI_ERA, checked there against SASBDB). Runs at 277 K, pH 7.4, in the two SAXS buffers of Du et al.:
  0.195 M   WT, S118D, pS118      (WT and S118D measured here)
  0.095 M   WT, pS118             (pS118 measured here; WT added so that phosphorylation and salt can be separated)
Variants A0 (frozen calibrated model) and S1 (charged-group forcing of Amendment 3), as make_pnt_calibrated_inputs.py,
and S1cap (added 2026-10-24): the charged-group form with both terms capped at c_sat = 2 contacts per residue and the
hydrophilic strength refitted on Wuttke (saturating_forcing.py; 2.893 CALVADOS-2, 0.927 Mpipi-Recharged). S1cap is
many-body, so it is a CustomGBForce: run_calvados_cal.py adds it to the CALVADOS-2 run (plain residues file), and
run_openmpipi_cal.py adds it when given a capped weight file (header w_kJmol, '# c_sat 2'):
  CALVADOS-2       lambda_i(277 K) = lambda_i + gamma_TEA (-dG_E,i) + gamma_h [du_h,i(277 K) - du_h,i(T0)] in a residues
                   file per form and variant; the phosphoserine (residue B) takes Ser's change, as the reweighting reads
                   the phosphosite as Ser (era_direct.py, page4_direct.py)
  Mpipi-Recharged  delta_i(277 K) per residue (kJ/mol) added by run_openmpipi_cal.py; the phosphoserine is a Glu bead with
                   its own particle type (run_openmpipi_page4.py --phospho 121:2) and takes Glu's weight, as the reweighting
                   reads the run sequence (Glu at the phosphosite)
Checked by era_cal_check.py on the plain 277 K ER-NTD frames (direct energy difference = gamma . E of the reweighting).
Run scripts, run length (time factor 3, as GAIVI_ERA) and submit logic are GAIVI_PAGE4's (copied in on GAIVI).
    python make_era_calibrated_inputs.py [--out $CS_DATA/Scripts/GAIVI_ERACAL]
Writes: calvados/{forms.fasta, residues_<form>_<variant>_277K.csv, residues_<form>_plain.csv, w_S1cap_277K.tsv, runs.tsv,
        run_calvados_cal.py}, mpipi/{forms.fasta, delta_<variant>_277K.tsv, w_S1cap_277K.tsv, runs.tsv,
        residues_CALVADOS2.csv, run_openmpipi_cal.py}, forcing.tsv, submit_eracal.sh, README.md
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402

import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import make_pnt_calibrated_inputs as MP  # noqa: E402
import saturating_forcing as SF  # noqa: E402

ERA = os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_ERA'))
T = 277
RUNS = [('WT', 0.195), ('S118D', 0.195), ('pS118', 0.195), ('WT', 0.095), ('pS118', 0.095)]
PH = 7.4
SITE = 121                      # construct position of UniProt Ser118
AA = MP.AA
CSAT = 2.0


def capped_gamma(model):
    t = pd.read_csv(os.path.join(SF.OUT, f'saturating_{model}.csv'))
    row = t[(t.variant == 'S1') & (t.c_sat == CSAT)]
    assert len(row) == 1
    return SF.HX.S1_TEA[model], float(row.gamma_h.iloc[0])


def capped_weights(model):
    """w_i = gamma_TEA w_TEA + gamma_h w_hyd at 277 K per residue letter (CALVADOS-2 lambda units; Mpipi-Recharged kJ/mol)."""
    gt, gh = capped_gamma(model)
    wt, wh = SF.weights(model, AA, float(T), 'DEKR')
    return {a: gt * wt[k] + gh * wh[k] for k, a in enumerate(AA)}, gt, gh


def read_fasta(path):
    out = {}
    for block in open(path).read().split('>')[1:]:
        h, *s = block.split('\n')
        out[h.strip().split()[0]] = ''.join(s).strip()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_ERACAL')))
    a = ap.parse_args()
    cdir, mdir = os.path.join(a.out, 'calvados'), os.path.join(a.out, 'mpipi')
    os.makedirs(cdir, exist_ok=True)
    os.makedirs(mdir, exist_ok=True)
    cfa, mfa = read_fasta(f'{ERA}/calvados/forms.fasta'), read_fasta(f'{ERA}/mpipi/forms.fasta')
    forms = sorted({f for f, _ in RUNS})
    for f in forms:
        assert len(cfa[f]) == len(mfa[f]) == 187
    assert cfa['pS118'][SITE - 1] == 'B' and mfa['pS118'][SITE - 1] == 'E' and cfa['WT'][SITE - 1] == 'S'
    eraruns = pd.read_csv(f'{ERA}/mpipi/runs.tsv', sep='\t')
    spec = eraruns[eraruns.form == 'pS118'].phospho.iloc[0]
    assert spec == f'{SITE}:2', spec
    with open(os.path.join(cdir, 'forms.fasta'), 'w') as fc, open(os.path.join(mdir, 'forms.fasta'), 'w') as fm:
        for f in forms:
            fc.write(f'>{f}\n{cfa[f]}\n')
            fm.write(f'>{f}\n{mfa[f]}\n')
    crun, mrun, rows = [], [], []
    for variant in ('A0', 'S1'):
        dl = MP.forcing('calvados', variant, float(T))
        dl['B'] = dl['S']                                     # phosphoserine keeps Ser's hydration change
        for f in forms:
            base = pd.read_csv(f'{ERA}/calvados/residues_{f}.csv')
            r = base.copy()
            missing = set(r['one']) - set(dl)
            assert not missing, missing
            r['lambdas'] = r['lambdas'] + r['one'].map(dl)
            r.to_csv(os.path.join(cdir, f'residues_{f}_{variant}_{T}K.csv'), index=False)
        de = {k: 4.184 * v for k, v in MP.forcing('recharged', variant, float(T)).items()}
        gt, gh = MP.gam('recharged', variant)
        with open(os.path.join(mdir, f'delta_{variant}_{T}K.tsv'), 'w') as fh:
            fh.write(f'# {variant} at {T} K: delta_i in kJ/mol on the attractive Wang-Frenkel shape (gamma_TEA {gt}, '
                     f'gamma_h {gh} on {MP.HYD[variant]}); phosphate-like Glu beads take E\nresidue\tdelta_kJmol\n')
            for k in AA:
                fh.write(f'{k}\t{de[k]:.10g}\n')
        for f, ionic in RUNS:
            tag = f'I{round(ionic * 1000):03d}'
            crun.append({'name': f'ERcal-{f}-{tag}-{variant}-pH74', 'form': f, 'set': variant, 'temp': T, 'ionic': ionic,
                         'ph': PH, 'residues': f'residues_{f}_{variant}_{T}K.csv', 'cap': '-'})
            mrun.append({'name': f'ERcal-{f}-{tag}-{variant}', 'form': f, 'set': variant, 'temp': T, 'ionic': ionic,
                         'phospho': spec if f == 'pS118' else '-', 'delta': f'delta_{variant}_{T}K.tsv'})
        for k in sorted(dl):
            rows.append({'variant': variant, 'T': T, 'residue': k, 'CALVADOS2_dlambda': dl[k],
                         'MpipiRecharged_delta_kJmol': de.get(k, de['E'] if k == 'B' else None)})
    # S1cap: plain residues files, capped weights per residue, runs
    wc, gtc, ghc = capped_weights('calvados')
    wc['B'] = wc['S']                                        # phosphoserine keeps Ser's hydration, as the reweighting
    with open(os.path.join(cdir, f'w_S1cap_{T}K.tsv'), 'w') as fh:
        fh.write(f'# S1cap at {T} K: w_i = gamma_TEA w_TEA + gamma_h w_hyd in lambda units (gamma_TEA {gtc}, gamma_h {ghc:.6g} on DEKR), '
                 f'phosphoserine B takes S\n# c_sat {CSAT:g}\nresidue\tw_lambda\n')
        for k in sorted(wc):
            fh.write(f'{k}\t{wc[k]:.10g}\n')
    wm, gtm, ghm = capped_weights('recharged')
    with open(os.path.join(mdir, f'w_S1cap_{T}K.tsv'), 'w') as fh:
        fh.write(f'# S1cap at {T} K: w_i = gamma_TEA w_TEA + gamma_h w_hyd in kJ/mol (gamma_TEA {gtm}, gamma_h {ghm:.6g} on DEKR), '
                 f'phosphate-like Glu beads take E\n# c_sat {CSAT:g}\nresidue\tw_kJmol\n')
        for k in AA:
            fh.write(f'{k}\t{wm[k]:.10g}\n')
    for f in forms:
        shutil.copy2(f'{ERA}/calvados/residues_{f}.csv', os.path.join(cdir, f'residues_{f}_plain.csv'))
    for f, ionic in RUNS:
        tag = f'I{round(ionic * 1000):03d}'
        crun.append({'name': f'ERcal-{f}-{tag}-S1cap-pH74', 'form': f, 'set': 'S1cap', 'temp': T, 'ionic': ionic, 'ph': PH,
                     'residues': f'residues_{f}_plain.csv', 'cap': f'w_S1cap_{T}K.tsv'})
        mrun.append({'name': f'ERcal-{f}-{tag}-S1cap', 'form': f, 'set': 'S1cap', 'temp': T, 'ionic': ionic,
                     'phospho': spec if f == 'pS118' else '-', 'delta': f'w_S1cap_{T}K.tsv'})
    for k in sorted(wc):
        rows.append({'variant': 'S1cap', 'T': T, 'residue': k, 'CALVADOS2_w_lambda': wc[k],
                     'MpipiRecharged_w_kJmol': wm.get(k, wm['E'] if k == 'B' else None), 'c_sat': CSAT})
    shutil.copy2(os.path.join(HERE, 'gaivi_templates', 'run_calvados_cal.py'), os.path.join(cdir, 'run_calvados_cal.py'))
    shutil.copy2(MP.MH.RESIDUES, os.path.join(mdir, 'residues_CALVADOS2.csv'))
    shutil.copy2(os.path.join(HERE, 'gaivi_templates', 'run_openmpipi_cal.py'), os.path.join(mdir, 'run_openmpipi_cal.py'))
    C, M = pd.DataFrame(crun), pd.DataFrame(mrun)
    assert not C.duplicated(['name', 'temp']).any() and not M.duplicated(['name', 'temp']).any()
    C.to_csv(os.path.join(cdir, 'runs.tsv'), sep='\t', index=False)
    M.to_csv(os.path.join(mdir, 'runs.tsv'), sep='\t', index=False)
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'forcing.tsv'), sep='\t', index=False)
    sub = open(os.path.join(HERE, 'gaivi_templates', 'submit_page4.sh')).read()
    rep = [
        ('submit_page4.sh — PAGE4 phosphoform runs on GAIVI (CALVADOS-2 and Mpipi-Recharged)',
         'submit_eracal.sh — ER-NTD Ser118 forms at 277 K under the calibrated models (A0, S1) on GAIVI;\n'
         '# generated by BENDER-MT/code/make_era_calibrated_inputs.py from GAIVI_PAGE4/submit_page4.sh'),
        ('submit_page4.sh', 'submit_eracal.sh'),
        ('SETS="heat,smFRET,SAXS"', 'SETS="A0,S1,S1cap"'),
        ('/data/taseefr/page4', '/data/taseefr/eracal'), ('/general/taseefr/page4', '/general/taseefr/eracal'),
        ('PAGE4_cal', 'ERACAL_cal'), ('PAGE4_mp', 'ERACAL_mp'), ('PAGE4_resub_', 'ERACAL_resub_'),
        ('# the 36 CALVADOS-2 runs in calvados/runs.tsv', f'# the {len(C)} CALVADOS-2 runs in calvados/runs.tsv'),
        ('# the 46 Mpipi-Recharged runs in mpipi/runs.tsv', f'# the {len(M)} Mpipi-Recharged runs in mpipi/runs.tsv'),
        ('#   options: --sets heat            only the 278-320 K series (or smFRET, SAXS; comma-separated)',
         '#   options: --sets A0              some variants only (A0 frozen model, S1 charged-group, S1cap capped; default all)'),
        ('(default 3 = 300 ns for PAGE4)', '(default 3, as GAIVI_ERA, about 420 ns)'),
        ('TEMP=\\$(echo "\\$LINE" | cut -f4); IONIC=\\$(echo "\\$LINE" | cut -f5); PHOS=\\$(echo "\\$LINE" | cut -f6)',
         'TEMP=\\$(echo "\\$LINE" | cut -f4); IONIC=\\$(echo "\\$LINE" | cut -f5); PHOS=\\$(echo "\\$LINE" | cut -f6); '
         'DELTA=\\$(echo "\\$LINE" | cut -f7)'),
        ('python "$SCRIPT_DIR/mpipi/run_openmpipi_page4.py" --name',
         'python "$SCRIPT_DIR/mpipi/run_openmpipi_cal.py" --delta "$SCRIPT_DIR/mpipi/\\$DELTA" --name'),
        ('PH=\\$(echo "\\$LINE" | cut -f6); RES=\\$(echo "\\$LINE" | cut -f7)',
         'PH=\\$(echo "\\$LINE" | cut -f6); RES=\\$(echo "\\$LINE" | cut -f7); CAP=\\$(echo "\\$LINE" | cut -f8)'),
        ('python "\\$RUNDIR/run.py" --path "\\$RUNDIR"',
         'grep -q "fresidues: .*\\$RES\\$" "\\$RUNDIR/components.yaml" || { echo "components.yaml does not use \\$RES"; exit 1; }\n'
         'if [[ "\\$CAP" == "-" ]]; then python "\\$RUNDIR/run.py" --path "\\$RUNDIR"; '
         'else python "$SCRIPT_DIR/calvados/run_calvados_cal.py" --cap "$SCRIPT_DIR/calvados/\\$CAP" --path "\\$RUNDIR"; fi'),
    ]
    for x, y in rep:
        assert x in sub, x
        sub = sub.replace(x, y)
    hook = 'submit() {   # $1 model dir, $2 outroot'
    assert hook in sub
    sub = sub.replace(hook, '# the run scripts are GAIVI_PAGE4\'s, unchanged; copy them in on first use\n'
                      'PAGE4_DIR="${PAGE4_DIR:-$HOME/GAIVI_PAGE4}"\n'
                      'for f in calvados/prepare_page4.py mpipi/run_openmpipi_page4.py; do\n'
                      '    [[ -f "$SCRIPT_DIR/$f" ]] || cp "$PAGE4_DIR/$f" "$SCRIPT_DIR/$f" || { echo "missing $PAGE4_DIR/$f"; exit 1; }\n'
                      'done\n\n' + hook)
    for bad in ('36 CALVADOS-2', '46 Mpipi', 'smFRET, SAXS', '/page4', 'PAGE4_cal', 'PAGE4_mp'):
        assert bad not in sub, bad
    open(os.path.join(a.out, 'submit_eracal.sh'), 'w').write(sub)
    os.chmod(os.path.join(a.out, 'submit_eracal.sh'), 0o755)
    open(os.path.join(a.out, 'README.md'), 'w').write(
        '# GAIVI_ERACAL: ER-NTD Ser118 forms at 277 K under the calibrated models\n\n'
        'Generated by `BENDER-MT/code/make_era_calibrated_inputs.py`. Construct SNA + ERalpha 1-184 (187 residues), the SAXS '
        'buffers of Du et al. 2025 (0.195 M: WT, S118D, pS118; 0.095 M: WT, pS118), with the hydration forcing built into '
        'the model instead of reweighting.\n\n'
        f'- {len(C)} CALVADOS-2 runs (`calvados/runs.tsv`): A0 and S1 by residues files with lambda_i(277 K), S1cap by the added force\n'
        f'- {len(M)} Mpipi-Recharged runs (`mpipi/runs.tsv`): `run_openmpipi_cal.py` adds the forcing (delta file per '
        'variant) to the unchanged `run_openmpipi_page4.py`; pS118 keeps `--phospho 121:2`\n'
        '- variants: A0 (frozen calibrated model), S1 (charged-group forcing) and S1cap (charged-group, capped at 2 contacts '
        'per residue; `run_calvados_cal.py` / capped weight file); values in `forcing.tsv`\n\n'
        'On GAIVI (activate the calvados env first):\n\n'
        '```bash\nscp -r GAIVI_ERACAL taseefr@<gaivi>:~/\ncd ~/GAIVI_ERACAL && bash submit_eracal.sh both\n```\n\n'
        'The run scripts are copied from `~/GAIVI_PAGE4` on first use.\n'
        'Outputs: `/data/taseefr/eracal` (CALVADOS-2) and `/general/taseefr/eracal` (Mpipi-Recharged).\n')
    print(f'CALVADOS-2: {len(C)} runs; Mpipi-Recharged: {len(M)} runs; written to {a.out}')
    print(C.to_string(index=False))
    print(M.to_string(index=False))


if __name__ == '__main__':
    main()
