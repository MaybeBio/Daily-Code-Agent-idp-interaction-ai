"""
make_page4_calibrated_inputs.py — GAIVI inputs for direct simulations of the PAGE4 phosphoforms under the calibrated
model (charged-group hydrophilic and TEA forcing, each residue saturating at c_sat = 2 contacts; the main model since
2026-10-25), replacing the reweighting of the plain heating series. In Mpipi-Recharged that reweighting kept 0-3%
effective sample at the ends of the range (page4_capped_fdt.py), so the heating curves of the calibrated model were only
approximate; CALVADOS-2 is included so that both force fields are checked directly.

Forms and phosphate models as make_page4_inputs.py (GAIVI_PAGE4, regenerated locally and checked against the sequences
of the runs): CALVADOS-2 WT, HIPK1occ, HIPK1T51, CLK2occ, CLK2full, Nmotif; Mpipi-Recharged the same plus the Glu
mimics HIPK1glu and CLK2glu. Temperatures 278, 288, 310 and 320 K at 150 mM (CALVADOS-2 pH 7.5): the forcing vanishes
at 300 K by construction, so the existing plain 300 K runs are already runs of the calibrated model. The forcing is the
saturating CustomGBForce of run_calvados_cal.py / run_openmpipi_cal.py, given a weight file per temperature with
  w_i(T) = gamma_TEA w_TEA,i(T) + gamma_h w_hyd,i(T)   (saturating_forcing.weights, charged-group set DEKR; strengths of
                                                         page4_capped_predictions.capped_gamma: 2.0 and 2.893 CALVADOS-2,
                                                         0.4 and 0.927 Mpipi-Recharged)
in lambda units (CALVADOS-2) or kJ/mol (Mpipi-Recharged), and '# c_sat 2'. As in the reweighting, CALVADOS-2
phosphoresidues (B, J and the partial-occupancy letters) take the weight of their Ser or Thr, and Mpipi-Recharged
phosphate-like Glu beads take Glu's.
Checked by page4_cal_check.py on frames of the plain runs at each temperature (direct energy difference = the saturating
energy of the analysis). Run scripts, run length (time factor 3, 300 ns) and submit logic are GAIVI_PAGE4's (copied in
on GAIVI from ~/GAIVI_PAGE4).
    python make_page4_calibrated_inputs.py [--out $CS_DATA/Scripts/GAIVI_PAGE4CAL]
Writes: calvados/{forms.fasta, residues_<form>_plain.csv, w_cal_<form>_<T>K.tsv (per form: the custom phosphoresidue
        letters are assigned form by form), runs.tsv, run_calvados_cal.py},
        mpipi/{forms.fasta, w_cal_<T>K.tsv, runs.tsv, residues_CALVADOS2.csv, run_openmpipi_cal.py}, forcing.tsv,
        submit_page4cal.sh, README.md
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
import make_era_calibrated_inputs as ME  # noqa: E402
import make_page4_inputs as MPI  # noqa: E402
import make_pnt_calibrated_inputs as MP  # noqa: E402
import saturating_forcing as SF  # noqa: E402

PAGE4 = os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_PAGE4'))
TEMPS = [278, 288, 310, 320]
CSAT = 2.0
AA = MP.AA


def weights(model, T):
    gt, gh = ME.capped_gamma(model)
    wt, wh = SF.weights(model, AA, float(T), 'DEKR')
    return {a: gt * wt[k] + gh * wh[k] for k, a in enumerate(AA)}, gt, gh


def phospho_base(residues_csv):
    """Custom CALVADOS letters of a residues file mapped to the Ser or Thr they phosphorylate (three-letter code SEP/TPO
    for a full phosphate, S<occ>/T<occ> for partial occupancy)."""
    r = pd.read_csv(residues_csv)
    out = {}
    for _, x in r[~r.one.isin(list(AA))].iterrows():
        three = str(x.three)
        base = 'S' if three.startswith(('SEP', 'S')) else 'T' if three.startswith(('TPO', 'T')) else None
        assert base is not None, (residues_csv, x.one, three)
        out[x.one] = base
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_PAGE4CAL')))
    a = ap.parse_args()
    cdir, mdir = os.path.join(a.out, 'calvados'), os.path.join(a.out, 'mpipi')
    os.makedirs(cdir, exist_ok=True)
    os.makedirs(mdir, exist_ok=True)
    cfa, mfa = ME.read_fasta(f'{PAGE4}/calvados/forms.fasta'), ME.read_fasta(f'{PAGE4}/mpipi/forms.fasta')
    assert set(cfa) == set(MPI.FORMS) and set(mfa) == set(MPI.FORMS) | set(MPI.GLU_MIMIC)
    mruns0 = pd.read_csv(f'{PAGE4}/mpipi/runs.tsv', sep='\t')
    spec = mruns0.drop_duplicates('form').set_index('form').phospho.to_dict()
    letters = {}                                         # per form: the custom letters are assigned form by form
    for f in cfa:
        letters[f] = phospho_base(f'{PAGE4}/calvados/residues_{f}.csv')
        assert set(cfa[f]) <= set(AA) | set(letters[f]), f
        shutil.copy2(f'{PAGE4}/calvados/residues_{f}.csv', os.path.join(cdir, f'residues_{f}_plain.csv'))
    with open(os.path.join(cdir, 'forms.fasta'), 'w') as fc:
        for f, s in cfa.items():
            fc.write(f'>{f}\n{s}\n')
    with open(os.path.join(mdir, 'forms.fasta'), 'w') as fm:
        for f, s in mfa.items():
            fm.write(f'>{f}\n{s}\n')
    crun, mrun, rows = [], [], []
    for T in TEMPS:
        w0, gtc, ghc = weights('calvados', T)
        for f in cfa:                                    # one file per form: letter meanings differ between forms
            wc = dict(w0)
            for k, base in letters[f].items():
                wc[k] = w0[base]                         # phosphoresidues keep the hydration of their Ser or Thr
            with open(os.path.join(cdir, f'w_cal_{f}_{T}K.tsv'), 'w') as fh:
                fh.write(f'# calibrated model at {T} K, form {f}: w_i = gamma_TEA w_TEA + gamma_h w_hyd in lambda units '
                         f'(gamma_TEA {gtc}, gamma_h {ghc:.6g} on DEKR); phosphoresidues take their Ser or Thr '
                         f'({", ".join(f"{k}={v}" for k, v in sorted(letters[f].items())) or "none"})\n'
                         f'# c_sat {CSAT:g}\nresidue\tw_lambda\n')
                for k in sorted(wc):
                    fh.write(f'{k}\t{wc[k]:.10g}\n')
            for k in sorted(wc):
                rows.append({'T': T, 'form': f, 'residue': k, 'CALVADOS2_w_lambda': wc[k], 'c_sat': CSAT})
        wm, gtm, ghm = weights('recharged', T)
        with open(os.path.join(mdir, f'w_cal_{T}K.tsv'), 'w') as fh:
            fh.write(f'# calibrated model at {T} K: w_i = gamma_TEA w_TEA + gamma_h w_hyd in kJ/mol (gamma_TEA {gtm}, '
                     f'gamma_h {ghm:.6g} on DEKR); phosphate-like Glu beads take E\n# c_sat {CSAT:g}\nresidue\tw_kJmol\n')
            for k in AA:
                fh.write(f'{k}\t{wm[k]:.10g}\n')
        for f in cfa:
            crun.append({'name': f'PAGE4cal-{f}-I150-pH75', 'form': f, 'set': 'cal', 'temp': T, 'ionic': 0.15, 'ph': 7.5,
                         'residues': f'residues_{f}_plain.csv', 'cap': f'w_cal_{f}_{T}K.tsv'})
        for f in mfa:
            mrun.append({'name': f'PAGE4cal-{f}-I150', 'form': f, 'set': 'cal', 'temp': T, 'ionic': 0.15,
                         'phospho': spec[f], 'delta': f'w_cal_{T}K.tsv'})
        for k in AA:
            rows.append({'T': T, 'form': 'all (Mpipi-Recharged)', 'residue': k, 'MpipiRecharged_w_kJmol': wm[k], 'c_sat': CSAT})
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
         'submit_page4cal.sh — PAGE4 phosphoforms at 278-320 K under the calibrated model (saturating forcing) on GAIVI;\n'
         '# generated by BENDER-MT/code/make_page4_calibrated_inputs.py from GAIVI_PAGE4/submit_page4.sh'),
        ('submit_page4.sh', 'submit_page4cal.sh'),
        ('SETS="heat,smFRET,SAXS"', 'SETS="cal"'),
        ('/data/taseefr/page4', '/data/taseefr/page4cal'), ('/general/taseefr/page4', '/general/taseefr/page4cal'),
        ('PAGE4_cal"', 'P4CAL_cal"'), ('PAGE4_mp"', 'P4CAL_mp"'), ('JN=PAGE4_cal', 'JN=P4CAL_cal'), ('JN=PAGE4_mp', 'JN=P4CAL_mp'),
        ('PAGE4_resub_', 'P4CAL_resub_'),
        ('# the 36 CALVADOS-2 runs in calvados/runs.tsv', f'# the {len(C)} CALVADOS-2 runs in calvados/runs.tsv'),
        ('# the 46 Mpipi-Recharged runs in mpipi/runs.tsv', f'# the {len(M)} Mpipi-Recharged runs in mpipi/runs.tsv'),
        ('#   options: --sets heat            only the 278-320 K series (or smFRET, SAXS; comma-separated)',
         '#   options: --sets cal             (one set: the calibrated model at 278, 288, 310 and 320 K)'),
        ('TEMP=\\$(echo "\\$LINE" | cut -f4); IONIC=\\$(echo "\\$LINE" | cut -f5); PHOS=\\$(echo "\\$LINE" | cut -f6)',
         'TEMP=\\$(echo "\\$LINE" | cut -f4); IONIC=\\$(echo "\\$LINE" | cut -f5); PHOS=\\$(echo "\\$LINE" | cut -f6); '
         'DELTA=\\$(echo "\\$LINE" | cut -f7)'),
        ('python "$SCRIPT_DIR/mpipi/run_openmpipi_page4.py" --name',
         'python "$SCRIPT_DIR/mpipi/run_openmpipi_cal.py" --delta "$SCRIPT_DIR/mpipi/\\$DELTA" --name'),
        ('PH=\\$(echo "\\$LINE" | cut -f6); RES=\\$(echo "\\$LINE" | cut -f7)',
         'PH=\\$(echo "\\$LINE" | cut -f6); RES=\\$(echo "\\$LINE" | cut -f7); CAP=\\$(echo "\\$LINE" | cut -f8)'),
        ('python "\\$RUNDIR/run.py" --path "\\$RUNDIR"',
         'grep -q "fresidues: .*\\$RES\\$" "\\$RUNDIR/components.yaml" || { echo "components.yaml does not use \\$RES"; exit 1; }\n'
         'python "$SCRIPT_DIR/calvados/run_calvados_cal.py" --cap "$SCRIPT_DIR/calvados/\\$CAP" --path "\\$RUNDIR"'),
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
    for bad in ('36 CALVADOS-2', '46 Mpipi', 'smFRET, SAXS', 'PAGE4_cal"', 'PAGE4_mp"', '/page4"'):
        assert bad not in sub, bad
    open(os.path.join(a.out, 'submit_page4cal.sh'), 'w').write(sub)
    os.chmod(os.path.join(a.out, 'submit_page4cal.sh'), 0o755)
    open(os.path.join(a.out, 'README.md'), 'w').write(
        '# GAIVI_PAGE4CAL: PAGE4 phosphoforms under the calibrated model\n\n'
        'Generated by `BENDER-MT/code/make_page4_calibrated_inputs.py`. Heating series at 278, 288, 310 and 320 K, 150 mM, '
        'with the saturating hydration forcing of the calibrated model built into the force field instead of reweighting. '
        'The forcing vanishes at 300 K, so the existing GAIVI_PAGE4 300 K runs complete each series.\n\n'
        f'- {len(C)} CALVADOS-2 runs (`calvados/runs.tsv`): plain residues files plus the saturating force added by '
        '`run_calvados_cal.py --cap w_cal_<form>_<T>K.tsv`\n'
        f'- {len(M)} Mpipi-Recharged runs (`mpipi/runs.tsv`): `run_openmpipi_cal.py --delta w_cal_<T>K.tsv` around the '
        'unchanged `run_openmpipi_page4.py`; phosphate-like forms keep their `--phospho` specification\n'
        '- weights per residue and temperature in `forcing.tsv` (c_sat 2)\n\n'
        'On GAIVI (activate the calvados env first). Mpipi-Recharged is the priority, because its reweighted heating curves '
        'were the unreliable ones:\n\n'
        '```bash\nscp -r GAIVI_PAGE4CAL taseefr@<gaivi>:~/\ncd ~/GAIVI_PAGE4CAL && bash submit_page4cal.sh mpipi    # or both\n```\n\n'
        'The run scripts are copied from `~/GAIVI_PAGE4` on first use.\n'
        'Outputs: `/data/taseefr/page4cal` (CALVADOS-2) and `/general/taseefr/page4cal` (Mpipi-Recharged).\n')
    print(f'CALVADOS-2: {len(C)} runs; Mpipi-Recharged: {len(M)} runs; written to {a.out}')
    print(C.groupby('temp').size().to_dict(), M.groupby('temp').size().to_dict())


if __name__ == '__main__':
    main()
