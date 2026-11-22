"""
make_pnt_calibrated_inputs.py — GAIVI inputs for direct simulations of PNt under the calibrated models, replacing the
reweighted held-out values whose effective sample fraction is 0.001-0.005 at the end temperatures (Table S1b). Written
2026-10-23 after residual_sequence_check.py; exploratory, nothing frozen is changed.

PNt (334 residues, make_heldout_inputs.PNT) at the held-out temperatures and buffers of Baxa et al. 2024:
  PNt-I025   0.025 M, 283 303 318 331 K;   PNt-I175   0.175 M, 283 293 303 313 323 331 K;   pH 7.4
Variants (one run per variant, salt and temperature; 2 x 10 runs per force field):
  A0   frozen calibrated model: TEA on every residue (gamma_TEA) + hydrophilic forcing on D E K R N Q S T H G (gamma_hyd),
       heldout_analysis.FROZEN
  S1   charged-group forcing of Amendment 3: TEA (same gamma_TEA) + hydrophilic forcing on D E K R only (gamma_chg),
       heldout_exploratory.S1_GAMMA
The forcing is the energy the reweighting adds, built into the model at each temperature:
  CALVADOS-2       lambda_i(T) = lambda_i + gamma_TEA * (-dG_E,i(T)) + gamma_h * [du_h,i(T) - du_h,i(T0)]  (h residues)
                   in a residues file per run (Ashbaugh-Hatch is affine in lambda_ij = mean of lambda_i, lambda_j, so
                   this equals H + sum_c gamma_c E_c of tea_reweight.class_energies / tea_calibrate.calvados_hyd_energies)
  Mpipi-Recharged  delta_i(T) in kJ/mol (the same sum x 4.184) on the attractive Wang-Frenkel shape of
                   tea_calibrate.recharged_class_energies, added by run_openmpipi_cal.py as one extra force
Checked by pnt_cal_check.py on the held-out PNt frames (direct energy difference = gamma . E of the reweighting).
Run scripts, run length (time factor 3, as the held-out runs) and submit logic are GAIVI_PAGE4's (copied in on GAIVI).
    python make_pnt_calibrated_inputs.py [--out $CS_DATA/Scripts/GAIVI_PNTCAL]
Writes: calvados/{forms.fasta, residues_<variant>_<T>K.csv, runs.tsv}, mpipi/{forms.fasta, delta_<variant>_<T>K.tsv,
        runs.tsv, residues_CALVADOS2.csv, run_openmpipi_cal.py}, forcing.tsv, submit_pntcal.sh, README.md
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
import heldout_analysis as HA  # noqa: E402
import heldout_exploratory as HX  # noqa: E402
import make_heldout_inputs as MH  # noqa: E402
import tea_reweight as W  # noqa: E402

SETS = {'PNt-I025': (0.025, [283, 303, 318, 331]), 'PNt-I175': (0.175, [283, 293, 303, 313, 323, 331])}
PH = 7.4
T0 = W.T0
HYD = {'A0': HA.TCal.HYD_RES, 'S1': 'DEKR'}
AA = 'ACDEFGHIKLMNPQRSTVWY'


def gam(model, variant):
    """(gamma_TEA, gamma_h) of the variant, mol/kcal (CALVADOS-2) or per unit gamma' (Mpipi-Recharged)."""
    if variant == 'A0':
        return HA.FROZEN[model]['tea'], HA.FROZEN[model]['hyd']
    return HX.S1_TEA[model], HX.S1_GAMMA[model]


def forcing(model, variant, T):
    """Per-residue change of the stickiness at T (CALVADOS-2: lambda units; Mpipi-Recharged: kcal/mol per unit eps)."""
    gt, gh = gam(model, variant)
    return {a: gt * (-W.dGE(a, T)) + (gh * (W.du_h(a, T) - W.du_h(a, T0)) if a in HYD[variant] else 0.0) for a in AA}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_PNTCAL')))
    a = ap.parse_args()
    cdir, mdir = os.path.join(a.out, 'calvados'), os.path.join(a.out, 'mpipi')
    os.makedirs(cdir, exist_ok=True)
    os.makedirs(mdir, exist_ok=True)
    seq = MH.PNT
    import mdtraj as md
    for m, d in (('CALVADOS2', 'HO-PNt-I175-pH74'), ('MPIPI_Recharged', 'HO-PNt-I175')):
        top = md.load_topology(f'{HA.BM}/heldout_{m}/303K/{d}/top.pdb')
        assert ''.join(r.code for r in top.residues) == seq, f'PNt differs from the held-out {m} run'
    assert len(seq) == 334
    base = pd.read_csv(MH.RESIDUES)
    for d in (cdir, mdir):
        open(os.path.join(d, 'forms.fasta'), 'w').write(f'>PNt\n{seq}\n')
    crun, mrun, rows = [], [], []
    for variant in ('A0', 'S1'):
        for st, (ionic, temps) in SETS.items():
            for T in temps:
                dl = forcing('calvados', variant, float(T))
                r = base.copy()
                r['lambdas'] = r['lambdas'] + r['one'].map(dl).fillna(0.0)
                fn = f'residues_{variant}_{T}K.csv'
                r.to_csv(os.path.join(cdir, fn), index=False)
                crun.append({'name': f'HC-{st}-{variant}-pH{round(PH * 10)}', 'form': 'PNt', 'set': variant, 'temp': T,
                             'ionic': ionic, 'ph': PH, 'residues': fn})
                de = {k: 4.184 * v for k, v in forcing('recharged', variant, float(T)).items()}
                fd = f'delta_{variant}_{T}K.tsv'
                with open(os.path.join(mdir, fd), 'w') as fh:
                    fh.write(f'# {variant} at {T} K: delta_i in kJ/mol on the attractive Wang-Frenkel shape '
                             f'(gamma_TEA {gam("recharged", variant)[0]}, gamma_h {gam("recharged", variant)[1]} on '
                             f'{HYD[variant]})\nresidue\tdelta_kJmol\n')
                    for k in AA:
                        fh.write(f'{k}\t{de[k]:.10g}\n')
                mrun.append({'name': f'HC-{st}-{variant}', 'form': 'PNt', 'set': variant, 'temp': T, 'ionic': ionic,
                             'phospho': '-', 'delta': fd})
                for k in AA:
                    rows.append({'variant': variant, 'T': T, 'residue': k, 'CALVADOS2_dlambda': dl[k],
                                 'CALVADOS2_lambda': float(base.set_index('one').loc[k, 'lambdas'] + dl[k]),
                                 'MpipiRecharged_delta_kJmol': de[k]})
    shutil.copy2(MH.RESIDUES, os.path.join(mdir, 'residues_CALVADOS2.csv'))
    shutil.copy2(os.path.join(HERE, 'gaivi_templates', 'run_openmpipi_cal.py'), os.path.join(mdir, 'run_openmpipi_cal.py'))
    C, M = pd.DataFrame(crun), pd.DataFrame(mrun)
    assert not C.duplicated(['name', 'temp']).any() and not M.duplicated(['name', 'temp']).any()
    C.to_csv(os.path.join(cdir, 'runs.tsv'), sep='\t', index=False)
    M.to_csv(os.path.join(mdir, 'runs.tsv'), sep='\t', index=False)
    F = pd.DataFrame(rows).drop_duplicates(['variant', 'T', 'residue'])
    F.to_csv(os.path.join(a.out, 'forcing.tsv'), sep='\t', index=False)
    # submit script: submit_page4.sh logic, renamed; Mpipi jobs go through run_openmpipi_cal.py with the run's delta file
    sub = open(os.path.join(HERE, 'gaivi_templates', 'submit_page4.sh')).read()
    rep = [
        ('submit_page4.sh — PAGE4 phosphoform runs on GAIVI (CALVADOS-2 and Mpipi-Recharged)',
         'submit_pntcal.sh — PNt under the calibrated models (A0, S1) on GAIVI (CALVADOS-2 and Mpipi-Recharged);\n'
         '# generated by BENDER-MT/code/make_pnt_calibrated_inputs.py from GAIVI_PAGE4/submit_page4.sh'),
        ('submit_page4.sh', 'submit_pntcal.sh'),
        ('SETS="heat,smFRET,SAXS"', 'SETS="A0,S1"'),
        ('/data/taseefr/page4', '/data/taseefr/pntcal'), ('/general/taseefr/page4', '/general/taseefr/pntcal'),
        ('PAGE4_cal', 'PNTCAL_cal'), ('PAGE4_mp', 'PNTCAL_mp'), ('PAGE4_resub_', 'PNTCAL_resub_'),
        ('# the 36 CALVADOS-2 runs in calvados/runs.tsv', f'# the {len(C)} CALVADOS-2 runs in calvados/runs.tsv'),
        ('# the 46 Mpipi-Recharged runs in mpipi/runs.tsv', f'# the {len(M)} Mpipi-Recharged runs in mpipi/runs.tsv'),
        ('#   options: --sets heat            only the 278-320 K series (or smFRET, SAXS; comma-separated)',
         '#   options: --sets A0              only one variant (A0 frozen model, S1 charged-group forcing; default both)'),
        ('(default 3 = 300 ns for PAGE4)', '(default 3, as the held-out PNt runs, about 1 us)'),
        ('TEMP=\\$(echo "\\$LINE" | cut -f4); IONIC=\\$(echo "\\$LINE" | cut -f5); PHOS=\\$(echo "\\$LINE" | cut -f6)',
         'TEMP=\\$(echo "\\$LINE" | cut -f4); IONIC=\\$(echo "\\$LINE" | cut -f5); PHOS=\\$(echo "\\$LINE" | cut -f6); '
         'DELTA=\\$(echo "\\$LINE" | cut -f7)'),
        ('python "$SCRIPT_DIR/mpipi/run_openmpipi_page4.py" --name',
         'python "$SCRIPT_DIR/mpipi/run_openmpipi_cal.py" --delta "$SCRIPT_DIR/mpipi/\\$DELTA" --name'),
        ('python "\\$RUNDIR/run.py" --path "\\$RUNDIR"',
         'grep -q "fresidues: .*\\$RES\\$" "\\$RUNDIR/components.yaml" || { echo "components.yaml does not use \\$RES"; exit 1; }\n'
         'python "\\$RUNDIR/run.py" --path "\\$RUNDIR"'),
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
    open(os.path.join(a.out, 'submit_pntcal.sh'), 'w').write(sub)
    os.chmod(os.path.join(a.out, 'submit_pntcal.sh'), 0o755)
    open(os.path.join(a.out, 'README.md'), 'w').write(
        '# GAIVI_PNTCAL: PNt under the calibrated models\n\n'
        'Generated by `BENDER-MT/code/make_pnt_calibrated_inputs.py`. PNt (334 residues) at the held-out temperatures and '
        'buffers (0.025 M: 283-331 K, 4 temperatures; 0.175 M: 283-331 K, 6 temperatures), with the hydration forcing '
        'built into the model instead of reweighting.\n\n'
        f'- {len(C)} CALVADOS-2 runs (`calvados/runs.tsv`): per-run residues file with lambda_i(T)\n'
        f'- {len(M)} Mpipi-Recharged runs (`mpipi/runs.tsv`): `run_openmpipi_cal.py` adds the forcing (per-run delta file) '
        'to the unchanged `run_openmpipi_page4.py`\n'
        '- variants: A0 (frozen calibrated model) and S1 (charged-group forcing, Amendment 3); values in `forcing.tsv`\n\n'
        'On GAIVI (activate the calvados env first):\n\n'
        '```bash\nscp -r GAIVI_PNTCAL taseefr@<gaivi>:~/\ncd ~/GAIVI_PNTCAL && bash submit_pntcal.sh both\n```\n\n'
        'Only A0: `bash submit_pntcal.sh both --sets A0`. The run scripts are copied from `~/GAIVI_PAGE4` on first use.\n'
        'Outputs: `/data/taseefr/pntcal` (CALVADOS-2) and `/general/taseefr/pntcal` (Mpipi-Recharged).\n')
    print(f'CALVADOS-2: {len(C)} runs; Mpipi-Recharged: {len(M)} runs; written to {a.out}')
    print(F[F.residue.isin(list('DKLG'))].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
