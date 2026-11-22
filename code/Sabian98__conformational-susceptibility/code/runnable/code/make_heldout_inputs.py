"""
make_heldout_inputs.py — GAIVI inputs for the held-out test of the calibrated heat response: independent experimental
temperature series of disordered proteins that were not used in the Wuttke calibration (tea_calibrate.py).

Data sets (downloaded 2026-10-11 into BENDER-MT/heldout_temperature/raw/, see datasets.tsv for the reference values):
  PNt     Baxa et al. Protein Sci 2024 33:e4986 (PMC11010952). N-terminal 334 residues of B. pertussis pertactin
          (sequence from their Table S1). SEC-SAXS Rg (molecular form factor) in 10 mM sodium phosphate pH 7.4 with
          0 / 150 / 500 mM KCl (ionic strength 0.025 / 0.175 / 0.525 M; 0.025 as stated in the paper); raw profiles in
          github.com/sosnicklab/SAXSonIDPs. Temperatures as measured: 0 KCl 10, 30, 45, 58 C; 150 mM 10, 20, 30, 40, 50,
          58 C; 500 mM 7, 20, 30, 45, 60 C.
  sfAFP   same paper: reduced snow flea antifreeze protein (Gly-rich control, flat with temperature), sequence 'sfAFP'
          of Table S1 (with the C-terminal SNAC and His6 tags as listed there; whether the SNAC tag was cleaved for the
          SAXS sample is not stated), 150 mM KCl, pH 7.5, 10, 30, 45, 50, 58 C.
  asyn    Woodard et al. Biophys J 2018 115:1190 (PMC6170596). Human alpha-synuclein (P37840) F4W/V26C and Y39W/V66C,
          Trp donor and coumarin (I-Cou) acceptor on the Cys, R0 ~2.4 nm; trFRET distance distributions (Table 1) at 5,
          10, 20, 30, 40 C in 20 mM Tris, 2 mM EDTA, pH 7.5 (ionic strength estimated 0.03 M).
  Aap     Yarawsky et al. J Mol Biol 2017 429:261 (PMC5363081). S. epidermidis RP62A Aap proline/glycine-rich region,
          residues 2225-2359 (the paper cites the deleted UniProt entry Q9L470; sequence from the current full-length
          entry Q5HKE8, which has the flanking EYGPT and LPDTG motifs the paper shows). DLS Rh at 5, 15, 25, 35, 45 C in
          20 mM KPO4 pH 7.4, 150 mM NaCl (ionic strength ~0.195 M). Per-temperature DLS values are only plotted (Fig. 4c);
          the text gives 38.9 to 37.8 A and a linear-fit slope of -0.015 A per C, and the supplement (Table S2, supplied
          by the user 2026-10-11) gives AUC frictional ratios of 2.14, 2.14 and 2.10 at 4, 20 and 37 C.
  p53     Human p53 N-terminal region 1-93 (P04637), DLS Rh in 10 mM sodium phosphate, 100 mM NaCl, pH 7 (ionic strength
          ~0.118 M; Perez et al. Proteins 2014 82:3373, PMC4237723, same lab and buffer as Langridge et al. Proteins 2014
          82:668, whose temperature series Yarawsky 2017 re-plots and fits with a slope of -0.076 A per C, against -0.015
          for Aap in the same panel). Rh 32.8 A at 25 C (Perez 2014). Simulated at 5, 15, 25, 35, 45 C, as Aap.
CALVADOS-2 runs use the standard residues table (pH sets the His charge); Mpipi-Recharged runs use the plain sequences
(no pH dependence). Run scripts, submit logic and run length are those of GAIVI_PAGE4 (copied unchanged), time factor 3.
    python make_heldout_inputs.py [--out $CS_DATA/Scripts/GAIVI_HELDOUT]
Writes: calvados/{forms.fasta, residues_CALVADOS2.csv, runs.tsv, prepare_page4.py}, mpipi/{forms.fasta, runs.tsv,
        run_openmpipi_page4.py, residues_CALVADOS2.csv}, submit_heldout.sh, datasets.tsv (the package README.md is kept)
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402

import pandas as pd  # noqa: E402

BM = (__import__('os').environ['CS_DATA'] + '/BENDER-MT')
RAW = f'{BM}/heldout_temperature/raw'
PAGE4 = os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_PAGE4'))
RESIDUES = os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv'))
SI_TXT = '/private/tmp/claude-501/-Volumes-User-Homes-Nu-Project-BENDER-MT/8a63b16b-e317-4515-bb69-cc052f438731/scratchpad/baxa_si.txt'

PNT = ('DWNNQSIVKTGERQHGIHIQGSDPGGVRTASGTTIKVSGRQAQGILLENPAAELQFRNGSVTSSGQLSDDGIRRFLGTVTVKAGKLVADHATLANVGDTWDDDGIALY'
       'VAGEQAQASIADSTLQGAGGVQIERGANVTVQRSAIVDGGLHIGALQSLQPEDLPPSRVVLRDTNVTAVPASGAPAAVSVLGASELTLDGGHITGGRAAGVAAMQGAVVH'
       'LQRATIRRGEALAGGAVPGGAVPGGAVPGGFGPGGFGPVLDGWYGVDVSGSSVELAQSIVEAPELGAAIRVGRGARVTVPGGSLSAPHGNVIETGGARRFAPQAAPLSIT'
       'LQAGAH')
SFAFP = 'CKGADGAHGVNGCPGTAGAAGSVGGPGCDGGHGGNGGNGNPGCAGGVGGAGGASGGTGVGGWGGKGGSGTPKGADGAPGAPGSHHWHHHHHH'
C2K = 273.15


def fasta(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>'))


def mutate(seq, muts):
    s = list(seq)
    for m in muts:
        wt, pos, new = m[0], int(m[1:-1]), m[-1]
        assert s[pos - 1] == wt, (m, s[pos - 1])
        s[pos - 1] = new
    return ''.join(s)


def forms():
    asyn = fasta(f'{RAW}/uniprot/P37840.fasta')
    aap = fasta(f'{RAW}/uniprot/Q5HKE8.fasta')
    assert len(asyn) == 140 and aap[2218:2223] == 'EYGPT' and aap[2359:2364] == 'LPDTG'
    p53 = fasta(f'{RAW}/uniprot/P04637.fasta')
    f = {'PNt': PNT, 'sfAFP': SFAFP, 'asyn4W26C': mutate(asyn, ['F4W', 'V26C']),
         'asyn39W66C': mutate(asyn, ['Y39W', 'V66C']), 'AapPGR': aap[2224:2359], 'p53N93': p53[:93]}
    if os.path.exists(SI_TXT):                       # cross-check against Table S1 as extracted from the Baxa SI
        si = open(SI_TXT).read()
        assert PNT in si and SFAFP in si, 'PNt or sfAFP differs from Table S1'
    assert len(f['PNt']) == 334 and len(f['AapPGR']) == 135 and f['p53N93'].startswith('MEEPQSDPSV') and len(f['p53N93']) == 93
    return f


# set -> (form, ionic M, pH, temperatures in C, readout)
SETS = {
    'PNt-I025': ('PNt', 0.025, 7.4, [10, 30, 45, 58], 'SAXS Rg'),
    'PNt-I175': ('PNt', 0.175, 7.4, [10, 20, 30, 40, 50, 58], 'SAXS Rg'),
    'PNt-I525': ('PNt', 0.525, 7.4, [7, 20, 30, 45, 60], 'SAXS Rg'),
    'sfAFP-I175': ('sfAFP', 0.175, 7.5, [10, 30, 45, 50, 58], 'SAXS Rg'),
    'asyn4W26C-I030': ('asyn4W26C', 0.030, 7.5, [5, 10, 20, 30, 40], 'trFRET P(r) 4-26'),
    'asyn39W66C-I030': ('asyn39W66C', 0.030, 7.5, [5, 10, 20, 30, 40], 'trFRET P(r) 39-66'),
    'AapPGR-I195': ('AapPGR', 0.195, 7.4, [5, 15, 25, 35, 45], 'DLS Rh'),
    'p53N93-I118': ('p53N93', 0.118, 7.0, [5, 15, 25, 35, 45], 'DLS Rh'),
}

# reference values from the papers (A); NaN where only plotted or in the inaccessible SI
REF = [
    # set, T (C), observable, value, low, high, source
    ('PNt-I175', 10, 'Rg', 52.72, 52.58, 52.86, 'Baxa 2024 text'), ('PNt-I175', 58, 'Rg', 46.98, 46.84, 47.12, 'Baxa 2024 text'),
    ('PNt-I025', 10, 'Rg', 47.07, 46.95, 47.19, 'Baxa 2024 text'), ('PNt-I025', 58, 'Rg', 41.83, 41.62, 42.04, 'Baxa 2024 text'),
    ('sfAFP-I175', 10, 'Rg', 23.69, 23.57, 23.81, 'Baxa 2024 text'), ('sfAFP-I175', 58, 'Rg', 24.29, 24.21, 24.37, 'Baxa 2024 text'),
    ('AapPGR-I195', 5, 'Rh', 38.9, None, None, 'Yarawsky 2017 text (range 38.9 to 37.8 over 5-45 C)'),
    ('AapPGR-I195', 45, 'Rh', 37.8, None, None, 'Yarawsky 2017 text'),
    ('AapPGR-I195', None, 'dRh/dT (A per C)', -0.015, None, None,
     'Yarawsky 2017 Fig. 4c linear fit of DLS 5-45 C (p53(1-93) in the same panel: -0.076)'),
    ('AapPGR-I195', 4, 'AUC f/f0', 2.14, None, None, 'Yarawsky 2017 Table S2 (s20,w 1.06 S; 20 mM KPO4, 150 mM NaCl)'),
    ('AapPGR-I195', 20, 'AUC f/f0', 2.14, None, None, 'Yarawsky 2017 Table S2 (s20,w 1.05 S)'),
    ('AapPGR-I195', 37, 'AUC f/f0', 2.10, None, None, 'Yarawsky 2017 Table S2 (s20,w 1.03 S)'),
    ('p53N93-I118', 25, 'Rh', 32.8, None, None, 'Perez 2014 text (DLS; SEC gives nearly the same, Table I)'),
    ('p53N93-I118', None, 'dRh/dT (A per C)', -0.076, None, None,
     'Yarawsky 2017 Fig. 4c linear fit of the DLS data of Langridge 2014 (Aap in the same panel: -0.015)'),
]
WOODARD = """4-26 5 33.8 32.6 35.6 25.1 23.9 26.9 25.4 3.5
4-26 10 33.8 33.0 34.5 25.3 24.4 26.2 25.1 3.4
4-26 20 32.3 30.4 34.8 26.4 24.2 29.1 22.0 3.0
4-26 30 30.9 27.6 34.7 26.7 23.9 30.6 19.4 2.8
4-26 40 29.3 23.5 35.4 28.9 24.5 34.4 13.1 2.1
39-66 5 38.9 36.1 42.8 31.6 28.6 35.8 26.7 2.1
39-66 10 37.6 35.2 40.8 30.3 27.7 33.9 26.1 2.3
39-66 20 34.8 32.8 37.1 27.6 25.5 30.1 24.5 2.8
39-66 30 33.2 27.6 40.4 28.4 23.8 36.0 21.2 2.5
39-66 40 32.9 22.6 43.7 29.4 22.3 41.0 18.5 2.2"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_HELDOUT')))
    a = ap.parse_args()
    F = forms()
    cdir, mdir = os.path.join(a.out, 'calvados'), os.path.join(a.out, 'mpipi')
    os.makedirs(cdir, exist_ok=True)
    os.makedirs(mdir, exist_ok=True)
    std = set(pd.read_csv(RESIDUES).one)
    for k, s in F.items():
        assert set(s) <= std, (k, set(s) - std)
    for d in (cdir, mdir):
        with open(os.path.join(d, 'forms.fasta'), 'w') as fh:
            for k, s in F.items():
                fh.write(f'>{k}\n{s}\n')
    shutil.copy2(RESIDUES, os.path.join(cdir, 'residues_CALVADOS2.csv'))
    shutil.copy2(RESIDUES, os.path.join(mdir, 'residues_CALVADOS2.csv'))   # run_openmpipi_page4 writes top.pdb from it
    shutil.copy2(os.path.join(PAGE4, 'calvados', 'prepare_page4.py'), cdir)
    shutil.copy2(os.path.join(PAGE4, 'mpipi', 'run_openmpipi_page4.py'), mdir)
    crun, mrun = [], []
    for st, (form, I, ph, temps, readout) in SETS.items():
        for tc in temps:
            T = int(round(tc + C2K))
            crun.append({'name': f'HO-{st}-pH{round(ph * 10)}', 'form': form, 'set': st, 'temp': T, 'ionic': I, 'ph': ph,
                         'residues': 'residues_CALVADOS2.csv'})
            mrun.append({'name': f'HO-{st}', 'form': form, 'set': st, 'temp': T, 'ionic': I, 'phospho': '-'})
    pd.DataFrame(crun).to_csv(os.path.join(cdir, 'runs.tsv'), sep='\t', index=False)
    pd.DataFrame(mrun).to_csv(os.path.join(mdir, 'runs.tsv'), sep='\t', index=False)
    rows = [{'set': s, 'T_C': t, 'observable': o, 'value_A': v, 'low_A': lo, 'high_A': hi, 'source': src}
            for s, t, o, v, lo, hi, src in REF]
    for line in WOODARD.splitlines():
        lp, t, pk, plo, phi, fw, flo, fhi, sa, sb = line.split()
        rows.append({'set': 'asyn4W26C-I030' if lp == '4-26' else 'asyn39W66C-I030',
                     'T_C': int(t), 'observable': f'P(r) peak {lp}', 'value_A': float(pk), 'low_A': float(plo),
                     'high_A': float(phi), 'source': f'Woodard 2018 Table 1 (FWHM {fw} A, a {sa}, b {sb}e-3; 66% CI)'})
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'datasets.tsv'), sep='\t', index=False)
    # submit script: GAIVI_PAGE4's, with its own job names, output roots and set handling
    sub = open(os.path.join(PAGE4, 'submit_page4.sh')).read()
    sub = sub.replace('submit_page4.sh — PAGE4 phosphoform runs on GAIVI (CALVADOS-2 and Mpipi-Recharged)',
                      'submit_heldout.sh — held-out temperature-series runs on GAIVI (CALVADOS-2 and Mpipi-Recharged);\n'
                      '# generated by BENDER-MT/code/make_heldout_inputs.py from GAIVI_PAGE4/submit_page4.sh')
    sub = sub.replace('submit_page4.sh', 'submit_heldout.sh')
    sub = sub.replace('SETS="heat,smFRET,SAXS"', 'SETS="all"')
    sub = sub.replace('[[ ",$SETS," == *",$SET,"* ]] || continue', '[[ "$SETS" == all || ",$SETS," == *",$SET,"* ]] || continue')
    sub = sub.replace('/data/taseefr/page4', '/data/taseefr/heldout').replace('/general/taseefr/page4', '/general/taseefr/heldout')
    sub = sub.replace('PAGE4_cal', 'HO_cal').replace('PAGE4_mp', 'HO_mp').replace('PAGE4_resub_', 'HO_resub_')
    sub = re.sub(r'#   bash submit_heldout\.sh calvados\s+# the 36 CALVADOS-2 runs in calvados/runs\.tsv',
                 f'#   bash submit_heldout.sh calvados                 # the {len(crun)} CALVADOS-2 runs in calvados/runs.tsv', sub)
    sub = re.sub(r'#   bash submit_heldout\.sh mpipi\s+# the 46 Mpipi-Recharged runs in mpipi/runs\.tsv',
                 f'#   bash submit_heldout.sh mpipi                    # the {len(mrun)} Mpipi-Recharged runs in mpipi/runs.tsv', sub)
    sub = sub.replace('#   options: --sets heat            only the 278-320 K series (or smFRET, SAXS; comma-separated)',
                      '#   options: --sets PNt-I175        only some sets (comma-separated; default all)')
    sub = sub.replace('(default 3 = 300 ns for PAGE4)', '(default 3; ~300 ns for chains up to 150 residues, ~1 us for PNt)')
    # the per-run check for queued or running jobs comes from submit_page4.sh (job names renamed above)
    for bad in ('36 CALVADOS-2', '46 Mpipi', 'smFRET, SAXS'):
        assert bad not in sub, bad
    open(os.path.join(a.out, 'submit_heldout.sh'), 'w').write(sub)
    os.chmod(os.path.join(a.out, 'submit_heldout.sh'), 0o755)
    print(f'CALVADOS-2: {len(crun)} runs; Mpipi-Recharged: {len(mrun)} runs; forms {", ".join(f"{k} ({len(v)})" for k, v in F.items())}; '
          f'written to {a.out}')


if __name__ == '__main__':
    main()
