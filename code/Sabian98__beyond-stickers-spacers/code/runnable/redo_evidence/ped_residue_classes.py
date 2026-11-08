"""
ped_residue_classes.py — in experiment-restrained ensembles of DISORDERED proteins (PED), which residue classes are
central? Tests the simulation result that aromatics are hubs and acidic residues never are (audit §3d), without
any simulation.

Ensembles: PED entries on the NAS (ZOIDBERG .../raw/pre/ped) whose ensemble is disordered by the PED statistics
(relative ASA >= 0.55 and DSSP entropy >= 0.15; folded NMR references score 0.42-0.52 and <= 0.12). The full-length
hnRNPA1 ensembles (two folded RRMs; 0.47 / 0.25) fail this and are excluded. PED00543 (all-atom MD) is a model,
reported separately.
Per model: Calpha contact graph (fixed 0.8 nm, |i-j| >= 3; or bead-size rule d < 1.37 * sigma_ij with CALVADOS-2
sigma), exact betweenness; mean over models; position-adjusted within-protein z. Class means: aromatic F/W/Y, His (separate: it does not behave like F/W/Y here),
acidic D/E, basic K/R, G/S/A, other.
Standard library only.  python3 ped_residue_classes.py  ->  ped_residue_classes_summary.txt, ped_residue_classes.csv
"""
import csv
import json
import math
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_anchor as ea  # noqa: E402
import hub_frames_tests as ht  # noqa: E402

PED = ea.PED
CAT = (__import__('os').environ['BSS_DATA'] + '/ZOIDBERG-ICML 2027/experimental_data/raw/ped/ped_all_entries.json')
ENSEMBLES = {  # entry: (ensemble id in the catalogue, file, kind)
    'PED00003': ('e001', 'PED00003/PED00003e001.pdb', 'PRE'),
    'PED00022': ('e001', 'PED00022/PED00022e001.pdb', 'NOE/PRE'),
    'PED00156': ('e001', 'PED00156/PED00156e002.pdb', 'NOE/PRE/smFRET'),
    'PED00511': ('e001', 'PED00511/pdbfile.pdb', 'NOE/PRE'),
    'PED00013': ('e001', 'PED00013/PED00013e001.pdb', 'NOE'),
    'PED00004': ('e001', 'PED00004/PED00004e001.pdb', 'NOE'),
    'PED00023': ('e001', 'PED00023/PED00023e001.pdb', 'PRE'),
    'PED00498': ('e002', 'PED00498/pdbfile.pdb', 'NOE/PRE'),
    'PED00499': ('e002', 'PED00499/pdbfile.pdb', 'NOE/PRE'),
    'PED00500': ('e002', 'PED00500/pdbfile.pdb', 'NOE/PRE'),
    'PED00494': ('e001', 'PED00494/pdbfile.pdb', 'DEER'),
    'PED00024': ('e001', 'PED00024/PED00024e001.pdb', 'PRE'),
    'PED00225': ('e001', 'PED00225/PED00225e000.pdb', 'NOE/PRE/RDC/SAXS'),
    'PED00212': ('e001', 'PED00212/PED00212e000.pdb', 'DEER/SAXS/PRE'),
    'PED00215': ('e001', 'PED00215/PED00215e000.pdb', 'PRE'),
    'PED00543': ('e001', 'PED00543/pdbfile.pdb', 'MD model'),
}
CLASS = {**dict.fromkeys('FWY', 'aromatic FWY'), 'H': 'His', **dict.fromkeys('DE', 'acidic DE'), **dict.fromkeys('KR', 'basic KR'),
         **dict.fromkeys('GSA', 'GSA')}
CLASSES = ['aromatic FWY', 'His', 'acidic DE', 'basic KR', 'GSA', 'other']
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def one(job):
    entry, rule, c = job
    ens, path, kind = ENSEMBLES[entry]
    seq, models = ea.ped_models(os.path.join(PED, path))
    n = len(seq)
    models = [m for m in models if len(m[0]) == n]
    if not models or n < 10:
        return entry, rule, None
    z, used = ea.ensemble_centrality(seq, models, rule, c)
    return entry, rule, (seq, z, used)


def main():
    cat = {e['entry_id']: e for e in json.load(open(CAT))}
    c_cal = json.load(open(ht.C_FILE))['calvados']
    meta = {}
    for e, (ens, path, kind) in ENSEMBLES.items():
        x = cat[e]
        m = [en for en in x['ensembles'] if en['ensemble_id'] == ens] or x['ensembles'][:1]
        ch = m[0]['chains'][0]
        rasa, ent = ch.get('relative_asa_mean') or 0, ch.get('entropy_dssp_mean') or 0
        meta[e] = (x['description']['title'], kind, rasa, ent, rasa >= 0.55 and ent >= 0.15)
    jobs = [(e, rule, c_cal) for e in ENSEMBLES if meta[e][4] for rule in ('fixed', 'size')]
    with Pool(9) as pool:
        res = {(e, r): v for e, r, v in pool.imap_unordered(one, jobs)}
    rows = []
    p('Residue-class centrality in experiment-restrained ensembles of disordered proteins (PED; mean position-adjusted z of '
      'mean exact betweenness over models)')
    p('Excluded as not disordered (rel ASA < 0.55 or DSSP entropy < 0.15): '
      + ', '.join(f'{e} ({meta[e][2]:.2f}/{meta[e][3]:.2f})' for e in ENSEMBLES if not meta[e][4]))
    pooled = {rule: defaultdict(list) for rule in ('fixed', 'size')}
    tally = {rule: [0, 0, 0] for rule in ('fixed', 'size')}     # aromatic > acidic, aromatic top class, n
    for rule in ('fixed', 'size'):
        p(f'\nRule: {rule}')
        p(f'  {"entry":9s} {"N":>4s} {"models":>6s}  ' + ''.join(f'{c:>16s}' for c in CLASSES) + '   restraints / protein')
        for e in ENSEMBLES:
            if not meta[e][4] or (e, rule) not in res or res[(e, rule)] is None:
                continue
            seq, z, used = res[(e, rule)]
            by = defaultdict(list)
            for aa, x in zip(seq, z):
                by[CLASS.get(aa, 'other')].append(x)
                if e != 'PED00543':
                    pooled[rule][aa].append(x)
            mean = {c: (sum(v) / len(v), len(v)) for c, v in by.items()}
            line = f'  {e:9s} {len(seq):4d} {used:6d}  '
            for c in CLASSES:
                line += f'{mean[c][0]:+9.2f} (n {mean[c][1]:3d})' if c in mean else f'{"-":>16s}'
            p(line + f'   {meta[e][1]} / {meta[e][0][:45]}')
            rows.append([e, rule, len(seq), used] + [f'{mean[c][0]:.4f}' if c in mean else '' for c in CLASSES] +
                        [mean[c][1] if c in mean else 0 for c in CLASSES])
            if e != 'PED00543' and 'aromatic FWY' in mean and 'acidic DE' in mean:
                t = tally[rule]
                t[2] += 1
                t[0] += mean['aromatic FWY'][0] > mean['acidic DE'][0]
                t[1] += max((c for c in mean if mean[c][1] >= 3), key=lambda c: mean[c][0]) == 'aromatic FWY'
        t = tally[rule]
        sign_p = sum(math.comb(t[2], k) for k in range(t[0], t[2] + 1)) / 2 ** t[2] if t[2] else float('nan')
        p(f'  experimental ensembles with both classes: {t[2]}; aromatic > acidic in {t[0]} (sign test p = {sign_p:.4f}); '
          f'aromatic is the most central class in {t[1]}')
        p('  pooled over experimental ensembles, by residue type (mean z, n): ' + '  '.join(
            f'{aa} {sum(v) / len(v):+.2f} ({len(v)})' for aa, v in sorted(pooled[rule].items(), key=lambda kv: -sum(kv[1]) / len(kv[1]))
            if len(v) >= 10))
    with open(os.path.join(HERE, 'ped_residue_classes.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['entry', 'rule', 'N', 'models'] + [f'z_{c}' for c in CLASSES] + [f'n_{c}' for c in CLASSES])
        w.writerows(rows)
    open(os.path.join(HERE, 'ped_residue_classes_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
