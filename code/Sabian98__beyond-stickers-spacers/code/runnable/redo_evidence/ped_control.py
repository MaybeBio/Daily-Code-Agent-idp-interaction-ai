"""
ped_control.py — is the "aromatics central, acidic residues not" pattern of the experiment-restrained PED ensembles
(ped_residue_classes.py) coming from the experimental data, or from the software that builds ensembles?

For the same proteins, PED also holds ensembles made WITHOUT experimental restraints:
  initial pool  IDPConformerGenerator with fractional secondary structure (local sampling, no long-range attraction):
                α-syn PED00425, Sic1 PED00423, drkN PED00427
  ML            idpGAN / cg2all (trained on simulations): β-syn PED00485, PaaA2 PED00482, Sic1 PED00487, drkN PED00479,
                PDEγ PED00471; AlphaFlex (AF2 + IDP conformer generator): α-syn PED02470, FUS PED01281
                (FUS here is the full-length protein; only residues aligned to the experimental FUS 1-267 construct are
                kept, before the graphs are built)
Same analysis as ped_residue_classes.py: Calpha graphs (fixed 0.8 nm; bead-size rule with CALVADOS-2 sigma), exact
betweenness, mean over models, position-adjusted z; class means F/W/Y, His, D/E, K/R, G/S/A, other; and the per-protein
contrast (F/W/Y) - (D/E).
Standard library only.  python3 ped_control.py  ->  ped_control_summary.txt
"""
import difflib
import json
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_anchor as ea  # noqa: E402
import hub_frames_tests as ht  # noqa: E402

PED = ea.PED
FILES = {
    'PED00024': 'PED00024/PED00024e001.pdb', 'PED00023': 'PED00023/PED00023e001.pdb',
    'PED00022': 'PED00022/PED00022e001.pdb', 'PED00156': 'PED00156/PED00156e002.pdb',
    'PED00003': 'PED00003/PED00003e001.pdb', 'PED00013': 'PED00013/PED00013e001.pdb',
    'PED00511': 'PED00511/pdbfile.pdb', 'PED00494': 'PED00494/pdbfile.pdb',
    'PED00425': 'PED00425/pdbfile.pdb', 'PED00423': 'PED00423/pdbfile.pdb', 'PED00427': 'PED00427/pdbfile.pdb',
    'PED00485': 'PED00485/pdbfile.pdb', 'PED00482': 'PED00482/pdbfile.pdb', 'PED00487': 'PED00487/pdbfile.pdb',
    'PED00479': 'PED00479/pdbfile.pdb', 'PED00471': 'PED00471/pdbfile.pdb', 'PED02470': 'PED02470/pdbfile.pdb',
    'PED01281': 'PED01281/pdbfile.pdb',
}
PROTEINS = [  # protein, experimental entries, [(control entry, kind)]
    ('α-synuclein', ['PED00024'], [('PED00425', 'initial pool'), ('PED02470', 'AlphaFlex')]),
    ('Sic1', ['PED00023'], [('PED00423', 'initial pool'), ('PED00487', 'idpGAN')]),
    ('drkN', ['PED00022', 'PED00156'], [('PED00427', 'initial pool'), ('PED00479', 'idpGAN')]),
    ('β-synuclein', ['PED00003'], [('PED00485', 'idpGAN')]),
    ('PaaA2', ['PED00013'], [('PED00482', 'idpGAN')]),
    ('PDEγ', ['PED00511'], [('PED00471', 'idpGAN/cg2all')]),
    ('FUS NTD', ['PED00494'], [('PED01281', 'AlphaFlex (aligned 1-267)')]),
]
CLASS = {**dict.fromkeys('FWY', 'FWY'), 'H': 'His', **dict.fromkeys('DE', 'DE'), **dict.fromkeys('KR', 'KR'),
         **dict.fromkeys('GSA', 'GSA')}
CLASSES = ['FWY', 'His', 'DE', 'KR', 'GSA', 'other']
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def load(entry, restrict_to=None):
    seq, models = ea.ped_models(os.path.join(PED, FILES[entry]))
    n = len(seq)
    models = [m for m in models if len(m[0]) == n]
    if restrict_to:
        sm = difflib.SequenceMatcher(None, restrict_to, seq, autojunk=False)
        keep = [blk.b + k for blk in sm.get_matching_blocks() for k in range(blk.size)]
        seq = ''.join(seq[i] for i in keep)
        models = [([m[0][i] for i in keep], [m[1][i] for i in keep], [m[2][i] for i in keep]) for m in models]
    return seq, models


def one(job):
    entry, rule, c, restrict_to = job
    seq, models = load(entry, restrict_to)
    z, used = ea.ensemble_centrality(seq, models, rule, c)
    by = defaultdict(list)
    for aa, x in zip(seq, z):
        by[CLASS.get(aa, 'other')].append(x)
    return entry, rule, {k: (sum(v) / len(v), len(v)) for k, v in by.items()}, used, len(seq)


def main():
    c = json.load(open(ht.C_FILE))['calvados']
    fus_seq = ea.ped_models(os.path.join(PED, FILES['PED00494']))[0]
    jobs = []
    for prot, exps, ctrls in PROTEINS:
        for e in exps:
            jobs += [(e, r, c, None) for r in ('fixed', 'size')]
        for e, kind in ctrls:
            jobs += [(e, r, c, fus_seq if e == 'PED01281' else None) for r in ('fixed', 'size')]
    with Pool(9) as pool:
        res = {(e, r): (m, used, n) for e, r, m, used, n in pool.imap_unordered(one, jobs)}

    def contrast(m):
        return m['FWY'][0] - m['DE'][0] if 'FWY' in m and 'DE' in m else None

    p('Control: experiment-restrained vs restraint-free / ML ensembles of the same disordered proteins '
      '(mean position-adjusted z by class; contrast = FWY - DE)')
    summary = {r: defaultdict(list) for r in ('fixed', 'size')}
    for rule in ('fixed', 'size'):
        p(f'\nRule: {rule}')
        p(f'  {"protein":12s} {"ensemble":10s} {"kind":26s} {"N":>4s} {"models":>6s} ' + ''.join(f'{k:>8s}' for k in CLASSES) + '  FWY-DE')
        for prot, exps, ctrls in PROTEINS:
            for e, kind in [(x, 'experimental') for x in exps] + ctrls:
                m, used, n = res[(e, rule)]
                cv = contrast(m)
                p(f'  {prot:12s} {e:10s} {kind:26s} {n:4d} {used:6d} ' + ''.join(
                    f'{m[k][0]:+8.2f}' if k in m else f'{"-":>8s}' for k in CLASSES) + (f'  {cv:+.2f}' if cv is not None else '     -'))
                if cv is not None:
                    summary[rule]['experimental' if kind == 'experimental' else
                                  ('initial pool' if kind == 'initial pool' else 'ML')].append((prot, cv))
        for k, v in summary[rule].items():
            pos = sum(x > 0 for _, x in v)
            p(f'  {k:13s}: FWY > DE in {pos}/{len(v)} ensembles; mean contrast {sum(x for _, x in v) / len(v):+.2f}')
    open(os.path.join(HERE, 'ped_control_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
