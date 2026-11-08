"""
ped_v2_build.py — sequences, model counts, AIUPred disorder and sequence novelty for the PED v2 extension (ped_v2_plan.md).
Reads only the first model of each deposited file (no labels). Writes ped_v2_set.csv (one row per entry, the input for
predict_frozen.py and for labelling) and ped_v2_build_summary.txt.
  python3 ped_v2_build.py AIUPRED_DIR      # AIUPRED_DIR holds <acc>.fasta, <acc>.json and entry_accs.json
"""
import csv
import difflib
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ped_null as pn  # noqa: E402

AI = sys.argv[1]
SETS = {
    'A': {'PED00229': 'p53-TAD1', 'PED00236': 'p53-TAD1', 'PED00231': 'RS-repeat peptide', 'PED00238': 'RS-repeat peptide',
          'PED00020': 'MeV NTAIL', 'PED00025': 'Sendai NTAIL', 'PED00016': 'p15PAF', 'PED00142': 'KISS-1',
          'PED00153': 'MKK7', 'PED00155': 'Histatin-5', 'PED00175': 'E1A 36-146', 'PED00181': 'N-WASP domain V',
          'PED00184': 'Tif2 NRID', 'PED00206': 'AR Tau-5 R2-R3', 'PED00210': 'Tir 388-550',
          'PED00522': 'SARS-CoV-2 N linker', 'PED07323': 'IAPP'},
    'B': {'PED00493': 'FUS NTD condensed', 'PED00454': 'pSic1', 'PED00192': 'Tau K18', 'PED00230': 'ACTR-AD1',
          'PED00237': 'ACTR-AD1'},
    'C': {'PED00531': 'Aβ40 C36m', 'PED00532': 'Aβ40 C22*', 'PED00533': 'Aβ40 a99SB-disp', 'PED00534': 'ACTR C36m',
          'PED00535': 'ACTR C22*', 'PED00536': 'ACTR a99SB-disp', 'PED00537': 'drkN C36m', 'PED00538': 'drkN C22*',
          'PED00539': 'drkN a99SB-disp', 'PED00540': 'PaaA2 C36m', 'PED00541': 'PaaA2 C22*', 'PED00542': 'PaaA2 a99SB-disp',
          'PED00543': 'α-synuclein C36m', 'PED00544': 'α-synuclein C22*', 'PED00545': 'α-synuclein a99SB-disp'},
}
SIZE_AS = {'x': 'S', 'y': 'T', 'z': 'Y'}


def files(entry):
    """All deposited PDB files of an entry, newest layout (<entry>/<ens>/*.pdb) first, then the older flat layout."""
    fs = sorted(f for f in glob.glob(os.path.join(pn.PED, entry, 'e*', '*.pdb')) if not os.path.basename(f).startswith('._'))
    if not fs:
        fs = sorted(f for f in glob.glob(os.path.join(pn.PED, entry, '*.pdb')) if not os.path.basename(f).startswith('._'))
    return fs


def first_seq_and_count(path):
    """Sequence of the first model and the number of models, with the same residue mapping as ped_null.load."""
    seq, n_models, first, cur_has = [], 0, True, False
    for line in open(path, errors='replace'):
        if line.startswith('MODEL'):
            if cur_has:
                n_models += 1
                first = False
            cur_has = False
        elif line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == 'CA' and line[16] in ' A':
            cur_has = True
            if first:
                rn = line[17:20]
                seq.append(pn.PHOS.get(rn) or pn.ea.AA3_VARIANTS.get(rn, 'X'))
        elif line.startswith('ENDMDL') and cur_has:
            n_models += 1
            first, cur_has = False, False
    if cur_has:
        n_models += 1
    return ''.join(seq), n_models


def training_kmers(k=12):
    seqs = {}
    for d in ('hub_frames_v3_cal', 'hub_frames_v3_mpipi', 'hub_frames_v3_recharged'):
        for r in csv.DictReader(open(os.path.join(HERE, d, 'per_residue.csv'))):
            seqs.setdefault(r['UniProt_ID'], []).append(r['aa'])
    km = set()
    for s in seqs.values():
        s = ''.join(s)
        km.update(s[i:i + k] for i in range(len(s) - k + 1))
    return km


def main():
    accs = json.load(open(os.path.join(AI, 'entry_accs.json')))
    km = training_kmers()
    rows, out = [], []
    for st, ents in SETS.items():
        for e, prot in ents.items():
            fs = files(e)
            if not fs:
                out.append(f'{st} {e} {prot}: no files')
                continue
            info = [(f, *first_seq_and_count(f)) for f in fs]
            seqs = {s for _, s, _ in info}
            keep = [(f, s, n) for f, s, n in info if s == info[0][1]]
            seq = keep[0][1]
            plain = ''.join(SIZE_AS.get(a, a) for a in seq)
            ai_mean = None
            for a in accs.get(e, []):
                fa = os.path.join(AI, f'{a}.fasta')
                fj = os.path.join(AI, f'{a}.json')
                if not (os.path.exists(fa) and os.path.exists(fj)):
                    continue
                useq = ''.join(open(fa).read().split('\n')[1:])
                sc = json.load(open(fj))['AIUPred']
                sm = difflib.SequenceMatcher(None, plain, useq, autojunk=False)
                idx = [blk.b + k for blk in sm.get_matching_blocks() for k in range(blk.size)]
                if len(idx) >= 0.8 * len(plain):
                    ai_mean = sum(sc[i] for i in idx) / len(idx)
                    break
            shared = sum(plain[i:i + 12] in km for i in range(len(plain) - 11))
            n_models = sum(n for _, _, n in keep)
            rel = [os.path.relpath(f, pn.PED) for f, _, _ in keep]
            rows.append({'set': st, 'id': e, 'protein': prot, 'files': ';'.join(rel), 'n': len(seq), 'models': n_models,
                         'aiupred_mean': f'{ai_mean:.3f}' if ai_mean is not None else '',
                         'shared_12mers': shared, 'sequence': seq})
            flag = []
            if len(seqs) > 1:
                flag.append(f'{len(info) - len(keep)} file(s) with a different sequence dropped')
            if st in 'AB' and ai_mean is not None and ai_mean < 0.5:
                flag.append('AIUPred < 0.5: excluded from primary')
            if shared:
                flag.append(f'shares {shared} 12-mers with training: related')
            out.append(f'{st} {e} {prot:22s} N={len(seq):4d} models={n_models:6d} files={len(keep)} '
                       f'AIUPred={ai_mean if ai_mean is None else round(ai_mean, 2)} ' + '; '.join(flag))
    with open(os.path.join(HERE, 'ped_v2_set.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    open(os.path.join(HERE, 'ped_v2_build_summary.txt'), 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
