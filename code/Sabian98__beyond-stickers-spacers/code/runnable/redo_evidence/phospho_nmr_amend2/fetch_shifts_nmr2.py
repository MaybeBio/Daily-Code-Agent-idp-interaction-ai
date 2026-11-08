"""
fetch_shifts_nmr2.py — downloads the assigned chemical shifts (H, N, CA) of the eight BMRB entries of
phospho_nmr_amend2_new_pairs.md from the BMRB API, in the column layout of the ZOIDBERG table used by phospho_nmr_test.py
(bmrb_id, seq_id, comp_id, aa_1letter, cs_H, cs_N, cs_CA). seq_id is the entity residue number (Comp_index_ID), the
numbering of the entity sequence the constructs were built from. Entity 1 and the first assigned shift list only.
Written and hashed before it was run.
  python3 fetch_shifts_nmr2.py   ->   shifts_nmr2.csv, shifts_nmr2_raw/<id>.json
"""
import csv
import json
import os
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
IDS = [27074, 27075, 27076, 27077, 19650, 19655, 50758, 50759]
AA = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
      'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
      'TYR': 'Y', 'VAL': 'V', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y'}


def main():
    raw = os.path.join(HERE, 'shifts_nmr2_raw')
    os.makedirs(raw, exist_ok=True)
    rows = []
    for i in IDS:
        url = f'https://api.bmrb.io/v2/entry/{i}?loop=Atom_chem_shift'
        txt = urllib.request.urlopen(url, timeout=60).read().decode()
        open(os.path.join(raw, f'{i}.json'), 'w').write(txt)
        loop = json.loads(txt)[str(i)]['Atom_chem_shift'][0]
        tags = [t.split('.')[-1] for t in loop['tags']]
        col = {t: k for k, t in enumerate(tags)}
        lists = sorted({r[col['Assigned_chem_shift_list_ID']] for r in loop['data']})
        res = {}
        for r in loop['data']:
            if r[col['Entity_ID']] != '1' or r[col['Assigned_chem_shift_list_ID']] != lists[0]:
                continue
            s = int(r[col['Comp_index_ID']])
            d = res.setdefault(s, {'bmrb_id': i, 'seq_id': s, 'comp_id': r[col['Comp_ID']],
                                   'aa_1letter': AA.get(r[col['Comp_ID']], 'X'), 'cs_H': '', 'cs_N': '', 'cs_CA': ''})
            a = r[col['Atom_ID']]
            if a in ('H', 'N', 'CA'):
                d[f'cs_{a}'] = r[col['Val']]
        rows += [res[s] for s in sorted(res)]
        print(f'{i}: {len(res)} residues, {sum(1 for v in res.values() if v["cs_H"] and v["cs_N"])} with H and N '
              f'(lists {lists})')
    with open(os.path.join(HERE, 'shifts_nmr2.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['bmrb_id', 'seq_id', 'comp_id', 'aa_1letter', 'cs_H', 'cs_N', 'cs_CA'])
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
