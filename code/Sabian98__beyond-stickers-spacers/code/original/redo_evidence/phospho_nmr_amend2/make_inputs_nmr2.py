"""
make_inputs_nmr2.py — inputs for phospho_nmr_amend2_new_pairs.md. Builds the unmodified and phosphorylated constructs
from the BMRB entity sequences in raw/ (unmodified entry), checks that the phosphorylated entry marks exactly the listed
sites (X in the sequence or SEP/PTR in the comp index) and that the residue there is S or Y, and writes one FASTA per pair
(temperature in the file name) to CALVADOS IDP/input_nmr2/ and OpenMPIPI-runs/input_nmr2/, plus the pTyr residues file.
  python3 make_inputs_nmr2.py
"""
import json
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
CAL = os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP')
MPI = os.path.expanduser('~/Downloads/Scripts/OpenMPIPI-runs')
PHOS = {'S': 'B', 'T': 'J', 'Y': 'Z'}
# name, unmodified entry, phosphorylated entry, sites, T (K), pH, ionic strength (M)
PAIRS = [('aSynS129', 27074, 27075, [129], 273, 7.4, 0.123),
         ('aSynY125', 27076, 27077, [125], 273, 7.4, 0.123),
         ('CD79b', 19650, 19655, [18, 29], 283, 7.0, 0.20),
         ('Tir', 50758, 50759, [69, 89, 98, 126], 283, 6.65, 0.150)]


def meta(i):
    return list(json.load(open(os.path.join(HERE, 'raw', f'meta_{i}.json'))).values())[0]


def seq(v):
    return v['_Entity.Polymer_seq_one_letter_code'][0].replace('\n', '').replace(' ', '').strip()


def main():
    for d in (os.path.join(CAL, 'input_nmr2'), os.path.join(MPI, 'input_nmr2')):
        os.makedirs(d, exist_ok=True)
    res = open(os.path.join(CAL, 'input_phospho', 'residues_phospho.csv')).read().rstrip('\n')
    assert '\nZ,' not in res
    res += '\nZ,PTR,243.16,0.977461144934,0.646,-2.0,0.38\n'
    open(os.path.join(CAL, 'input_nmr2', 'residues_phospho_ptyr.csv'), 'w').write(res)
    rows = []
    for name, u, p, sites, T, ph, ionic in PAIRS:
        vu, vp = meta(u), meta(p)
        wt, ps = seq(vu), seq(vp)
        assert len(wt) == len(ps), name
        assert set(wt) <= set('ACDEFGHIKLMNPQRSTVWY'), (name, set(wt) - set('ACDEFGHIKLMNPQRSTVWY'))
        marked = sorted({k + 1 for k, c in enumerate(ps) if c == 'X'} |
                        {int(i) for i, c in zip(vp.get('_Entity_comp_index.ID', []), vp.get('_Entity_comp_index.Comp_ID', []))
                         if c in ('SEP', 'TPO', 'PTR')})
        assert marked == sites, (name, marked, sites)
        assert all(a == b for k, (a, b) in enumerate(zip(wt, ps)) if k + 1 not in sites), name
        assert all(wt[s - 1] in PHOS for s in sites), name
        ph_seq = list(wt)
        for s in sites:
            ph_seq[s - 1] = PHOS[wt[s - 1]]
        ph_seq = ''.join(ph_seq)
        fa = f'nmr2_{name}_{T}K.fasta'
        txt = ''
        for form, sq in (('WT', wt), ('P', ph_seq)):
            for r in (1, 2, 3):
                txt += f'>NMR2-{name}-{form}-r{r}|{name}\n' + '\n'.join(sq[k:k + 60] for k in range(0, len(sq), 60)) + '\n'
        for d in (os.path.join(CAL, 'input_nmr2'), os.path.join(MPI, 'input_nmr2')):
            open(os.path.join(d, fa), 'w').write(txt)
        rows.append((name, u, p, len(wt), ' '.join(f'{wt[s - 1]}{s}' for s in sites), T, ph, ionic, fa, wt, ph_seq))
        print(f'{name:9s} N {len(wt):3d}  {" ".join(f"{wt[s - 1]}{s}" for s in sites):20s} {T} K pH {ph} I {ionic} -> {fa}')
    with open(os.path.join(HERE, 'pairs_nmr2.csv'), 'w') as fh:
        fh.write('protein,bmrb_unmod,bmrb_phos,N,sites,T_K,pH,ionic_M,fasta,wt_seq,phos_seq\n')
        for r in rows:
            fh.write(','.join(str(x) for x in r) + '\n')
    cal_cmds = [f'bash submit_calvados.sh input_nmr2/{r[8]} --ph {r[6]} --ionic {r[7]} --outroot /data/taseefr/phospho_nmr2 '
                f'--residues input_nmr2/residues_phospho_ptyr.csv' for r in rows]
    mpi_cmds = [f'bash submit_openmpipi.sh input_nmr2/{r[8]} --phospho --ionic {r[7]} --no-chain '
                f'--outroot /general/taseefr/phospho_nmr2_recharged' for r in rows]
    common = ['', '| pair | BMRB unmod / phospho | N | sites | T (K) | pH | I (M) |', '|---|---|---|---|---|---|---|'] + \
             [f'| {r[0]} | {r[1]} / {r[2]} | {r[3]} | {r[4]} | {r[5]} | {r[6]} | {r[7]} |' for r in rows]
    open(os.path.join(CAL, 'input_nmr2', 'README.md'), 'w').write('\n'.join([
        '# NMR case series, Amendment 2 (CALVADOS-2)', '',
        'Inputs for `phospho_nmr_amend2_new_pairs.md` (hashed 2026-10-15). B = pSer, J = pThr, Z = pTyr (Tyr lambda and',
        'sigma, charge -2). Always pass `--residues input_nmr2/residues_phospho_ptyr.csv`; without it the phosphosites',
        'are not phosphorylated. 24 runs.', '', '```bash', 'source ~/calvados_env/bin/activate', 'cd "<CALVADOS IDP dir>"']
        + cal_cmds + ['```'] + common + ['', 'Copy `/data/taseefr/phospho_nmr2/<T>K/NMR2-*` (dcd + top.pdb) back to '
                                          '`BENDER-DYN/phospho_nmr2/<T>K/`.']) + '\n')
    open(os.path.join(MPI, 'input_nmr2', 'README.md'), 'w').write('\n'.join([
        '# NMR case series, Amendment 2 (Mpipi-Recharged)', '',
        'Inputs for `phospho_nmr_amend2_new_pairs.md` (hashed 2026-10-15). Needs run_openmpipi_idp.py 34378cab or later',
        '(pTyr support: Z = Tyr bead with phosphate Yukawa). Always pass `--phospho`. Mpipi-Recharged has no pH. 24 runs.',
        '', '```bash', 'source ~/calvados_env/bin/activate', 'cd "<OpenMPIPI-runs dir>"'] + mpi_cmds + ['```'] + common +
        ['', 'Copy `/general/taseefr/phospho_nmr2_recharged/<T>K/NMR2-*` (dcd + top.pdb) back to '
         '`BENDER-DYN/phospho_nmr2_recharged/<T>K/`.']) + '\n')


if __name__ == '__main__':
    main()
