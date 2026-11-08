"""
hub_dataset.py — build the training/test data for a sequence -> per-residue hub predictor.

Sources (all in ../, produced by hub_frames.py and friends; v2 contact rule, |i-j| >= 3, 0.8 nm, exact betweenness):
  CALVADOS-2       hub_frames/, hub_frames_paired_calvados/, hub_frames_all/, hub_frames_rest/
                   (all 10,282 BENDER sequences; hub_frames_rest/ = the 558 compact chains that fail only the nu-fit
                   r^2 filter; 1,000 frames; 500 / 250 for N > 300 / 600)
  MPIPI-GG         hub_frames_mpipi/                                            (1,000 frames)
  Mpipi-Recharged  hub_frames_recharged/ (289 pilot), hub_frames_recharged_all/ (full BENDER run, labelled on GAIVI
                   with ../gaivi/hub_labels_gaivi.py; identical pipeline and frame counts)
  chemistry-only   sim_null_per_residue_all.csv  (simulation z - Rg-matched residue-neutral null z, where available)
  experiment       ped_null_per_residue.csv      (PED ensembles of disordered proteins: z of experimental ensemble and null)
Writes (to this folder):
  proteins.csv      protein_id, force_field, group, N, fold, frames_used, splithalf_sb (time-averaged betweenness),
                    splithalf_sb_occ (hub occupancy; Spearman-Brown corrected even/odd split), source, sequence
  residues.csv      protein_id, force_field, pos, aa, z_btw (position-adjusted), z_btw_raw, mean_btw, occupancy,
                    burstiness, deg_mean, excess_fixed (blank if no null)
  experimental.csv  entry, protein, restraints, pos, aa, z_exp, z_null, z_exp_size, z_null_size
Folds (0-4) are set by a hash of the protein ID, so a protein is in the same fold in every force field.
Standard library only.  python3 hub_dataset.py [OUT_DIR]
"""
import csv
import os
import sys
import zlib
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
OUT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE   # output folder (default: this folder)
SOURCES = [('CALVADOS-2', 'hub_frames'), ('CALVADOS-2', 'hub_frames_paired_calvados'), ('CALVADOS-2', 'hub_frames_all'),
           ('CALVADOS-2', 'hub_frames_rest'), ('MPIPI-GG', 'hub_frames_mpipi'), ('Mpipi-Recharged', 'hub_frames_recharged'),
           ('Mpipi-Recharged', 'hub_frames_recharged_all')]
if os.environ.get('HUB_LABELS') == 'v3':        # labels recomputed with the code-review fixes of 2026-10-03
    SOURCES = [('CALVADOS-2', 'hub_frames_v3_cal'), ('MPIPI-GG', 'hub_frames_v3_mpipi'),
               ('Mpipi-Recharged', 'hub_frames_v3_recharged'), ('Mpipi-Recharged', 'hub_frames_recharged_all')]


def fold_of(uid):
    return zlib.crc32(uid.encode()) % 5


def main():
    excess = {}
    f = os.path.join(EV, 'sim_null_per_residue_all.csv')
    if os.path.exists(f):
        for r in csv.DictReader(open(f)):
            if r['z_sim_fixed'] and r['z_null_fixed']:
                excess[(r['force_field'], r['UniProt_ID'], int(r['pos']))] = float(r['z_sim_fixed']) - float(r['z_null_fixed'])
    seen = set()
    n_res = defaultdict(int)
    with open(os.path.join(OUT, 'proteins.csv'), 'w', newline='') as fp, \
            open(os.path.join(OUT, 'residues.csv'), 'w', newline='') as fr:
        wp, wr = csv.writer(fp), csv.writer(fr)
        wp.writerow(['protein_id', 'force_field', 'group', 'N', 'fold', 'frames_used', 'splithalf_sb', 'splithalf_sb_occ',
                     'source', 'sequence'])
        wr.writerow(['protein_id', 'force_field', 'pos', 'aa', 'z_btw', 'z_btw_raw', 'mean_btw', 'occupancy', 'burstiness',
                     'deg_mean', 'excess_fixed'])
        for ff, d in SOURCES:
            pf, rf = os.path.join(EV, d, 'per_protein.csv'), os.path.join(EV, d, 'per_residue.csv')
            if not (os.path.exists(pf) and os.path.exists(rf)):
                print(f'skip {d}: not found')
                continue
            meta = {r['UniProt_ID']: r for r in csv.DictReader(open(pf))}
            res = defaultdict(dict)
            for r in csv.DictReader(open(rf)):
                if (ff, r['UniProt_ID']) in seen:
                    continue
                res[r['UniProt_ID']][int(r['pos'])] = r
            for uid, rows in res.items():
                if (ff, uid) in seen or uid not in meta:
                    continue
                seen.add((ff, uid))
                n = len(rows)
                seq = ''.join(rows[i]['aa'] for i in range(n))
                # raw within-protein z of mean betweenness
                b = [float(rows[i]['mean_btw']) for i in range(n)]
                m = sum(b) / n
                sd = (sum((x - m) ** 2 for x in b) / n) ** 0.5 or 1.0
                m_ = meta[uid]
                ho = m_.get('splithalf_r_occupancy', '')
                sb_occ = f'{2 * float(ho) / (1 + float(ho)):.4f}' if ho not in ('', 'nan') and float(ho) > -1 else ''
                wp.writerow([uid, ff, rows[0]['group'], n, fold_of(uid), m_['frames_used'], m_['splithalf_spearman_brown_btw'],
                             sb_occ, d, seq])
                for i in range(n):
                    r = rows[i]
                    ex = excess.get((ff, uid, i))
                    wr.writerow([uid, ff, i, r['aa'], r['z_posadj'], f'{(b[i] - m) / sd:.4f}', r['mean_btw'], r['occupancy'],
                                 r['burstiness'], r['deg_mean'], f'{ex:.4f}' if ex is not None else ''])
                n_res[ff] += n
    f = os.path.join(EV, 'ped_null_per_residue.csv')
    if os.path.exists(f):
        with open(os.path.join(OUT, 'experimental.csv'), 'w', newline='') as fe:
            w = csv.writer(fe)
            w.writerow(['entry', 'protein', 'restraints', 'pos', 'aa', 'z_exp', 'z_null', 'z_exp_size', 'z_null_size'])
            for r in csv.DictReader(open(f)):
                w.writerow([r['entry'], r['protein'], r['restraints'], r['pos'], r['aa'], r['z_exp_fixed'], r['z_null_fixed'],
                            r['z_exp_size'], r['z_null_size']])
    by = defaultdict(int)
    for ff, uid in seen:
        by[ff] += 1
    print('proteins:', dict(by), '| residues:', dict(n_res), '| chemistry-only labels:', len(excess))


if __name__ == '__main__':
    main()
