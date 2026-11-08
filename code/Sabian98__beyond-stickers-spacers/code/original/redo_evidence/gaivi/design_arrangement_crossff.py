"""
design_arrangement_crossff.py — sequences for repeating the arrangement test (Fig. 3A) in Mpipi-GG and Mpipi-Recharged.
Plan: arrangement_crossff_plan.md (hash both before any of these runs exists).

Two arms of 150 proteins each, disjoint, chosen deterministically (seed 0):
  shuffle arm  proteins of the CALVADOS-2 shuffle cohort (597); one shuffle each (_global_01), the exact sequence that
               CALVADOS-2 simulated (hub_frames_v3_shuffle)
  Y-to-F arm   proteins of the CALVADOS-2 Y-to-F cohort (800); the exact simulated Y-to-F sequence (hub_frames_v3_yf)
Preference order: wild type already labelled in both Mpipi models, then in Mpipi-Recharged, then the rest. A wild type is
simulated again only where the force field has no label for it. All IDs carry the prefix XFF_ so that no run directory
collides with a BENDER or design-batch run.
Writes (standard library only):
  arrangement_crossff_recharged_300K.fasta   Mpipi-Recharged runs (variants + missing wild types)
  arrangement_crossff_mpipigg_300K.fasta     Mpipi-GG runs (variants + missing wild types)
  arrangement_crossff_table.csv              one row per protein and arm
  python3 design_arrangement_crossff.py
"""
import csv
import os
import random
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
N_ARM = 150


def seqs(dirname):
    s, g = defaultdict(list), {}
    for r in csv.DictReader(open(os.path.join(EV, dirname, 'per_residue.csv'))):
        s[r['UniProt_ID']].append(r['aa'])
        g[r['UniProt_ID']] = r['group']
    return {u: ''.join(v) for u, v in s.items()}, g


def pick(cands, GG, RC, rng, exclude=()):
    cands = sorted(c for c in cands if c not in exclude)
    tiers = [[c for c in cands if c in GG and c in RC], [c for c in cands if c in RC and c not in GG],
             [c for c in cands if c not in RC]]
    out = []
    for t in tiers:
        rng.shuffle(t)
        out += t[:N_ARM - len(out)]
        if len(out) == N_ARM:
            break
    return out


def main():
    rng = random.Random(0)
    CAL, grp = seqs('hub_frames_v3_cal')
    SH, _ = seqs('hub_frames_v3_shuffle')
    YF, _ = seqs('hub_frames_v3_yf')
    GG, _ = seqs('hub_frames_v3_mpipi')
    RC, _ = seqs('hub_frames_v3_recharged')
    sh_c = {u.rsplit('_global_', 1)[0] for u in SH if u.endswith('_global_01')}
    yf_c = {u[:-3] for u in YF}
    sh_arm = pick([u for u in sh_c if u in CAL], GG, RC, rng)
    yf_arm = pick([u for u in yf_c if u in CAL], GG, RC, rng, exclude=set(sh_arm))
    rows, fa = [], {'recharged': [], 'mpipigg': []}
    for arm, prots in (('shuffle', sh_arm), ('YF', yf_arm)):
        for u in prots:
            var_id = f'{u}_global_01' if arm == 'shuffle' else f'{u}_YF'
            var_seq = SH[var_id] if arm == 'shuffle' else YF[var_id]
            wt = CAL[u]
            assert len(var_seq) == len(wt), u
            for ff, have in (('recharged', RC), ('mpipigg', GG)):
                fa[ff].append((f'XFF_{var_id}', f'{arm}_variant', var_seq, u))
                if u not in have:
                    fa[ff].append((f'XFF_{u}_WT', 'wild_type', wt, u))
            rows.append(dict(arm=arm, parent=u, group=grp.get(u, ''), N=len(wt), variant_id=f'XFF_{var_id}',
                             wt_in_recharged=u in RC, wt_in_mpipigg=u in GG,
                             n_changed=sum(a != b for a, b in zip(var_seq, wt))))
    for ff, recs in fa.items():
        with open(os.path.join(HERE, f'arrangement_crossff_{ff}_300K.fasta'), 'w') as fh:
            for i, cat, s, parent in recs:
                fh.write(f'>{i}|{cat}|N={len(s)};parent={parent}\n{s}\n')
    with open(os.path.join(HERE, 'arrangement_crossff_table.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for ff, recs in fa.items():
        nv = sum(c != 'wild_type' for _, c, _, _ in recs)
        print(f'{ff:10s} runs {len(recs):4d} (variants {nv}, wild types {len(recs) - nv}), '
              f'{sum(len(s) for _, _, s, _ in recs):,} residues')
    for arm in ('shuffle', 'YF'):
        a = [r for r in rows if r['arm'] == arm]
        print(f'{arm:8s} {len(a)} proteins; WT in both Mpipi {sum(r["wt_in_recharged"] and r["wt_in_mpipigg"] for r in a)}, '
              f'in Recharged {sum(r["wt_in_recharged"] for r in a)}; median N {sorted(r["N"] for r in a)[len(a) // 2]}, '
              f'max N {max(r["N"] for r in a)}')


if __name__ == '__main__':
    main()
