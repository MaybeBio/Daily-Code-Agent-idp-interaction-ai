"""
design_glycine.py — the glycine sub-batch: does a glycine-rich neighbourhood make a residue a hub, causally?

Bases are proteins with complete labels in CALVADOS-2 and Mpipi-Recharged (data_v3) and 60-160 residues. A base needs
a shared hub h, i.e. a residue in the top 10% of position-adjusted occupancy in both force fields (and in Mpipi-GG when
the protein has Mpipi-GG labels) and at least 8 residues from either end, whose +-5 window holds k >= 3 glycines besides h, and a far window
of 11 residues whose centre is at least 15 residues from h and 6 from either end and which holds at least k glycines
and k exchange partners (A S T N Q D E K). Selection uses the simulation labels only, never a model prediction. One base
per protein (its hub with the most window glycines, then the best worst-case rank); bases are taken in order of k, then
of that rank, then of ID, with at most 4 per kingdom, up to 24. Bases without Mpipi-GG labels need one extra
unmodified (WT) run in Mpipi-GG, listed in hub_design_glycine_wt_mpipigg_300K.fasta.
Variants per base (the natural sequence is already simulated in all three force fields with the same protocol):
  G2S_near   the k glycines of the hub window -> S
  G2S_far    the k glycines of the far window closest to its centre -> S (same composition change, away from the hub)
  G2A_near   the k glycines of the hub window -> A
  Gmove      the k glycines of the hub window exchange places with the k partner residues of the far window closest to
             its centre, so composition is unchanged and the glycines move away from the hub
  P_near     the residues at h-2 and h+2 -> P (h-3 / h+3 when those are G, P or aromatic; skipped if neither works)
The hub residue itself is never changed. Positions are 0-based, as in hub_design_table.csv.
  ~/.venvs/bender/bin/python design_glycine.py  ->  hub_design_glycine_300K.fasta, hub_design_glycine_table.csv,
                                                     hub_design_glycine_seqs.csv (WT + variants, for predict_frozen.py),
                                                     hub_design_glycine_wt_mpipigg_300K.fasta
"""
import csv
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
sys.path.insert(0, EV)
sys.path.insert(0, os.path.join(EV, 'hub_predictor'))
import arrangement_test as AT  # noqa: E402
import train_hub_model as T  # noqa: E402

FFS = ('CALVADOS-2', 'MPIPI-GG', 'Mpipi-Recharged')
PARTNERS = set('ASTNQDEK')
NOT_FOR_P = set('GPFWY')
MAX_BASES, MAX_PER_GROUP = 24, 4


def shared_hubs(seq, occ):
    n = len(seq)
    k10 = max(1, int(round(0.1 * n)))
    ranks = []
    for ff in occ:
        r = AT.posadj(np.asarray(occ[ff], float))
        order = np.argsort(-r)
        rk = np.empty(n, int)
        rk[order] = np.arange(n)
        ranks.append(rk)
    worst = np.max(ranks, axis=0)
    return [(h, int(worst[h])) for h in range(8, n - 8) if worst[h] < k10]


def design_for(seq, h):
    n = len(seq)
    near = [i for i in range(h - 5, h + 6) if i != h and seq[i] == 'G']
    k = len(near)
    if k < 3:
        return None
    best = None
    for c in range(6, n - 6):
        if abs(c - h) < 15:
            continue
        win = [i for i in range(c - 5, c + 6) if 0 <= i < n]
        g = sorted((i for i in win if seq[i] == 'G'), key=lambda i: (abs(i - c), i))
        part = sorted((i for i in win if seq[i] in PARTNERS), key=lambda i: (abs(i - c), i))
        if len(g) >= k and len(part) >= k:
            key = (-len(g), -abs(c - h), c)          # most glycines, then farthest from the hub
            if best is None or key < best[0]:
                best = (key, c, g[:k], part[:k])
    if best is None:
        return None
    _, c, far_g, partners = best
    p_pos = []
    for d in (2, 3):
        cand = [h - d, h + d]
        if all(seq[i] not in NOT_FOR_P for i in cand):
            p_pos = cand
            break
    return {'k': k, 'near': near, 'far_centre': c, 'far': far_g, 'partners': partners, 'p_pos': p_pos}


def mutate(seq, changes):
    s = list(seq)
    for i, a in changes.items():
        s[i] = a
    return ''.join(s)


def main():
    prots, y = T.load(os.path.join(EV, 'hub_predictor', 'data_v3'), 'occupancy')
    group = {k[1]: r['group'] for k, r in prots.items() if k[0] == 'CALVADOS-2'}      # kingdom of every protein
    complete = {ff: {k[1] for k in prots if k[0] == ff and len(y.get(k, {})) == len(prots[k]['sequence'])} for ff in FFS}
    ids = sorted(complete['CALVADOS-2'] & complete['Mpipi-Recharged'])
    cands = []
    for u in ids:
        have = [ff for ff in FFS if u in complete[ff]]
        seqs = {prots[(ff, u)]['sequence'] for ff in have}
        if len(seqs) != 1:
            continue
        seq = seqs.pop()
        if not 60 <= len(seq) <= 160:
            continue
        occ = {ff: [y[(ff, u)][i] for i in range(len(seq))] for ff in have}
        best = None
        for h, worst in shared_hubs(seq, occ):
            d = design_for(seq, h)
            if d is not None:
                key = (-d['k'], worst, h)
                if best is None or key < best[0]:
                    best = (key, h, worst, d)
        if best is not None:
            cands.append((best[0][0], best[2], u, best[1], seq, best[3]))
    cands.sort()
    chosen, per_group = [], Counter()
    for negk, worst, u, h, seq, d in cands:
        g = group.get(u, '?')
        if per_group[g] >= MAX_PER_GROUP:
            continue
        chosen.append((u, g, h, seq, d, worst))
        per_group[g] += 1
        if len(chosen) == MAX_BASES:
            break
    print(f'{len(ids)} proteins with CALVADOS-2 and Mpipi-Recharged labels; {len(cands)} qualify; {len(chosen)} chosen; '
          f'per kingdom {dict(per_group)}; with Mpipi-GG labels {sum(u in complete["MPIPI-GG"] for u, *_ in chosen)}')
    rows, fasta, seqs_csv, wt_gg = [], [], [], []
    for u, g, h, seq, d, worst in chosen:
        k = d['k']
        variants = {
            'G2S_near': ({i: 'S' for i in d['near']}, f'{k} glycines within 5 of hub {h} -> S'),
            'G2S_far': ({i: 'S' for i in d['far']}, f'{k} glycines of the far window centred at {d["far_centre"]} -> S'),
            'G2A_near': ({i: 'A' for i in d['near']}, f'{k} glycines within 5 of hub {h} -> A'),
            'Gmove': ({**{i: seq[j] for i, j in zip(d['near'], d['partners'])}, **{j: 'G' for j in d['partners']}},
                      f'{k} glycines near hub {h} exchanged with residues {" ".join(seq[j] for j in d["partners"])} '
                      f'of the far window centred at {d["far_centre"]}'),
        }
        if d['p_pos']:
            variants['P_near'] = ({i: 'P' for i in d['p_pos']}, f'residues {" ".join(seq[i] for i in d["p_pos"])} at '
                                  f'{d["p_pos"][0]} and {d["p_pos"][1]} flanking hub {h} -> P')
        seqs_csv.append({'id': f'GLY_{u}_WT', 'sequence': seq})
        if u not in complete['MPIPI-GG']:
            wt_gg.append(f'>GLY_{u}_WT|glycine_wt|N={len(seq)};parent={u}\n{seq}')
        for name, (ch, note) in variants.items():
            vid = f'GLY_{u}_{name}'
            vs = mutate(seq, ch)
            assert len(vs) == len(seq) and vs[h] == seq[h]
            if name == 'Gmove':
                assert sorted(vs) == sorted(seq)
            rows.append({'id': vid, 'category': 'glycine', 'n': len(seq), 'parent': u, 'host': '', 'guest': name,
                         'positions': ' '.join(str(i) for i in sorted(ch)),
                         'notes': f'{note}; hub {seq[h]}{h}, worst rank {worst} of {len(seq)}; kingdom {g}', 'sequence': vs})
            fasta.append(f'>{vid}|glycine|N={len(seq)};parent={u}\n{vs}')
            seqs_csv.append({'id': vid, 'sequence': vs})
        print(f'  {u:12s} {g:9s} n={len(seq):3d} hub {seq[h]}{h:<3d} k={k} far centre {d["far_centre"]:3d} '
              f'variants {len(variants)}  window {seq[h - 5:h + 6]}')
    open(os.path.join(HERE, 'hub_design_glycine_300K.fasta'), 'w').write('\n'.join(fasta) + '\n')
    with open(os.path.join(HERE, 'hub_design_glycine_table.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, 'hub_design_glycine_seqs.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['id', 'sequence'])
        w.writeheader()
        w.writerows(seqs_csv)
    open(os.path.join(HERE, 'hub_design_glycine_wt_mpipigg_300K.fasta'), 'w').write('\n'.join(wt_gg) + '\n')
    print(f'{len(rows)} designed sequences, {sum(r["n"] for r in rows)} residues -> hub_design_glycine_300K.fasta; '
          f'{len(wt_gg)} WT runs for Mpipi-GG only -> hub_design_glycine_wt_mpipigg_300K.fasta')


if __name__ == '__main__':
    main()
