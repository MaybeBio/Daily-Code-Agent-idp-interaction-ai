"""
hub_frames_summary.py — summary of <dir>/per_residue.csv and per_protein.csv (see hub_frames.py).
Standard library only.  python3 hub_frames_summary.py [dir] [label]  ->  prints and writes <dir>/summary.txt
(dir defaults to hub_frames, label to CALVADOS-2)
"""
import csv
import os
import statistics as st
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, sys.argv[1] if len(sys.argv) > 1 else 'hub_frames')
LABEL = sys.argv[2] if len(sys.argv) > 2 else 'CALVADOS-2'
CLASSES = ['aro', 'pos', 'neg', 'GSA', 'ILVM', 'P', 'other']
AAS = 'GSACYFWNTPQHDRELIVMK'
MIN_HUB_FRAMES = 5
out = []


def p(s=''):
    out.append(s)
    print(s)


def f(x):
    return float(x) if x not in ('', None) else None


def med(v):
    v = [x for x in v if x is not None]
    return st.median(v) if v else float('nan')


def mean(v):
    v = [x for x in v if x is not None]
    return sum(v) / len(v) if v else float('nan')


def q(v, a):
    v = sorted(x for x in v if x is not None)
    return v[int(a * (len(v) - 1))] if v else float('nan')


prot = list(csv.DictReader(open(f'{D}/per_protein.csv')))
res = list(csv.DictReader(open(f'{D}/per_residue.csv')))
fails = list(csv.DictReader(open(f'{D}/failed.csv'))) if os.path.exists(f'{D}/failed.csv') else []
groups = sorted({r['group'] for r in prot})
counts = ', '.join('%s %d' % (g, sum(r['group'] == g for r in prot)) for g in groups)
n_unwrapped = sum(int(r['bonds_unwrapped']) for r in prot)
p(f'Per-frame hubs, {LABEL} 300 K: {len(prot)} proteins ({counts}), '
  f'{len(res)} residues, {prot[0]["frames_used"]} frames each from the production frames; {len(fails)} failed.')
p(f'Chains rebuilt across the periodic box: {n_unwrapped} bond corrections in total.')

p('\n1. Do hubs stay hubs from one snapshot to the next?')
jo = [f(r['jaccard_successive']) for r in prot]
jr = [f(r['jaccard_random_approx']) for r in prot]
p(f'   Overlap (Jaccard) of the top-10% hub sets of successive sampled frames: median {med(jo):.3f} '
  f'(IQR {q(jo, .25):.3f}-{q(jo, .75):.3f}); random sets of the same size: {med(jr):.3f}; '
  f'median ratio {med([a / b for a, b in zip(jo, jr) if a and b]):.2f}')
sh = [f(r['splithalf_r_btw']) for r in prot]
sb = [f(r['splithalf_spearman_brown_btw']) for r in prot]
so = [f(r['splithalf_r_occupancy']) for r in prot]
rv = [f(r['r_vs_v2_npz_mean']) for r in prot]
p(f'   Reliability of a residue\'s time-averaged betweenness, half vs half of the frames: median r {med(sh):.2f} '
  f'(IQR {q(sh, .25):.2f}-{q(sh, .75):.2f}); Spearman-Brown for all frames {med(sb):.2f}; '
  f'occupancy half vs half {med(so):.2f}')
p(f'   Agreement with the stored v2 means (50 frames, k = 30 pivots): median r {med(rv):.2f} '
  f'(IQR {q(rv, .25):.2f}-{q(rv, .75):.2f})')

p('\n2. Cross-validated: class from even frames, hub occupancy measured on odd frames (null = 0.10)')
by = defaultdict(list)
for r in res:
    if r['class_evenframes']:
        by[r['class_evenframes']].append(f(r['occupancy_oddframes']))
for c in ('hub', 'middle', 'anti'):
    v = by[c]
    p(f'   {c:7s} n {len(v):6d}  mean occupancy {mean(v):.3f}  median {med(v):.3f}  '
      f'>25% of frames: {sum(x > 0.25 for x in v) / len(v):.1%}  >50%: {sum(x > 0.5 for x in v) / len(v):.1%}')

p('\n3. How a residue\'s centrality is spread over frames (class from all frames)')
p('   class    participation (frames with betweenness > 0)   burstiness (share from its top 10% of frames; 0.10 = even)   CV')
byc = defaultdict(list)
for r in res:
    byc[r['class']].append(r)
for c in ('hub', 'middle', 'anti'):
    v = byc[c]
    p(f'   {c:7s}  {med([f(r["participation"]) for r in v]):.3f}{"":40s}{med([f(r["burstiness"]) for r in v]):.3f}'
      f'{"":52s}{med([f(r["cv_btw"]) for r in v]):.2f}')

p(f'\n4. Contacts in the frames where a residue is a hub vs the frames where it is not '
  f'(residues with >= {MIN_HUB_FRAMES} hub frames; medians over residues)')
sel = [r for r in res if int(r['n_hubframes']) >= MIN_HUB_FRAMES and r['deg_nonhubframes']]
p(f'   {len(sel)} residues. Contacts per frame: {med([f(r["deg_hubframes"]) for r in sel]):.2f} vs '
  f'{med([f(r["deg_nonhubframes"]) for r in sel]):.2f}; long-range share (|i-j| >= 6): '
  f'{med([f(r["lr_share_hubframes"]) for r in sel]):.2f} vs {med([f(r["lr_share_nonhubframes"]) for r in sel]):.2f}; '
  f'mean |i-j| of contacts: {med([f(r["mean_sep_hubframes"]) for r in sel]):.1f} vs '
  f'{med([f(r["mean_sep_nonhubframes"]) for r in sel]):.1f}')
tot_h = defaultdict(float)
tot_n = defaultdict(float)
for r in res:
    for c in CLASSES:
        tot_h[c] += int(r[f'hub_{c}'])
        tot_n[c] += int(r[f'non_{c}'])
sh_, sn_ = sum(tot_h.values()), sum(tot_n.values())
p('   Partner classes of all contacts made in hub frames vs non-hub frames (share of contacts):')
p('   ' + '  '.join(f'{c} {tot_h[c] / sh_:.3f}/{tot_n[c] / sn_:.3f}' for c in CLASSES))

p('\n5. By residue type (all residues): how often a hub, and what its contacts are when it is one')
p('   aa   n_res  occupancy  partic.  burst.  frac_hub_class  deg hub/non  LR share hub/non  mean|i-j| hub/non  '
  'top partner shift in hub frames')
aa_rows = defaultdict(list)
for r in res:
    aa_rows[r['aa']].append(r)
for aa in AAS:
    v = aa_rows.get(aa, [])
    if not v:
        continue
    s = [r for r in v if int(r['n_hubframes']) >= MIN_HUB_FRAMES and r['deg_nonhubframes']]
    ph = defaultdict(float)
    pn = defaultdict(float)
    for r in v:
        for c in CLASSES:
            ph[c] += int(r[f'hub_{c}'])
            pn[c] += int(r[f'non_{c}'])
    th, tn = sum(ph.values()) or 1, sum(pn.values()) or 1
    shift = sorted(CLASSES, key=lambda c: -(ph[c] / th - pn[c] / tn))[:2]
    p(f'   {aa}  {len(v):6d}   {mean([f(r["occupancy"]) for r in v]):.3f}    {mean([f(r["participation"]) for r in v]):.3f}'
      f'   {med([f(r["burstiness"]) for r in v]):.3f}   {sum(r["class"] == "hub" for r in v) / len(v):.3f}        '
      f'{med([f(r["deg_hubframes"]) for r in s]):.2f}/{med([f(r["deg_nonhubframes"]) for r in s]):.2f}    '
      f'{med([f(r["lr_share_hubframes"]) for r in s]):.2f}/{med([f(r["lr_share_nonhubframes"]) for r in s]):.2f}         '
      f'{med([f(r["mean_sep_hubframes"]) for r in s]):.1f}/{med([f(r["mean_sep_nonhubframes"]) for r in s]):.1f}        '
      + ', '.join(f'{c} {ph[c] / th - pn[c] / tn:+.3f}' for c in shift))

p('\n6. Charged residues: share of contacts with opposite-charge partners, hub frames vs non-hub frames')
for aa, opp in (('E', 'pos'), ('D', 'pos'), ('K', 'neg'), ('R', 'neg')):
    v = aa_rows.get(aa, [])
    h = sum(int(r[f'hub_{opp}']) for r in v)
    n_ = sum(int(r[f'non_{opp}']) for r in v)
    th = sum(int(r['contacts_hubframes']) for r in v) or 1
    tn = sum(int(r['contacts_nonhubframes']) for r in v) or 1
    p(f'   {aa}: {h / th:.3f} vs {n_ / tn:.3f}')

p('\n7. Hub occupancy by group (mean over residues; cross-validated hubs: occupancy on odd frames)')
for g in groups:
    v = [r for r in res if r['group'] == g]
    hv = [f(r['occupancy_oddframes']) for r in v if r['class_evenframes'] == 'hub']
    p(f'   {g:9s} all residues {mean([f(r["occupancy"]) for r in v]):.3f}   hubs (cross-validated) {mean(hv):.3f}')

open(f'{D}/summary.txt', 'w').write('\n'.join(out) + '\n')
