"""
exp_anchor.py — do the simulations' hubs (CALVADOS-2, MPIPI-GG; consensus vs single-force-field) match experiment?

Simulations: the SAXS-benchmark constructs, 1,000 frames each (hub_frames_saxs_*; and the bead-size rule, run here
for these constructs): α-synuclein SASDQJ7, eIF4G1 1-249 P39935, hnRNPA1 LCD V5 SASDTN2, tau 2N4R SASDLU4.

A. Experimental ensembles (PED, multi-model PDB; Calpha coordinates), hubs computed exactly like the simulations
   (per-model contact graphs, exact betweenness, mean over models, position-adjusted z):
     α-syn  PED00024 PRE-restrained (CALVADOS-2 was fitted partly to α-syn PRE: independent only for MPIPI-GG)
            PED00543 all-atom C36m MD, unweighted (an independent model, not experiment)
     eIF4G1 PED00225 PRE/NOE/RDC/SAXS-restrained (independent of both force fields)
     hnRNPA1 PED00212 (DDR/SAXS/PRE) and PED00215 (PRE), full length; the simulated construct is the LCD V5 variant,
            so only exactly aligned residues are compared and the full-length context differs.
   Metrics: within-protein Spearman of residue centrality (simulation vs ensemble); hit rate of simulation hub classes
   (consensus = top 10% in both force fields, single-force-field hubs) among the ensemble's top 20% (chance 20%);
   mean ensemble centrality of aromatic (F W Y H), acidic (D E) and G/S/A residues.
B. Raw PRE intensity ratios (no ensemble model): α-syn (Ahmed 12 labels; Martini set 5 labels) and tau 2N4R (htau40,
   10 labels), from Pilot_hub_pre_anchor/pre_profiles.csv. Per label and residue with |i - label| >= 10:
   effect = 1 - I_para/I_dia, minus its mean at the same sequence separation (all labels pooled), averaged per residue
   -> experimental long-range contact score. Compared with simulation centrality, hub classes and residue types
   (permutation p-values, 10,000 shuffles). All these proteins were used to fit or test CALVADOS-2.

  python3 exp_anchor.py run       # bead-size per-frame runs for the four constructs (both force fields)
  python3 exp_anchor.py           # analysis -> exp_anchor_summary.txt
"""
import csv
import difflib
import json
import math
import os
import random
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hub_frames as hf  # noqa: E402
import hub_frames_mpipi as hm  # noqa: E402
import hub_frames_saxs as hs  # noqa: E402
import hub_frames_tests as ht  # noqa: E402
import xtc_stdlib as xtc  # noqa: E402
from hub_z_v2_stdlib import zs, quad_resid  # noqa: E402

PED = '/Volumes/User Homes/ZOIDBERG-ICML 2027/experimental_data/raw/pre/ped'
PRE_CSV = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Pilot_hub_pre_anchor/pre_profiles.csv'
CONSTRUCTS = {'SASDQJ7': 'α-synuclein', 'P39935': 'eIF4G1 1-249', 'SASDTN2': 'hnRNPA1 LCD V5', 'SASDLU4': 'tau 2N4R'}
ENSEMBLES = [
    ('SASDQJ7', 'PED00024 α-syn PRE-restrained', 'PED00024/PED00024e001.pdb', 'CALVADOS-2 partly fitted to α-syn PRE'),
    ('SASDQJ7', 'PED00543 α-syn all-atom MD (model)', 'PED00543/pdbfile.pdb', 'independent model, not experiment'),
    ('P39935', 'PED00225 eIF4G1 PRE/NOE/RDC/SAXS', 'PED00225/PED00225e000.pdb', 'independent of both force fields'),
    ('SASDTN2', 'PED00212 hnRNPA1 DDR/SAXS/PRE', 'PED00212/PED00212e000.pdb', 'full length vs simulated LCD variant'),
    ('SASDTN2', 'PED00215 hnRNPA1 PRE', 'PED00215/PED00215e000.pdb', 'full length vs simulated LCD variant'),
]
PRE_SETS = [('α-synuclein (Ahmed set, 12 labels)', 'SASDQJ7'), ('α-synuclein (Martini set, 5 labels)', 'SASDQJ7'),
            ('tau 2N4R (htau40, 10 labels)', 'SASDLU4')]
SIM = {('CALVADOS-2', 'fixed'): 'hub_frames_saxs_calvados', ('MPIPI-GG', 'fixed'): 'hub_frames_saxs_mpipi',
       ('CALVADOS-2', 'size'): 'hub_frames_saxs_sizeaware_calvados', ('MPIPI-GG', 'size'): 'hub_frames_saxs_sizeaware_mpipi'}
AROM, ACID, SMALL = set('FWYH'), set('DE'), set('GSA')
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


# ── bead-size runs for the four constructs ───────────────────────────────────
def mpi_saxs_sigma(uid):
    d = os.path.join(hs.MPI_SAXS, uid)
    data = open(os.path.join(d, 'seq.dat')).read()
    types = []
    for line in data.split('Atoms')[1].strip().splitlines():
        t = line.split()
        if len(t) >= 7 and t[0].isdigit():
            types.append((int(t[0]), int(t[2])))
        elif types:
            break
    types = [ty for _, ty in sorted(types)]
    pair = {}
    for m in re.finditer(r'^pair_coeff\s+(\d+)\s+(\d+)\s+wf/cut\s+[\d.eE+-]+\s+([\d.eE+-]+)',
                         open(os.path.join(d, 'run.in')).read(), re.M):
        i, j, s = int(m.group(1)), int(m.group(2)), float(m.group(3))
        pair[(i, j)] = pair[(j, i)] = s
    return [[pair[(a, b)] for b in types] for a in types]


def one_cal_size(job):
    uid, group, c = job
    d = os.path.join(hs.CAL_SAXS, uid)
    try:
        seq = ''.join(hf.AA3.get(line[17:20], 'X') for line in open(os.path.join(d, 'top.pdb'))
                      if line.startswith('ATOM') and line[12:16].strip() == 'CA')
        fh, off, frame, natom, has_cell, n_frames = hf.dcd_open(os.path.join(d, f'{uid}.dcd'))
        start = n_frames // 2
        idx = [start + (k * (n_frames - start)) // hf.FRAMES for k in range(hf.FRAMES)]
        C2, maxcut = ht.c2_matrix(ht.cal_sigma_matrix(seq), c)
        try:
            return hf.analyse(uid, group, seq, lambda k: hf.dcd_frame(fh, off, frame, natom, has_cell, k), n_frames,
                              idx, ref=False, C2=C2, maxcut=maxcut)
        finally:
            fh.close()
    except Exception as e:
        return ('fail', uid, group, repr(e))


def one_mpi_size(job):
    uid, group, c = job
    d = os.path.join(hs.MPI_SAXS, uid)
    try:
        seq = hm.read_sequence(os.path.join(d, 'seq.dat'))
        data = open(os.path.join(d, 'traj.xtc'), 'rb').read()
        offs = xtc.index(data)
        n_frames = len(offs)
        idx = [(k * n_frames) // hf.FRAMES for k in range(hf.FRAMES)]
        C2, maxcut = ht.c2_matrix(mpi_saxs_sigma(uid), c)

        def get_frame(k):
            natoms, step, time, box, X, Y, Z = xtc.decode(data, offs[k])
            return [[v * 10.0 for v in X], [v * 10.0 for v in Y], [v * 10.0 for v in Z]], \
                (box[0] * 10.0, box[4] * 10.0, box[8] * 10.0)
        return hf.analyse(uid, group, seq, get_frame, n_frames, idx, ref=False, C2=C2, maxcut=maxcut)
    except Exception as e:
        return ('fail', uid, group, repr(e))


def run():
    cs = json.load(open(ht.C_FILE))
    jobs = [(u, 'SAXS') for u in CONSTRUCTS]
    hf.main(worker=one_cal_size, jobs=[(u, g, cs['calvados']) for u, g in jobs],
            out=os.path.join(HERE, SIM[('CALVADOS-2', 'size')]))
    hf.main(worker=one_mpi_size, jobs=[(u, g, cs['mpipi']) for u, g in jobs],
            out=os.path.join(HERE, SIM[('MPIPI-GG', 'size')]))


# ── helpers ──────────────────────────────────────────────────────────────────
def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2.0
        i = j + 1
    return r


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb) if sa and sb else float('nan')


def spearman(a, b):
    return pearson(ranks(a), ranks(b))


def sim_results(dirname, uid):
    path = os.path.join(HERE, dirname, 'per_residue.csv')
    if not os.path.exists(path):
        return None
    rows = [r for r in csv.DictReader(open(path)) if r['UniProt_ID'] == uid]
    if not rows:
        return None
    rows.sort(key=lambda r: int(r['pos']))
    return ''.join(r['aa'] for r in rows), [float(r['z_posadj']) for r in rows], [r['class'] for r in rows]


AA3_VARIANTS = {**hf.AA3, 'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'HSD': 'H', 'HSE': 'H', 'HSP': 'H', 'SEP': 'S',
                'TPO': 'T', 'PTR': 'Y', 'MSE': 'M', 'CYX': 'C', 'CYM': 'C', 'ASH': 'D', 'GLH': 'E', 'LYN': 'K'}


def ped_models(path):
    seq, models, cur, first = [], [], None, True
    for line in open(path):
        if line.startswith('MODEL'):
            if cur is not None and cur[0]:      # previous model had no ENDMDL (seen in MMMx files, e.g. PED00494)
                models.append(cur)
                first = False
            cur = ([], [], [])
        elif line.startswith('ATOM') and line[12:16].strip() == 'CA' and line[16] in ' A':
            if cur is None:
                cur = ([], [], [])
            cur[0].append(float(line[30:38]))
            cur[1].append(float(line[38:46]))
            cur[2].append(float(line[46:54]))
            if first:
                seq.append(AA3_VARIANTS.get(line[17:20], 'X'))
        elif line.startswith('ENDMDL') and cur is not None:
            models.append(cur)
            cur, first = None, False
    if cur and cur[0]:
        models.append(cur)
    return ''.join(seq), models


def ensemble_centrality(seq, models, rule, c_cal):
    """Position-adjusted z of mean exact betweenness over models (fixed 0.8 nm rule or CALVADOS-2 bead-size rule)."""
    n = len(seq)
    C2, maxcut = (None, 8.0)
    if rule == 'size':
        C2, maxcut = ht.c2_matrix(ht.cal_sigma_matrix(seq), c_cal)
    tot = [0.0] * n
    used = 0
    for m in models:
        if len(m[0]) != n:
            continue
        adj = hf.contact_lists([list(m[0]), list(m[1]), list(m[2])], n, C2, maxcut)
        b = hf.betweenness(adj, n)
        tot = [x + y for x, y in zip(tot, b)]
        used += 1
    z = zs([x / used for x in tot])
    return quad_resid(z), used


def align(a, b):
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    return [(blk.a + k, blk.b + k) for blk in sm.get_matching_blocks() for k in range(blk.size)]


def top_set(vals, frac):
    k = max(1, round(frac * len(vals)))
    return set(sorted(range(len(vals)), key=lambda i: -vals[i])[:k])


def perm_p(scores, group_idx, n_perm=10000, seed=0):
    """One-sided p that the mean score of group_idx exceeds a random set of the same size."""
    obs = sum(scores[i] for i in group_idx) / len(group_idx)
    rng = random.Random(seed)
    idx = list(range(len(scores)))
    k = len(group_idx)
    ge = sum(1 for _ in range(n_perm) if sum(scores[i] for i in rng.sample(idx, k)) / k >= obs)
    return obs, (ge + 1) / (n_perm + 1)


def classes_for(uid, rule):
    """Simulation residue classes for both force fields under one rule: consensus / single / anti (sim indices)."""
    a = sim_results(SIM[('CALVADOS-2', rule)], uid)
    b = sim_results(SIM[('MPIPI-GG', rule)], uid)
    if not a or not b or a[0] != b[0]:
        return None
    seq = a[0]
    hc = {i for i, c in enumerate(a[2]) if c == 'hub'}
    hm_ = {i for i, c in enumerate(b[2]) if c == 'hub'}
    ac = {i for i, c in enumerate(a[2]) if c == 'anti'}
    am = {i for i, c in enumerate(b[2]) if c == 'anti'}
    return {'seq': seq, 'zc': a[1], 'zm': b[1], 'consensus hub': hc & hm_, 'CALVADOS-2-only hub': hc - hm_,
            'MPIPI-GG-only hub': hm_ - hc, 'consensus anti-hub': ac & am}


# ── analysis ─────────────────────────────────────────────────────────────────
def main():
    c_cal = json.load(open(ht.C_FILE))['calvados']
    p('Experimental anchor for simulation hubs (1,000 frames per force field; hubs = top 10% by position-adjusted betweenness)')
    p('\nA. Experimental (PED) ensembles, hubs computed the same way; ensemble hubs = top 20% of aligned residues (chance 20%)')
    for uid, name, path, note in ENSEMBLES:
        seq_e, models = ped_models(os.path.join(PED, path))
        p(f'\n  {name}: {len(models)} models, {len(seq_e)} residues  [{note}]')
        for rule in ('fixed', 'size'):
            cl = classes_for(uid, rule)
            if not cl:
                p(f'    rule {rule}: simulation results missing')
                continue
            ze, used = ensemble_centrality(seq_e, models, rule, c_cal)
            pairs = align(cl['seq'], seq_e)
            si = [i for i, j in pairs]
            ej = {i: j for i, j in pairs}
            ev = [ze[ej[i]] for i in si]
            hub_e = {si[k] for k in top_set(ev, 0.20)}
            r_c = spearman([cl['zc'][i] for i in si], ev)
            r_m = spearman([cl['zm'][i] for i in si], ev)
            p(f'    rule {rule:5s} ({len(pairs)} aligned residues; ensemble rule matches): Spearman with ensemble '
              f'CALVADOS-2 {r_c:+.2f}, MPIPI-GG {r_m:+.2f}')
            line = '      hits in the ensemble top 20%: '
            for cname in ('consensus hub', 'CALVADOS-2-only hub', 'MPIPI-GG-only hub', 'consensus anti-hub'):
                s = [i for i in cl[cname] if i in ej]
                line += f'{cname} {sum(i in hub_e for i in s)}/{len(s)}   '
            p(line)
            grp = defaultdict(list)
            for i in si:
                aa = cl['seq'][i]
                g = 'aromatic FWYH' if aa in AROM else 'acidic DE' if aa in ACID else 'GSA' if aa in SMALL else 'other'
                grp[g].append(ze[ej[i]])
            p('      ensemble centrality by residue class (mean z, n): ' + '  '.join(
                f'{g} {sum(v) / len(v):+.2f} ({len(v)})' for g, v in sorted(grp.items())))

    p('\nB. Raw PRE (1 - I_para/I_dia, |i - label| >= 10, minus the mean at the same separation; mean per residue)')
    rows = [r for r in csv.DictReader(open(PRE_CSV)) if r['ff'] == 'CALVADOS-2']
    for ds, uid in PRE_SETS:
        pts = [(int(r['label']), int(r['residue']), float(r['exp'])) for r in rows
               if r['data_set'] == ds and r['exp'] not in ('', 'nan') and not math.isnan(float(r['exp']))]
        pts = [(L, res, e) for L, res, e in pts if abs(res - L) >= 10]
        trend = defaultdict(list)
        for L, res, e in pts:
            trend[abs(res - L) // 5].append(1 - e)
        tmean = {k: sum(v) / len(v) for k, v in trend.items()}
        per = defaultdict(list)
        for L, res, e in pts:
            per[res - 1].append((1 - e) - tmean[abs(res - L) // 5])   # 1-based residue -> 0-based position
        score = {i: sum(v) / len(v) for i, v in per.items()}
        p(f'\n  {ds}: {len(pts)} long-range points, {len(score)} residues scored')
        for rule in ('fixed', 'size'):
            cl = classes_for(uid, rule)
            if not cl:
                p(f'    rule {rule}: simulation results missing')
                continue
            ids = [i for i in sorted(score) if i < len(cl['seq'])]
            sv = [score[i] for i in ids]
            p(f'    rule {rule:5s}: Spearman(PRE score, centrality) CALVADOS-2 {spearman(sv, [cl["zc"][i] for i in ids]):+.2f}, '
              f'MPIPI-GG {spearman(sv, [cl["zm"][i] for i in ids]):+.2f}')
            pos = {i: k for k, i in enumerate(ids)}
            line = '      mean PRE score (permutation p, higher = more long-range contact):'
            for cname in ('consensus hub', 'CALVADOS-2-only hub', 'MPIPI-GG-only hub', 'consensus anti-hub'):
                s = [pos[i] for i in cl[cname] if i in pos]
                if len(s) >= 3:
                    obs, pv = perm_p(sv, s)
                    line += f'  {cname} {obs:+.3f} (n {len(s)}, p {pv:.3f})'
            p(line)
        cl = classes_for(uid, 'fixed')
        if cl:
            ids = [i for i in sorted(score) if i < len(cl['seq'])]
            grp = defaultdict(list)
            for i in ids:
                aa = cl['seq'][i]
                g = 'aromatic FWYH' if aa in AROM else 'acidic DE' if aa in ACID else 'GSA' if aa in SMALL else 'other'
                grp[g].append(score[i])
            p('      PRE score by residue class (mean, n): ' + '  '.join(
                f'{g} {sum(v) / len(v):+.3f} ({len(v)})' for g, v in sorted(grp.items())))
    open(os.path.join(HERE, 'exp_anchor_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'run':
        run()
    else:
        main()
