"""
consensus_hubs.py — residues that are hubs (or anti-hubs) in BOTH CALVADOS-2 and MPIPI-GG: how many, which
residue types, and in what sequence context; and whether bead size (shared by both models) could explain them.

Input: hub_frames_paired_calvados/ and hub_frames_mpipi/ per_residue.csv (908 sequences, 1,000 frames each,
class = top / bottom 10% of the protein by position-adjusted time-averaged betweenness).
Enrichment is observed / expected with the expectation from each protein's own composition.
Standard library only.
  python3 consensus_hubs.py                                   # CALVADOS-2 vs MPIPI-GG, 0.8 nm rule
  python3 consensus_hubs.py DIR_A DIR_B LABEL_A LABEL_B OUTNAME   # any two hub_frames outputs
-> prints and writes OUTNAME (default consensus_hubs_summary.txt)
"""
import csv
import math
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ARGS = sys.argv[1:]
CAL = os.path.join(HERE, ARGS[0] if len(ARGS) > 0 else 'hub_frames_paired_calvados')
MPI = os.path.join(HERE, ARGS[1] if len(ARGS) > 1 else 'hub_frames_mpipi')
LA = ARGS[2] if len(ARGS) > 2 else 'CALVADOS-2'
LB = ARGS[3] if len(ARGS) > 3 else 'MPIPI-GG'
OUTNAME = ARGS[4] if len(ARGS) > 4 else 'consensus_hubs_summary.txt'
CAL_PARAMS = '/Users/taseef/Downloads/Scripts/CALVADOS-Slab/input/residues_CALVADOS2.csv'
MPI_RUNS = '/Volumes/User Homes/Nu Project/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi 300k '
AAS = 'GSACYFWNTPQHDRELIVMK'
AROM, CHG = set('FWY'), set('DEKR')
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb)


def load(d):
    prot = defaultdict(dict)
    for r in csv.DictReader(open(f'{d}/per_residue.csv')):
        prot[r['UniProt_ID']][int(r['pos'])] = (r['aa'], r['class'], float(r['z_posadj']), r['group'])
    return prot


def mpipi_sigmas(n_runs=60):
    """sigma_ii (Å) per residue type from the wf/cut pair_coeff lines of a few run.in files."""
    sig = {}
    for u in sorted(os.listdir(MPI_RUNS))[:n_runs]:
        d = os.path.join(MPI_RUNS, u)
        try:
            masses = open(os.path.join(d, 'seq.dat')).read().split('Masses')[1].split('Atoms')[0]
            t2aa = {int(m.group(1)): m.group(2) for m in re.finditer(r'^\s*(\d+)\s+[\d.]+\s+#\s*([A-Z])', masses, re.M)}
            for m in re.finditer(r'^pair_coeff\s+(\d+)\s+(\d+)\s+wf/cut\s+[\d.]+\s+([\d.]+)', open(os.path.join(d, 'run.in')).read(), re.M):
                i, j, s = int(m.group(1)), int(m.group(2)), float(m.group(3))
                if i == j and i in t2aa:
                    sig[t2aa[i]] = s
        except Exception:
            continue
        if len(sig) == 20:
            break
    return sig


def main():
    pc, pm = load(CAL), load(MPI)
    ids = [u for u in pc if u in pm and len(pc[u]) == len(pm[u])]
    cats = ('both_hub', 'cal_only_hub', 'mpi_only_hub', 'both_anti', 'cal_only_anti', 'mpi_only_anti', 'all')
    obs = defaultdict(lambda: defaultdict(int))
    exp = defaultdict(lambda: defaultdict(float))
    ctx = defaultdict(lambda: defaultdict(list))
    n_both, n_exp_both, n_both_anti, n_exp_anti, n_hub = 0, 0.0, 0, 0.0, 0
    for u in ids:
        a, b = pc[u], pm[u]
        n = len(a)
        seq = ''.join(a[i][0] for i in range(n))
        comp = defaultdict(int)
        for ch in seq:
            comp[ch] += 1
        hub_c = {i for i in range(n) if a[i][1] == 'hub'}
        hub_m = {i for i in range(n) if b[i][1] == 'hub'}
        anti_c = {i for i in range(n) if a[i][1] == 'anti'}
        anti_m = {i for i in range(n) if b[i][1] == 'anti'}
        m = len(hub_c)
        n_hub += m
        n_both += len(hub_c & hub_m)
        n_exp_both += m * len(hub_m) / n
        n_both_anti += len(anti_c & anti_m)
        n_exp_anti += len(anti_c) * len(anti_m) / n
        sets = {'both_hub': hub_c & hub_m, 'cal_only_hub': hub_c - hub_m, 'mpi_only_hub': hub_m - hub_c,
                'both_anti': anti_c & anti_m, 'cal_only_anti': anti_c - anti_m, 'mpi_only_anti': anti_m - anti_c,
                'all': set(range(n))}
        for c, s in sets.items():
            for i in s:
                obs[c][seq[i]] += 1
                w = seq[max(0, i - 3):i] + seq[i + 1:i + 4]
                w5 = seq[max(0, i - 5):i + 6]
                ctx[c]['end_dist'].append(min(i, n - 1 - i) / n)
                ctx[c]['aro_w3'].append(sum(ch in AROM for ch in w) / max(1, len(w)))
                ctx[c]['chg_w3'].append(sum(ch in CHG for ch in w) / max(1, len(w)))
                ctx[c]['GS_w3'].append(sum(ch in 'GS' for ch in w) / max(1, len(w)))
                ctx[c]['same_w3'].append(sum(ch == seq[i] for ch in w))
                ctx[c]['ncpr_w5'].append((sum(ch in 'KR' for ch in w5) - sum(ch in 'DE' for ch in w5)) / len(w5))
                ctx[c]['N'].append(n)
            for ch, k in comp.items():
                exp[c][ch] += len(s) * k / n
    p(f'Consensus hubs, {LA} ({os.path.basename(CAL)}) and {LB} ({os.path.basename(MPI)}), same {len(ids)} sequences, '
      f'1,000 frames each (hub / anti-hub = top / bottom 10% by position-adjusted time-averaged betweenness); '
      f'columns "cal" = {LA}, "mpi" = {LB}')
    p(f'  hubs in both: {n_both} of {n_hub} hubs ({n_both / n_hub:.0%}); expected by chance {n_exp_both:.0f} '
      f'({n_both / n_exp_both:.1f}x chance)')
    p(f'  anti-hubs in both: {n_both_anti} ({n_both_anti / n_exp_anti:.1f}x chance)')

    # residue-type profiles (mean position-adjusted z) of each model and their agreement
    za, zb = defaultdict(list), defaultdict(list)
    for u in ids:
        for i, (aa, cl, z, g) in pc[u].items():
            za[aa].append(z)
        for i, (aa, cl, z, g) in pm[u].items():
            zb[aa].append(z)
    ks0 = [aa for aa in AAS if za[aa] and zb[aa]]
    ma = [sum(za[aa]) / len(za[aa]) for aa in ks0]
    mb = [sum(zb[aa]) / len(zb[aa]) for aa in ks0]
    p('  residue-type profile (mean position-adjusted z) ' + LA + ' / ' + LB + ':  '
      + '  '.join(f'{aa} {x:+.2f}/{y:+.2f}' for aa, x, y in zip(ks0, ma, mb)))
    ks1 = [i for i, aa in enumerate(ks0) if aa != 'P']
    p(f'  profile agreement across types r = {pearson(ma, mb):+.2f} ({pearson([ma[i] for i in ks1], [mb[i] for i in ks1]):+.2f} without P)')

    p('\nResidue-type enrichment, observed / expected from each protein\'s composition (1 = neutral)')
    p('  aa   ' + ''.join(f'{c:>14s}' for c in cats[:-1]))
    for aa in AAS:
        p(f'  {aa}    ' + ''.join(f'{obs[c][aa] / exp[c][aa]:14.2f}' if exp[c][aa] > 0 else f'{"-":>14s}' for c in cats[:-1]))

    p('\nSequence context (means): relative distance to the nearer chain end; within +-3: aromatic, charged, G/S share, '
      'same-type neighbours; NCPR within +-5; chain length')
    keys = ['end_dist', 'aro_w3', 'chg_w3', 'GS_w3', 'same_w3', 'ncpr_w5', 'N']
    p('  ' + f'{"":15s}' + ''.join(f'{k:>10s}' for k in keys))
    for c in cats:
        p(f'  {c:15s}' + ''.join(f'{sum(ctx[c][k]) / len(ctx[c][k]):10.3f}' for k in keys))

    # bead size: shared by both models?
    cp = {r['one']: r for r in csv.DictReader(open(CAL_PARAMS))}
    sig_c = {aa: float(cp[aa]['sigmas']) * 10 for aa in AAS}
    lam_c = {aa: float(cp[aa]['lambdas']) for aa in AAS}
    sig_m = mpipi_sigmas()
    ks = [aa for aa in AAS if aa in sig_m]
    enr_both = [obs['both_hub'][aa] / exp['both_hub'][aa] for aa in ks]
    enr_anti = [obs['both_anti'][aa] / exp['both_anti'][aa] for aa in ks]
    p(f'\nBead size (Å): CALVADOS-2 sigma / MPIPI-GG sigma_ii   ' + '  '.join(f'{aa} {sig_c[aa]:.1f}/{sig_m[aa]:.1f}' for aa in ks))
    p(f'  r(CALVADOS-2 sigma, MPIPI-GG sigma) across types = {pearson([sig_c[a] for a in ks], [sig_m[a] for a in ks]):+.2f}')
    p(f'  consensus-hub enrichment vs CALVADOS-2 sigma r = {pearson(enr_both, [sig_c[a] for a in ks]):+.2f}, '
      f'vs MPIPI-GG sigma r = {pearson(enr_both, [sig_m[a] for a in ks]):+.2f}, vs CALVADOS-2 lambda r = '
      f'{pearson(enr_both, [lam_c[a] for a in ks]):+.2f}')
    p(f'  consensus-anti-hub enrichment vs CALVADOS-2 sigma r = {pearson(enr_anti, [sig_c[a] for a in ks]):+.2f}, '
      f'vs MPIPI-GG sigma r = {pearson(enr_anti, [sig_m[a] for a in ks]):+.2f}')
    open(os.path.join(HERE, OUTNAME), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
