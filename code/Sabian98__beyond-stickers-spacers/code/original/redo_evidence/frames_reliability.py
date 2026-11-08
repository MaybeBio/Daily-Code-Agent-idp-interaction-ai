"""
frames_reliability.py — how many frames does a residue's time-averaged betweenness need? (Gap G2 of methods_paper/PLAN.md.)

100 CALVADOS-2 proteins (random from the 1,198 of hub_frames/), 1,000 production frames each, |i-j| >= 3, 0.8 nm.
Per frame: exact betweenness, and the k = 30 random-pivot approximation used by the original pipeline (networkx-style
scaling n/k). For F in 25 ... 500 frames:
  split-half  r between the means of two disjoint random sets of F frames (exact)
  vs full     r between the mean of F frames and the mean of all 1,000 exact frames
and the original setting (50 frames, k = 30 pivots) against the full exact mean. Medians over proteins.
Standard library only.  python3 frames_reliability.py  ->  frames_reliability_summary.txt
"""
import csv
import math
import os
import random
import sys
from collections import deque
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hub_frames as hf  # noqa: E402
import hub_frames_tests as ht  # noqa: E402

N_PROT, K_PIVOTS, SEED = 100, 30, 11
FS = [25, 50, 100, 200, 500]


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb) if sa and sb else float('nan')


def pivot_betweenness(adj, n, k, rng):
    """Brandes from k random sources, scaled like networkx (normalised, times n / k)."""
    bc = [0.0] * n
    for s in rng.sample(range(n), min(k, n)):
        if not adj[s]:
            continue
        S, preds, sigma, dist = [], {s: []}, {s: 1.0}, {s: 0}
        Q = deque([s])
        while Q:
            v = Q.popleft()
            S.append(v)
            dv = dist[v] + 1
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dv
                    sigma[w] = 0.0
                    preds[w] = []
                    Q.append(w)
                if dist[w] == dv:
                    sigma[w] += sigma[v]
                    preds[w].append(v)
        delta = dict.fromkeys(S, 0.0)
        while S:
            w = S.pop()
            dw = (1.0 + delta[w]) / sigma[w]
            for v in preds[w]:
                delta[v] += sigma[v] * dw
            if w != s:
                bc[w] += delta[w]
    scale = (1.0 / ((n - 1) * (n - 2))) * (n / min(k, n)) if n > 2 else 0.0
    return [x * scale for x in bc]


def one(job):
    uid, group, seed = job
    rng = random.Random(seed)
    try:
        seq, get_frame, n_frames, idx, close = ht.cal_frames(uid, group)
    except Exception:
        return None
    n = len(seq)
    exact, piv = [], []
    for k in idx:
        xyz, box = get_frame(k)
        hf.unwrap(xyz, box)
        adj = hf.contact_lists(xyz, n)
        exact.append(hf.betweenness(adj, n))
        piv.append(pivot_betweenness(adj, n, K_PIVOTS, rng))
    close()
    F = len(exact)
    full = [sum(exact[f][i] for f in range(F)) / F for i in range(n)]
    out = {'uid': uid, 'N': n}
    for nf in FS:
        perm = rng.sample(range(F), 2 * nf)
        a, b = perm[:nf], perm[nf:]
        ma = [sum(exact[f][i] for f in a) / nf for i in range(n)]
        mb = [sum(exact[f][i] for f in b) / nf for i in range(n)]
        out[f'split_{nf}'] = pearson(ma, mb)
        out[f'full_{nf}'] = pearson(ma, full)
    sub = [round(i * F / 50) for i in range(50)]                 # 50 evenly spaced frames, as the original pipeline
    m50p = [sum(piv[f][i] for f in sub) / 50 for i in range(n)]
    m50e = [sum(exact[f][i] for f in sub) / 50 for i in range(n)]
    out['orig_50_pivot_vs_full'] = pearson(m50p, full)
    out['exact_50_vs_full'] = pearson(m50e, full)
    return out


def med(v):
    v = sorted(x for x in v if x == x)
    return v[len(v) // 2] if v else float('nan')


def main():
    prot = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'hub_frames', 'per_protein.csv'))):
        prot[r['UniProt_ID']] = r['group']
    ids = random.Random(SEED).sample(sorted(prot.items()), N_PROT)
    jobs = [(u, g, SEED + k) for k, (u, g) in enumerate(ids)]
    with Pool(9) as pool:
        res = [r for r in pool.imap_unordered(one, jobs) if r]
    lines = [f'Reliability of time-averaged betweenness vs number of frames: {len(res)} CALVADOS-2 proteins '
             f'(median N {med([r["N"] for r in res])}), 1,000 production frames, |i-j| >= 3, 0.8 nm. Medians over proteins.',
             f'  {"frames":>7s} {"split-half r":>13s} {"r vs 1,000-frame mean":>23s}']
    for nf in FS:
        lines.append(f'  {nf:7d} {med([r[f"split_{nf}"] for r in res]):13.2f} {med([r[f"full_{nf}"] for r in res]):23.2f}')
    lines.append(f'\n  original setting (50 evenly spaced frames, k = {K_PIVOTS} pivots) vs 1,000-frame exact mean: '
                 f'median r {med([r["orig_50_pivot_vs_full"] for r in res]):.2f}; '
                 f'50 frames exact: {med([r["exact_50_vs_full"] for r in res]):.2f}')
    open(os.path.join(HERE, 'frames_reliability_summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
