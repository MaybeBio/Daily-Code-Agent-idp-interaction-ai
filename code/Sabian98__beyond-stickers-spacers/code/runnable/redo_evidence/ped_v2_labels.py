"""
ped_v2_labels.py — hub occupancy labels and the five-pool residue-neutral null for the PED v2 entries (ped_v2_set.csv),
exactly as exp_occupancy.py and null_multi.py, with the files of an entry pooled before thinning (ped_v2_plan.md).
Refuses to label set A or B unless the frozen v2 predictions and the plan are registered in gaivi/FROZEN.sha256.
  python3 ped_v2_labels.py [SETS=ABC] [K=5]   ->  hub_predictor/ped_v2_occupancy.csv
"""
import csv
import hashlib
import math
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_occupancy as eo  # noqa: E402
import ped_null as pn  # noqa: E402

SETS = sys.argv[1] if len(sys.argv) > 1 else 'ABC'
K = int(sys.argv[2]) if len(sys.argv) > 2 else 5
FROZEN = os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02_ped_v2.csv')
PLAN = os.path.join(HERE, 'ped_v2_plan.md')


def registered(path):
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    return h in open(os.path.join(HERE, 'gaivi', 'FROZEN.sha256')).read()


def load_pooled(files):
    seq, models = None, []
    for f in files.split(';'):
        s, m = pn.load(f)
        if seq is None:
            seq = s
        if s == seq:
            models += m
    if len(models) > 1000:
        models = models[::math.ceil(len(models) / 1000)]
    return seq, models


def one(job):
    j, e, files = job
    seq, models = load_pooled(files)
    oe, used = eo.occupancy(seq, models)
    nulls = []
    for k in range(K):
        nm = pn.null_models(seq, models, 70000 + 97 * k + j)[0]
        nulls.append(eo.occupancy(seq, nm)[0])
    return e, seq, oe, nulls, used


def main():
    rows = [r for r in csv.DictReader(open(os.path.join(HERE, 'ped_v2_set.csv'))) if r['set'] in SETS]
    if any(r['set'] in 'AB' for r in rows):
        if not (os.path.exists(FROZEN) and registered(FROZEN) and registered(PLAN)):
            raise SystemExit('v2 predictions or plan not registered in FROZEN.sha256; refusing to label sets A/B')
    meta = {r['id']: r for r in rows}
    jobs = sorted(((j, r['id'], r['files']) for j, r in enumerate(rows)), key=lambda x: -int(meta[x[1]]['n']))
    res = []
    with Pool(int(os.environ.get('WORKERS', '8'))) as p:
        for n, r in enumerate(p.imap_unordered(one, jobs), 1):
            res.append(r)
            print(f'{n}/{len(jobs)} {r[0]} N={len(r[1])} models={r[4]}', flush=True)
    out = os.path.join(HERE, 'hub_predictor', 'ped_v2_occupancy.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['set', 'entry', 'protein', 'pos', 'aa', 'occ_exp', 'occ_null_mean', 'models']
                   + [f'occ_null_{k}' for k in range(K)])
        for e, seq, oe, nulls, used in sorted(res):
            m = meta[e]
            for i, a in enumerate(seq):
                vals = [o[i] for o in nulls]
                w.writerow([m['set'], e, m['protein'], i, a, f'{oe[i]:.4f}', f'{sum(vals) / len(vals):.4f}', used]
                           + [f'{v:.4f}' for v in vals])
    print(f'{len(res)} entries -> {out}')


if __name__ == '__main__':
    main()
