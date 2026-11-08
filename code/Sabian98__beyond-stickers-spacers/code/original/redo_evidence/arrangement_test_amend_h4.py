"""
arrangement_test_amend_h4.py — amendment to H4 of arrangement_test.py (2026-10-04, after H1-H3 had been seen).
The frozen predictions for the shuffles in the pre-registered run were made from shuffle_sequences.csv of the earlier
BENDER-DYN version, which turned out not to be the shuffles that were simulated (identical in 3 of 1,791; the simulated
sequences are those of each trajectory's top.pdb, which the labels use). This recomputes H4 for the shuffles with the
same hash-verified frozen models on the simulated sequences
(gaivi/predictions_frozen_2026-10-02_arrangement_simseq_2026-10-04.csv). The models have no free parameters, so
predicting after H1-H3 were seen cannot be tuned to the labels. H1-H3 and the Y-to-F parts were unaffected.
  ~/.venvs/bender/bin/python arrangement_test_amend_h4.py  ->  arrangement_test_amend_h4.txt
"""
import csv
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

import arrangement_test as at

HERE = os.path.dirname(os.path.abspath(__file__))
PRED = os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02_arrangement_simseq_2026-10-04.csv')


def main():
    SH = at.read('hub_frames_v3_shuffle')
    CAL = at.read('hub_frames_v3_cal')
    pred = defaultdict(list)
    for r in csv.DictReader(open(PRED)):
        pred[r['id']].append((float(r['consensus']), float(r['cal'])))
    lines = ['H4 (amended): frozen models on the simulated shuffle sequences (median per-protein Spearman, raw / position-partial)']
    for j, model in enumerate(('consensus', 'CALVADOS-2')):
        raw, par, flat = [], [], []
        for u, (s, o) in SH.items():
            if u not in pred or len(pred[u]) != len(o):
                continue
            x = np.array([q[j] for q in pred[u]])
            raw.append(spearmanr(x, o).correlation)
            par.append(at.partial(x, o))
        lines.append(f'  shuffles  {model:11s} n={len(raw):5d}  {np.nanmedian(raw):.2f} / {np.nanmedian(par):.2f}')
    # how variable are the shuffled profiles compared with the natural ones? (signal available to predict)
    sd_s = np.median([at.posadj(o).std() for s, o in SH.values()])
    wts = {u.rsplit('_global_', 1)[0] for u in SH}
    sd_n = np.median([at.posadj(CAL[u][1]).std() for u in wts if u in CAL])
    lines.append(f'  median SD of position-adjusted occupancy: shuffles {sd_s:.4f}, their natural sequences {sd_n:.4f}')
    open(os.path.join(HERE, 'arrangement_test_amend_h4.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
