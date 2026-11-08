"""
freeze_degree_consensus.py — degree_models_plan.md (1a1deeb9), analysis 4: the consensus degree model, trained on the mean
degree over the force fields that simulated a protein (proteins with labels in >= 2 force fields, same construction as the
consensus model of train_hub_multiff.py), all of data_v3, same features and booster. Pickled; no prediction is scored here.
  ~/.venvs/bender/bin/python freeze_degree_consensus.py   ->   degree_consensus_2026-10-13.pkl
"""
import os
import pickle
from collections import defaultdict

import numpy as np

import train_hub_model as T
from train_hub_multiff import fit_w

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'degree_consensus_2026-10-13.pkl')


def main():
    assert not os.path.exists(OUT), 'exists; refusing to overwrite'
    prots, y = T.load(os.path.join(HERE, 'data_v3'), 'deg_mean')
    seq_of = {}
    for (ff, uid), r in prots.items():
        seq_of.setdefault(uid, r['sequence'])
    by_uid = defaultdict(dict)
    for (ff, uid), lab in y.items():
        if len(lab) == len(seq_of[uid]) and prots[(ff, uid)]['sequence'] == seq_of[uid]:
            by_uid[uid][ff] = np.array([lab[i] for i in range(len(lab))])
    multi = sorted(u for u, d in by_uid.items() if len(d) >= 2)
    X = np.vstack([T.features(seq_of[u], 10) for u in multi])
    Y = np.concatenate([np.mean(list(by_uid[u].values()), axis=0) for u in multi])
    m = fit_w('hgb', X, Y, np.ones(len(Y)))
    pickle.dump(m, open(OUT, 'wb'))
    print(f'consensus degree model: {len(multi)} proteins, {len(Y)} residues -> {OUT}')


if __name__ == '__main__':
    main()
