#!/usr/bin/env python3
"""Quantify pairwise label permutation gains without redefining physical states."""
import csv
import itertools
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent
source=OUT.parent/"state_sampling/TDP43_condensate_300K_cross_state_joint.csv"
with source.open() as handle:
    rows=list(csv.DictReader(handle))
matrices={}
for row in rows:
    key=(int(row["attempt_a"]),int(row["attempt_b"]))
    mat=matrices.setdefault(key,np.zeros((5,5)))
    mat[int(row["state_a"]),int(row["state_b"])]=float(row["hard_joint"])
records=[]
for (a,b),mat in matrices.items():
    assert np.isclose(mat.sum(),1)
    candidates=[(sum(mat[i,p[i]] for i in range(5)),p) for p in itertools.permutations(range(5))]
    score,perm=max(candidates,key=lambda x:x[0])
    records.append(dict(attempt_a=a,attempt_b=b,existing_agreement=np.trace(mat),
        optimal_pairwise_permutation_agreement=score,gain=score-np.trace(mat),
        a_to_b_permutation=" ".join(map(str,perm))))
with (OUT/"tdp43_label_permutation_diagnostic.csv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
old=np.array([r["existing_agreement"] for r in records])
new=np.array([r["optimal_pairwise_permutation_agreement"] for r in records])
report=f"""# TDP-43 condensate at high temperature: label-permutation diagnostic

Input: same-frame hard joint assignments from state_sampling, sampled every 10
saved frames. All {len(records)} pairs of 20 models are included.

- Mean hard agreement under the existing structural mappings: {old.mean():.6f}.
- Mean agreement after an optimal pairwise 5-label permutation: {new.mean():.6f}.
- Median gain: {np.median(new-old):.6f}.
- Pairs with a gain greater than 0.1: {int(((new-old)>.1).sum())}/{len(records)}.
- Mean disagreement remaining even after permutation: {(1-new).mean():.6f}.

This is a diagnostic, not a replacement alignment. It optimizes and evaluates on
the same sampled frames; rare-state correspondences can be ambiguous. Pairwise
optimal permutations need not define a cycle-consistent set of 20 mappings and do
not establish correspondence across temperatures. No stored state map, matrix,
population, or inference has been overwritten. Structural contact-map checks,
held-out agreement, and a common representation would be needed before accepting
a new alignment. Agreement gains indicate a label-correspondence contribution to
the disagreement; residual disagreement indicates that a permutation alone cannot
make the partitions identical. Neither result proves splitting of one physical state.
"""
(OUT/"label_permutation_diagnostic.md").write_text(report)
print(report)
