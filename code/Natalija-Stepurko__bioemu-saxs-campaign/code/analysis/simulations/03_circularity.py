"""Circularity check (design §7): each simulated protein against every protein of the CALVADOS training and
validation sets. Local alignment (BLOSUM62, gap open -10, extend -0.5); identity = identical aligned
positions over the length of the shorter sequence. A protein is flagged at identity >= 0.90 to any listed
protein. Writes results/simulations/circularity.csv (one row per simulated protein: best hit overall, flag,
and the sets in which a hit of >= 0.90 occurs)."""

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from common import OUT

THRESHOLD = 0.90
al = Align.PairwiseAligner(
    mode="local",
    substitution_matrix=substitution_matrices.load("BLOSUM62"),
    open_gap_score=-10,
    extend_gap_score=-0.5,
)
CLEAN = str.maketrans({c: "X" for c in "BJOUZ"})


def identity(a: str, b: str) -> float:
    aln = al.align(a.translate(CLEAN), b.translate(CLEAN))[0]
    return aln.counts().identities / min(len(a), len(b))


inp = pd.read_csv(OUT / "inputs.csv")
sets = pd.read_csv(OUT / "calvados_sets.csv")
uniq = sets.drop_duplicates("sequence")[["sequence"]].reset_index(drop=True)
rows = []
for _, r in inp.iterrows():
    ids = uniq.assign(identity=[identity(r["sequence"], s) for s in uniq["sequence"]])
    hit = sets.merge(ids, on="sequence")
    best = hit.sort_values("identity", ascending=False).iloc[0]
    close = hit[hit["identity"] >= THRESHOLD]
    rows.append(
        {
            "label": r["label"],
            "length": r["length"],
            "best_identity": round(float(best["identity"]), 4),
            "best_name": best["name"],
            "best_set": f"{best['set']} {best['role']}",
            "flagged": bool(len(close)),
            "flag_sets": ";".join(
                sorted({f"{s} {ro}" for s, ro in zip(close["set"], close["role"], strict=True)})
            ),
            "flag_names": ";".join(sorted(set(close["name"]))),
        }
    )
out = pd.DataFrame(rows)
out.to_csv(OUT / "circularity.csv", index=False)
print(f"{out['flagged'].sum()} of {len(out)} flagged at identity >= {THRESHOLD}")
print(out[out["flagged"] | (out["best_identity"] > 0.5)].to_string())
