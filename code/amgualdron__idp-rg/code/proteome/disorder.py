"""Disorder of every protein in the human reference proteome, with metapredict V3 (Lotthammer et al. 2024).

    conda run -n disorder python proteome/disorder.py

Needs metapredict (and torch), which live in their own environment, not in idprg's. Writes to proteome/results/:

    disorder_proteins.csv   one row per protein: accession, gene, length, mean disorder score, fraction of residues
                            scoring >= 0.5, fraction of residues inside predicted IDRs, number of IDRs
    disorder_idrs.csv       one row per predicted IDR (metapredict defaults: threshold 0.5, IDRs >= 12 residues):
                            accession, start (0-based), end (exclusive), length, mean score, sequence
    disorder_scores.npz     per-residue scores of all proteins, concatenated: acc, offset (len + 1), score
"""
import pathlib
import re
import sys
import time

import metapredict as meta
import numpy as np
import pandas as pd
from Bio import SeqIO

HERE = pathlib.Path(__file__).resolve().parent
FASTA = HERE.parent / "data" / "UP000005640_9606.fasta"
OUT = HERE / "results"
THRESHOLD = 0.5
BATCH = 32        # small batches: batches are padded to their longest chain (titin, ~35k residues) and must fit the GPU


def main(fasta=FASTA):
    OUT.mkdir(exist_ok=True)
    recs = list(SeqIO.parse(fasta, "fasta"))
    # metapredict accepts the 20 standard residues only: U (selenocysteine) is scored as C and X as G. The IDR
    # sequences are cut from the original sequence, so chains containing U or X are dropped by idprg downstream.
    seqs = {r.id.split("|")[1]: str(r.seq).replace("U", "C").replace("X", "G") for r in recs}
    t0 = time.time()
    res = meta.predict_disorder(seqs, version="V3", return_domains=True, show_progress_bar=True,
                               batch_size=BATCH)
    print(f"metapredict V3: {len(seqs)} proteins in {time.time() - t0:.0f} s", file=sys.stderr)

    prot, idrs, scores = [], [], []
    for r in recs:
        db, acc, entry = r.id.split("|")
        gene = re.search(r" GN=(\S+)", r.description)
        o = res[acc]
        d = np.asarray(o.disorder, np.float32)
        bounds = o.disordered_domain_boundaries
        in_idr = sum(e - s for s, e in bounds)
        prot.append(dict(accession=acc, entry=entry, reviewed=db == "sp", gene=gene.group(1) if gene else "",
                         description=r.description.split(" ", 1)[1].split(" OS=")[0], length=len(d),
                         mean_disorder=float(d.mean()), frac_disordered=float(np.mean(d >= THRESHOLD)),
                         frac_in_idr=in_idr / len(d), n_idr=len(bounds)))
        for k, (s, e) in enumerate(bounds):
            idrs.append(dict(accession=acc, idr=k + 1, start=s, end=e, length=e - s,
                             mean_disorder=float(d[s:e].mean()), sequence=str(r.seq)[s:e]))
        scores.append(d)

    pd.DataFrame(prot).to_csv(OUT / "disorder_proteins.csv", index=False, float_format="%.4g")
    pd.DataFrame(idrs).to_csv(OUT / "disorder_idrs.csv", index=False, float_format="%.4g")
    np.savez_compressed(OUT / "disorder_scores.npz", acc=np.array([p["accession"] for p in prot]),
                        offset=np.concatenate([[0], np.cumsum([len(s) for s in scores])]), score=np.concatenate(scores))
    print(f"{len(prot)} proteins, {len(idrs)} IDRs -> {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:])
