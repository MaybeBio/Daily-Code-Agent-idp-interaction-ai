"""Sequence-only inference: FASTA (or precomputed .npz) → residue-level probabilities.

No structure and no pre-annotated IDRs are required: FlexBind jointly predicts IDRs and
protein-/RNA-/DNA-binding residues.

    python scripts/predict.py --checkpoint checkpoints/DP81/seed42/best.pt \
        --fasta proteins.fasta --out predictions.tsv
Output TSV columns: id, position (1-based), residue, p_<task> ×4, pred_<task> ×4
(binary predictions use the validation-selected thresholds stored in the checkpoint).
"""
import argparse
import os
import sys

import numpy as np
import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.dataset import TASKS  # noqa: E402
from scripts.evaluate import load_flexbind  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--checkpoint", required=True)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--fasta")
    src.add_argument("--npz", help="precomputed ProtT5 embeddings (ids, per_residue)")
    p.add_argument("--out", required=True)
    p.add_argument("--batch_size", type=int, default=8)
    p.add_argument("--prott5_state_dict", default=None)
    args = p.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    seqs = {}
    if args.fasta:
        from scripts.extract_embeddings import embed, load_prott5, read_fasta
        records = read_fasta(args.fasta)
        tok, enc = load_prott5("Rostlab/prot_t5_xl_uniref50", args.prott5_state_dict, device, half=True)
        emb = embed(records, tok, enc, device)
        del enc
        ids = [sid for sid, _ in records]
        seqs = dict(records)
    else:
        npz = np.load(args.npz, allow_pickle=True)
        ids = [str(x) for x in npz["ids"]]
        emb = {sid: np.asarray(e, dtype=np.float32) for sid, e in zip(ids, npz["per_residue"])}
        if "seqs" in npz.files:
            seqs = dict(zip(ids, map(str, npz["seqs"])))

    model, thresholds = load_flexbind(args.checkpoint, device)
    thresholds = thresholds or {t: 0.5 for t in TASKS}

    with open(args.out, "w") as f:
        f.write("id\tposition\tresidue\t" + "\t".join(f"p_{t}" for t in TASKS) + "\t"
                + "\t".join(f"pred_{t}" for t in TASKS) + "\n")
        for s in range(0, len(ids), args.batch_size):
            chunk = ids[s:s + args.batch_size]
            Ls = [emb[i].shape[0] for i in chunk]
            H = torch.zeros(len(chunk), max(Ls), emb[chunk[0]].shape[1])
            M = torch.zeros(len(chunk), max(Ls), dtype=torch.bool)
            for k, i in enumerate(chunk):
                H[k, :Ls[k]] = torch.from_numpy(emb[i])
                M[k, :Ls[k]] = True
            probs = {t: v.float().cpu().numpy() for t, v in model.predict_proba(H.to(device), M.to(device)).items()}
            for k, i in enumerate(chunk):
                seq = seqs.get(i, "")
                for pos in range(Ls[k]):
                    ps = [float(probs[t][k, pos]) for t in TASKS]
                    preds = [int(v >= thresholds[t]) for v, t in zip(ps, TASKS)]
                    res = seq[pos] if pos < len(seq) else "X"
                    f.write(f"{i}\t{pos + 1}\t{res}\t" + "\t".join(f"{v:.4f}" for v in ps) + "\t"
                            + "\t".join(map(str, preds)) + "\n")
    print(f"Wrote predictions for {len(ids)} proteins to {args.out}")


if __name__ == "__main__":
    main()
