"""Residue-level ProtT5 embeddings (frozen encoder, d1 = 1024) → .npz in the format used by FlexBind.

Full-length sequences are embedded without truncation or chunking.

    python scripts/extract_embeddings.py --fasta proteins.fasta --out data/processed/my_set.npz
    # optional: use the weights downloaded by weights/download_weights.sh
    python scripts/extract_embeddings.py --fasta proteins.fasta --out my_set.npz \
        --state_dict weights/prott5_weights.bin
"""
import argparse
import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch


def read_fasta(path: str) -> List[Tuple[str, str]]:
    records, sid, seq = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None:
                    records.append((sid, "".join(seq)))
                sid, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if sid is not None:
        records.append((sid, "".join(seq)))
    return records


def load_prott5(model_name: str, state_dict: Optional[str], device: str, half: bool):
    from transformers import T5EncoderModel, T5Tokenizer
    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False, legacy=True)
    model = T5EncoderModel.from_pretrained(model_name)
    if state_dict:
        sd = torch.load(state_dict, map_location="cpu")
        missing, _ = model.load_state_dict(sd, strict=False)   # decoder keys are ignored
        if missing:
            print(f"[WARN] {len(missing)} encoder keys missing from {state_dict}")
    model = model.to(device).eval()
    if half and device == "cuda":
        model = model.half()
    for p in model.parameters():
        p.requires_grad_(False)                                  # frozen encoder
    return tokenizer, model


@torch.no_grad()
def embed(records: List[Tuple[str, str]], tokenizer, model, device: str,
          max_residues_per_batch: int = 4000) -> Dict[str, np.ndarray]:
    out = {}
    order = sorted(range(len(records)), key=lambda i: -len(records[i][1]))
    batch: List[int] = []

    def flush(b):
        seqs = [" ".join(re.sub(r"[UZOB]", "X", records[i][1].upper())) for i in b]
        enc = tokenizer(seqs, add_special_tokens=True, padding="longest", return_tensors="pt").to(device)
        h = model(input_ids=enc["input_ids"], attention_mask=enc["attention_mask"]).last_hidden_state
        for k, i in enumerate(b):
            L = len(records[i][1])
            out[records[i][0]] = h[k, :L].float().cpu().numpy()   # drop </s> and padding

    for i in order:
        L = len(records[i][1])
        if batch and (len(batch) + 1) * max(L, len(records[batch[0]][1])) > max_residues_per_batch:
            flush(batch)
            batch = []
        batch.append(i)
    if batch:
        flush(batch)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--model_name", default="Rostlab/prot_t5_xl_uniref50")
    p.add_argument("--state_dict", default=None, help="optional local ProtT5 weights (e.g. weights/prott5_weights.bin)")
    p.add_argument("--max_residues_per_batch", type=int, default=4000)
    p.add_argument("--fp32", action="store_true", help="disable half precision on GPU")
    args = p.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    records = read_fasta(args.fasta)
    tokenizer, model = load_prott5(args.model_name, args.state_dict, device, half=not args.fp32)
    emb = embed(records, tokenizer, model, device, args.max_residues_per_batch)
    ids = [sid for sid, _ in records]
    per_residue = np.empty(len(ids), dtype=object)
    for k, sid in enumerate(ids):
        per_residue[k] = emb[sid]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    np.savez(args.out, ids=np.array(ids), per_residue=per_residue, dim=1024,
             seqs=np.array([s for _, s in records]))
    print(f"Saved {len(ids)} embeddings to {args.out}")


if __name__ == "__main__":
    main()
