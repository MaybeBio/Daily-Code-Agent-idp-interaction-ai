"""Generate ESM-2 reference activations for validating the Rust port.

The Rust implementation in `crates/opendna-ml` is checked against the
reference PyTorch implementation. Rather than depend on Python at test time,
this script is run once and its outputs are committed to
`crates/opendna-ml/tests/data`.

    pip install torch transformers
    python scripts/gen_esm2_reference.py crates/opendna-ml/tests/data

Re-run it only when the set of test sequences or the reference model changes.
"""
import json, os, sys
import numpy as np
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

MODEL = "facebook/esm2_t6_8M_UR50D"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

tok = AutoTokenizer.from_pretrained(MODEL)
model = EsmForMaskedLM.from_pretrained(MODEL)
model.eval()

# Short, chemically varied sequences. The second repeats a residue and the
# third is a single residue, to exercise padding and edge cases.
SEQS = [
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY",
    "M",
]

cfg = model.config
meta = {
    "model": MODEL,
    "vocab_size": cfg.vocab_size,
    "hidden_size": cfg.hidden_size,
    "num_hidden_layers": cfg.num_hidden_layers,
    "num_attention_heads": cfg.num_attention_heads,
    "intermediate_size": cfg.intermediate_size,
    "layer_norm_eps": cfg.layer_norm_eps,
    "position_embedding_type": cfg.position_embedding_type,
    "pad_token_id": cfg.pad_token_id,
    "mask_token_id": cfg.mask_token_id,
    "token_dropout": cfg.token_dropout,
    "emb_layer_norm_before": getattr(cfg, "emb_layer_norm_before", None),
    "vocab": tok.get_vocab(),
    "cases": [],
}

for i, s in enumerate(SEQS):
    enc = tok(s, return_tensors="pt")
    ids = enc["input_ids"]
    with torch.no_grad():
        out = model(**enc, output_hidden_states=True)
    hidden = out.hidden_states[-1][0].numpy().astype(np.float32)
    logits = out.logits[0].numpy().astype(np.float32)
    # Mean-pooled embedding over real residues (excludes <cls> and <eos>),
    # which is what the downstream heads consume.
    pooled = hidden[1:-1].mean(axis=0)
    hidden.tofile(os.path.join(OUT, f"case{i}_hidden.f32"))
    logits.tofile(os.path.join(OUT, f"case{i}_logits.f32"))
    pooled.astype(np.float32).tofile(os.path.join(OUT, f"case{i}_pooled.f32"))
    meta["cases"].append({
        "sequence": s,
        "token_ids": ids[0].tolist(),
        "n_tokens": int(ids.shape[1]),
        "hidden_shape": list(hidden.shape),
        "logits_shape": list(logits.shape),
        "hidden_abs_max": float(np.abs(hidden).max()),
        "logits_abs_max": float(np.abs(logits).max()),
    })
    print(f"case{i}: {len(s)} aa -> {ids.shape[1]} tokens, hidden {hidden.shape}")

with open(os.path.join(OUT, "reference.json"), "w") as f:
    json.dump(meta, f, indent=2, sort_keys=True)
print("wrote", OUT)
