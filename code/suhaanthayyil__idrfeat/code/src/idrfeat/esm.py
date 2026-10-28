from __future__ import annotations

import numpy as np


def attention_position_features(att) -> dict[str, np.ndarray]:
    att = np.asarray(att, dtype=float)
    received = att.sum(axis=2).mean(axis=(0, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(att > 0, np.log(att), 0.0)
    row_entropy = -(att * logp).sum(axis=3)
    entropy = row_entropy.mean(axis=(0, 1))
    return {"received": received, "entropy": entropy}


def segment_attention_features(att, start: int, end: int) -> dict[str, float]:
    feats = attention_position_features(att)
    sl = slice(start - 1, end)
    return {
        "attn_received_mean": float(feats["received"][sl].mean()),
        "attn_entropy_mean": float(feats["entropy"][sl].mean()),
    }


def segment_embedding_mean(emb, start: int, end: int) -> np.ndarray:
    return np.asarray(emb, dtype=float)[start - 1 : end].mean(axis=0)


def embed_and_attend(
    seq: str,
    model_name: str = "esm2_t33_650M_UR50D",
    repr_layer: int = 33,
    device: str = "cpu",
):
    import esm
    import torch

    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.to(device).eval()
    _, _, tokens = alphabet.get_batch_converter()([("query", seq)])
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[repr_layer], need_head_weights=True, return_contacts=False)
    end = len(seq) + 1
    emb = out["representations"][repr_layer][0, 1:end].cpu().numpy()
    att = out["attentions"][0, :, :, 1:end, 1:end].cpu().numpy()
    return emb, att
