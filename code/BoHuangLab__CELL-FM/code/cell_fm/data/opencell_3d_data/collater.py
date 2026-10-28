# -*- coding: utf-8 -*-
from typing import List, Union
import torch
from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer


def pad_1d_unsqueeze(x: torch.Tensor, padlen: int, start: int, pad_num: Union[int, float]):
    xlen = x.size(0)
    assert start + xlen <= padlen, f"padlen {padlen} is too small for xlen {xlen} and start {start}"
    new_x = x.new_full([padlen], pad_num, dtype=x.dtype)
    new_x[start: start + xlen] = x
    return new_x.unsqueeze(0)


def collate_fn(
    samples: List[dict],
    vocab: EsmSequenceTokenizer,
    max_protein_sequence_len: int,
    min_protein_sequence_len: int = 0,
):
    samples = [
        s for s in samples
        if s["protein_seq"].size(0) <= (max_protein_sequence_len + 2)
        and s["protein_seq"].size(0) > (min_protein_sequence_len + 2)
    ]

    max_tokens = max(len(s["protein_seq"]) for s in samples)

    batch = {}

    batch["protein_seq"] = torch.cat([
        pad_1d_unsqueeze(s["protein_seq"], max_tokens, 0, vocab.pad_token_id)
        for s in samples
    ])

    batch["protein_seq_masked"] = torch.cat([
        pad_1d_unsqueeze(s["protein_seq_masked"], max_tokens, 0, vocab.pad_token_id)
        for s in samples
    ])

    batch["protein_seq_mask"] = torch.cat([
        pad_1d_unsqueeze(s["protein_seq_mask"], max_tokens, 0, 0)
        for s in samples
    ])

    # images: (1, D, H, W) → unsqueeze → (1, 1, D, H, W) → cat → (B, 1, D, H, W)
    batch["protein_img"] = torch.cat([s["protein_img"].unsqueeze(0) for s in samples])
    batch["nucleus_img"] = torch.cat([s["nucleus_img"].unsqueeze(0) for s in samples])

    batch["zm_label"] = torch.cat([s["zm_label"].unsqueeze(0) for s in samples])
    batch["location_label"] = torch.cat([s["location_label"].unsqueeze(0) for s in samples])
    batch["protein_idx"] = torch.cat([s["protein_idx"].unsqueeze(0) for s in samples])

    return {"batched_data": batch}
