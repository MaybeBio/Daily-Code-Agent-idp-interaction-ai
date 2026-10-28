# -*- coding: utf-8 -*-
import random
from typing import List, Union
import torch
from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer

# How over-long sequences are kept within max_protein_sequence_len.
SEQ_LENGTH_CONTROLS = ('filter', 'crop', 'random_crop')


def pad_1d_unsqueeze(x: torch.Tensor, padlen: int, start: int, pad_num: Union[int, float]):
    xlen = x.size(0)
    assert start + xlen <= padlen, f"padlen {padlen} is too small for xlen {xlen} and start {start}"
    new_x = x.new_full([padlen], pad_num, dtype=x.dtype)
    new_x[start: start + xlen] = x
    return new_x.unsqueeze(0)


def crop_1d(x: torch.Tensor, start: int, length: int):
    """Keep <cls> + x[start : start + length] residues + <eos> from a (<cls>, res..., <eos>) tensor."""
    return torch.cat([x[:1], x[1 + start: 1 + start + length], x[-1:]])


def crop_sample(sample: dict, max_protein_sequence_len: int, random_start: bool):
    """Crop the sequence tensors of one sample to at most max_protein_sequence_len residues.

    All three sequence tensors share the same (<cls>, res..., <eos>) layout, so the same window
    is applied to each to keep residues, masked residues and the mask aligned.
    """
    num_res = sample["protein_seq"].size(0) - 2
    if num_res <= max_protein_sequence_len:
        return sample

    start = random.randint(0, num_res - max_protein_sequence_len) if random_start else 0

    cropped = dict(sample)
    for key in ("protein_seq", "protein_seq_masked", "protein_seq_mask"):
        cropped[key] = crop_1d(sample[key], start, max_protein_sequence_len)
    return cropped


def collate_fn(
    samples: List[dict],
    vocab: EsmSequenceTokenizer,
    max_protein_sequence_len: int,
    min_protein_sequence_len: int = 0,
    seq_length_control: str = 'filter',
):
    """Batch samples, keeping sequence length in (min_protein_sequence_len, max_protein_sequence_len].

    seq_length_control decides what happens to sequences longer than max_protein_sequence_len:
        'filter'       drop them (default, previous behaviour)
        'crop'         keep the N-terminal max_protein_sequence_len residues
        'random_crop'  keep a random window of max_protein_sequence_len residues

    Sequences of at most min_protein_sequence_len residues are dropped under every mode --
    cropping cannot lengthen a sequence.
    """
    if seq_length_control not in SEQ_LENGTH_CONTROLS:
        raise ValueError(
            f"Unknown seq_length_control '{seq_length_control}'; expected one of {SEQ_LENGTH_CONTROLS}"
        )

    samples = [
        s for s in samples
        if s["protein_seq"].size(0) > (min_protein_sequence_len + 2)
    ]

    if seq_length_control == 'filter':
        samples = [
            s for s in samples
            if s["protein_seq"].size(0) <= (max_protein_sequence_len + 2)
        ]
    else:
        samples = [
            crop_sample(s, max_protein_sequence_len, seq_length_control == 'random_crop')
            for s in samples
        ]

    if not samples:
        raise ValueError(
            f"collate_fn dropped every sample in the batch under seq_length_control="
            f"'{seq_length_control}' with min_protein_sequence_len={min_protein_sequence_len} and "
            f"max_protein_sequence_len={max_protein_sequence_len}"
        )

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
