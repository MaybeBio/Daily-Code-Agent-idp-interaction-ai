"""
BiLSTM-based multi-task regression model for protein sequence properties.

Predicts 5 properties per sequence:
    - Asphericity
    - Prefactor
    - Scaled_RE   (end-to-end distance)
    - Scaled_RG   (radius of gyration)
    - Scaling_Exponent

Architecture
------------
1. A shared Bidirectional LSTM encoder reads the (CLS + residues) sequence.
   - Tokens are represented as FIXED ONE-HOT vectors (not learned
     embeddings). The LSTM's input_size is therefore forced to equal
     vocab_size.
   - PAD tokens are excluded from computation entirely via
     pack_padded_sequence / pad_packed_sequence (using real sequence
     lengths), so the LSTM never "sees" padding.
   - The PAD row of the one-hot table is zeroed out, so even if padding
     ever leaked into a computation it would contribute no signal (and,
     being a fixed buffer, it can't receive gradient anyway).
2. The sequence representation ("pooled_repr") is built from the final
   forward/backward hidden states of the top LSTM layer, concatenated.
   This plays the same architectural role the CLS-token embedding played
   in the Transformer version.
3. Five separate "property branches" (small MLPs) each take the pooled
   representation and produce:
       (a) a property-specific "label embedding"
       (b) an auxiliary scalar prediction for that property (supervised
           directly)
4. The pooled representation + all property label embeddings are
   concatenated and passed through a final fusion network that predicts
   all properties jointly.
5. Total loss = final joint loss + weighted sum of the auxiliary losses.
"""

import math
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence
from torch.utils.data import Dataset, DataLoader


# ---------------------------------------------------------------------------
# 1. Vocabulary / tokenizer
# ---------------------------------------------------------------------------

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"     # 20 standard residues
EXTRA_TOKENS = "XBZUO"                    # non-standard / ambiguous residues
PAD_TOKEN = "<PAD>"
CLS_TOKEN = "<CLS>"

VOCAB = [PAD_TOKEN, CLS_TOKEN] + list(AMINO_ACIDS) + list(EXTRA_TOKENS)
TOKEN_TO_IDX = {tok: i for i, tok in enumerate(VOCAB)}
PAD_IDX = TOKEN_TO_IDX[PAD_TOKEN]
CLS_IDX = TOKEN_TO_IDX[CLS_TOKEN]
UNK_IDX = TOKEN_TO_IDX["X"]
VOCAB_SIZE = len(VOCAB)


def encode_sequence(seq: str, max_len: int):
    """
    Encode a protein sequence into token ids with a leading CLS token
    and trailing PAD tokens up to max_len.

    Returns:
        input_ids: LongTensor [max_len]
        attention_mask: LongTensor [max_len]  (1 = real token incl. CLS, 0 = PAD)
        length: int (number of real, non-PAD tokens incl. CLS)
    """
    ids = [CLS_IDX] + [TOKEN_TO_IDX.get(aa, UNK_IDX) for aa in seq]
    ids = ids[:max_len]                       # truncate if longer than max_len
    length = len(ids)
    attn = [1] * length

    pad_len = max_len - length
    ids = ids + [PAD_IDX] * pad_len
    attn = attn + [0] * pad_len

    return (
        torch.tensor(ids, dtype=torch.long),
        torch.tensor(attn, dtype=torch.long),
        length,
    )


# ---------------------------------------------------------------------------
# 2. Dataset
# ---------------------------------------------------------------------------

class ProteinPropertyDataset(Dataset):
    """
    Expects a list of dicts (or a DataFrame-like object) with:
        'Sequence': str
        'Asphericity': float
        'Prefactor': float
        'Scaled_RE': float
        'Scaled_RG': float
        'Scaling_Exponent': float
    """

    PROPERTY_NAMES = ["Asphericity", "Prefactor", "Scaled_RE", "Scaled_RG", "Scaling_Exponent"]

    def __init__(self, records, max_len=None):
        self.records = records
        if max_len is None:
            # +1 for the CLS token
            max_len = max(len(r["Sequence"]) for r in records) + 1
        self.max_len = max_len

    def __len__(self):
        return len(self.records)

    def __getitem__(self, idx):
        rec = self.records[idx]
        input_ids, attention_mask, length = encode_sequence(rec["Sequence"], self.max_len)
        labels = torch.tensor(
            [rec[name] for name in self.PROPERTY_NAMES], dtype=torch.float32
        )
        return {
            "input_ids": input_ids,
            "attention_mask": attention_mask,
            "length": length,
            "labels": labels,
        }


def collate_fn(batch):
    """Simple stacking collate (sequences already padded to dataset max_len)."""
    input_ids = torch.stack([b["input_ids"] for b in batch])
    attention_mask = torch.stack([b["attention_mask"] for b in batch])
    lengths = torch.tensor([b["length"] for b in batch], dtype=torch.long)
    labels = torch.stack([b["labels"] for b in batch])
    return {
        "input_ids": input_ids,
        "attention_mask": attention_mask,
        "lengths": lengths,
        "labels": labels,
    }


# ---------------------------------------------------------------------------
# 3. Shared BiLSTM encoder
# ---------------------------------------------------------------------------

class ProteinBiLSTMEncoder(nn.Module):
    """
    Bidirectional multi-layer LSTM encoder.

    Tokens are represented as FIXED one-hot vectors instead of a learned
    nn.Embedding table. Because of this, the LSTM's input_size is forced
    to equal vocab_size -- there's no separate embed_dim hyperparameter
    anymore (projecting the one-hot vectors to an arbitrary embed_dim
    would just reintroduce a learned linear layer, defeating the point).

    Position information in an LSTM is implicit in the recurrence, so no
    separate positional encoding module is needed (same as before).
    """

    def __init__(
        self,
        vocab_size=VOCAB_SIZE,
        hidden_size=256,
        num_layers=2,
        dropout=0.1,
        pad_idx=PAD_IDX,
    ):
        super().__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers

        # input_dim is derived, not a free hyperparameter, since each
        # token is represented by a fixed vocab_size-dim one-hot vector.
        self.input_dim = vocab_size

        one_hot = torch.eye(vocab_size)
        one_hot[pad_idx] = 0.0  # PAD contributes no signal (fixed buffer, no gradient)
        self.register_buffer("one_hot_table", one_hot)  # [vocab_size, vocab_size]

        self.input_dropout = nn.Dropout(dropout)

        self.lstm = nn.LSTM(
            input_size=self.input_dim,
            hidden_size=hidden_size,
            num_layers=num_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if num_layers > 1 else 0.0,
        )

        # output representation dim = hidden_size * 2 (fwd + bwd)
        self.output_dim = hidden_size * 2

    def forward(self, input_ids, attention_mask, lengths):
        """
        input_ids:      [B, L]
        attention_mask: [B, L]  (1 = real token, 0 = PAD)
        lengths:        [B]     number of real (non-PAD) tokens per sample

        Returns:
            pooled_repr: [B, hidden_size * 2]  sequence-level representation
            encoded:     [B, L, hidden_size * 2]  per-token representations
                         (PAD positions are zero-filled)
        """
        x = self.one_hot_table[input_ids]  # [B, L, vocab_size], fixed one-hot vectors
        x = self.input_dropout(x)

        # Pack so the LSTM skips PAD positions entirely (no wasted compute,
        # and PAD contributes nothing to gradients or hidden states).
        packed = pack_padded_sequence(
            x, lengths.cpu(), batch_first=True, enforce_sorted=False
        )
        packed_out, (h_n, c_n) = self.lstm(packed)
        encoded, _ = pad_packed_sequence(
            packed_out, batch_first=True, total_length=input_ids.size(1)
        )  # [B, L, hidden_size * 2]

        # h_n: [num_layers * 2, B, hidden_size] — last layer's final hidden
        # states for the forward and backward directions.
        h_n = h_n.view(self.num_layers, 2, -1, self.hidden_size)
        last_layer_fwd = h_n[-1, 0]  # [B, hidden_size]
        last_layer_bwd = h_n[-1, 1]  # [B, hidden_size]
        pooled_repr = torch.cat([last_layer_fwd, last_layer_bwd], dim=-1)  # [B, hidden*2]

        return pooled_repr, encoded


# ---------------------------------------------------------------------------
# 4. Per-property branch
# ---------------------------------------------------------------------------

class PropertyBranch(nn.Module):
    """Takes the shared pooled representation and produces a property-
    specific 'label embedding' plus an auxiliary scalar prediction."""

    def __init__(self, in_dim, hidden_dim=128, embed_dim=64, dropout=0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, embed_dim),
            nn.GELU(),
        )
        self.head = nn.Linear(embed_dim, 1)

    def forward(self, pooled_repr):
        embed = self.net(pooled_repr)        # [B, embed_dim]
        pred = self.head(embed).squeeze(-1)  # [B]
        return embed, pred


# ---------------------------------------------------------------------------
# 5. Full multi-task model
# ---------------------------------------------------------------------------

class ProteinMultiTaskRegressor(nn.Module):

    PROPERTY_NAMES = ["Asphericity", "Prefactor", "Scaled_RE", "Scaled_RG", "Scaling_Exponent"]

    def __init__(
        self,
        vocab_size=VOCAB_SIZE,
        hidden_size=256,
        num_layers=2,
        dropout=0.1,
        branch_hidden=128,
        branch_embed=256,
        fusion_hidden=256,
        pad_idx=PAD_IDX,
    ):
        super().__init__()
        self.encoder = ProteinBiLSTMEncoder(
            vocab_size=vocab_size,
            hidden_size=hidden_size,
            num_layers=num_layers,
            dropout=dropout,
            pad_idx=pad_idx,
        )
        repr_dim = self.encoder.output_dim  # hidden_size * 2

        self.branches = nn.ModuleDict(
            {
                name: PropertyBranch(repr_dim, branch_hidden, branch_embed, dropout)
                for name in self.PROPERTY_NAMES
            }
        )

        fusion_input_dim = repr_dim + branch_embed * len(self.PROPERTY_NAMES)
        self.fusion = nn.Sequential(
            nn.Linear(fusion_input_dim, fusion_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(fusion_hidden, fusion_hidden // 2),
            nn.GELU(),
            nn.Linear(fusion_hidden // 2, len(self.PROPERTY_NAMES)),
        )

    def forward(self, input_ids, attention_mask, lengths):
        pooled_repr, _ = self.encoder(input_ids, attention_mask, lengths)  # [B, repr_dim]

        embeds = []
        aux_preds = {}
        for name in self.PROPERTY_NAMES:
            embed, pred = self.branches[name](pooled_repr)
            embeds.append(embed)
            aux_preds[name] = pred

        fusion_input = torch.cat([pooled_repr] + embeds, dim=-1)
        final_pred = self.fusion(fusion_input)  # [B, 5] joint prediction

        return {
            "final_pred": final_pred,      # order matches PROPERTY_NAMES
            "aux_preds": aux_preds,        # dict[name] -> [B]
        }


# ---------------------------------------------------------------------------
# 6. Loss function
# ---------------------------------------------------------------------------

def compute_loss(outputs, labels, aux_weight=0.1):
    """
    outputs: dict from model.forward()
    labels:  [B, 5] in the order ProteinMultiTaskRegressor.PROPERTY_NAMES
    """
    mse = nn.MSELoss()

    final_loss = mse(outputs["final_pred"], labels)

    aux_loss = 0.0
    for i, name in enumerate(ProteinMultiTaskRegressor.PROPERTY_NAMES):
        aux_loss = aux_loss + mse(outputs["aux_preds"][name], labels[:, i])
    aux_loss = aux_loss / len(ProteinMultiTaskRegressor.PROPERTY_NAMES)

    total_loss = final_loss + aux_weight * aux_loss
    return total_loss, final_loss, aux_loss
