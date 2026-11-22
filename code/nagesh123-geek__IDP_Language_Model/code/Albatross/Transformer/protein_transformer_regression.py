"""
Transformer-based multi-task regression model for protein sequence properties.

Predicts 5 properties per sequence:
    - asphericity
    - prefactor
    - scaled RE  (end-to-end distance)
    - scaled RG  (radius of gyration)
    - scaling exponent

Architecture
------------
1. A shared Transformer encoder reads the (CLS + residues + PAD) sequence.
   - PAD tokens are excluded from attention via a key_padding_mask.
   - PAD embeddings never receive gradient (padding_idx in nn.Embedding).
2. The CLS token's output embedding is extracted as the sequence representation.
3. Four separate "property branches" (small MLPs) each take the CLS embedding
   and produce:
       (a) a property-specific "label embedding"
       (b) an auxiliary scalar prediction for that property (supervised directly)
4. The CLS embedding + all 4 label embeddings are concatenated and passed
   through a final fusion network that predicts all 4 properties jointly.
5. Total loss = final joint loss + weighted sum of the 4 auxiliary losses.

"""

import math
import torch
import torch.nn as nn
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
    """
    ids = [CLS_IDX] + [TOKEN_TO_IDX.get(aa, UNK_IDX) for aa in seq]
    ids = ids[:max_len]                       # truncate if longer than max_len
    attn = [1] * len(ids)

    pad_len = max_len - len(ids)
    ids = ids + [PAD_IDX] * pad_len
    attn = attn + [0] * pad_len

    return torch.tensor(ids, dtype=torch.long), torch.tensor(attn, dtype=torch.long)


# ---------------------------------------------------------------------------
# 2. Dataset
# ---------------------------------------------------------------------------

class ProteinPropertyDataset(Dataset):
    """
    Expects a list of dicts (or a DataFrame-like object) with:
        'sequence': str
        'asphericity': float
        'scaled_re': float
        'scaled_rg': float
        'scaling_exponent': float
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
        input_ids, attention_mask = encode_sequence(rec["Sequence"], self.max_len)
        labels = torch.tensor(
            [rec[name] for name in self.PROPERTY_NAMES], dtype=torch.float32
        )
        return {
            "input_ids": input_ids,
            "attention_mask": attention_mask,
            "labels": labels,
        }


def collate_fn(batch):
    """Simple stacking collate (sequences already padded to dataset max_len).
    If you prefer dynamic per-batch padding for efficiency, pad here instead
    of in the Dataset and re-derive attention_mask per batch."""
    input_ids = torch.stack([b["input_ids"] for b in batch])
    attention_mask = torch.stack([b["attention_mask"] for b in batch])
    labels = torch.stack([b["labels"] for b in batch])
    return {"input_ids": input_ids, "attention_mask": attention_mask, "labels": labels}


# ---------------------------------------------------------------------------
# 3. Positional encoding
# ---------------------------------------------------------------------------

class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=1024, dropout=0.1):
        super().__init__()
        self.dropout = nn.Dropout(dropout)

        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float32).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2, dtype=torch.float32) * (-math.log(10000.0) / d_model)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer("pe", pe.unsqueeze(0))  # [1, max_len, d_model]

    def forward(self, x):
        x = x + self.pe[:, : x.size(1), :]
        return self.dropout(x)


# ---------------------------------------------------------------------------
# 4. Shared Transformer encoder
# ---------------------------------------------------------------------------

class ProteinTransformerEncoder(nn.Module):
    def __init__(
        self,
        vocab_size=VOCAB_SIZE,
        d_model=256,
        nhead=8,
        num_layers=4,
        dim_feedforward=512,
        dropout=0.1,
        max_len=1024,
        pad_idx=PAD_IDX,
    ):
        super().__init__()
        self.d_model = d_model
        # padding_idx=pad_idx -> embedding row for PAD is fixed at 0 and
        # never receives a gradient update.
        self.embedding = nn.Embedding(vocab_size, d_model, padding_idx=pad_idx)
        self.pos_encoding = PositionalEncoding(d_model, max_len, dropout)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            batch_first=True,
            activation="gelu",
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

    def forward(self, input_ids, attention_mask):
        """
        input_ids:      [B, L]
        attention_mask: [B, L]  (1 = real token, 0 = PAD)
        """
        # key_padding_mask expects True where positions should be IGNORED
        key_padding_mask = attention_mask == 0

        x = self.embedding(input_ids) * math.sqrt(self.d_model)
        x = self.pos_encoding(x)
        encoded = self.encoder(x, src_key_padding_mask=key_padding_mask)  # [B, L, d_model]

        cls_repr = encoded[:, 0, :]  # CLS token is always position 0
        return cls_repr, encoded


# ---------------------------------------------------------------------------
# 5. Per-property branch
# ---------------------------------------------------------------------------

class PropertyBranch(nn.Module):
    """Takes the shared CLS embedding and produces a property-specific
    'label embedding' plus an auxiliary scalar prediction for that property."""

    def __init__(self, d_model, hidden_dim=128, embed_dim=64, dropout=0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(d_model, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, embed_dim),
            nn.GELU(),
        )
        self.head = nn.Linear(embed_dim, 1)

    def forward(self, cls_repr):
        embed = self.net(cls_repr)          # [B, embed_dim]
        pred = self.head(embed).squeeze(-1)  # [B]
        return embed, pred


# ---------------------------------------------------------------------------
# 6. Full multi-task model
# ---------------------------------------------------------------------------

class ProteinMultiTaskRegressor(nn.Module):

    PROPERTY_NAMES = ["Asphericity", "Prefactor", "Scaled_RE", "Scaled_RG", "Scaling_Exponent"]

    def __init__(
        self,
        vocab_size=VOCAB_SIZE,
        d_model=256,
        nhead=8,
        num_layers=4,
        dim_feedforward=512,
        dropout=0.1,
        max_len=1024,
        branch_hidden=128,
        branch_embed=256,
        fusion_hidden=256,
        pad_idx=PAD_IDX,
    ):
        super().__init__()
        self.encoder = ProteinTransformerEncoder(
            vocab_size=vocab_size,
            d_model=d_model,
            nhead=nhead,
            num_layers=num_layers,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            max_len=max_len,
            pad_idx=pad_idx,
        )

        self.branches = nn.ModuleDict(
            {
                name: PropertyBranch(d_model, branch_hidden, branch_embed, dropout)
                for name in self.PROPERTY_NAMES
            }
        )

        fusion_input_dim = d_model + branch_embed * len(self.PROPERTY_NAMES)
        self.fusion = nn.Sequential(
            nn.Linear(fusion_input_dim, fusion_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(fusion_hidden, fusion_hidden // 2),
            nn.GELU(),
            nn.Linear(fusion_hidden // 2, len(self.PROPERTY_NAMES)),
        )

    def forward(self, input_ids, attention_mask):
        cls_repr, _ = self.encoder(input_ids, attention_mask)  # [B, d_model]

        embeds = []
        aux_preds = {}
        for name in self.PROPERTY_NAMES:
            embed, pred = self.branches[name](cls_repr)
            embeds.append(embed)
            aux_preds[name] = pred

        fusion_input = torch.cat([cls_repr] + embeds, dim=-1)
        final_pred = self.fusion(fusion_input)  # [B, 4] joint prediction

        return {
            "final_pred": final_pred,      # order matches PROPERTY_NAMES
            "aux_preds": aux_preds,        # dict[name] -> [B]
        }


# ---------------------------------------------------------------------------
# 7. Loss function
# ---------------------------------------------------------------------------

def compute_loss(outputs, labels, aux_weight=0.1):
    """
    outputs: dict from model.forward()
    labels:  [B, 4] in the order ProteinMultiTaskRegressor.PROPERTY_NAMES
    """
    mse = nn.MSELoss()

    final_loss = mse(outputs["final_pred"], labels)

    aux_loss = 0.0
    for i, name in enumerate(ProteinMultiTaskRegressor.PROPERTY_NAMES):
        aux_loss = aux_loss + mse(outputs["aux_preds"][name], labels[:, i])
    aux_loss = aux_loss / len(ProteinMultiTaskRegressor.PROPERTY_NAMES)

    total_loss = final_loss + aux_weight * aux_loss
    return total_loss, final_loss, aux_loss




