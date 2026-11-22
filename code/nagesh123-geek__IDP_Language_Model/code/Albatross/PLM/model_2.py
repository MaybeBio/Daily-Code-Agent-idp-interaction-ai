"""
Multi-task regression model that takes a PRE-POOLED protein language model (PLM)
embedding (one fixed-length vector per sequence) and predicts 5 properties:
    - asphericity
    - prefactor
    - scaled RE
    - scaled RG
    - scaling exponent

No tokenizer, no CLS/PAD tokens, no transformer encoder is needed here because
the PLM embedding already summarizes the whole sequence into one vector.

Architecture
------------
1. LayerNorm on the raw PLM embedding for stability.
2. Four property branches (small MLPs), each taking the embedding and producing:
       (a) a property-specific "label embedding"
       (b) an auxiliary scalar prediction for that property
3. The (normalized) embedding + all 4 label embeddings are concatenated and
   passed through a fusion network that predicts all 4 properties jointly.
4. Loss = final joint loss + weighted sum of the 4 auxiliary losses.
"""

import torch
import torch.nn as nn



class ResidualBlock(nn.Module):
    def __init__(self, dim, hidden_dim, dropout=0.1):
        super().__init__()

        self.net = nn.Sequential(
            nn.LayerNorm(dim),
            nn.Linear(dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, dim),
            nn.Dropout(dropout),
        )

    def forward(self, x):
        return x + self.net(x)

# ---------------------------------------------------------------------------
# Per-property branch 
# ---------------------------------------------------------------------------
class PropertyBranch(nn.Module):
    def __init__(self, input_dim, hidden_dim, embed_dim, dropout=0.1):
        super().__init__()

        self.encoder = nn.Sequential(
            nn.LayerNorm(input_dim),

            nn.Linear(input_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),

            nn.Linear(hidden_dim, embed_dim),
        )

        self.residual = ResidualBlock(embed_dim, embed_dim * 2, dropout)

        self.predictor = nn.Sequential(
            nn.LayerNorm(embed_dim),
            nn.Linear(embed_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, x):
        emb = self.encoder(x)
        emb = self.residual(emb)
        pred = self.predictor(emb).squeeze(-1)
        return emb, pred


# ---------------------------------------------------------------------------
# Full model
# ---------------------------------------------------------------------------

class EmbeddingMultiTaskRegressor(nn.Module):

    PROPERTY_NAMES = [
        "asphericity",
        "prefactor",
        "scaled_re",
        "scaled_rg",
        "scaling_exponent",
    ]

    def __init__(
        self,
        embedding_dim,
        branch_hidden=128,
        branch_embed=256,
        fusion_hidden=256,
        dropout=0.1,
        normalize_input=True,
    ):
        super().__init__()

        self.embedding_dim = embedding_dim
        self.normalize_input = normalize_input

        if normalize_input:
            self.input_norm = nn.LayerNorm(embedding_dim)

        self.branches = nn.ModuleDict(
            {
                name: PropertyBranch(
                    embedding_dim,
                    branch_hidden,
                    branch_embed,
                    dropout,
                )
                for name in self.PROPERTY_NAMES
            }
        )

        fusion_input_dim = embedding_dim + branch_embed * len(self.PROPERTY_NAMES)

        self.fusion = nn.Sequential(
            nn.LayerNorm(fusion_input_dim),

            nn.Linear(fusion_input_dim, fusion_hidden),
            nn.GELU(),
            nn.Dropout(dropout),

            ResidualBlock(fusion_hidden, fusion_hidden * 2, dropout),
            ResidualBlock(fusion_hidden, fusion_hidden * 2, dropout),

            nn.LayerNorm(fusion_hidden),

            nn.Linear(fusion_hidden, fusion_hidden // 2),
            nn.GELU(),
            nn.Dropout(dropout),

            nn.LayerNorm(fusion_hidden // 2),

            nn.Linear(fusion_hidden // 2, len(self.PROPERTY_NAMES)),
        )

    def forward(self, embedding):

        x = self.input_norm(embedding) if self.normalize_input else embedding

        branch_embeddings = []
        aux_preds = {}

        for name in self.PROPERTY_NAMES:
            emb, pred = self.branches[name](x)
            branch_embeddings.append(emb)
            aux_preds[name] = pred

        fusion_input = torch.cat([x] + branch_embeddings, dim=-1)

        final_pred = self.fusion(fusion_input)

        return {
            "final_pred": final_pred,
            "aux_preds": aux_preds,
        }

# ---------------------------------------------------------------------------
# Loss function (identical logic to the sequence-based version)
# ---------------------------------------------------------------------------

def compute_loss(outputs, labels, aux_weight=0.1):
    mse = nn.MSELoss()

    final_loss = mse(outputs["final_pred"], labels)

    aux_loss = 0.0
    for i, name in enumerate(EmbeddingMultiTaskRegressor.PROPERTY_NAMES):
        aux_loss = aux_loss + mse(outputs["aux_preds"][name], labels[:, i])
    aux_loss = aux_loss / len(EmbeddingMultiTaskRegressor.PROPERTY_NAMES)

    total_loss = final_loss + aux_weight * aux_loss
    return total_loss, final_loss, aux_loss
