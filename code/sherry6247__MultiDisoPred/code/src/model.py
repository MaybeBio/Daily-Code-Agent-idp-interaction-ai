'''
Date: 2026-03-02 21:02:21
LastEditors: Sicen Liu
LastEditTime: 2026-11-01 14:43:39
FilePath: /liusicen/mygithub/MultiDisoPred/src/model.py
Description:  
Copyright: © 2025 liusicen@smbu.edu.cn. All rights reserved.
'''
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence


# -------------------------
# Utils
# -------------------------
def masked_mean(x: torch.Tensor, mask: Optional[torch.Tensor]) -> torch.Tensor:
    """
    x: (B, L, D)
    mask: (B, L) True=valid, False=PAD
    """
    if mask is None:
        return x.mean(dim=1)
    m = mask.float().unsqueeze(-1)  # (B,L,1)
    denom = m.sum(dim=1).clamp_min(1.0)
    return (x * m).sum(dim=1) / denom


def rand_like(x: torch.Tensor, scale: float = 1.0) -> torch.Tensor:
    return torch.randn_like(x) * scale


# -------------------------
# Modality Encoder (simple BiLSTM)
# -------------------------
class BiLSTMEncoder(nn.Module):
    def __init__(self, in_dim: int, d_model: int, dropout: float = 0.1):
        super().__init__()
        self.proj = nn.Linear(in_dim, d_model)
        self.ln = nn.LayerNorm(d_model)
        self.lstm = nn.LSTM(
            input_size=d_model,
            hidden_size=d_model // 2,
            num_layers=1,
            batch_first=True,
            bidirectional=True,
        )
        self.drop = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, seq_mask: Optional[torch.Tensor] = None) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        x: (B, L, in_dim)
        return:
          h: (B, L, D)
          g: (B, D) pooled seq embedding
        """
        # padding 位置的原始特征在 collate 中被填 0，据此得到每条序列的真实长度，
        # 避免 padding 进入 LSTM 污染双向（尤其反向）隐藏状态。
        lengths = (x.abs().sum(dim=-1) > 0).sum(dim=1).clamp(min=1).long()

        h = self.drop(self.ln(self.proj(x)))
        packed = pack_padded_sequence(h, lengths.cpu(), batch_first=True, enforce_sorted=False)
        packed_h, _ = self.lstm(packed)
        h, _ = pad_packed_sequence(packed_h, batch_first=True, total_length=h.size(1))
        g = masked_mean(h, seq_mask)
        return h, g


# -------------------------
# Re-weighting model W (sample-level weights)
# -------------------------
class ReWeightMLP(nn.Module):
    def __init__(self, in_dim: int, num_experts: int, hidden: int = 256, depth: int = 2, dropout: float = 0.1, tau: float = 1.0):
        super().__init__()
        self.tau = tau
        layers = []
        d = in_dim
        for _ in range(depth - 1):
            layers += [nn.Linear(d, hidden), nn.ReLU(), nn.Dropout(dropout)]
            d = hidden
        layers += [nn.Linear(d, num_experts)]
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: (B, in_dim)
        return weights: (B, E) softmax
        """
        logits = self.net(x) / max(self.tau, 1e-6)
        return F.softmax(logits, dim=-1)


# -------------------------
# Fusion inside each interaction expert
# (multi-modal attention over modalities per residue)
# -------------------------
class MultiModalFusion(nn.Module):
    def __init__(self, d_model: int, dropout: float = 0.1):
        super().__init__()
        self.Wq = nn.Linear(d_model, d_model, bias=False)
        self.Wk = nn.Linear(d_model, d_model, bias=False)
        self.ffn = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, d_model * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model * 4, d_model),
        )
        self.dropout = nn.Dropout(dropout)

    def forward(self, *hs: torch.Tensor, seq_ctx: torch.Tensor, seq_mask: Optional[torch.Tensor] = None) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        n_mod = len(hs)
        B, L, D = hs[0].shape
        if n_mod == 1:
            z = masked_mean(hs[0], seq_mask)
            attn_w = torch.ones(B, L, 1, device=hs[0].device, dtype=hs[0].dtype)
            return hs[0], z, attn_w
        
        H = torch.stack(hs, dim=2)  # (B, L, n_mod, D)
        q = self.Wq(seq_ctx).unsqueeze(1).unsqueeze(2)  # (B, 1, 1, D)
        k = self.Wk(H)  # (B, L, n_mod, D)
        scores = torch.matmul(q, k.transpose(-2, -1)) / math.sqrt(D)  # (B, L, 1, n_mod)
        attn_w = F.softmax(scores, dim=-1)  # (B, L, 1, n_mod)
        out = torch.matmul(attn_w, H).squeeze(2)  # (B, L, D)
        out = out + self.dropout(self.ffn(out))
        z = masked_mean(out, seq_mask)
        return out, z, attn_w


# -------------------------
# Interaction Expert
# -------------------------
class InteractionExpert(nn.Module):
    def __init__(self, d_model: int, dropout: float = 0.1, proj_dim: int = 128):
        super().__init__()
        self.fusion = MultiModalFusion(d_model, dropout)
        self.head = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, 1),
        )
        self.proj = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Linear(d_model, proj_dim),
        )

    def forward(self, *hs, seq_ctx=None, seq_mask=None):
        fused, z, attn_w = self.fusion(*hs, seq_ctx=seq_ctx, seq_mask=seq_mask)
        logits = self.head(fused).squeeze(-1)
        emb = self.proj(z)
        return {"logits": logits, "emb": emb, "attn_w": attn_w}


# -------------------------
# Residue-I2MoE main model (n modalities -> 5 interaction experts)
# -------------------------
class MultiDisoPred(nn.Module):
    """
    Experts:
      uni1, uni2, uni3, syn, red
    Training uses (1 + n) forward passes per expert by masking one modality embedding with random vectors.
    Supports dynamic number of modalities via feature_dims and feature_indices.
    """
    def __init__(
        self,
        feature_dims,
        d_model: int = 256,
        dropout: float = 0.1,
        proj_dim: int = 128,
        # lambda_int: float = 0.1,
        triplet_margin: float = 0.2,
        rw_hidden: int = 256,
        rw_depth: int = 2,
        rw_tau: float = 1.0,
        rand_scale: float = 1.0,
        active_expert_names: Optional[List[str]] = None,
        feature_indices: Optional[List[int]] = None,
        no_synergy_loss: bool = False,
        no_redundancy_loss: bool = False,
    ):
        super().__init__()
        self.feature_dims = feature_dims if isinstance(feature_dims, list) else list(feature_dims)
        self.n_mod = len(self.feature_dims)
        self.d_model = d_model
        # self.lambda_int = lambda_int
        self.rand_scale = rand_scale
        self.no_synergy_loss = no_synergy_loss
        self.no_redundancy_loss = no_redundancy_loss

        if feature_indices is None:
            feature_indices = list(range(self.n_mod))
        self.feature_indices = feature_indices

        # modality encoders (dynamic)
        self.encoders = nn.ModuleList([
            BiLSTMEncoder(fd, d_model, dropout) for fd in self.feature_dims
        ])
        # backward compatibility
        for i in range(3):
            if i < self.n_mod:
                setattr(self, f'enc{i+1}', self.encoders[i])
            else:
                setattr(self, f'enc{i+1}', None)

        # 5 interaction experts
        self.expert_uni1 = InteractionExpert(d_model, dropout, proj_dim)
        self.expert_uni2 = InteractionExpert(d_model, dropout, proj_dim)
        self.expert_uni3 = InteractionExpert(d_model, dropout, proj_dim)
        self.expert_syn  = InteractionExpert(d_model, dropout, proj_dim)
        self.expert_red  = InteractionExpert(d_model, dropout, proj_dim)

        self.all_experts = [
            self.expert_uni1,
            self.expert_uni2,
            self.expert_uni3,
            self.expert_syn,
            self.expert_red,
        ]
        self.all_expert_names = ["uni1", "uni2", "uni3", "syn", "red"]

        # active experts for ablation study
        if active_expert_names is None:
            active_expert_names = self.all_expert_names[:]
        self.active_expert_names = active_expert_names
        self.expert_names = active_expert_names
        self.experts = [self.all_experts[self.all_expert_names.index(n)] for n in active_expert_names]

        # register buffer for expert mask (used during forward)
        mask = torch.tensor([1.0 if n in active_expert_names else 0.0 for n in self.all_expert_names])
        self.register_buffer("expert_mask", mask)

        # reweighting W: input uses concatenated gs
        self.reweighter = ReWeightMLP(
            in_dim=self.n_mod * d_model,
            num_experts=len(self.all_experts),
            hidden=rw_hidden,
            depth=rw_depth,
            dropout=dropout,
            tau=rw_tau,
        )

        self.triplet = nn.TripletMarginLoss(margin=triplet_margin, p=2)

    def _mask_modality(
        self,
        *hs: torch.Tensor,
        which: int
    ) -> Tuple[torch.Tensor, ...]:
        """
        which: 0-based index -> replace that modality with random vectors
        """
        mhs = list(hs)
        if 0 <= which < len(mhs):
            mhs[which] = rand_like(mhs[which], self.rand_scale)
        return tuple(mhs)

    def forward(
        self,
        x1: torch.Tensor, x2: torch.Tensor, x3: torch.Tensor,
        seq_mask: Optional[torch.Tensor] = None,
        compute_interaction_loss: bool = True,
        return_features: bool = False,
    ):
        """
        Returns:
          final_logits: (B,L)
          aux: dict with weights, per-expert logits, interaction loss, etc.
          features: dict with intermediate features (if return_features=True)
        """
        xs = [x1, x2, x3]

        # encode selected modalities
        hs = []
        gs = []
        for feat_idx, enc in zip(self.feature_indices, self.encoders):
            h, g = enc(xs[feat_idx], seq_mask)
            hs.append(h)
            gs.append(g)

        # re-weighting weights (sample-level)
        rw_in = torch.cat(gs, dim=-1)                # (B, n_mod*D)
        w = self.reweighter(rw_in)                   # (B, 5)

        # apply expert mask for ablation (zero out disabled experts, then renormalize)
        w = w * self.expert_mask.unsqueeze(0)        # (B, 5)
        w = w / (w.sum(dim=-1, keepdim=True) + 1e-8)  # renormalize

        # forward full + masked for each expert (always compute all 5 for interaction loss consistency)
        out_full = {}
        out_mask = {m_idx: {} for m_idx in range(self.n_mod)}

        seq_ctx = sum(gs) / len(gs)
        for name, expert in zip(self.all_expert_names, self.all_experts):
            # full
            o0 = expert(*hs, seq_ctx=seq_ctx, seq_mask=seq_mask)
            out_full[name] = o0

            # masked
            for m_idx in range(self.n_mod):
                mhs = self._mask_modality(*hs, which=m_idx)
                om = expert(*mhs, seq_ctx=seq_ctx, seq_mask=seq_mask)
                out_mask[m_idx][name] = om

        # weighted sum of full expert logits (use all 5 experts so dims match w)
        expert_logits = torch.stack([out_full[n]["logits"] for n in self.all_expert_names], dim=-1)  # (B,L,5)
        final_logits = (expert_logits * w.unsqueeze(1)).sum(dim=-1)                              # (B,L)

        # Prepare features dict if requested
        features = None
        if return_features:
            features = {
                # Encoded features (residue-level)
                'h1': hs[0].detach() if self.n_mod > 0 else None,
                'h2': hs[1].detach() if self.n_mod > 1 else None,
                'h3': hs[2].detach() if self.n_mod > 2 else None,
                # Global features (sequence-level)
                'g1': gs[0].detach() if self.n_mod > 0 else None,
                'g2': gs[1].detach() if self.n_mod > 1 else None,
                'g3': gs[2].detach() if self.n_mod > 2 else None,
                # Expert embeddings
                'expert_emb_uni1': out_full['uni1']['emb'].detach(),
                'expert_emb_uni2': out_full['uni2']['emb'].detach(),
                'expert_emb_uni3': out_full['uni3']['emb'].detach(),
                'expert_emb_syn': out_full['syn']['emb'].detach(),
                'expert_emb_red': out_full['red']['emb'].detach(),
                # Attention weights (from fusion)
                'attn_w_uni1': out_full['uni1']['attn_w'].detach(),
                'attn_w_uni2': out_full['uni2']['attn_w'].detach(),
                'attn_w_uni3': out_full['uni3']['attn_w'].detach(),
                'attn_w_syn': out_full['syn']['attn_w'].detach(),
                'attn_w_red': out_full['red']['attn_w'].detach(),
                'expert_logits': expert_logits.detach(),
            }

        # interaction loss (Appendix C style)
        lint = final_logits.new_zeros(())
        lint_dict = {}

        if compute_interaction_loss:
            # normalize embeddings
            def norm(x): return F.normalize(x, dim=-1)

            # uniqueness experts (uni1..uni3): triplet(anchor, positive(mask others), negative(mask self))
            # 特征名固定绑定: uni{i+1} 恒对应原始特征索引 i (PSSM=0 / ESM2=1 / DR-BERT=2)
            for feat_idx, uni_name in [(0, "uni1"), (1, "uni2"), (2, "uni3")]:
                if feat_idx not in self.feature_indices:
                    continue
                pos = self.feature_indices.index(feat_idx)   # 该特征在当前 hs 列表中的位置
                anchor = norm(out_full[uni_name]["emb"])
                negative = norm(out_mask[pos][uni_name]["emb"])
                trip_losses = []
                for j in range(self.n_mod):
                    if j == pos:
                        continue
                    positive = norm(out_mask[j][uni_name]["emb"])
                    trip_losses.append(self.triplet(anchor, positive, negative))
                if trip_losses:
                    li = sum(trip_losses) / len(trip_losses)
                else:
                    li = final_logits.new_zeros(())
                lint_dict[f"Lint_{uni_name}"] = li

            if self.n_mod > 1:
                if not self.no_synergy_loss:
                    # synergy expert: mean CosSim(anchor, masked)  (minimize similarity)
                    syn_anchor = norm(out_full["syn"]["emb"])
                    syn_cs = []
                    for j in range(self.n_mod):
                        syn_cs.append(F.cosine_similarity(syn_anchor, norm(out_mask[j]["syn"]["emb"]), dim=-1).mean())
                    lint_dict["Lint_syn"] = sum(syn_cs) / self.n_mod
                else:
                    lint_dict["Lint_syn"] = final_logits.new_zeros(())

                if not self.no_redundancy_loss:
                    # redundancy expert: mean (1 - CosSim(anchor, masked)) (minimize -> enforce similarity)
                    red_anchor = norm(out_full["red"]["emb"])
                    red_cs = []
                    for j in range(self.n_mod):
                        red_cs.append(1.0 - F.cosine_similarity(red_anchor, norm(out_mask[j]["red"]["emb"]), dim=-1).mean())
                    lint_dict["Lint_red"] = sum(red_cs) / self.n_mod
                else:
                    lint_dict["Lint_red"] = final_logits.new_zeros(())
            else:
                lint_dict["Lint_syn"] = final_logits.new_zeros(())
                lint_dict["Lint_red"] = final_logits.new_zeros(())

            # average over experts
            lint = sum(lint_dict.values()) / len(lint_dict)

        aux = {
            "weights": w,                       # (B,E) local interpretation
            "expert_logits": expert_logits,      # (B,L,E)
            "interaction_loss": lint,
            "interaction_loss_terms": lint_dict,
        }
        
        if return_features:
            return final_logits, aux, features
        return final_logits, aux


