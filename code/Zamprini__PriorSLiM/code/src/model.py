from __future__ import annotations
import random
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F
WINDOW_LENGTH = 29
LATENT_DIM = 40960
MOTIF_INDICES = (54, 55, 56)
CHEMISTRY_INDICES = (*range(54), 57)
CHEMISTRY_GROUPS = (tuple(range(0, 48, 6)), tuple(range(2, 48, 6)), tuple(range(4, 48, 6)), (48, 49, 50, 51, 52, 53, 54), tuple(range(1, 48, 2)))

def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

def _center_mask(length: int, width: int, device: torch.device) -> torch.Tensor:
    center = length // 2
    radius = width // 2
    positions = torch.arange(length, device=device)
    return (positions >= center - radius) & (positions <= center + radius)

def _masked_mean(values: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    denominator = mask.sum(dim=1).clamp_min(1).to(values.dtype)
    return values.masked_fill(~mask, 0.0).sum(dim=1) / denominator

def _masked_max(values: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    maxima = torch.where(mask, values, torch.full_like(values, -1000000000.0)).max(dim=1).values
    return torch.where(maxima < -100000000.0, torch.zeros_like(maxima), maxima)

def _masked_topk_mean(values: torch.Tensor, mask: torch.Tensor, k: int) -> torch.Tensor:
    masked = torch.where(mask, values, torch.full_like(values, -1000000000.0))
    top = torch.topk(masked, k=min(k, values.shape[1]), dim=1).values
    valid = top > -100000000.0
    top = torch.where(valid, top, torch.zeros_like(top))
    return top.sum(dim=1) / valid.sum(dim=1).clamp_min(1).to(values.dtype)

def _masked_peak_fraction(values: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    peaks = (values > 0.0) & mask
    denominator = mask.sum(dim=1).clamp_min(1).to(values.dtype)
    return peaks.sum(dim=1).to(values.dtype) / denominator

def compute_pssm_center_summaries(x_pssm: torch.Tensor, valid_mask: torch.Tensor) -> torch.Tensor:
    length = x_pssm.shape[1]
    center7 = _center_mask(length, 7, x_pssm.device).unsqueeze(0) & valid_mask
    center15 = _center_mask(length, 15, x_pssm.device).unsqueeze(0) & valid_mask
    flank = (~_center_mask(length, 15, x_pssm.device)).unsqueeze(0) & valid_mask
    summaries = []
    for values in (x_pssm[..., 4], x_pssm[..., 5]):
        flank_mean = _masked_mean(values, flank)
        center7_mean = _masked_mean(values, center7)
        center15_mean = _masked_mean(values, center15)
        summaries.extend((center7_mean, _masked_max(values, center7), _masked_topk_mean(values, center7, 3), center15_mean, _masked_max(values, center15), _masked_topk_mean(values, center15, 5), flank_mean, center7_mean - flank_mean, center15_mean - flank_mean, _masked_peak_fraction(values, center7), _masked_peak_fraction(values, center15)))
    return torch.stack(summaries, dim=1)

class ResidualConvBlock(nn.Module):

    def __init__(self, channels: int, kernel_size: int, dropout: float=0.15):
        super().__init__()
        self.conv = nn.Conv1d(channels, channels, kernel_size=kernel_size, padding=kernel_size // 2)
        self.norm = nn.LayerNorm(channels)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = x
        x = self.conv(x).permute(0, 2, 1)
        x = self.norm(x).permute(0, 2, 1)
        return residual + self.dropout(F.gelu(x))

class MaskedAttentionPooling(nn.Module):

    def __init__(self, dimension: int, hidden_dimension: int=32):
        super().__init__()
        self.score = nn.Sequential(nn.LayerNorm(dimension), nn.Linear(dimension, hidden_dimension), nn.GELU(), nn.Dropout(0.1), nn.Linear(hidden_dimension, 1))

    def forward(self, x: torch.Tensor, valid_mask: torch.Tensor) -> torch.Tensor:
        scores = self.score(x).squeeze(-1).masked_fill(~valid_mask, -1000000000.0)
        weights = torch.softmax(scores, dim=1)
        weights = torch.where(valid_mask, weights, torch.zeros_like(weights))
        weights = weights / weights.sum(dim=1, keepdim=True).clamp_min(1e-08)
        return torch.sum(x * weights.unsqueeze(-1), dim=1)

class PSSMSummaryBranch(nn.Module):

    def __init__(self):
        super().__init__()
        self.cons_channels = 7
        self.seq_len = WINDOW_LENGTH
        self.valid_mask_index = 6
        self.input_proj = nn.Linear(7, 32)
        self.pos_embedding = nn.Parameter(torch.empty(1, WINDOW_LENGTH, 32))
        nn.init.normal_(self.pos_embedding, mean=0.0, std=0.02)
        self.conv_blocks = nn.Sequential(ResidualConvBlock(32, kernel_size=5, dropout=0.15), ResidualConvBlock(32, kernel_size=9, dropout=0.15))
        self.attn_pool = MaskedAttentionPooling(32, hidden_dimension=32)
        self.pooled_out = nn.Sequential(nn.Linear(32, 24), nn.LayerNorm(24), nn.GELU(), nn.Dropout(0.2))
        self.summary_out = nn.Sequential(nn.LayerNorm(22), nn.Linear(22, 32), nn.GELU(), nn.Dropout(0.15), nn.Linear(32, 16), nn.LayerNorm(16), nn.GELU())

    def forward(self, x_pssm: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        valid_mask = x_pssm[..., self.valid_mask_index] > -1.0
        valid_mask = torch.where(valid_mask.any(dim=1, keepdim=True), valid_mask, torch.ones_like(valid_mask))
        summaries = compute_pssm_center_summaries(x_pssm, valid_mask)
        summary_input = summaries
        if summaries.device.type == 'mps' and torch.is_grad_enabled() and (not summaries.requires_grad):
            summary_input = summaries.detach().requires_grad_(True)
        encoded = self.input_proj(x_pssm) + self.pos_embedding[:, :x_pssm.shape[1], :]
        encoded = self.conv_blocks(encoded.permute(0, 2, 1)).permute(0, 2, 1)
        pooled = self.attn_pool(encoded, valid_mask)
        embedding = torch.cat((self.pooled_out(pooled), self.summary_out(summary_input)), dim=1)
        return (embedding, summaries)

class ConservationEncoder(nn.Module):

    def __init__(self):
        super().__init__()
        self.pssm_branch = PSSMSummaryBranch()
        self.fusion = nn.Sequential(nn.Linear(41, 32), nn.LayerNorm(32), nn.GELU(), nn.Dropout(0.15))

    def forward(self, pssm: torch.Tensor, pssm_window: torch.Tensor) -> torch.Tensor:
        embedding, _summaries = self.pssm_branch(pssm)
        return self.fusion(torch.cat((embedding, pssm_window), dim=1))

class MotifEncoder(nn.Module):

    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(3, 16), nn.GELU(), nn.Dropout(0.1), nn.Linear(16, 32), nn.LayerNorm(32), nn.GELU())

    def forward(self, motif: torch.Tensor) -> torch.Tensor:
        return self.net(motif)

class GatedGroupEncoder(nn.Module):

    def __init__(self, input_dimension: int):
        super().__init__()
        self.value = nn.Linear(input_dimension, 16)
        self.gate = nn.Linear(input_dimension, 16)
        self.dropout = nn.Dropout(0.15)
        self.out = nn.Sequential(nn.Linear(16, 8), nn.LayerNorm(8), nn.GELU())

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        gated = F.silu(self.value(x)) * torch.sigmoid(self.gate(x))
        return self.out(self.dropout(gated))

class ChemistryEncoder(nn.Module):

    def __init__(self):
        super().__init__()
        flattened = sorted((index for group in CHEMISTRY_GROUPS for index in group))
        if flattened != list(range(55)):
            raise RuntimeError('Chemical groups do not cover all 55 non-motif features')
        self.group_indices = [torch.tensor(group, dtype=torch.long) for group in CHEMISTRY_GROUPS]
        self.encoders = nn.ModuleList((GatedGroupEncoder(len(group)) for group in CHEMISTRY_GROUPS))
        self.fusion = nn.Sequential(nn.Linear(40, 48), nn.LayerNorm(48), nn.GELU(), nn.Dropout(0.2), nn.Linear(48, 32), nn.LayerNorm(32), nn.GELU(), nn.Dropout(0.2))

    def forward(self, chemistry: torch.Tensor) -> torch.Tensor:
        embeddings = [encoder(chemistry.index_select(1, indices.to(chemistry.device))) for indices, encoder in zip(self.group_indices, self.encoders)]
        return self.fusion(torch.cat(embeddings, dim=1))

def _interaction_fusion() -> nn.Sequential:
    return nn.Sequential(nn.Linear(96, 32), nn.LayerNorm(32), nn.GELU(), nn.Dropout(0.15))

class AccessibilityResidual(nn.Module):

    def __init__(self) -> None:
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(7, 32), nn.LayerNorm(32), nn.GELU())
        self.head = nn.Sequential(nn.Linear(32, 32), nn.GELU(), nn.Dropout(0.1), nn.Linear(32, 1))
        self.b_A = nn.Parameter(torch.zeros(()))
        nn.init.zeros_(self.head[-1].weight)
        nn.init.zeros_(self.head[-1].bias)

    def forward(self, accessibility: torch.Tensor) -> torch.Tensor:
        if accessibility.ndim != 2 or accessibility.shape[1] != 7:
            raise ValueError('Accessibility must have shape [B,7]')
        return self.head(self.encoder(accessibility)).squeeze(-1) + self.b_A

class IntrinsicEncoder(nn.Module):

    def __init__(self, scalers, motif_dropout=0.4):
        super().__init__()
        flags = dict(conservation=True, motif=True, chemistry=True)
        self.flags = flags
        self.motif_dropout = motif_dropout
        self.motif_indices = MOTIF_INDICES
        self.chemistry_indices = CHEMISTRY_INDICES
        for key, values in scalers.items():
            for field in ('mean', 'scale'):
                self.register_buffer(f'{key}_{field}', values[field].clone())
        self.conservation = ConservationEncoder() if flags['conservation'] else None
        self.motif = MotifEncoder() if flags['motif'] else None
        self.chemistry = ChemistryEncoder() if flags['chemistry'] else None
        self.cm = _interaction_fusion() if self.conservation is not None and self.motif is not None else None
        self.qm = _interaction_fusion() if self.chemistry is not None and self.motif is not None else None
        count = sum((m is not None for m in (self.conservation, self.motif, self.chemistry, self.cm, self.qm)))
        if not count:
            raise ValueError('A uses an intercept, not IntrinsicAblation')
        self.fusion = nn.Sequential(nn.Linear(count * 32, 96), nn.LayerNorm(96), nn.GELU(), nn.Dropout(0.2), nn.Linear(96, 64), nn.LayerNorm(64), nn.GELU(), nn.Dropout(0.15))
        self.head = nn.Linear(64, 1)

    def forward(self, residue, mask, peak, chemical):
        pieces = []
        cons = motif = chem = None
        if self.conservation is not None:
            p = torch.cat((residue, mask.unsqueeze(-1)), dim=-1).reshape(len(residue), -1)
            p = ((p - self.pssm_mean) / self.pssm_scale).reshape(-1, 29, 7)
            cons = self.conservation(p, (peak - self.pssm_window_mean) / self.pssm_window_scale)
            pieces.append(cons)
        if self.motif is not None:
            idx = list(self.motif_indices)
            m = (chemical[:, idx] - self.chemical_mean[idx]) / self.chemical_scale[idx]
            if self.training and self.motif_dropout:
                m = m * (torch.rand((len(m), 1), device=m.device) >= self.motif_dropout)
            motif = self.motif(m)
            pieces.append(motif)
        if self.chemistry is not None:
            idx = list(self.chemistry_indices)
            chem = self.chemistry((chemical[:, idx] - self.chemical_mean[idx]) / self.chemical_scale[idx])
            pieces.append(chem)
        if self.cm is not None:
            pieces.append(self.cm(torch.cat((cons, motif, cons * motif), dim=1)))
        if self.qm is not None:
            pieces.append(self.qm(torch.cat((chem, motif, chem * motif), dim=1)))
        return self.head(self.fusion(torch.cat(pieces, dim=1))).squeeze(-1)

class PositionalReadout(nn.Module):

    def __init__(self, feature_count, weights=None):
        super().__init__()
        self.weights = nn.Parameter(torch.zeros(feature_count, 29) if weights is None else torch.as_tensor(weights, dtype=torch.float32).clone())

    def forward(self, events):
        col, val, owner, n = events
        contributions = self.weights.reshape(-1)[col] * val
        result = val.new_zeros(n)
        return result.index_add(0, owner, contributions)

class JointModel(nn.Module):

    def __init__(self, scalers, feature_count, seed, motif_dropout=0.4):
        super().__init__()
        seed_all(seed)
        self.intrinsic = IntrinsicEncoder(scalers, motif_dropout)
        seed_all(seed + 100)
        self.access = AccessibilityResidual()
        seed_all(seed + 200)
        self.profile = PositionalReadout(feature_count)
        self.register_buffer('access_mean', scalers['accessibility']['mean'].clone())
        self.register_buffer('access_scale', scalers['accessibility']['scale'].clone())

    def forward(self, modalities, events, return_parts=False):
        intrinsic = self.intrinsic(*modalities[:4])
        access = self.access((modalities[4] - self.access_mean) / self.access_scale)
        explicit = intrinsic + access
        motifae = self.profile(events)
        final = explicit + motifae
        return (final, explicit, motifae) if return_parts else final

    def load_explicit(self, base_state, access_state):
        self.intrinsic.load_state_dict(base_state, strict=True)
        self.access.load_state_dict(access_state, strict=True)

    def parameter_groups(self, training):
        return [dict(params=self.intrinsic.parameters(), lr=training['base_lr'], weight_decay=training['base_wd']), dict(params=self.access.parameters(), lr=training['access_lr'], weight_decay=training['base_wd']), dict(params=self.profile.parameters(), lr=training['profile_lr'], weight_decay=training['profile_wd'])]

class SelectedRefinement(torch.nn.Module):

    def __init__(self, selected, W0, candidate_count=5845):
        super().__init__()
        lookup = torch.full((candidate_count,), -1, dtype=torch.long)
        lookup[selected] = torch.arange(len(selected))
        self.register_buffer('lookup', lookup)
        self.register_buffer('W0', W0.detach().clone().float())
        self.DeltaW = torch.nn.Parameter(torch.zeros_like(self.W0))

    def forward(self, events):
        col, value, owner, n = events
        j = self.lookup[col // 29]
        keep = j >= 0
        contribution = (self.W0 + self.DeltaW)[j[keep], col[keep] % 29] * value[keep]
        return torch.zeros(n, device=value.device, dtype=value.dtype).index_add(0, owner[keep], contribution)

    def drift_loss(self):
        return self.DeltaW.square().sum() / (self.W0.square().sum() + 1e-12)

class PriorSLiMPredictor(nn.Module):

    def __init__(self, bundle):
        super().__init__()
        if bundle.get('format') != 'jwarm_bundle_v1':
            raise ValueError('Expected a standalone J_WARM bundle')
        self.feature_manifest = bundle['feature_manifest']
        self.threshold = float(bundle['threshold'])
        self.network = JointModel(bundle['preprocessing'], len(bundle['feature_ids']), bundle['seed'], bundle['motif_dropout'])
        state = dict(bundle['explicit_state'], **{'profile.weights': bundle['W_final']})
        self.network.load_state_dict(state, strict=True)
        self.register_buffer('feature_ids', bundle['feature_ids'].long().clone())
        self.register_buffer('scales', bundle['train_p95_scales'].float().clone())
        lookup = torch.full((LATENT_DIM,), -1, dtype=torch.long)
        lookup[self.feature_ids] = torch.arange(len(self.feature_ids))
        self.register_buffer('lookup', lookup)
        self.requires_grad_(False).eval()

    def events(self, sparse):
        ids = sparse['feature_indices']
        residues = sparse['residue_indices']
        values = sparse['feature_values']
        j = self.lookup[ids]
        keep = j >= 0
        cols = j[keep] * 29 + residues[keep] % 29
        return (cols, values[keep] / (self.scales[j[keep]] + 1e-08), residues[keep] // 29, len(sparse['valid_mask']))

    def forward(self, modalities, sparse, return_parts=False):
        return self.network(modalities, self.events(sparse), return_parts=return_parts)
