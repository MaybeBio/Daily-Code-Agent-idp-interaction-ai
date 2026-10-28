from __future__ import annotations
import math
import torch
from torch import nn
from torch.nn import functional as F
WINDOWS = (33, 129, 257)

def pool_windows(embeddings: torch.Tensor, valid_length: int | None=None) -> torch.Tensor:
    """Original PSTP boundary-corrected means, without padded residues."""
    if embeddings.ndim != 2 or embeddings.shape[1] != 650:
        raise ValueError('Expected [sequence_length, 650] embeddings')
    length = embeddings.shape[0] if valid_length is None else valid_length
    if not 0 < length <= embeddings.shape[0]:
        raise ValueError('Invalid effective sequence length')
    x = embeddings[:length].T.unsqueeze(0)
    positions = torch.arange(length, device=x.device)
    result = []
    for window in WINDOWS:
        radius = window // 2
        counts = (positions + radius + 1).clamp(max=length) - (positions - radius).clamp(min=0)
        mean = F.avg_pool1d(x, window, stride=1, padding=radius)[0].T
        result.append(mean * (window / counts.to(x.dtype)).unsqueeze(1))
    return torch.stack(result)

class BaselineHead(nn.Module):
    """Shared head across the three original scales; encoders stay fixed."""

    def __init__(self) -> None:
        super().__init__()
        self.layer1 = nn.Linear(650, 20)
        self.layer2 = nn.Linear(20, 5)
        self.layer3 = nn.Linear(5, 1)
        for layer in (self.layer1, self.layer2, self.layer3):
            nn.init.kaiming_normal_(layer.weight, nonlinearity='leaky_relu')
            nn.init.zeros_(layer.bias)

    def position_logits(self, pooled: torch.Tensor) -> torch.Tensor:
        x = F.leaky_relu(self.layer1(pooled))
        x = self.layer2(x)
        return self.layer3(x).squeeze(-1)

    def aggregate(self, logits: torch.Tensor) -> torch.Tensor:
        return torch.sigmoid(logits).max(dim=-1).values

    def forward(self, pooled: torch.Tensor) -> torch.Tensor:
        if pooled.ndim != 3 or pooled.shape[0] != 3 or pooled.shape[2] != 650 or (pooled.shape[1] == 0):
            raise ValueError('Expected three scales: [3, sequence_length, 650]')
        return self.aggregate(self.position_logits(pooled)).mean()

    def forward_batch(self, pooled: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        if pooled.ndim != 4 or pooled.shape[1] != 3 or pooled.shape[3] != 650:
            raise ValueError('Expected batched scales [batch, 3, length, 650]')
        if mask.shape != (pooled.shape[0], pooled.shape[2]) or not mask.any(dim=1).all():
            raise ValueError('Invalid residue mask')
        logits = self.position_logits(pooled).masked_fill(~mask[:, None, :], -torch.inf)
        return torch.sigmoid(logits).max(dim=-1).values.mean(dim=1)

def paper_regularizer(model: BaselineHead) -> torch.Tensor:
    """Original gfs=True: 0.003 group norms plus L1 of remaining parameters."""
    weight = model.layer1.weight
    groups = weight[:, :320].square().sum().sqrt() + weight[:, 320:].square().sum().sqrt()
    return 0.003 * groups + model.layer1.bias.abs().sum() + sum((parameter.abs().sum() for layer in (model.layer2, model.layer3) for parameter in layer.parameters()))

def _pool_channels(values: torch.Tensor) -> torch.Tensor:
    if values.ndim != 2 or values.shape[0] == 0:
        raise ValueError('Expected non-empty [length, channels] tensor')
    length = values.shape[0]
    x = values.T.unsqueeze(0)
    positions = torch.arange(length, device=values.device)
    result = []
    for window in WINDOWS:
        radius = window // 2
        counts = (positions + radius + 1).clamp(max=length) - (positions - radius).clamp(min=0)
        mean = F.avg_pool1d(x, window, stride=1, padding=radius)[0].T
        result.append(mean * (window / counts.to(values.dtype)).unsqueeze(1))
    return torch.stack(result)

def _pool_channels_batch(values: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if values.ndim != 3 or mask.shape != values.shape[:2] or (not mask.any(dim=1).all()):
        raise ValueError('Expected batched values and a non-empty residue mask')
    batch, padded_length, _ = values.shape
    lengths = mask.sum(dim=1)
    x = (values * mask[..., None]).transpose(1, 2)
    positions = torch.arange(padded_length, device=values.device)[None, :]
    result = []
    for window in WINDOWS:
        radius = window // 2
        left = (positions - radius).clamp(min=0)
        right = torch.minimum(positions + radius + 1, lengths[:, None])
        counts = (right - left).clamp(min=1).to(values.dtype)
        sums = F.avg_pool1d(x, window, stride=1, padding=radius) * window
        means = sums.transpose(1, 2) / counts[..., None]
        result.append(means * mask[..., None])
    return torch.stack(result, dim=1)

class DepthwiseSeparableResidual(nn.Module):

    def __init__(self, channels: int, kernel_size: int, dilation: int) -> None:
        super().__init__()
        padding = dilation * (kernel_size // 2)
        self.depthwise = nn.Conv1d(channels, channels, kernel_size, padding=padding, dilation=dilation, groups=channels)
        self.pointwise = nn.Conv1d(channels, channels, 1)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        if values.ndim not in (2, 3):
            raise ValueError('Expected [length, channels] or [batch, length, channels]')
        single = values.ndim == 2
        encoded = (values.unsqueeze(0) if single else values).transpose(1, 2)
        update = self.pointwise(self.depthwise(encoded))
        result = (values.unsqueeze(0) if single else values) + F.leaky_relu(update.transpose(1, 2))
        return result[0] if single else result

class PhaseLACE(nn.Module):
    """Final local-global residual member; the pretrained baseline is frozen."""

    def __init__(self, baseline: BaselineHead, residual_scale: float=0.3) -> None:
        super().__init__()
        if not math.isfinite(residual_scale) or residual_scale <= 0:
            raise ValueError('residual_scale must be finite and positive')
        self.baseline = baseline
        self.residual_scale = residual_scale
        self.esm_projection = nn.Linear(320, 64)
        self.local_blocks = nn.ModuleList([DepthwiseSeparableResidual(64, 9, 1), DepthwiseSeparableResidual(64, 9, 4)])
        self.global_projection = nn.Linear(320, 32)
        self.delta_head = nn.Sequential(nn.Linear(96, 20), nn.LeakyReLU(), nn.Linear(20, 1))
        nn.init.zeros_(self.delta_head[-1].weight)
        nn.init.zeros_(self.delta_head[-1].bias)
        for parameter in self.baseline.parameters():
            parameter.requires_grad_(False)
        self.baseline.eval()

    def train(self, mode: bool=True) -> 'PhaseLACE':
        super().train(mode)
        self.baseline.eval()
        return self

    def forward_components(self, raw: torch.Tensor) -> dict[str, torch.Tensor]:
        if raw.ndim != 2 or raw.shape[1] != 650 or raw.shape[0] == 0:
            raise ValueError('Expected raw residue embeddings [length, 650]')
        pooled_raw = pool_windows(raw)
        with torch.no_grad():
            base_position = self.baseline.position_logits(pooled_raw)
        local = F.leaky_relu(self.esm_projection(raw[:, :320]))
        for block in self.local_blocks:
            local = block(local)
        local_scales = _pool_channels(local)
        global_input = raw[:, :320]
        context = F.leaky_relu(self.global_projection(global_input.mean(dim=0)))
        broadcast = context.view(1, 1, -1).expand(3, raw.shape[0], -1)
        delta_position = self.delta_head(torch.cat((local_scales, broadcast), dim=-1)).squeeze(-1)
        final_position = base_position + self.residual_scale * delta_position
        scale_logits = final_position.max(dim=1).values
        probability = torch.sigmoid(scale_logits).mean()
        return {'probability': probability, 'base_position_logits': base_position, 'delta_position_logits': delta_position, 'scale_logits': scale_logits}

    def forward(self, raw: torch.Tensor) -> torch.Tensor:
        return self.forward_components(raw)['probability']

    def forward_batch(self, raw: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        if raw.ndim != 3 or raw.shape[2] != 650 or mask.shape != raw.shape[:2] or (not mask.any(dim=1).all()):
            raise ValueError('Expected batched raw embeddings and a non-empty mask')
        pooled_raw = _pool_channels_batch(raw, mask)
        with torch.no_grad():
            base_position = self.baseline.position_logits(pooled_raw)
        local = F.leaky_relu(self.esm_projection(raw[:, :, :320])) * mask[..., None]
        for block in self.local_blocks:
            local = block(local) * mask[..., None]
        local_scales = _pool_channels_batch(local, mask)
        global_input = raw[:, :, :320]
        context_input = (global_input * mask[..., None]).sum(dim=1) / mask.sum(dim=1, keepdim=True)
        context = F.leaky_relu(self.global_projection(context_input))
        broadcast = context[:, None, None, :].expand(-1, 3, raw.shape[1], -1)
        delta_position = self.delta_head(torch.cat((local_scales, broadcast), dim=-1)).squeeze(-1)
        final_position = (base_position + self.residual_scale * delta_position).masked_fill(~mask[:, None, :], -torch.inf)
        return torch.sigmoid(final_position.max(dim=2).values).mean(dim=1)

def residual_regularizer(model: PhaseLACE) -> torch.Tensor:
    return sum((parameter.abs().sum() for parameter in model.parameters() if parameter.requires_grad))
