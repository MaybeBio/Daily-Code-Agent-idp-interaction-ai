import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class ResBlock3D(nn.Module):
    """3-D residual block conditioned on SD3's embedding (timestep + pooled sequence).

    The embedding scale-shifts the second GroupNorm, as in ADM. `conv2` starts at zero, so a new
    block is the identity.
    """

    def __init__(self, channels: int, cond_dim: int):
        super().__init__()
        groups = math.gcd(32, channels)
        self.norm1 = nn.GroupNorm(groups, channels, eps=1e-6)
        self.conv1 = nn.Conv3d(channels, channels, kernel_size=3, padding=1)
        self.cond_proj = nn.Linear(cond_dim, 2 * channels)
        self.norm2 = nn.GroupNorm(groups, channels, eps=1e-6)
        self.conv2 = nn.Conv3d(channels, channels, kernel_size=3, padding=1)
        nn.init.zeros_(self.conv2.weight)
        nn.init.zeros_(self.conv2.bias)

    def forward(self, x: torch.Tensor, temb: torch.Tensor) -> torch.Tensor:
        h = self.conv1(F.silu(self.norm1(x)))
        scale, shift = self.cond_proj(F.silu(temb))[:, :, None, None, None].chunk(2, dim=1)
        h = self.norm2(h) * (1 + scale) + shift
        return x + self.conv2(F.silu(h))
