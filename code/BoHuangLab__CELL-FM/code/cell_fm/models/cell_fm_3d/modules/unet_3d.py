import math

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.checkpoint import checkpoint

from .res_block_3d import ResBlock3D


class UNet3D(nn.Module):
    """Encoder and decoder halves of a 3-D UNet; the caller runs the bottleneck between them.

    `block_out_channels` is the width of each resolution level, so it has one more entry than there
    are 2x downsamples. Every level runs `layers_per_block` ResBlocks conditioned on `temb` on the
    way down and again on the way up. The decoder doubles the resolution with a transposed conv,
    concatenates the encoder's features from that level and projects them back with a 1x1 conv.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        block_out_channels: list,
        layers_per_block: int,
        cond_dim: int,
        checkpointing: bool = True,
    ):
        super().__init__()
        chans = list(block_out_channels)
        pairs = list(zip(chans[:-1], chans[1:]))

        def res_blocks(channels):
            return nn.ModuleList([ResBlock3D(channels, cond_dim) for _ in range(layers_per_block)])

        self.conv_in = nn.Conv3d(in_channels, chans[0], kernel_size=3, padding=1)
        self.down_blocks = nn.ModuleList([res_blocks(c) for c in chans])
        self.downsamples = nn.ModuleList(
            [nn.Conv3d(c, c_next, kernel_size=2, stride=2) for c, c_next in pairs]
        )
        self.upsamples = nn.ModuleList(
            [nn.ConvTranspose3d(c_next, c, kernel_size=2, stride=2) for c, c_next in pairs]
        )
        self.skip_projs = nn.ModuleList([nn.Conv3d(2 * c, c, kernel_size=1) for c in chans[:-1]])
        self.up_blocks = nn.ModuleList([res_blocks(c) for c in chans])
        self.norm_out = nn.GroupNorm(math.gcd(32, chans[0]), chans[0], eps=1e-6)
        self.conv_out = nn.Conv3d(chans[0], out_channels, kernel_size=3, padding=1)
        # Zero-init as in ADM, so a new model starts by predicting zeros.
        nn.init.zeros_(self.conv_out.weight)
        nn.init.zeros_(self.conv_out.bias)

        # Not named gradient_checkpointing: HF Trainer reads that attribute as model-wide
        # checkpointing and would switch DDP's find_unused_parameters off.
        self.checkpointing = checkpointing

    def encode(self, x: torch.Tensor, temb: torch.Tensor):
        """Returns the lowest-resolution features and the per-level skips for `decode`."""
        x = self.conv_in(x)
        skips = []
        for level, blocks in enumerate(self.down_blocks):
            x = self._run(blocks, x, temb)
            if level < len(self.downsamples):
                skips.append(x)
                x = self.downsamples[level](x)
        return x, skips

    def decode(self, x: torch.Tensor, skips: list, temb: torch.Tensor) -> torch.Tensor:
        for level in reversed(range(len(self.up_blocks))):
            if level < len(self.upsamples):
                x = torch.cat([self.upsamples[level](x), skips[level]], dim=1)
                x = self.skip_projs[level](x)
            x = self._run(self.up_blocks[level], x, temb)
        return self.conv_out(F.silu(self.norm_out(x)))

    def _run(self, blocks: nn.ModuleList, x: torch.Tensor, temb: torch.Tensor) -> torch.Tensor:
        for block in blocks:
            # Full-resolution 3-D activations dominate memory; recomputing a block in backward
            # costs one extra forward of it.
            if self.checkpointing and self.training and torch.is_grad_enabled():
                x = checkpoint(block, x, temb, use_reentrant=False)
            else:
                x = block(x, temb)
        return x
