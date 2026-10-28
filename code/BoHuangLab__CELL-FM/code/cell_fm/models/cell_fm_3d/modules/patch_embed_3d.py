import numpy as np
import torch
import torch.nn as nn

from cell_fm.models.cell_fm.modules.positional_embedding import get_1d_sincos_pos_embed_from_grid


def get_3d_sincos_pos_embed(embed_dim: int, nd: int, nh: int, nw: int) -> np.ndarray:
    """Return (nd*nh*nw, embed_dim) 3-D sincos positional embedding.

    embed_dim is split roughly equally across D, H, W axes (each dim even).
    """
    base = (embed_dim // 6) * 2   # largest even number ≤ embed_dim / 3
    dim_d = base
    dim_h = base
    dim_w = embed_dim - 2 * base  # also even: embed_dim even, base even

    gd, gh, gw = np.meshgrid(
        np.arange(nd, dtype=np.float32),
        np.arange(nh, dtype=np.float32),
        np.arange(nw, dtype=np.float32),
        indexing='ij',
    )
    return np.concatenate([
        get_1d_sincos_pos_embed_from_grid(dim_d, gd.flatten()),
        get_1d_sincos_pos_embed_from_grid(dim_h, gh.flatten()),
        get_1d_sincos_pos_embed_from_grid(dim_w, gw.flatten()),
    ], axis=1)  # (nd*nh*nw, embed_dim)


class PatchEmbed3D(nn.Module):
    """Conv3d patch embedding with 3-D sincos positional embedding."""

    def __init__(
        self,
        depth: int,
        height: int,
        width: int,
        patch_d: int,
        patch_hw: int,
        in_channels: int,
        embed_dim: int,
    ):
        super().__init__()
        self.proj = nn.Conv3d(
            in_channels, embed_dim,
            kernel_size=(patch_d, patch_hw, patch_hw),
            stride=(patch_d, patch_hw, patch_hw),
        )
        nd, nh, nw = depth // patch_d, height // patch_hw, width // patch_hw
        pos_embed = get_3d_sincos_pos_embed(embed_dim, nd, nh, nw)
        self.register_buffer(
            "pos_embed",
            torch.from_numpy(pos_embed).float().unsqueeze(0),  # (1, N, D)
        )
        self.num_patches = nd * nh * nw

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: (B, C, D, H, W)
        x = self.proj(x)                  # (B, embed_dim, nd, nh, nw)
        x = x.flatten(2).transpose(1, 2)  # (B, N, embed_dim)
        return x + self.pos_embed
