import torch.nn as nn

from .patch_embed_3d import PatchEmbed3D
from .res_block_3d import ResBlock3D
from .transformer_3d import SD3Transformer3DModel
from .unet_3d import UNet3D


class CondConvNet3D(nn.Module):
    """Downsamples a 3-D conditioning volume to match the VAE latent spatial size."""

    def __init__(self, in_channels, cond_out_channels):
        super().__init__()
        layers = []
        layers.append(
            nn.Conv3d(in_channels, cond_out_channels[0], kernel_size=2, stride=2, padding=0)
        )
        for i in range(len(cond_out_channels) - 1):
            layers.append(nn.SiLU())
            layers.append(
                nn.Conv3d(cond_out_channels[i], cond_out_channels[i + 1], kernel_size=2, stride=2, padding=0)
            )
        self.model = nn.Sequential(*layers)

    def forward(self, x):
        return self.model(x)
