# -*- coding: utf-8 -*-
from typing import Tuple

import torch
import torch.nn as nn

from .blocks import UNetMidBlock3D, UpDecoderBlock3D
from diffusers.utils import is_torch_version


class Decoder(nn.Module):
    def __init__(
        self,
        in_channels: int = 4,
        out_channels: int = 1,
        num_up_blocks: int = 4,
        block_out_channels: Tuple[int, ...] = (32, 64, 128, 256),
        layers_per_block: int = 2,
        norm_num_groups: int = 32,
        act_fn: str = "silu",
    ):
        super().__init__()
        self.layers_per_block = layers_per_block

        self.conv_in = nn.Conv3d(in_channels, block_out_channels[-1], kernel_size=3, stride=1, padding=1)

        self.mid_block = UNetMidBlock3D(
            in_channels=block_out_channels[-1],
            resnet_eps=1e-6,
            resnet_act_fn=act_fn,
            output_scale_factor=1,
            resnet_groups=norm_num_groups,
        )

        self.up_blocks = nn.ModuleList([])
        reversed_block_out_channels = list(reversed(block_out_channels))
        output_channel = reversed_block_out_channels[0]
        for i in range(num_up_blocks):
            prev_output_channel = output_channel
            output_channel = reversed_block_out_channels[i]
            is_final_block = i == num_up_blocks - 1
            self.up_blocks.append(
                UpDecoderBlock3D(
                    in_channels=prev_output_channel,
                    out_channels=output_channel,
                    num_layers=layers_per_block + 1,
                    resnet_eps=1e-6,
                    resnet_act_fn=act_fn,
                    resnet_groups=norm_num_groups,
                    add_upsample=not is_final_block,
                )
            )

        self.conv_norm_out = nn.GroupNorm(num_channels=block_out_channels[0], num_groups=norm_num_groups, eps=1e-6)
        self.conv_act = nn.SiLU()
        self.conv_out = nn.Conv3d(block_out_channels[0], out_channels, kernel_size=3, padding=1)

        self.gradient_checkpointing = False

    def forward(self, sample: torch.Tensor) -> torch.Tensor:
        sample = self.conv_in(sample)
        upscale_dtype = next(iter(self.up_blocks.parameters())).dtype

        if self.training and self.gradient_checkpointing:
            def create_custom_forward(module):
                def custom_forward(*inputs):
                    return module(*inputs)
                return custom_forward

            if is_torch_version(">=", "1.11.0"):
                sample = torch.utils.checkpoint.checkpoint(
                    create_custom_forward(self.mid_block), sample, use_reentrant=False
                )
                sample = sample.to(upscale_dtype)
                for up_block in self.up_blocks:
                    sample = torch.utils.checkpoint.checkpoint(
                        create_custom_forward(up_block), sample, use_reentrant=False
                    )
            else:
                sample = torch.utils.checkpoint.checkpoint(create_custom_forward(self.mid_block), sample)
                sample = sample.to(upscale_dtype)
                for up_block in self.up_blocks:
                    sample = torch.utils.checkpoint.checkpoint(create_custom_forward(up_block), sample)
        else:
            sample = self.mid_block(sample)
            sample = sample.to(upscale_dtype)
            for up_block in self.up_blocks:
                sample = up_block(sample)

        sample = self.conv_norm_out(sample)
        sample = self.conv_act(sample)
        sample = self.conv_out(sample)
        return sample
