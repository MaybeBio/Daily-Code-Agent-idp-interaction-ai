# -*- coding: utf-8 -*-
from typing import Tuple

import torch
import torch.nn as nn

from .blocks import DownEncoderBlock3D, UNetMidBlock3D
from diffusers.utils import is_torch_version


class Encoder(nn.Module):
    def __init__(
        self,
        in_channels: int = 1,
        out_channels: int = 4,
        num_down_blocks: int = 4,
        block_out_channels: Tuple[int, ...] = (32, 64, 128, 256),
        layers_per_block: int = 2,
        norm_num_groups: int = 32,
        act_fn: str = "silu",
        double_z: bool = True,
    ):
        super().__init__()
        self.layers_per_block = layers_per_block

        self.conv_in = nn.Conv3d(
            in_channels, block_out_channels[0], kernel_size=3, stride=1, padding=1, padding_mode='reflect'
        )

        self.down_blocks = nn.ModuleList([])
        output_channel = block_out_channels[0]
        for i in range(num_down_blocks):
            input_channel = output_channel
            output_channel = block_out_channels[i]
            is_final_block = i == num_down_blocks - 1
            self.down_blocks.append(
                DownEncoderBlock3D(
                    in_channels=input_channel,
                    out_channels=output_channel,
                    num_layers=layers_per_block,
                    add_downsample=not is_final_block,
                    resnet_eps=1e-6,
                    resnet_act_fn=act_fn,
                    resnet_groups=norm_num_groups,
                    downsample_padding=0,
                )
            )

        self.mid_block = UNetMidBlock3D(
            in_channels=block_out_channels[-1],
            resnet_eps=1e-6,
            resnet_act_fn=act_fn,
            output_scale_factor=1,
            resnet_groups=norm_num_groups,
        )

        self.conv_norm_out = nn.GroupNorm(num_channels=block_out_channels[-1], num_groups=norm_num_groups, eps=1e-6)
        self.conv_act = nn.SiLU()
        conv_out_channels = 2 * out_channels if double_z else out_channels
        self.conv_out = nn.Conv3d(block_out_channels[-1], conv_out_channels, kernel_size=3, padding=1)

        self.gradient_checkpointing = False

    def forward(self, sample: torch.Tensor) -> torch.Tensor:
        sample = self.conv_in(sample)

        if self.training and self.gradient_checkpointing:
            def create_custom_forward(module):
                def custom_forward(*inputs):
                    return module(*inputs)
                return custom_forward

            if is_torch_version(">=", "1.11.0"):
                for down_block in self.down_blocks:
                    sample = torch.utils.checkpoint.checkpoint(
                        create_custom_forward(down_block), sample, use_reentrant=False
                    )
                sample = torch.utils.checkpoint.checkpoint(
                    create_custom_forward(self.mid_block), sample, use_reentrant=False
                )
            else:
                for down_block in self.down_blocks:
                    sample = torch.utils.checkpoint.checkpoint(create_custom_forward(down_block), sample)
                sample = torch.utils.checkpoint.checkpoint(create_custom_forward(self.mid_block), sample)
        else:
            for down_block in self.down_blocks:
                sample = down_block(sample)
            sample = self.mid_block(sample)

        sample = self.conv_norm_out(sample)
        sample = self.conv_act(sample)
        sample = self.conv_out(sample)
        return sample
