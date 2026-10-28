# -*- coding: utf-8 -*-
from typing import Optional

import torch
import torch.nn.functional as F
from torch import nn
from diffusers.models.normalization import RMSNorm
from diffusers.utils import is_torch_version
from diffusers.models.activations import get_activation


class ResnetBlock3D(nn.Module):
    def __init__(
        self,
        *,
        in_channels: int,
        out_channels: Optional[int] = None,
        conv_shortcut: bool = False,
        dropout: float = 0.0,
        groups: int = 32,
        groups_out: Optional[int] = None,
        eps: float = 1e-6,
        non_linearity: str = "swish",
        output_scale_factor: float = 1.0,
        use_in_shortcut: Optional[bool] = None,
        conv_shortcut_bias: bool = True,
    ):
        super().__init__()

        self.in_channels = in_channels
        out_channels = in_channels if out_channels is None else out_channels
        self.out_channels = out_channels
        self.use_conv_shortcut = conv_shortcut
        self.output_scale_factor = output_scale_factor

        if groups_out is None:
            groups_out = groups

        self.norm1 = nn.GroupNorm(num_groups=groups, num_channels=in_channels, eps=eps, affine=True)
        self.conv1 = nn.Conv3d(in_channels, out_channels, kernel_size=3, stride=1, padding=1)
        self.norm2 = nn.GroupNorm(num_groups=groups_out, num_channels=out_channels, eps=eps, affine=True)
        self.dropout = nn.Dropout(dropout)
        self.conv2 = nn.Conv3d(out_channels, out_channels, kernel_size=3, stride=1, padding=1)
        self.nonlinearity = get_activation(non_linearity)

        self.use_in_shortcut = (in_channels != out_channels) if use_in_shortcut is None else use_in_shortcut
        self.conv_shortcut = None
        if self.use_in_shortcut:
            self.conv_shortcut = nn.Conv3d(
                in_channels, out_channels, kernel_size=1, stride=1, padding=0, bias=conv_shortcut_bias,
            )

    def forward(self, input_tensor: torch.Tensor) -> torch.Tensor:
        hidden_states = input_tensor
        hidden_states = self.norm1(hidden_states)
        hidden_states = self.nonlinearity(hidden_states)
        hidden_states = self.conv1(hidden_states)
        hidden_states = self.norm2(hidden_states)
        hidden_states = self.nonlinearity(hidden_states)
        hidden_states = self.dropout(hidden_states)
        hidden_states = self.conv2(hidden_states)

        if self.conv_shortcut is not None:
            input_tensor = self.conv_shortcut(input_tensor)

        return (input_tensor + hidden_states) / self.output_scale_factor


class DownEncoderBlock3D(nn.Module):
    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        dropout: float = 0.0,
        num_layers: int = 1,
        resnet_eps: float = 1e-6,
        resnet_act_fn: str = "swish",
        resnet_groups: int = 32,
        output_scale_factor: float = 1.0,
        add_downsample: bool = True,
        downsample_padding: int = 1,
    ):
        super().__init__()
        resnets = []
        for i in range(num_layers):
            resnets.append(
                ResnetBlock3D(
                    in_channels=in_channels if i == 0 else out_channels,
                    out_channels=out_channels,
                    eps=resnet_eps,
                    groups=resnet_groups,
                    dropout=dropout,
                    non_linearity=resnet_act_fn,
                    output_scale_factor=output_scale_factor,
                )
            )
        self.resnets = nn.ModuleList(resnets)
        self.downsamplers = (
            nn.ModuleList([Downsample3D(out_channels, use_conv=True, out_channels=out_channels, padding=downsample_padding)])
            if add_downsample else None
        )

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        for resnet in self.resnets:
            hidden_states = resnet(hidden_states)
        if self.downsamplers is not None:
            for downsampler in self.downsamplers:
                hidden_states = downsampler(hidden_states)
        return hidden_states


class UpDecoderBlock3D(nn.Module):
    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        dropout: float = 0.0,
        num_layers: int = 1,
        resnet_eps: float = 1e-6,
        resnet_act_fn: str = "swish",
        resnet_groups: int = 32,
        output_scale_factor: float = 1.0,
        add_upsample: bool = True,
    ):
        super().__init__()
        resnets = []
        for i in range(num_layers):
            resnets.append(
                ResnetBlock3D(
                    in_channels=in_channels if i == 0 else out_channels,
                    out_channels=out_channels,
                    eps=resnet_eps,
                    groups=resnet_groups,
                    dropout=dropout,
                    non_linearity=resnet_act_fn,
                    output_scale_factor=output_scale_factor,
                )
            )
        self.resnets = nn.ModuleList(resnets)
        self.upsamplers = (
            nn.ModuleList([Upsample3D(out_channels, out_channels=out_channels)])
            if add_upsample else None
        )

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        for resnet in self.resnets:
            hidden_states = resnet(hidden_states)
        if self.upsamplers is not None:
            for upsampler in self.upsamplers:
                hidden_states = upsampler(hidden_states)
        return hidden_states


class Downsample3D(nn.Module):
    def __init__(
        self,
        channels: int,
        use_conv: bool = False,
        out_channels: Optional[int] = None,
        padding: int = 1,
        kernel_size: int = 3,
        norm_type: Optional[str] = None,
        eps: Optional[float] = 1e-5,
        elementwise_affine: Optional[bool] = True,
        bias: bool = True,
    ):
        super().__init__()
        self.channels = channels
        self.out_channels = out_channels or channels
        self.use_conv = use_conv
        self.padding = padding
        stride = 2

        if norm_type == "ln_norm":
            self.norm = nn.LayerNorm(self.channels, eps=eps, elementwise_affine=elementwise_affine)
        elif norm_type == "rms_norm":
            self.norm = RMSNorm(channels, eps, elementwise_affine)
        elif norm_type is None:
            self.norm = None
        else:
            raise ValueError(f"Unknown norm_type: {norm_type}")

        if use_conv:
            self.conv = nn.Conv3d(self.channels, self.out_channels, kernel_size=kernel_size, stride=stride, padding=padding, bias=bias)
        else:
            assert self.channels == self.out_channels
            self.conv = nn.AvgPool3d(kernel_size=stride, stride=stride)

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        assert hidden_states.shape[1] == self.channels

        if self.norm is not None:
            hidden_states = self.norm(hidden_states.permute(0, 2, 3, 4, 1)).permute(0, 4, 1, 2, 3)

        if self.use_conv and self.padding == 0:
            hidden_states = F.pad(hidden_states, (0, 1, 0, 1, 0, 1), mode="constant", value=0)

        return self.conv(hidden_states)


class Upsample3D(nn.Module):
    def __init__(
        self,
        channels: int,
        out_channels: Optional[int] = None,
        kernel_size: int = 3,
        padding: int = 1,
        norm_type: Optional[str] = None,
        eps: Optional[float] = None,
        elementwise_affine: Optional[bool] = None,
        bias: bool = True,
    ):
        super().__init__()
        self.channels = channels
        self.out_channels = out_channels or channels

        if norm_type == "ln_norm":
            self.norm = nn.LayerNorm(self.channels, eps=eps, elementwise_affine=elementwise_affine)
        elif norm_type == "rms_norm":
            self.norm = RMSNorm(channels, eps, elementwise_affine)
        elif norm_type is None:
            self.norm = None
        else:
            raise ValueError(f"Unknown norm_type: {norm_type}")

        self.conv = nn.Conv3d(self.channels, self.out_channels, kernel_size=kernel_size, padding=padding, bias=bias)

    def forward(self, hidden_states: torch.Tensor, output_size: Optional[int] = None) -> torch.Tensor:
        assert hidden_states.shape[1] == self.channels

        if self.norm is not None:
            hidden_states = self.norm(hidden_states.permute(0, 2, 3, 4, 1)).permute(0, 4, 1, 2, 3)

        dtype = hidden_states.dtype
        if dtype == torch.bfloat16 and is_torch_version("<", "2.1"):
            hidden_states = hidden_states.to(torch.float32)

        if output_size is None:
            D = hidden_states.shape[2]
            if D >= 2:
                # Split depth to avoid OOM on large 3D volumes
                mid = D // 2
                x1 = F.interpolate(hidden_states[:, :, :mid], scale_factor=2.0, mode="nearest")
                x2 = F.interpolate(hidden_states[:, :, mid:], scale_factor=2.0, mode="nearest")
                hidden_states = torch.cat([x1, x2], dim=2)
            else:
                hidden_states = F.interpolate(hidden_states, scale_factor=2.0, mode="nearest")
        else:
            hidden_states = F.interpolate(hidden_states, size=output_size, mode="nearest")

        if dtype == torch.bfloat16 and is_torch_version("<", "2.1"):
            hidden_states = hidden_states.to(dtype)

        return self.conv(hidden_states)


class UNetMidBlock3D(nn.Module):
    def __init__(
        self,
        in_channels: int,
        dropout: float = 0.0,
        num_layers: int = 1,
        resnet_eps: float = 1e-6,
        resnet_act_fn: str = "swish",
        resnet_groups: int = 32,
        output_scale_factor: float = 1.0,
    ):
        super().__init__()
        resnet_groups = resnet_groups if resnet_groups is not None else min(in_channels // 4, 32)

        # always at least one resnet before and one after
        resnets = [
            ResnetBlock3D(
                in_channels=in_channels,
                out_channels=in_channels,
                eps=resnet_eps,
                groups=resnet_groups,
                dropout=dropout,
                non_linearity=resnet_act_fn,
                output_scale_factor=output_scale_factor,
            )
        ]
        for _ in range(num_layers):
            resnets.append(
                ResnetBlock3D(
                    in_channels=in_channels,
                    out_channels=in_channels,
                    eps=resnet_eps,
                    groups=resnet_groups,
                    dropout=dropout,
                    non_linearity=resnet_act_fn,
                    output_scale_factor=output_scale_factor,
                )
            )
        self.resnets = nn.ModuleList(resnets)

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        for resnet in self.resnets:
            hidden_states = resnet(hidden_states)
        return hidden_states
