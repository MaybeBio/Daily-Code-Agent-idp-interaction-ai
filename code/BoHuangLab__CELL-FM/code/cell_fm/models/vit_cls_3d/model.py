# -*- coding: utf-8 -*-
import os
import json
import torch
import torch.nn as nn
from functools import partial

from transformers import PreTrainedModel
from torchvision.models.vision_transformer import VisionTransformer

from cell_fm.logging import logger
from cell_fm.pipeline.utils import VanillaOutput
from .config import ViT3DConfig


class ViT3DModel(PreTrainedModel):
    config_class = ViT3DConfig

    def __init__(self, config: ViT3DConfig):
        super().__init__(config)
        self.config = config
        self.loss = nn.CrossEntropyLoss()
        self.net = ViT3D(config)
        self.load_pretrained_weights(config, checkpoint_path=config.loadcheck_path)

    def load_pretrained_weights(self, config, checkpoint_path):
        if not (config.ft or config.infer):
            return

        if config.ft:
            logger.info(f"Finetune from checkpoint: {checkpoint_path}")
        else:
            logger.info(f"Infer from checkpoint: {checkpoint_path}")

        if os.path.isdir(checkpoint_path):
            with open(os.path.join(checkpoint_path, "pytorch_model.bin.index.json"), "r") as f:
                index_data = json.load(f)
            weight_map = index_data["weight_map"]
            shard_files = set(weight_map.values())
            shard_weights = {}
            for shard_file in shard_files:
                shard_weights[shard_file] = torch.load(
                    os.path.join(checkpoint_path, shard_file), map_location="cpu"
                )
            checkpoints_state = {
                k: shard_weights[v][k] for k, v in weight_map.items()
            }
        elif os.path.splitext(checkpoint_path)[1] == '.safetensors':
            from safetensors.torch import load_file
            checkpoints_state = load_file(checkpoint_path)
        else:
            checkpoints_state = torch.load(checkpoint_path, map_location="cpu")

        if "model" in checkpoints_state:
            checkpoints_state = checkpoints_state["model"]
        elif "module" in checkpoints_state:
            checkpoints_state = checkpoints_state["module"]

        model_state_dict = self.state_dict()
        filtered_state_dict = {
            k: v for k, v in checkpoints_state.items()
            if k in model_state_dict and v.size() == model_state_dict[k].size()
        }

        IncompatibleKeys = self.load_state_dict(filtered_state_dict, strict=False)._asdict()

        missing_keys = [k for k in IncompatibleKeys["missing_keys"] if "dummy" not in k]
        unexpected_keys = [k for k in IncompatibleKeys["unexpected_keys"] if "dummy" not in k]

        if missing_keys:
            logger.info(f"Missing keys in {checkpoint_path}: {missing_keys}")
        if unexpected_keys:
            logger.info(f"Unexpected keys in {checkpoint_path}: {unexpected_keys}")

    @torch.no_grad()
    def prepare_data(self, batched_data):
        # Concat nucleus + protein along channel dim: (B, 2, D, H, W)
        return torch.cat([batched_data['nucleus_img'], batched_data['protein_img']], dim=1)

    def forward(self, batched_data, **kwargs):
        input_img = self.prepare_data(batched_data)
        labels = batched_data['protein_idx']
        logits = self.net(input_img)          # (B, num_classes)
        loss = self.loss(logits, labels.view(-1))
        loss = loss.mean()
        return VanillaOutput(loss=loss)

    def embed(self, x: torch.Tensor):
        with torch.no_grad():
            return self.net.embed(x)


class ViT3D(VisionTransformer):
    def __init__(self, config: ViT3DConfig):
        D, H, W = config.input_spatial_size

        # Build the parent with dummy 2D size so encoder/heads/class_token are initialized
        super().__init__(
            image_size=H,
            patch_size=config.patch_size,
            num_layers=config.num_layers,
            num_heads=config.num_heads,
            hidden_dim=config.hidden_dim,
            mlp_dim=config.mlp_dim,
            dropout=0.0,
            attention_dropout=0.0,
            num_classes=config.num_classes,
            representation_size=None,
            norm_layer=partial(nn.LayerNorm, eps=1e-6),
            conv_stem_configs=None,
        )

        self.patch_d = config.patch_d
        self.patch_hw = config.patch_size
        self.hidden_dim = config.hidden_dim

        nd = D // config.patch_d
        nh = H // config.patch_size
        nw = W // config.patch_size
        n_patches = nd * nh * nw

        # Replace Conv2d patch projector with Conv3d
        self.conv_proj = nn.Conv3d(
            in_channels=config.in_channels,
            out_channels=config.hidden_dim,
            kernel_size=(config.patch_d, config.patch_size, config.patch_size),
            stride=(config.patch_d, config.patch_size, config.patch_size),
        )

        # Resize positional embedding to match 3D token count (n_patches + class token)
        self.encoder.pos_embedding = nn.Parameter(
            torch.zeros(1, n_patches + 1, config.hidden_dim)
        )
        nn.init.trunc_normal_(self.encoder.pos_embedding, std=0.02)

    def _process_input(self, x: torch.Tensor) -> torch.Tensor:
        # x: (B, C, D, H, W)
        n = x.shape[0]
        x = self.conv_proj(x)                                  # (B, hidden_dim, nd, nh, nw)
        x = x.reshape(n, self.hidden_dim, -1)                  # (B, hidden_dim, n_patches)
        x = x.permute(0, 2, 1)                                 # (B, n_patches, hidden_dim)
        return x

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self._process_input(x)
        n = x.shape[0]
        batch_class_token = self.class_token.expand(n, -1, -1)
        x = torch.cat([batch_class_token, x], dim=1)
        x = self.encoder(x)
        x = x.mean(dim=1)   # global average pooling
        x = self.heads(x)
        return x

    def embed(self, x: torch.Tensor) -> torch.Tensor:
        x = self._process_input(x)
        n = x.shape[0]
        batch_class_token = self.class_token.expand(n, -1, -1)
        x = torch.cat([batch_class_token, x], dim=1)
        x = self.encoder(x)
        return x.mean(dim=1)
