import os
import torch
import torch.nn as nn
from cell_fm.logging import logger
from transformers import PreTrainedModel
from .config import ViTConfig

from torchvision.models.vision_transformer import VisionTransformer
from functools import partial

from cell_fm.pipeline.utils import VanillaOutput
import json


class ViTModel(PreTrainedModel):
    config_class = ViTConfig

    def __init__(self, config: ViTConfig):
        super().__init__(config)
        self.config = config

        self.loss = nn.CrossEntropyLoss()
        self.net = ViT(config)

        self.load_pretrained_weights(config, checkpoint_path=config.loadcheck_path)

    def load_pretrained_weights(self, config, checkpoint_path):
        """
        Load pretrained weights from a given state_dict.
        """
        if config.ft or config.infer:
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
                    shard_weights[shard_file] = torch.load(os.path.join(checkpoint_path, shard_file), map_location="cpu")

                checkpoints_state = {}
                for param_name, shard_file in weight_map.items():
                    checkpoints_state[param_name] = shard_weights[shard_file][param_name]
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
            filtered_state_dict = {k: v for k, v in checkpoints_state.items() if k in model_state_dict and v.size() == model_state_dict[k].size()}

            IncompatibleKeys = self.load_state_dict(filtered_state_dict, strict=False)
            IncompatibleKeys = IncompatibleKeys._asdict()

            missing_keys = []
            for keys in IncompatibleKeys["missing_keys"]:
                if keys.find("dummy") == -1:
                    missing_keys.append(keys)

            unexpected_keys = []
            for keys in IncompatibleKeys["unexpected_keys"]:
                if keys.find("dummy") == -1:
                    unexpected_keys.append(keys)

            if len(missing_keys) > 0:
                logger.info(
                    "Missing keys in {}: {}".format(
                        checkpoint_path,
                        missing_keys,
                    )
                )

            if len(unexpected_keys) > 0:
                logger.info(
                    "Unexpected keys {}: {}".format(
                        checkpoint_path,
                        unexpected_keys,
                    )
                )

    def forward(self, batched_data, **kwargs):
        nucleus_img = batched_data['nucleus_img']
        protein_img = batched_data['protein_img']
        label = batched_data['label']

        input_img = torch.cat([nucleus_img, protein_img], dim=1)  # [B, 2, H, W]

        logits = self.net(input_img)  # [B, num_classes]

        loss = self.loss(logits, label.view(-1))
        loss = loss.mean()

        return VanillaOutput(loss=loss)

    def predict(self, x):
        with torch.no_grad():
            input_img = x  # [B, 2, H, W]
            logits = self.net(input_img)  # [B, num_classes]

        return logits

class ViT(VisionTransformer):
    def __init__(self, config: ViTConfig):
        super().__init__(
            image_size=config.image_size,
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

        old_conv = self.conv_proj
        in_channels = config.in_channels

        self.conv_proj = nn.Conv2d(
            in_channels=in_channels,
            out_channels=old_conv.out_channels,
            kernel_size=old_conv.kernel_size,
            stride=old_conv.stride,
            padding=old_conv.padding,
            bias=old_conv.bias is not None,
        )

    def forward(self, x: torch.Tensor):
        # Reshape and permute the input tensor
        x = self._process_input(x)
        n = x.shape[0]

        # Expand the class token to the full batch
        batch_class_token = self.class_token.expand(n, -1, -1)
        x = torch.cat([batch_class_token, x], dim=1)

        x = self.encoder(x)
        x = x.mean(dim=1)  # Global average pooling over the sequence length
        x = self.heads(x)

        return x

    def embed(self, x: torch.Tensor):
        # Reshape and permute the input tensor
        x = self._process_input(x)
        n = x.shape[0]

        # Expand the class token to the full batch
        batch_class_token = self.class_token.expand(n, -1, -1)
        x = torch.cat([batch_class_token, x], dim=1)

        x = self.encoder(x)
        embed = x.mean(dim=1)  # Global average pooling over the sequence length

        return embed