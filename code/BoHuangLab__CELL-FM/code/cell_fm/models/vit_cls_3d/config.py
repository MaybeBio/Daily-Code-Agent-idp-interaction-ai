# -*- coding: utf-8 -*-
from dataclasses import dataclass, field
from typing import List
from transformers import PretrainedConfig


@dataclass
class ViT3DConfig(PretrainedConfig):
    model_type: str = 'vit_3d'

    # Dataset parameters
    data_path: str = ""
    split_key: str = 'all'
    phase: str = 'train'

    data_aug: bool = False

    # Input volume shape — stored as comma-separated string, parsed to int list in __init__
    input_spatial_size: str = '48,192,192'

    # Model parameters
    in_channels: int = 2          # nucleus + protein
    patch_d: int = 8              # depth patch size
    patch_size: int = 16          # H/W patch size
    num_layers: int = 12
    num_heads: int = 8
    hidden_dim: int = 768
    mlp_dim: int = 3072
    num_classes: int = 1311

    # Training parameters
    loadcheck_path: str = '.'
    ft: bool = False
    infer: bool = False
    ifresume: bool = False

    output_dir: str = ""
    learning_rate: float = 1e-4
    weight_decay: float = 0.0
    gradient_accumulation_steps: int = 1
    per_device_train_batch_size: int = 4
    per_device_eval_batch_size: int = 4

    num_train_epochs: int = 10
    fp16: bool = False
    bf16: bool = False
    logging_dir: str = ""
    logging_steps: int = 100
    max_steps: int = -1
    warmup_steps: int = 1000
    save_steps: int = 2000

    dataloader_num_workers: int = 8
    seed: int = 6

    # evaluation
    num_steps: int = 100

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Parse input_spatial_size from comma-separated string to int list
        if isinstance(self.input_spatial_size, str):
            self.input_spatial_size = [int(x) for x in self.input_spatial_size.split(',')]
