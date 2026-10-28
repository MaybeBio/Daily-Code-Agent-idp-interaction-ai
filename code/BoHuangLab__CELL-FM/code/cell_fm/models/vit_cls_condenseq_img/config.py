# -*- coding: utf-8 -*-
from dataclasses import dataclass
from transformers import PretrainedConfig

@dataclass
class ViTConfig(PretrainedConfig):
    model_type: str = 'vit'

    # Model parameters
    image_size: int = 160
    in_channels: int = 2
    patch_size: int = 4
    num_layers: int = 8
    num_heads: int = 8
    hidden_dim: int = 768
    mlp_dim: int = 3072
    num_classes: int = 2

    # Training parameters
    loadcheck_path: str = '.'
    ft: bool = False
    infer: bool = False
    ifresume: bool = False

    output_dir: str = ""
    learning_rate: float = 1e-4
    weight_decay: float = 0.0
    gradient_accumulation_steps: int = 1
    per_device_train_batch_size: int = 2
    per_device_eval_batch_size: int = 2

    num_train_epochs: int = 10
    fp16: bool = False
    bf16: bool = False
    logging_dir: str = ""
    logging_steps: int = 10
    max_steps: int = -1
    warmup_steps: int = 1000
    save_steps: int = 1000

    dataloader_num_workers: int = 8
    seed: int = 6

    def __init__(self, **kwargs):
        # Use `super().__init__` to handle arguments from PretrainedConfig
        super().__init__(**kwargs)