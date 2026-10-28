# -*- coding: utf-8 -*-
from dataclasses import dataclass
from transformers import PretrainedConfig

@dataclass
class ViTConfig(PretrainedConfig):
    model_type: str = 'vit'

    # Dataset parameters
    data_path: str = ""
    split_key: str = 'train'
    phase: str = 'train'  # 'train', 'val', 'test'

    data_aug: bool = False
    normalize: bool = False
    img_resize: int = 256
    img_crop_size: int = 1024

    cell_image: str = 'nucl' # 'nucl', 'nucl,er', 'nucl,mt', 'nucl,er,mt'

    # Model parameters
    patch_size: int = 4
    num_layers: int = 12
    num_heads: int = 8
    hidden_dim: int = 768
    mlp_dim: int = 3072
    num_classes: int = 1000

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
    
    # evaluation
    num_steps: int = 100

    def __init__(self, **kwargs):
        # Use `super().__init__` to handle arguments from PretrainedConfig
        super().__init__(**kwargs)