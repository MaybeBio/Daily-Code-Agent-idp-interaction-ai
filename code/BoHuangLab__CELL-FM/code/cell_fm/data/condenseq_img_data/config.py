# -*- coding: utf-8 -*-
from dataclasses import dataclass
from transformers import PretrainedConfig

@dataclass
class CondenSeqImageDatasetConfig(PretrainedConfig):
    # Dataset parameters
    data_path: str = ""
    split_key: str = 'train'

    data_aug: bool = False
    img_resize: int = 160

    def __init__(self, **kwargs):
        # Use `super().__init__` to handle arguments from PretrainedConfig
        super().__init__(**kwargs)