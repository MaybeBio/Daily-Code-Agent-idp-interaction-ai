# -*- coding: utf-8 -*-
from dataclasses import dataclass
from transformers import PretrainedConfig

@dataclass
class CondenSeqDatasetConfig(PretrainedConfig):
    # Dataset parameters
    data_path: str = ""
    split_key: str = 'train'
    phase: str = 'train'  # 'train', 'val', 'test'

    data_aug: bool = False
    img_resize: int = 512
    img_crop_size: int = 1024
    seq_zero_mask_ratio: float = 0.0
    seq_full_mask_ratio: float = 0.0
    img_type: str = 'GFP'  # 'GFP', 'SNAP'
    
    cell_image: str = 'nucl' # 'nucl', 'nucl,er', 'nucl,mt', 'nucl,er,mt'
    test_cell_image: str = 'nucl' # 'nucl', 'nucl,er', 'nucl,mt', 'nucl,er,mt'

    max_protein_sequence_len: int = 2048

    def __init__(self, **kwargs):
        # Use `super().__init__` to handle arguments from PretrainedConfig
        super().__init__(**kwargs)