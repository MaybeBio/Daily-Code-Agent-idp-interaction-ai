# -*- coding: utf-8 -*-
from dataclasses import dataclass
from transformers import PretrainedConfig

@dataclass
class CELLFMConfig(PretrainedConfig):
    model_type: str = 'flow_matching'

    # Dataset parameters
    data_path: str = ""
    split_key: str = 'train'
    phase: str = 'train'  # 'train', 'val', 'test'

    data_aug: bool = False
    img_resize: int = 512
    img_crop_size: int = 1024
    seq_zero_mask_ratio: float = 0.0
    seq_full_mask_ratio: float = 0.0
    pre_pad_seq: bool = False

    cell_image: str = 'nucl' # 'nucl', 'nucl,er', 'nucl,mt', 'nucl,er,mt'
    test_cell_image: str = 'nucl' # 'nucl', 'nucl,er', 'nucl,mt', 'nucl,er,mt'

    # Transport parameters
    path_type: str = 'Linear'
    prediction: str = 'velocity'
    loss_weight: str = None
    train_eps: float = None
    sample_eps: float = None

    # Model parameters
    ## VAE
    in_channels: int = 1
    out_channels: int = 1
    num_down_blocks: int = 3
    latent_channels: int = 4
    vae_block_out_channels: str = '128,256,512'

    ## CELL-Diff
    img_mask_ratio: float = 0.5
    cond_out_channels: str = '32,64'
    sample_size: int = 64
    esm_embedding: str = 'esmc_300m'
    esm_fixed_embedding: bool = True
    encoder_hidden_size: int = 320
    max_protein_sequence_len: int = 2048
    encoder_num_hidden_layers: int = 16
    num_heads: int = 8
    dim_head: int = 64
    vocab_size: int = 33
    dropout: float = 0.0
    final_dropout: float = 0.0
    encoder_patch_size: int = 8

    ### Image generator
    img_generator_num_layers: int = 18
    img_generator_patch_size: int = 4
    attention_head_dim: int = 64
    num_attention_heads: int = 16

    ### Image decoder
    img_decoder_num_hidden_layers: int = 8
    img_decoder_hidden_size: int = 512
    img_decoder_num_heads: int = 8
    img_decoder_dim_head: int = 64
    img_decoder_dropout: float = 0.0
    img_decoder_final_dropout: float = 0.0

    cell_image_ratio: float = 1.0

    # Loss parameters
    seq_loss_coeff: float = 1.0
    img_diff_loss_coeff: float = 1.0
    img_recon_loss_coeff: float = 1.0
    location_loss_coeff: float = 1.0

    # Training parameters
    vae_loadcheck_path: str = '.'
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

        self.vae_block_out_channels = kwargs.get("vae_block_out_channels", self.vae_block_out_channels)
        if not isinstance(self.vae_block_out_channels, list):
            self.vae_block_out_channels = [int(a) for a in self.vae_block_out_channels.split(',')]

        self.cond_out_channels = kwargs.get("cond_out_channels", self.cond_out_channels)
        if not isinstance(self.cond_out_channels, list):
            self.cond_out_channels = [int(a) for a in self.cond_out_channels.split(',')]