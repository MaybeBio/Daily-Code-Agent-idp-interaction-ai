# -*- coding: utf-8 -*-
from dataclasses import dataclass
from transformers import PretrainedConfig


@dataclass
class VAE3DConfig(PretrainedConfig):
    model_type: str = 'vae_3d'

    # Dataset parameters
    data_path: str = ""
    split_key: str = 'train'

    data_aug: bool = False
    img_resize: int = 1024
    img_crop_size: int = 256

    # Loss parameters
    recon_loss_coeff: float = 1.0
    kl_loss_coeff: float = 1.0

    # Model parameters
    in_channels: int = 1
    out_channels: int = 1
    num_down_blocks: int = 4
    latent_channels: int = 4
    vae_block_out_channels: str = '32,64,128,256'
    input_spatial_size: str = '32,512,512'   # D,H,W
    norm_num_groups: int = 32
    layers_per_block: int = 2

    # Training parameters
    vae_loadcheck_path: str = ""
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
        super().__init__(**kwargs)

        self.vae_block_out_channels = kwargs.get("vae_block_out_channels", self.vae_block_out_channels)
        if not isinstance(self.vae_block_out_channels, list):
            self.vae_block_out_channels = [int(c) for c in self.vae_block_out_channels.split(',')]

        self.input_spatial_size = kwargs.get("input_spatial_size", self.input_spatial_size)
        if not isinstance(self.input_spatial_size, list):
            self.input_spatial_size = [int(s) for s in self.input_spatial_size.split(',')]
