# -*- coding: utf-8 -*-
import os
import torch
import torch.nn as nn
from .vae_3d_config import VAE3DConfig
from .modules.autoencoders import Autoencoder3DKL
from transformers import PreTrainedModel
from cell_fm.pipeline.utils import VAEOutput
from cell_fm.logging import logger


class VAE3DModel(PreTrainedModel):
    config_class = VAE3DConfig

    def __init__(self, config: VAE3DConfig):
        super().__init__(config)
        self.config = config

        self.num_down_blocks = config.num_down_blocks
        self.num_up_blocks = config.num_down_blocks

        self.vae = Autoencoder3DKL(
            in_channels=config.in_channels,
            out_channels=config.out_channels,
            num_down_blocks=config.num_down_blocks,
            block_out_channels=tuple(config.vae_block_out_channels),
            layers_per_block=config.layers_per_block,
            latent_channels=config.latent_channels,
            norm_num_groups=config.norm_num_groups,
        )

        self.load_pretrained_weights(config, checkpoint_path=config.vae_loadcheck_path)

    def load_pretrained_weights(self, config, checkpoint_path):
        if config.ft or config.infer:
            if config.ft:
                logger.info(f"Finetune from checkpoint: {checkpoint_path}")
            else:
                logger.info(f"Infer from checkpoint: {checkpoint_path}")

            if os.path.splitext(checkpoint_path)[1] == '.safetensors':
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

    def encode(self, x: torch.Tensor):
        """Encodes input volume into latent distribution."""
        return self.vae.encode(x).latent_dist

    def decode(self, latents: torch.Tensor):
        """Decodes latents back to volume."""
        return self.vae.decode(latents)

    def forward(self, batched_data):
        x = batched_data['protein_img']   # (B, C, D, H, W)

        latent_dist = self.encode(x)
        latents = latent_dist.sample()    # reparameterization trick
        recon_x = self.decode(latents).sample
        total_loss, recon_loss, kl_loss = self.compute_loss(x, recon_x, latent_dist)

        log_loss = {
            "total_loss": total_loss.item(),
            "recon_loss": recon_loss.item(),
            "kl_loss":    kl_loss.item(),
        }
        return VAEOutput(total_loss, log_loss)

    def compute_loss(self, x, recon_x, latent_dist):
        recon_loss = nn.MSELoss()(recon_x, x)
        kl_loss = -0.5 * torch.mean(
            1 + latent_dist.logvar - latent_dist.mean.pow(2) - latent_dist.logvar.exp()
        )
        total_loss = self.config.recon_loss_coeff * recon_loss + self.config.kl_loss_coeff * kl_loss
        return total_loss, recon_loss, kl_loss

    def sample(self, num_samples: int = 1, device: str = "cpu") -> torch.Tensor:
        """Generate samples by decoding random latents.

        Latent spatial size is derived from input_spatial_size and num_down_blocks.
        The encoder performs (num_down_blocks - 1) actual downsamples (the last block
        skips downsampling), so the spatial reduction factor is 2^(num_down_blocks-1).
        """
        n_ds = self.config.num_down_blocks - 1
        d_lat  = self.config.input_spatial_size[0] // (2 ** n_ds)
        hw_lat = self.config.input_spatial_size[1] // (2 ** n_ds)
        latents = torch.randn(
            (num_samples, self.config.latent_channels, d_lat, hw_lat, hw_lat), device=device
        )
        with torch.no_grad():
            return self.decode(latents).sample

    def reconstruct(self, x: torch.Tensor) -> torch.Tensor:
        latent_dist = self.encode(x)
        latents = latent_dist.sample()
        return self.decode(latents).sample
