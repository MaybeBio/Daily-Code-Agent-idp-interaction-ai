# -*- coding: utf-8 -*-
import os
import json
from copy import deepcopy

import torch
import torch.nn as nn

from transformers import PreTrainedModel

from cell_fm.models.cell_fm.modules.protein_sequence_embedding import ESMEmbed
from cell_fm.models.cell_fm.modules.transport import create_transport, Sampler
from cell_fm.models.vae_3d.vae_3d_model import VAE3DModel
from cell_fm.models.vae_3d.vae_3d_config import VAE3DConfig
from cell_fm.pipeline.utils import CELLFMOutput
from cell_fm.logging import logger

from .cell_fm_3d_config import CELLFM3DConfig
from .modules import SD3Transformer3DModel, CondConvNet3D, UNet3D


# The transport runs on t in [0, 1], but the sinusoidal timestep embedding SD3 uses internally
# (diffusers' Timesteps) inherits ADM's frequency schedule, which is built for t in [0, 1000].
# Fed t in [0, 1] its lowest frequencies see arguments of ~1e-4 rad, so cos rounds to exactly 1
# and the code collapses: at t=0.3 vs t=0.7, 94 of 256 features are bit-identical in bf16 (2
# after scaling). Units, not a hyperparameter -- deliberately not a config field, since a
# pretrain/eval mismatch here would silently corrupt inference.
TIMESTEP_SCALE = 1000.0


class CELLFM3DModel(PreTrainedModel):
    config_class = CELLFM3DConfig

    def __init__(self, config: CELLFM3DConfig):
        super().__init__(config)
        self.config = config

        self.net = CELLFM3D(config)

        self.transport = create_transport(
            config.path_type,
            config.prediction,
            config.loss_weight,
            config.train_eps,
            config.sample_eps,
            timestep_sampler=config.timestep_sampler,
            logit_mean=config.logit_mean,
            logit_std=config.logit_std,
        )
        self.transport_sampler = Sampler(self.transport)

        self.vae = self.initialize_vae(config)
        self.load_pretrained_weights(config, checkpoint_path=config.loadcheck_path)
        self.vae = self.prepare_vae(self.vae, config.vae_loadcheck_path)

    def initialize_vae(self, config: CELLFM3DConfig) -> VAE3DModel:
        vae_config = VAE3DConfig(
            in_channels=config.in_channels,
            out_channels=config.out_channels,
            num_down_blocks=config.num_down_blocks,
            latent_channels=config.latent_channels,
            vae_block_out_channels=config.vae_block_out_channels,
            input_spatial_size=config.input_spatial_size,
            norm_num_groups=config.norm_num_groups,
            layers_per_block=config.layers_per_block,
            ft=False,
            infer=False,
            vae_loadcheck_path="",
        )
        return VAE3DModel(vae_config)

    def prepare_vae(self, vae: VAE3DModel, checkpoint_path: str) -> VAE3DModel:
        if checkpoint_path and checkpoint_path != '.':
            vae_config = deepcopy(vae.config)
            vae_config.ft = False
            vae_config.infer = True
            vae.load_pretrained_weights(vae_config, checkpoint_path=checkpoint_path)
        for param in vae.parameters():
            param.requires_grad = False
        vae.eval()
        return vae

    def load_pretrained_weights(self, config: CELLFM3DConfig, checkpoint_path: str):
        if not (config.ft or config.infer):
            return

        if config.ft:
            logger.info(f"Finetune from checkpoint: {checkpoint_path}")
        else:
            logger.info(f"Infer from checkpoint: {checkpoint_path}")

        if os.path.isdir(checkpoint_path):
            with open(os.path.join(checkpoint_path, "pytorch_model.bin.index.json")) as f:
                weight_map = json.load(f)["weight_map"]
            shard_weights = {}
            for shard_file in set(weight_map.values()):
                shard_weights[shard_file] = torch.load(
                    os.path.join(checkpoint_path, shard_file), map_location="cpu"
                )
            checkpoints_state = {k: shard_weights[v][k] for k, v in weight_map.items()}
        elif checkpoint_path.endswith('.safetensors'):
            from safetensors.torch import load_file
            checkpoints_state = load_file(checkpoint_path)
        else:
            checkpoints_state = torch.load(checkpoint_path, map_location="cpu")

        if "model" in checkpoints_state:
            checkpoints_state = checkpoints_state["model"]
        elif "module" in checkpoints_state:
            checkpoints_state = checkpoints_state["module"]

        model_state_dict = self.state_dict()
        filtered = {
            k: v for k, v in checkpoints_state.items()
            if k in model_state_dict and v.size() == model_state_dict[k].size()
        }
        self.load_state_dict(filtered, strict=False)

        # Anything dropped here silently keeps its random init and the model still runs, so a
        # config that disagrees with the checkpoint only ever surfaces as bad samples. Warn, and
        # keep the shape mismatches separate: those mean a geometry field (patch_size,
        # unet_block_out_channels, num_attention_heads, ...) differs from the run that wrote the
        # checkpoint, which is a mistake, whereas a plain absence is the ordinary case for a module
        # the checkpoint predates.
        mismatched = sorted(
            k for k, v in checkpoints_state.items()
            if k in model_state_dict and v.size() != model_state_dict[k].size()
        )
        absent = sorted(set(model_state_dict) - set(filtered) - set(mismatched))
        unknown = sorted(set(checkpoints_state) - set(model_state_dict))

        logger.info(
            f"Loaded {len(filtered)}/{len(model_state_dict)} tensors from {checkpoint_path}"
        )
        if mismatched:
            shapes = [
                (k, tuple(checkpoints_state[k].shape), tuple(model_state_dict[k].shape))
                for k in mismatched[:3]
            ]
            logger.warning(
                f"{len(mismatched)} tensors in {checkpoint_path} have a different shape and were "
                f"NOT loaded; they keep their random init. Check the geometry config against the "
                f"run that wrote this checkpoint. (checkpoint vs model): {shapes}"
            )
        if absent:
            logger.warning(
                f"{len(absent)} model tensors are absent from {checkpoint_path} and keep their "
                f"random init (e.g. {absent[:5]})"
            )
        if unknown:
            logger.info(
                f"{len(unknown)} tensors in {checkpoint_path} are not part of this model "
                f"(e.g. {unknown[:5]})"
            )

    @torch.no_grad()
    def prepare_data(self, batched_data):
        protein_img = batched_data['protein_img']               # (B, 1, D, H, W)
        protein_img_latent = self.vae.encode(protein_img).sample()

        nucleus_img = batched_data['nucleus_img']               # (B, 1, D, H, W)
        if self.config.cell_image == 'nucl':
            cell_img = nucleus_img
        else:
            raise ValueError(f"Cell image type '{self.config.cell_image}' is not supported")

        return protein_img_latent, cell_img

    def forward(self, batched_data, **kwargs):
        protein_seq = batched_data['protein_seq']
        protein_img_latent, cell_img = self.prepare_data(batched_data)

        # Upcast before sampling t: under bf16 autocast the VAE latent is bf16, and
        # `Transport.sample` ends with `t.to(x1)`, so t would inherit bf16's 8-bit mantissa. After
        # TIMESTEP_SCALE that leaves ~1378 distinct t_cond values over [0, 1000] (fp32 gives ~2e5)
        # with errors up to 4 -- enough to alias the top 13% of the 256 sinusoidal timestep
        # channels, and to put inference (which runs in fp32, off that grid) on a schedule training
        # never saw. Same units-vs-precision trap as TIMESTEP_SCALE, at the other end of the
        # frequency range. Also keeps the regression target ut out of bf16.
        t, x0, x1 = self.transport.sample(protein_img_latent.float())
        t, xt, ut = self.transport.path_sampler.plan(t, x0, x1)

        img_output = self.net(xt, cell_img, protein_seq, t)
        loss = self.transport.training_losses(img_output, x0, xt, ut, t)["loss"].mean()

        log_output = {
            "train_loss": loss.item(),
        }
        return CELLFMOutput(loss=loss, log_output=log_output)

    @torch.no_grad()
    def sequence_to_image(self, protein_seq: torch.Tensor, cell_img: torch.Tensor, num_steps: int = 100) -> torch.Tensor:
        device = protein_seq.device
        B = protein_seq.shape[0]

        # Derive latent spatial size: (num_down_blocks - 1) actual downsamples
        n_ds = self.config.num_down_blocks - 1
        D, H, W = self.config.input_spatial_size
        d_lat = D // (2 ** n_ds)
        h_lat = H // (2 ** n_ds)
        w_lat = W // (2 ** n_ds)

        noise = torch.randn(
            B, self.config.latent_channels, d_lat, h_lat, w_lat, device=device
        )
        sample_fn = self.transport_sampler.sample_ode(num_steps=num_steps)

        def model_fn(xt, t):
            return self.net(xt, cell_img, protein_seq, t)

        latent = sample_fn(noise, model_fn)[-1]
        return self.vae.decode(latent).sample


class CELLFM3D(nn.Module):

    def __init__(self, config: CELLFM3DConfig):
        super().__init__()

        D, H, W = config.input_spatial_size
        n_ds = config.num_down_blocks - 1
        d_lat = D // (2 ** n_ds)
        h_lat = H // (2 ** n_ds)
        w_lat = W // (2 ** n_ds)

        # SD3 runs at the UNet's lowest resolution: the latent halved once per extra UNet level
        unet_channels = config.unet_block_out_channels
        scale = 2 ** (len(unet_channels) - 1)
        d_sd, h_sd, w_sd = d_lat // scale, h_lat // scale, w_lat // scale
        pd, phw = config.patch_d, config.patch_size
        if d_lat % (scale * pd) or h_lat % (scale * phw) or w_lat % (scale * phw):
            raise ValueError(
                f"Latent {d_lat}x{h_lat}x{w_lat} must divide by the UNet's {scale}x downsample "
                f"times the SD3 patch {pd}x{phw}x{phw}"
            )

        # Cell image conditioning
        self.cond_conv = CondConvNet3D(
            config.in_channels * len(config.cell_image.split(',')),
            config.cond_out_channels,
        )

        # UNet on the noisy latent plus cell image features, conditioned on SD3's timestep +
        # pooled-sequence embedding. SD3 is its bottleneck.
        self.unet = UNet3D(
            in_channels=config.latent_channels + config.cond_out_channels[-1],
            out_channels=config.latent_channels,
            block_out_channels=unet_channels,
            layers_per_block=config.unet_layers_per_block,
            cond_dim=config.num_attention_heads * config.attention_head_dim,
            checkpointing=config.unet_checkpointing,
        )

        # Image generator (3-D SD3 transformer) at the UNet's lowest resolution
        self.img_generator = SD3Transformer3DModel(
            depth=d_sd,
            height=h_sd,
            width=w_sd,
            patch_d=config.patch_d,
            patch_hw=config.patch_size,
            latent_channels=unet_channels[-1],
            in_channels=unet_channels[-1],
            num_layers=config.img_generator_num_layers,
            attention_head_dim=config.attention_head_dim,
            num_attention_heads=config.num_attention_heads,
            joint_attention_dim=config.encoder_hidden_size,
            pooled_projection_dim=config.encoder_hidden_size,
            qk_norm=config.qk_norm,
        )

        # ESM-C embeddings go straight to the transformer. No adapter in between: SD3's
        # context_embedder is already an affine map on this tensor, so a Linear here would collapse
        # into it, and its AdaLayerNormZero re-normalizes the context at every block entry, so a
        # norm here buys nothing either. Position is likewise left to ESM-C's RoPE -- an additive
        # sincos term would carry norm sqrt(hidden/2) against an ESM signal of ~1.5 and bury it.
        self.initialize_protein_sequence_embedding(config)

        # Same id ESM-C uses to build its own internal mask, so the two agree by construction.
        self.seq_pad_token_id = self.protein_sequence_embedding.model.tokenizer.pad_token_id

    def initialize_protein_sequence_embedding(self, config: CELLFM3DConfig):
        self.protein_sequence_embedding = ESMEmbed(
            config.esm_embedding,
            config.encoder_hidden_size,
            config.esm_fixed_embedding,
        )

    def forward(
        self,
        protein_img_latent: torch.Tensor,
        cell_img: torch.Tensor,
        protein_seq: torch.Tensor,
        time: torch.Tensor,
    ) -> torch.Tensor:
        # Scaled here rather than at the call sites so training and sequence_to_image cannot diverge.
        t_cond = time * TIMESTEP_SCALE

        seq_embeds = self.protein_sequence_embedding(protein_seq)               # (B, T, hidden)
        seq_mask = protein_seq != self.seq_pad_token_id                         # (B, T), True = real

        # Global sequence summary for AdaLN, averaged over real tokens only. Dividing by T instead
        # would shrink a 127-residue protein's summary ~10x in a batch padded to 1367 -- the exact
        # failure that turned POLR2F into noise. Masked explicitly rather than leaning on ESM-C's
        # zero pad rows: ESMEmbed.scale_layer is an Identity only while encoder_hidden_size matches
        # the ESM width, and the biased Linear it becomes otherwise would make pad rows nonzero.
        seq_weight = seq_mask.unsqueeze(-1).to(seq_embeds.dtype)                # (B, T, 1)
        pooled_seq = (seq_embeds * seq_weight).sum(dim=1) / seq_weight.sum(dim=1).clamp(min=1.0)

        # SD3's timestep + pooled-sequence embedding, shared by the transformer and the UNet.
        temb = self.img_generator.time_text_embed(t_cond, pooled_seq)           # (B, inner_dim)

        cell_img_conv = self.cond_conv(cell_img)                                # (B, cond_C, d_lat, h_lat, w_lat)
        x = torch.cat([protein_img_latent, cell_img_conv], dim=1)               # (B, latent+cond_C, d_lat, h_lat, w_lat)

        x, skips = self.unet.encode(x, temb)                                    # (B, unet_C[-1], d_sd, h_sd, w_sd)
        x = self.img_generator(
            hidden_states=x,
            encoder_hidden_states=seq_embeds,
            temb=temb,
            encoder_attention_mask=seq_mask,
        )                                                                       # (B, unet_C[-1], d_sd, h_sd, w_sd)
        return self.unet.decode(x, skips, temb)                                 # (B, latent, d_lat, h_lat, w_lat)
