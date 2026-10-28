import os
import torch
import torch.nn as nn
import numpy as np
from cell_fm.logging import logger

from diffusers.models.embeddings import PatchEmbed
from diffusers.models.transformers import SD3Transformer2DModel

from .modules import TimestepEmbedder, CondConvNet
from .modules.positional_embedding import get_1d_sincos_pos_embed, get_2d_sincos_pos_embed
from .modules.transformer import TransformerBlock, MLMHead
from .modules.protein_sequence_embedding import ESMEmbed
from .modules.transport import create_transport, Sampler

from transformers import PreTrainedModel
from .cell_fm_config import CELLFMConfig

from cell_fm.pipeline.utils import CELLFMOutput
from cell_fm.models.vae.vae_model import VAEModel
import json
from copy import deepcopy


class CELLFMModel(PreTrainedModel):
    config_class = CELLFMConfig

    def __init__(self, config: CELLFMConfig, loss_fn=None):
        super().__init__(config)
        self.config = config

        self.loss = loss_fn(config)
        self.net = self.initialize_net(config)

        self.transport = create_transport(
            config.path_type, 
            config.prediction, 
            config.loss_weight, 
            config.train_eps, 
            config.sample_eps, 
        )
        self.transport_sampler = Sampler(self.transport)

        self.vae = self.initialize_vae()
        self.load_pretrained_weights(config, checkpoint_path=config.loadcheck_path)
        self.vae = self.prepare_vae(self.vae, self.config.vae_loadcheck_path)

    def initialize_net(self, config: CELLFMConfig):
        return CELLFM(config)

    def initialize_vae(self):
        vae_config = deepcopy(self.config)
        vae_config.ft = False
        vae_config.infer = False
        vae_config.vae_loadcheck_path = None
        vae = VAEModel(vae_config)

        return vae

    def prepare_vae(self, vae: VAEModel, checkpoint_path):
        if checkpoint_path is not None:
            vae_config = deepcopy(self.config)
            vae_config.ft = False
            vae_config.infer = True
            vae.load_pretrained_weights(vae_config, checkpoint_path=checkpoint_path)

        for param in vae.parameters():
            param.requires_grad = False
        vae.eval()

        return vae

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

    @torch.no_grad()
    def prepare_data(self, batched_data, vae: VAEModel):

        # protein_img = batched_data['protein_img']
        # nucleus_img = batched_data['nucleus_img']
        # microtubules_img = batched_data['microtubules_img']
        # protein_seq = batched_data['protein_seq']
        # protein_seq_masked = batched_data['protein_seq_masked']
        # protein_seq_mask = batched_data['protein_seq_mask']

        protein_img = batched_data['protein_img']
        nucleus_img = batched_data['nucleus_img']
        protein_img_latent = vae.encode(protein_img).sample()

        if self.config.cell_image == 'nucl':
            cell_img = nucleus_img

        elif self.config.cell_image == 'nucl,er':
            if 'ER_img' in batched_data:
                ER_img = batched_data['ER_img']
                mask = (torch.rand(ER_img.shape[0]) <= self.config.cell_image_ratio).bool()
                ER_img[~mask] = torch.full(size = [ER_img.shape[2], ER_img.shape[3]], fill_value=-2.0, device=ER_img.device)
            else:
                ER_img = torch.full_like(protein_img, fill_value=-2.0)
            cell_img = torch.cat([nucleus_img, ER_img], dim=1)

        elif self.config.cell_image == 'nucl,mt':
            if 'microtubules_img' in batched_data:
                microtubules_img = batched_data['microtubules_img']
                mask = (torch.rand(microtubules_img.shape[0]) <= self.config.cell_image_ratio).bool()
                microtubules_img[~mask] = torch.full(size = [microtubules_img.shape[2], microtubules_img.shape[3]], fill_value=-2.0, device=microtubules_img.device)
            else:
                microtubules_img = torch.full_like(protein_img, fill_value=-2.0)
            cell_img = torch.cat([nucleus_img, microtubules_img], dim=1)

        elif self.config.cell_image == 'nucl,er,mt':
            if 'ER_img' in batched_data:
                ER_img = batched_data['ER_img']
                mask = (torch.rand(ER_img.shape[0]) <= self.config.cell_image_ratio).bool()
                ER_img[~mask] = torch.full(size = [ER_img.shape[2], ER_img.shape[3]], fill_value=-2.0, device=ER_img.device)
            else:
                ER_img = torch.full_like(protein_img, fill_value=-2.0)

            if 'microtubules_img' in batched_data:
                microtubules_img = batched_data['microtubules_img']
                mask = (torch.rand(microtubules_img.shape[0]) <= self.config.cell_image_ratio).bool()
                microtubules_img[~mask] = torch.full(size = [microtubules_img.shape[2], microtubules_img.shape[3]], fill_value=-2.0, device=microtubules_img.device)
            else:
                microtubules_img = torch.full_like(protein_img, fill_value=-2.0)
            cell_img = torch.cat([nucleus_img, ER_img, microtubules_img], dim=1)
        else:
            raise ValueError(f"Cell image type: {self.config.cell_image} is not supported")
        
        return protein_img_latent, cell_img

    def forward(self, batched_data, **kwargs):
        protein_seq_masked = batched_data['protein_seq_masked']
        protein_img_latent, cell_img = self.prepare_data(batched_data, self.vae)

        t, x0, x1 = self.transport.sample(protein_img_latent)
        t, xt, ut = self.transport.path_sampler.plan(t, x0, x1)
        img_output, seq_output, img_recon = self.net(xt, cell_img, protein_seq_masked, t, self.config.img_mask_ratio)
        B, *_, C = xt.shape
        assert img_output.size() == (B, *xt.size()[1:-1], C)

        loss_dict = self.transport.training_losses(img_output, x0, xt, ut, t)

        img_recon_loss = torch.mean((img_recon - protein_img_latent)**2)
        
        model_output = (seq_output, loss_dict["loss"].mean(), img_recon_loss)
        
        loss, log_loss = self.loss(batched_data, model_output)

        return CELLFMOutput(loss=loss, log_output=log_loss)

    def sequence_to_image(self, protein_seq, cell_img, num_steps=100):
        protein_img_latent = torch.randn(cell_img.shape[0], self.config.latent_channels, self.config.sample_size, self.config.sample_size).to(cell_img.device)

        with torch.no_grad():
            sample_fn = self.transport_sampler.sample_ode(num_steps=num_steps) # default to ode sampling
            def fn(img_latent, t):
                return self.net.seq2img(img_latent, cell_img, protein_seq, t, 0)

            protein_img_latent = sample_fn(protein_img_latent, fn)[-1]

        return self.vae.decode(protein_img_latent).sample

    def image_to_sequence(self, protein_seq, protein_seq_mask, protein_img, cell_img, progress=True, sampling_strategy='oaardm', order='l2r', temperature=1.0):
        if sampling_strategy == "oaardm":
            return self.oaardm_sample(protein_seq, protein_seq_mask, protein_img, cell_img, progress, order, temperature)
        raise ValueError(f"Sampling strategy {sampling_strategy} is not supported")

    def oaardm_sample(self, protein_seq, protein_seq_mask, protein_img, cell_img, progress=True, order='l2r', temperature=1.0):
        # One sequence per call. The fill order below is a single list of positions shared by
        # the whole batch, which only means anything if there is one mask to follow. Batched,
        # torch.where on the 2-D mask returns row indices rather than positions: a batch of 2
        # would overwrite positions 0 and 1 -- the <cls> token and the first residue -- and
        # leave the masked tail untouched, decoding as '_' with no error raised.
        if protein_seq.shape[0] != 1:
            raise ValueError(
                f"oaardm_sample expects batch size 1, got {protein_seq.shape[0]}. "
                "Call it once per sequence."
            )

        sampled_seq = protein_seq.clone()

        loc = torch.where(protein_seq_mask.bool().reshape(-1))[0].cpu().numpy()

        if order == 'l2r':
            loc = np.sort(loc)
        elif order == 'random':
            np.random.shuffle(loc)

        if progress:
            # Lazy import so that we don't depend on tqdm.
            from tqdm.auto import tqdm
            loc = tqdm(loc)

        with torch.no_grad():
            protein_img_latent = self.vae.encode(protein_img).sample()
            for i in loc:
                t0, t1 = self.transport.check_interval(self.transport.train_eps, self.transport.sample_eps)
                # t1 = 1.0, t0 = 0.0
                t = torch.tensor([t1] * sampled_seq.shape[0], device=sampled_seq.device)
                seq_output = self.net(protein_img_latent, cell_img, sampled_seq, t, img_mask_ratio=0)[1]
                p = seq_output[:, i, 4:4+20] # sample at location i (random), dont let it predict non-standard AA
                p = torch.nn.functional.softmax(p / temperature, dim=1) # softmax over categorical probs
                p_sample = torch.multinomial(p, num_samples=1)
                sampled_seq[:, i] = p_sample.squeeze() + 4

        return sampled_seq

    def embed(self, protein_seq, protein_img, cell_img, img_mask_ratio=0):
        with torch.no_grad():
            protein_img_latent = self.vae.encode(protein_img).sample()

            t0, t1 = self.transport.check_interval(self.transport.train_eps, self.transport.sample_eps)
            t = torch.tensor([t1] * protein_seq.shape[0], device=protein_seq.device)
            img_feat_embed, seq_feat_embed = self.net.embed(protein_img_latent, cell_img, protein_seq, t, img_mask_ratio=img_mask_ratio)

        return img_feat_embed, seq_feat_embed

    def recon(self, protein_seq, protein_img, cell_img):
        with torch.no_grad():
            protein_img_latent = self.vae.encode(protein_img).sample()

            t0, t1 = self.transport.check_interval(self.transport.train_eps, self.transport.sample_eps)
            t = torch.tensor([t1] * protein_seq.shape[0], device=protein_seq.device)
            img_recon_latent = self.net(protein_img_latent, cell_img, protein_seq, t, img_mask_ratio=0.5)[2]

        img_recon = self.vae.decode(img_recon_latent).sample
        return img_recon


class CELLFM(nn.Module):
    def __init__(self, config: CELLFMConfig):
        super().__init__()
        self.config = config

        self.cond_conv = CondConvNet(
            config.in_channels * len(config.cell_image.split(',')), 
            config.cond_out_channels, 
        )

        # Time step embedding
        self.time_embedding = TimestepEmbedder(hidden_size=config.encoder_hidden_size)

        # Token type embeddings
        self.token_type_embeddings = nn.Embedding(2, config.encoder_hidden_size)

        # Image embedding without positional embedding
        self.img_embedding = PatchEmbed(
            height=config.sample_size, 
            width=config.sample_size, 
            patch_size=config.encoder_patch_size, 
            in_channels=config.latent_channels + config.cond_out_channels[-1], 
            embed_dim=config.encoder_hidden_size, 
            pos_embed_type=None, 
        )

        # Image positional embedding
        # Initialize pos_embed by sin-cos embedding:
        img_encoder_pos_embed = get_2d_sincos_pos_embed(
            config.encoder_hidden_size, self.img_embedding.base_size, cls_token=True, extra_tokens=1, 
        )
        # Use fixed sin-cos embedding:
        self.img_encoder_pos_embed = nn.Parameter(
            torch.from_numpy(img_encoder_pos_embed).float().unsqueeze(0), 
            requires_grad=False, 
        )

        # Sequence positional embedding
        # Initialize pos_embed by sin-cos embedding:
        protein_sequence_pos_embed = get_1d_sincos_pos_embed(
            config.encoder_hidden_size, self.config.max_protein_sequence_len+2, 
        )
        # Use fixed sin-cos embedding:
        self.protein_sequence_pos_embed = nn.Parameter(
            torch.from_numpy(protein_sequence_pos_embed).float().unsqueeze(0), 
            requires_grad=False, 
        )

        # Sequence projection
        self.seq_proj_in = nn.Linear(config.encoder_hidden_size, config.encoder_hidden_size, bias=True)

        # Unified encoder
        self.unified_encoder = nn.ModuleList([
            TransformerBlock(
                hidden_size=config.encoder_hidden_size, 
                num_heads=config.num_heads, 
                dim_head=config.dim_head, 
                dropout=config.dropout, 
                final_dropout=config.final_dropout, 
            ) for _ in range(config.encoder_num_hidden_layers)
        ])

        # Image generator
        self.img_generator = SD3Transformer2DModel(
            sample_size=config.sample_size, 
            patch_size=config.img_generator_patch_size, 
            in_channels=config.latent_channels + config.cond_out_channels[-1], 
            num_layers=config.img_generator_num_layers, 
            attention_head_dim=config.attention_head_dim, 
            num_attention_heads=config.num_attention_heads, 
            joint_attention_dim=config.encoder_hidden_size, 
            pooled_projection_dim=config.encoder_hidden_size, 
            caption_projection_dim=config.num_attention_heads * config.attention_head_dim, 
            out_channels=config.latent_channels, 
        )

        # Image decoder
        self.img_decoder_embed = nn.Linear(config.encoder_hidden_size, config.img_decoder_hidden_size, bias=True)
        self.mask_token = nn.Parameter(torch.zeros(1, 1, config.img_decoder_hidden_size))

        # Initialize pos_embed by sin-cos embedding:
        img_decoder_pos_embed = get_2d_sincos_pos_embed(
            config.img_decoder_hidden_size, self.img_embedding.base_size, cls_token=True, extra_tokens=1, 
        )
        # Use fixed sin-cos embedding:
        self.img_decoder_pos_embed = nn.Parameter(
            torch.from_numpy(img_decoder_pos_embed).float().unsqueeze(0), 
            requires_grad=False, 
        )

        self.img_decoder = nn.ModuleList([
            TransformerBlock(
                hidden_size=config.img_decoder_hidden_size, 
                num_heads=config.img_decoder_num_heads, 
                dim_head=config.img_decoder_dim_head, 
                dropout=config.img_decoder_dropout, 
                final_dropout=config.img_decoder_final_dropout, 
            ) for _ in range(config.img_decoder_num_hidden_layers)
        ])
        self.img_decoder_norm = nn.LayerNorm(config.img_decoder_hidden_size)
        self.img_decoder_pred = nn.Linear(config.img_decoder_hidden_size, config.encoder_patch_size**2 * config.latent_channels, bias=True) # decoder to patch

        # Sequence decoder
        self.seq_decoder = MLMHead(
            hidden_size=config.encoder_hidden_size, 
            vocab_size=config.vocab_size, 
        )

        self.initialize_weights()
        self.initialize_protein_sequence_embedding()

    def initialize_weights(self):
        # initialization
        torch.nn.init.normal_(self.mask_token, std=.02)

        # initialize nn.Linear and nn.LayerNorm
        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            # we use xavier_uniform following official JAX ViT:
            torch.nn.init.xavier_uniform_(m.weight)
            if isinstance(m, nn.Linear) and m.bias is not None:
                nn.init.constant_(m.bias, 0)
        elif isinstance(m, nn.LayerNorm):
            if m.bias is not None:
                nn.init.constant_(m.bias, 0)
            if m.weight is not None:
                nn.init.constant_(m.weight, 1.0)

    def initialize_protein_sequence_embedding(self):
        self.protein_sequence_embedding = ESMEmbed(self.config.esm_embedding, self.config.encoder_hidden_size, self.config.esm_fixed_embedding)

    def forward(self, protein_img, cell_img, protein_seq, time, img_mask_ratio=0.5):
        time_embeds = self.time_embedding(time)

        # Tokenize protein sequence
        # Size: B x T x C
        seq_embeds = self.protein_sequence_embedding(protein_seq)
        seq_embeds = self.seq_proj_in(seq_embeds)
        seq_embeds = seq_embeds + self.protein_sequence_pos_embed[:, :protein_seq.shape[1]]

        seq_token_type = torch.full(
            size=(seq_embeds.shape[0], 1), fill_value=1, device=seq_embeds.device
        ).long()

        seq_embeds = seq_embeds + self.token_type_embeddings(seq_token_type)

        cell_img_conv = self.cond_conv(cell_img)
        concat_img = torch.cat([protein_img, cell_img_conv], dim=1)
        img_embeds = self.img_embedding(concat_img)

        # add pos embed w/o time step token
        img_embeds = img_embeds + self.img_encoder_pos_embed[:, 1:, :]

        # masking: length -> length * mask_ratio
        img_embeds, img_mask, img_ids_restore = self.random_masking(img_embeds, img_mask_ratio)

        # append time step token
        time_embeds_with_pos_embed = time_embeds.unsqueeze(1) + self.img_encoder_pos_embed[:, :1, :]
        img_embeds = torch.cat([time_embeds_with_pos_embed, img_embeds], dim=1)        

        img_token_type = torch.full(
            size=(img_embeds.shape[0], 1), fill_value=0, device=img_embeds.device, 
        ).long()
        
        img_embeds = img_embeds + self.token_type_embeddings(img_token_type)
        co_embeds = torch.cat([img_embeds, seq_embeds], dim=1)

        x = co_embeds

        for i, block in enumerate(self.unified_encoder):
            x = block(x)

        co_feat = x

        img_output = self.img_generator(
            hidden_states=concat_img, timestep=time, encoder_hidden_states=co_feat, 
            pooled_projections=time_embeds, 
        ).sample

        seq_feat = co_feat[:, img_embeds.shape[1]:]
        img_feat = co_feat[:, :img_embeds.shape[1]]

        seq_output = self.seq_decoder(seq_feat)

        img_recon = self.forward_img_decoder(img_feat, img_ids_restore)
        img_recon = self.unpatchify(img_recon, self.config.latent_channels)

        return img_output, seq_output, img_recon

    def random_masking(self, x, mask_ratio):
        """
        Perform per-sample random masking by per-sample shuffling.
        Per-sample shuffling is done by argsort random noise.
        x: [N, L, D], sequence
        """
        N, L, D = x.shape  # batch, length, dim
        len_keep = int(L * (1 - mask_ratio))
        
        noise = torch.rand(N, L, device=x.device)  # noise in [0, 1]
        
        # sort noise for each sample
        ids_shuffle = torch.argsort(noise, dim=1)  # ascend: small is keep, large is remove
        ids_restore = torch.argsort(ids_shuffle, dim=1)

        # keep the first subset
        ids_keep = ids_shuffle[:, :len_keep]
        x_masked = torch.gather(x, dim=1, index=ids_keep.unsqueeze(-1).repeat(1, 1, D))

        # generate the binary mask: 0 is keep, 1 is remove
        mask = torch.ones([N, L], device=x.device)
        mask[:, :len_keep] = 0
        # unshuffle to get the binary mask
        mask = torch.gather(mask, dim=1, index=ids_restore)

        return x_masked, mask, ids_restore

    def forward_img_decoder(self, x, ids_restore):
        # embed tokens
        x = self.img_decoder_embed(x)

        # append mask tokens to sequence
        mask_tokens = self.mask_token.repeat(x.shape[0], ids_restore.shape[1] + 1 - x.shape[1], 1)
        x_ = torch.cat([x[:, 1:, :], mask_tokens], dim=1)  # no time step token
        x_ = torch.gather(x_, dim=1, index=ids_restore.unsqueeze(-1).repeat(1, 1, x.shape[2]))  # unshuffle
        x = torch.cat([x[:, :1, :], x_], dim=1)  # append time step token

        # add pos embed
        x = x + self.img_decoder_pos_embed

        # apply Transformer blocks
        for block in self.img_decoder:
            x = block(x)
        x = self.img_decoder_norm(x)

        # predictor projection
        x = self.img_decoder_pred(x)

        # remove time step token
        x = x[:, 1:, :]

        return x

    def unpatchify(self, x, num_channels):
        """
        x: (N, L, patch_size**2 *C)
        imgs: (N, C, H, W)
        """
        p = self.config.encoder_patch_size
        h = w = int(x.shape[1]**.5)
        assert h * w == x.shape[1]
        
        x = x.reshape(shape=(x.shape[0], h, w, p, p, num_channels))
        x = torch.einsum('nhwpqc->nchpwq', x)
        imgs = x.reshape(shape=(x.shape[0], num_channels, h * p, h * p))
        return imgs

    def embed(self, protein_img, cell_img, protein_seq, time, img_mask_ratio=0):
        time_embeds = self.time_embedding(time)

        # Tokenize protein sequence
        # Size: B x T x C
        seq_embeds = self.protein_sequence_embedding(protein_seq)
        seq_embeds = self.seq_proj_in(seq_embeds)
        seq_embeds = seq_embeds + self.protein_sequence_pos_embed[:, :protein_seq.shape[1]]

        seq_token_type = torch.full(
            size=(seq_embeds.shape[0], 1), fill_value=1, device=seq_embeds.device
        ).long()

        seq_embeds = seq_embeds + self.token_type_embeddings(seq_token_type)

        cell_img_conv = self.cond_conv(cell_img)
        concat_img = torch.cat([protein_img, cell_img_conv], dim=1)
        img_embeds = self.img_embedding(concat_img)

        # add pos embed w/o time step token
        img_embeds = img_embeds + self.img_encoder_pos_embed[:, 1:, :]

        # masking: length -> length * mask_ratio
        img_embeds, img_mask, img_ids_restore = self.random_masking(img_embeds, img_mask_ratio)

        # append time step token
        time_embeds_with_pos_embed = time_embeds.unsqueeze(1) + self.img_encoder_pos_embed[:, :1, :]
        img_embeds = torch.cat([time_embeds_with_pos_embed, img_embeds], dim=1)        

        img_token_type = torch.full(
            size=(img_embeds.shape[0], 1), fill_value=0, device=img_embeds.device, 
        ).long()
        
        img_embeds = img_embeds + self.token_type_embeddings(img_token_type)
        co_embeds = torch.cat([img_embeds, seq_embeds], dim=1)

        x = co_embeds

        for i, block in enumerate(self.unified_encoder):
            x = block(x)

        co_feat = x

        seq_feat = co_feat[:, img_embeds.shape[1]:]
        img_feat = co_feat[:, :img_embeds.shape[1]]

        img_feat_embed = img_feat
        seq_feat_embed = seq_feat

        return img_feat_embed, seq_feat_embed

    def seq2img(self, protein_img, cell_img, protein_seq, time, img_mask_ratio=0):
        time_embeds = self.time_embedding(time)

        # Tokenize protein sequence
        # Size: B x T x C
        seq_embeds = self.protein_sequence_embedding(protein_seq)
        seq_embeds = self.seq_proj_in(seq_embeds)
        seq_embeds = seq_embeds + self.protein_sequence_pos_embed[:, :protein_seq.shape[1]]

        seq_token_type = torch.full(
            size=(seq_embeds.shape[0], 1), fill_value=1, device=seq_embeds.device
        ).long()

        seq_embeds = seq_embeds + self.token_type_embeddings(seq_token_type)

        cell_img_conv = self.cond_conv(cell_img)
        concat_img = torch.cat([protein_img, cell_img_conv], dim=1)
        img_embeds = self.img_embedding(concat_img)

        # add pos embed w/o time step token
        img_embeds = img_embeds + self.img_encoder_pos_embed[:, 1:, :]

        # masking: length -> length * mask_ratio
        img_embeds, img_mask, img_ids_restore = self.random_masking(img_embeds, img_mask_ratio)

        # append time step token
        time_embeds_with_pos_embed = time_embeds.unsqueeze(1) + self.img_encoder_pos_embed[:, :1, :]
        img_embeds = torch.cat([time_embeds_with_pos_embed, img_embeds], dim=1)        

        img_token_type = torch.full(
            size=(img_embeds.shape[0], 1), fill_value=0, device=img_embeds.device, 
        ).long()
        
        img_embeds = img_embeds + self.token_type_embeddings(img_token_type)
        co_embeds = torch.cat([img_embeds, seq_embeds], dim=1)

        x = co_embeds

        for i, block in enumerate(self.unified_encoder):
            x = block(x)

        co_feat = x

        img_output = self.img_generator(
            hidden_states=concat_img, timestep=time, encoder_hidden_states=co_feat, 
            pooled_projections=time_embeds, 
        ).sample

        return img_output
