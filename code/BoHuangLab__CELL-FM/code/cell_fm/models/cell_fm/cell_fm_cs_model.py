import numpy as np
import torch
import torch.nn as nn

from cell_fm.models.cell_fm.cell_fm_model import CELLFM as BaseCELLFM
from cell_fm.models.cell_fm.cell_fm_model import CELLFMModel as BaseCELLFMModel
from cell_fm.models.cell_fm.modules.positional_embedding import get_2d_sincos_pos_embed
from cell_fm.pipeline.utils import CELLFMOutput

from .cell_fm_cs_config import CELLFMConfig


class CELLFMCSModel(BaseCELLFMModel):
    config_class = CELLFMConfig

    def initialize_net(self, config: CELLFMConfig):
        return CELLFMCS(config)

    def forward(self, batched_data, **kwargs):
        protein_seq_masked = batched_data["protein_seq_masked"]
        protein_img_latent, cell_img = self.prepare_data(batched_data, self.vae)
        protein_intensity_level = batched_data["protein_intensity_level"]

        t, x0, x1 = self.transport.sample(protein_img_latent)
        t, xt, ut = self.transport.path_sampler.plan(t, x0, x1)
        img_output, seq_output, img_recon = self.net(
            xt,
            cell_img,
            protein_seq_masked,
            protein_intensity_level,
            t,
            self.config.img_mask_ratio,
        )
        B, *_, C = xt.shape
        assert img_output.size() == (B, *xt.size()[1:-1], C)

        loss_dict = self.transport.training_losses(img_output, x0, xt, ut, t)
        img_recon_loss = torch.mean((img_recon - protein_img_latent) ** 2)
        model_output = (seq_output, loss_dict["loss"].mean(), img_recon_loss)
        loss, log_loss = self.loss(batched_data, model_output)

        return CELLFMOutput(loss=loss, log_output=log_loss)

    def sequence_to_image(self, protein_seq, cell_img, protein_intensity_level, num_steps=100):
        protein_img_latent = torch.randn(
            cell_img.shape[0],
            self.config.latent_channels,
            self.config.sample_size,
            self.config.sample_size,
            device=cell_img.device,
        )

        with torch.no_grad():
            sample_fn = self.transport_sampler.sample_ode(num_steps=num_steps)

            def fn(img_latent, t):
                return self.net.inference(
                    img_latent,
                    cell_img,
                    protein_seq,
                    protein_intensity_level,
                    t,
                    0,
                )[0]

            protein_img_latent = sample_fn(protein_img_latent, fn)[-1]

        return self.vae.decode(protein_img_latent).sample

    def image_to_sequence(
        self,
        protein_seq,
        protein_seq_mask,
        protein_img,
        cell_img,
        protein_intensity_level,
        progress=True,
        sampling_strategy="oaardm",
        order="l2r",
        temperature=1.0,
    ):
        if sampling_strategy == "oaardm":
            return self.oaardm_sample(
                protein_seq,
                protein_seq_mask,
                protein_img,
                cell_img,
                protein_intensity_level,
                progress,
                order,
                temperature,
            )
        raise ValueError(f"Sampling strategy {sampling_strategy} is not supported")

    def oaardm_sample(
        self,
        protein_seq,
        protein_seq_mask,
        protein_img,
        cell_img,
        protein_intensity_level,
        progress=True,
        order="l2r",
        temperature=1.0,
    ):
        from esm.utils.constants.esm3 import SEQUENCE_VOCAB

        device = protein_seq.device
        vocab_to_idx = {aa: idx for idx, aa in enumerate(SEQUENCE_VOCAB)}
        standard_aas = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        allowed_token_ids = torch.tensor(
            [vocab_to_idx[aa] for aa in standard_aas if aa in vocab_to_idx],
            dtype=torch.long,
            device=device,
        )

        # One sequence per call -- see the note in cell_fm_model.oaardm_sample. Batched,
        # torch.where on the 2-D mask returns row indices rather than positions, and the
        # sampler silently rewrites the wrong tokens.
        if protein_seq.shape[0] != 1:
            raise ValueError(
                f"oaardm_sample expects batch size 1, got {protein_seq.shape[0]}. "
                "Call it once per sequence."
            )

        # Fill a copy, not the caller's tensor. Writing through protein_seq leaves it with
        # no <mask> tokens left, so a caller looping over num_gen with the same prompt
        # resamples the previous draw instead of starting over and its draws come out as a
        # chain -- near-duplicates, and correlated samples under any test applied to them.
        sampled_seq = protein_seq.clone()

        loc = torch.where(protein_seq_mask.bool().reshape(-1))[0].cpu().numpy()
        if order == "l2r":
            loc = np.sort(loc)
        elif order == "random":
            np.random.shuffle(loc)

        if progress:
            from tqdm.auto import tqdm

            loc = tqdm(loc)

        with torch.no_grad():
            protein_img_latent = self.vae.encode(protein_img).sample()
            batch_size = protein_seq.shape[0]

            for i in loc:
                _, t1 = self.transport.check_interval(
                    self.transport.train_eps,
                    self.transport.sample_eps,
                )
                t = torch.full((batch_size,), t1, device=device)
                seq_output = self.net.inference(
                    protein_img_latent,
                    cell_img,
                    sampled_seq,
                    protein_intensity_level,
                    t,
                )[1]
                logits_i = seq_output[:, i, :]
                aa_logits = logits_i[:, allowed_token_ids]
                probs = torch.nn.functional.softmax(aa_logits / temperature, dim=-1)
                sample_idx = torch.multinomial(probs, num_samples=1).squeeze(-1)
                sampled_seq[:, i] = allowed_token_ids[sample_idx]

        return sampled_seq

    def recon(self, protein_seq, protein_img, cell_img, protein_intensity_level):
        with torch.no_grad():
            protein_img_latent = self.vae.encode(protein_img).sample()
            _, t1 = self.transport.check_interval(self.transport.train_eps, self.transport.sample_eps)
            t = torch.tensor([t1] * protein_seq.shape[0], device=protein_seq.device)
            img_recon_latent = self.net(
                protein_img_latent,
                cell_img,
                protein_seq,
                protein_intensity_level,
                t,
                img_mask_ratio=0.5,
            )[2]

        return self.vae.decode(img_recon_latent).sample


class CELLFMCS(BaseCELLFM):
    def __init__(self, config: CELLFMConfig):
        super().__init__(config)

        self.protein_intensity_level_embedding = nn.Linear(
            in_features=1,
            out_features=config.encoder_hidden_size,
            bias=True,
        )
        self.time_step_pos_embed = nn.Parameter(
            torch.randn(1, 1, config.encoder_hidden_size) * 0.02,
            requires_grad=True,
        )
        self.protein_intensity_level_pos_embed = nn.Parameter(
            torch.randn(1, 1, config.encoder_hidden_size) * 0.02,
            requires_grad=True,
        )

        img_encoder_pos_embed = get_2d_sincos_pos_embed(
            config.encoder_hidden_size,
            self.img_embedding.base_size,
        )
        self.img_encoder_pos_embed = nn.Parameter(
            torch.from_numpy(img_encoder_pos_embed).float().unsqueeze(0),
            requires_grad=False,
        )

        img_decoder_pos_embed = get_2d_sincos_pos_embed(
            config.img_decoder_hidden_size,
            self.img_embedding.base_size,
        )
        self.img_decoder_pos_embed = nn.Parameter(
            torch.from_numpy(img_decoder_pos_embed).float().unsqueeze(0),
            requires_grad=False,
        )

        self._init_weights(self.protein_intensity_level_embedding)

    def _encode_conditioned(
        self,
        protein_img,
        cell_img,
        protein_seq,
        protein_intensity_level,
        time,
        img_mask_ratio,
    ):
        time_embeds = self.time_embedding(time)
        intensity_embeds = self.protein_intensity_level_embedding(protein_intensity_level)

        seq_embeds = self.protein_sequence_embedding(protein_seq)
        seq_embeds = self.seq_proj_in(seq_embeds)
        seq_embeds = seq_embeds + self.protein_sequence_pos_embed[:, : protein_seq.shape[1]]
        seq_token_type = torch.full(
            size=(seq_embeds.shape[0], 1),
            fill_value=1,
            device=seq_embeds.device,
        ).long()
        seq_embeds = seq_embeds + self.token_type_embeddings(seq_token_type)

        cell_img_conv = self.cond_conv(cell_img)
        concat_img = torch.cat([protein_img, cell_img_conv], dim=1)
        img_embeds = self.img_embedding(concat_img)
        img_embeds = img_embeds + self.img_encoder_pos_embed
        img_embeds, _, img_ids_restore = self.random_masking(img_embeds, img_mask_ratio)

        time_embeds_with_pos = time_embeds.unsqueeze(1) + self.time_step_pos_embed
        intensity_embeds_with_pos = intensity_embeds.unsqueeze(1) + self.protein_intensity_level_pos_embed
        img_embeds = torch.cat([time_embeds_with_pos, intensity_embeds_with_pos, img_embeds], dim=1)
        img_token_type = torch.full(
            size=(img_embeds.shape[0], 1),
            fill_value=0,
            device=img_embeds.device,
        ).long()
        img_embeds = img_embeds + self.token_type_embeddings(img_token_type)

        co_embeds = torch.cat([img_embeds, seq_embeds], dim=1)
        co_feat = co_embeds
        for block in self.unified_encoder:
            co_feat = block(co_feat)

        return concat_img, img_embeds, img_ids_restore, co_feat, time_embeds

    def forward(self, protein_img, cell_img, protein_seq, protein_intensity_level, time, img_mask_ratio=0.5):
        concat_img, img_embeds, img_ids_restore, co_feat, time_embeds = self._encode_conditioned(
            protein_img,
            cell_img,
            protein_seq,
            protein_intensity_level,
            time,
            img_mask_ratio,
        )
        img_output = self.img_generator(
            hidden_states=concat_img,
            timestep=time,
            encoder_hidden_states=co_feat,
            pooled_projections=time_embeds,
        ).sample
        seq_feat = co_feat[:, img_embeds.shape[1] :]
        img_feat = co_feat[:, : img_embeds.shape[1]]
        seq_output = self.seq_decoder(seq_feat)
        img_recon = self.forward_img_decoder(img_feat, img_ids_restore)
        img_recon = self.unpatchify(img_recon, self.config.latent_channels)

        return img_output, seq_output, img_recon

    def forward_img_decoder(self, x, ids_restore):
        x = self.img_decoder_embed(x)
        mask_tokens = self.mask_token.repeat(x.shape[0], ids_restore.shape[1] + 2 - x.shape[1], 1)
        x_ = torch.cat([x[:, 2:, :], mask_tokens], dim=1)
        x_ = torch.gather(x_, dim=1, index=ids_restore.unsqueeze(-1).repeat(1, 1, x.shape[2]))
        x_ = x_ + self.img_decoder_pos_embed
        x = torch.cat([x[:, :2, :], x_], dim=1)

        for block in self.img_decoder:
            x = block(x)
        x = self.img_decoder_norm(x)
        x = self.img_decoder_pred(x)

        return x[:, 2:, :]

    def inference(self, protein_img, cell_img, protein_seq, protein_intensity_level, time, img_mask_ratio=0):
        concat_img, img_embeds, _, co_feat, time_embeds = self._encode_conditioned(
            protein_img,
            cell_img,
            protein_seq,
            protein_intensity_level,
            time,
            img_mask_ratio,
        )
        img_output = self.img_generator(
            hidden_states=concat_img,
            timestep=time,
            encoder_hidden_states=co_feat,
            pooled_projections=time_embeds,
        ).sample
        seq_feat = co_feat[:, img_embeds.shape[1] :]
        seq_output = self.seq_decoder(seq_feat)

        return img_output, seq_output
