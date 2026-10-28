import torch
import torch.nn as nn
import torch.nn.functional as F

from diffusers.models.transformers import SD3Transformer2DModel

from .patch_embed_3d import PatchEmbed3D


class MaskedJointAttnProcessor:
    """`JointAttnProcessor2_0` that honours a key padding mask on the context stream.

    Upstream SD3 cannot mask: `JointAttnProcessor2_0` accepts `attention_mask` and then calls
    `scaled_dot_product_attention` without it, and `JointTransformerBlock.forward` has no way to
    pass one down. Sequence padding therefore acts as real conditioning context. This processor
    reads the current mask off `owner` -- rewritten on every `SD3Transformer3DModel.forward`, so
    it can never be stale -- and masks the padded keys/values out of joint attention.

    Only keys are masked. Padded *queries* still attend to the image tokens, so no attention row
    is fully masked (which would produce NaNs); their outputs land in `encoder_hidden_states`,
    which this model discards.
    """

    def __init__(self, owner: nn.Module):
        # Plain attribute, not a submodule: nothing here enters the state dict.
        self.owner = owner

    def __call__(
        self,
        attn,
        hidden_states: torch.FloatTensor,
        encoder_hidden_states: torch.FloatTensor = None,
        attention_mask=None,
        *args,
        **kwargs,
    ) -> torch.FloatTensor:
        residual = hidden_states
        batch_size = hidden_states.shape[0]

        query = attn.to_q(hidden_states)
        key = attn.to_k(hidden_states)
        value = attn.to_v(hidden_states)

        inner_dim = key.shape[-1]
        head_dim = inner_dim // attn.heads

        query = query.view(batch_size, -1, attn.heads, head_dim).transpose(1, 2)
        key = key.view(batch_size, -1, attn.heads, head_dim).transpose(1, 2)
        value = value.view(batch_size, -1, attn.heads, head_dim).transpose(1, 2)

        if attn.norm_q is not None:
            query = attn.norm_q(query)
        if attn.norm_k is not None:
            key = attn.norm_k(key)

        attn_mask = None
        if encoder_hidden_states is not None:
            encoder_query = attn.add_q_proj(encoder_hidden_states)
            encoder_key = attn.add_k_proj(encoder_hidden_states)
            encoder_value = attn.add_v_proj(encoder_hidden_states)

            encoder_query = encoder_query.view(batch_size, -1, attn.heads, head_dim).transpose(1, 2)
            encoder_key = encoder_key.view(batch_size, -1, attn.heads, head_dim).transpose(1, 2)
            encoder_value = encoder_value.view(batch_size, -1, attn.heads, head_dim).transpose(1, 2)

            if attn.norm_added_q is not None:
                encoder_query = attn.norm_added_q(encoder_query)
            if attn.norm_added_k is not None:
                encoder_key = attn.norm_added_k(encoder_key)

            query = torch.cat([query, encoder_query], dim=2)
            key = torch.cat([key, encoder_key], dim=2)
            value = torch.cat([value, encoder_value], dim=2)

            # Keys run [image tokens ; sequence tokens]; image tokens are always kept.
            seq_mask = self.owner.encoder_attention_mask
            if seq_mask is not None:
                img_keep = seq_mask.new_ones((batch_size, residual.shape[1]))
                attn_mask = torch.cat([img_keep, seq_mask], dim=1)[:, None, None, :]

        hidden_states = F.scaled_dot_product_attention(
            query, key, value, attn_mask=attn_mask, dropout_p=0.0, is_causal=False
        )
        hidden_states = hidden_states.transpose(1, 2).reshape(batch_size, -1, attn.heads * head_dim)
        hidden_states = hidden_states.to(query.dtype)

        if encoder_hidden_states is None:
            return attn.to_out[1](attn.to_out[0](hidden_states))

        hidden_states, encoder_hidden_states = (
            hidden_states[:, : residual.shape[1]],
            hidden_states[:, residual.shape[1]:],
        )
        if not attn.context_pre_only:
            encoder_hidden_states = attn.to_add_out(encoder_hidden_states)

        hidden_states = attn.to_out[0](hidden_states)
        hidden_states = attn.to_out[1](hidden_states)

        return hidden_states, encoder_hidden_states


class SD3Transformer3DModel(nn.Module):
    """SD3 transformer adapted for 3-D volumetric latents.

    Builds an SD3Transformer2DModel to obtain its transformer blocks and
    conditioning modules, then replaces the 2-D patch embedding with a 3-D
    variant and handles 3-D unpatchification in forward().
    """

    def __init__(
        self,
        depth: int,
        height: int,
        width: int,
        patch_d: int,
        patch_hw: int,
        latent_channels: int,
        num_layers: int,
        attention_head_dim: int,
        num_attention_heads: int,
        joint_attention_dim: int,
        pooled_projection_dim: int,
        in_channels: int = None,
        qk_norm: str = 'rms_norm',
    ):
        super().__init__()
        if in_channels is None:
            in_channels = latent_channels
        inner_dim = num_attention_heads * attention_head_dim

        # Query/key normalization inside joint attention (SD3.5 default). MMDiT attention sets
        # added_kv_proj_dim for the context stream, and diffusers only builds norm_added_q/k for
        # 'rms_norm' and 'fp32_layer_norm' -- 'layer_norm' raises there, so reject it up front.
        if qk_norm in (None, '', 'none', 'None'):
            qk_norm = None
        elif qk_norm not in ('rms_norm', 'fp32_layer_norm'):
            raise ValueError(
                f"Unknown qk_norm '{qk_norm}'; expected one of ['fp32_layer_norm', 'none', 'rms_norm']"
            )

        # Build SD3 with a minimal 2-D config to instantiate its transformer blocks.
        _sd3 = SD3Transformer2DModel(
            sample_size=4,
            patch_size=4,
            in_channels=latent_channels,
            out_channels=latent_channels,
            num_layers=num_layers,
            attention_head_dim=attention_head_dim,
            num_attention_heads=num_attention_heads,
            joint_attention_dim=joint_attention_dim,
            pooled_projection_dim=pooled_projection_dim,
            caption_projection_dim=inner_dim,
            pos_embed_max_size=4,
            qk_norm=qk_norm,
        )

        # Extract and register reusable components; frees the dummy pos_embed.
        self.time_text_embed = _sd3.time_text_embed
        self.context_embedder = _sd3.context_embedder
        self.transformer_blocks = _sd3.transformer_blocks
        self.norm_out = _sd3.norm_out
        del _sd3

        # 3-D patch embedding and output projection
        self.pos_embed = PatchEmbed3D(
            depth, height, width, patch_d, patch_hw, in_channels, inner_dim
        )
        self.proj_out = nn.Linear(
            inner_dim, patch_d * patch_hw * patch_hw * latent_channels
        )

        nd, nh, nw = depth // patch_d, height // patch_hw, width // patch_hw
        self._patch_d = patch_d
        self._patch_hw = patch_hw
        self._nd, self._nh, self._nw = nd, nh, nw
        self._out_channels = latent_channels

        # Swap in the masking processor so sequence padding stays out of joint attention.
        self.encoder_attention_mask = None
        for block in self.transformer_blocks:
            block.attn.set_processor(MaskedJointAttnProcessor(self))

    def forward(
        self,
        hidden_states: torch.Tensor,
        encoder_hidden_states: torch.Tensor,
        temb: torch.Tensor,
        encoder_attention_mask: torch.Tensor = None,
    ) -> torch.Tensor:
        # temb is self.time_text_embed(timestep, pooled_projections); the caller computes it so the
        # ResBlocks around this transformer can share it.
        # Published to the attention processors for this forward; always rewritten, never stale.
        self.encoder_attention_mask = (
            encoder_attention_mask.bool() if encoder_attention_mask is not None else None
        )

        # hidden_states: (B, C, D, H, W)
        hidden_states = self.pos_embed(hidden_states)                            # (B, N, D_model)
        encoder_hidden_states = self.context_embedder(encoder_hidden_states)

        for block in self.transformer_blocks:
            encoder_hidden_states, hidden_states = block(
                hidden_states, encoder_hidden_states, temb
            )

        hidden_states = self.norm_out(hidden_states, temb)
        hidden_states = self.proj_out(hidden_states)                             # (B, N, pd*ph*pw*C)

        # 3-D unpatchify
        pd, ph, pw = self._patch_d, self._patch_hw, self._patch_hw
        nd, nh, nw = self._nd, self._nh, self._nw
        B, C = hidden_states.shape[0], self._out_channels

        hidden_states = hidden_states.reshape(B, nd, nh, nw, pd, ph, pw, C)
        hidden_states = hidden_states.permute(0, 7, 1, 4, 2, 5, 3, 6)          # (B, C, nd, pd, nh, ph, nw, pw)
        return hidden_states.reshape(B, C, nd * pd, nh * ph, nw * pw)           # (B, C, D, H, W)
