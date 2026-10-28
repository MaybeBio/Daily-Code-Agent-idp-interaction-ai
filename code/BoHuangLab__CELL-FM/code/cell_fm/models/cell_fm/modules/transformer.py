import torch
import torch.nn as nn
from diffusers.models.attention_processor import Attention
from diffusers.models.attention import FeedForward
from transformers.activations import ACT2FN

class BertPredictionHeadTransform(nn.Module):
    def __init__(
            self, 
            hidden_size, 
            hidden_act: str = 'gelu', 
            norm_eps: float = 1e-5,
        ):
        super().__init__()
        self.dense = nn.Linear(hidden_size, hidden_size)
        self.transform_act_fn = ACT2FN[hidden_act]
        self.LayerNorm = nn.LayerNorm(hidden_size, norm_eps)

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        hidden_states = self.dense(hidden_states)
        hidden_states = self.transform_act_fn(hidden_states)
        hidden_states = self.LayerNorm(hidden_states)
        return hidden_states

class MLMHead(nn.Module):
    def __init__(
            self,
            hidden_size, 
            vocab_size, 
            hidden_act: str = 'gelu', 
            norm_eps: float = 1e-5, 
            weight=None
        ):
        super().__init__()
        self.transform = BertPredictionHeadTransform(hidden_size, hidden_act, norm_eps)
        self.decoder = nn.Linear(hidden_size, vocab_size, bias=False)
        self.bias = nn.Parameter(torch.zeros(vocab_size))
        if weight is not None:
            self.decoder.weight = weight

    def forward(self, x):
        x = self.transform(x)
        x = self.decoder(x) + self.bias
        return x

class TransformerBlock(nn.Module):
    def __init__(
            self, 
            hidden_size, 
            num_heads, 
            dim_head, 
            dropout: float = 0.0, 
            final_dropout: float = 0.0, 
            norm_elementwise_affine: bool = True, 
            norm_eps: float = 1e-5,
            attention_bias: bool = False, 
            attention_out_bias: bool = True, 
            activation_fn: str = "geglu",
            ff_inner_dim = None, 
            ff_bias: bool = True,
        ):
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden_size, elementwise_affine=norm_elementwise_affine, eps=norm_eps)

        self.attn = Attention(
            query_dim=hidden_size, 
            heads=num_heads, 
            dim_head=dim_head, 
            dropout=dropout, 
            bias=attention_bias, 
            cross_attention_dim=None, 
            upcast_attention=False, 
            out_bias=attention_out_bias, 
        )

        self.ff = FeedForward(
            hidden_size, 
            dropout=dropout, 
            activation_fn=activation_fn, 
            final_dropout=final_dropout, 
            inner_dim=ff_inner_dim, 
            bias=ff_bias, 
        )

        self.norm2 = nn.LayerNorm(hidden_size, norm_eps, norm_elementwise_affine)

    def forward(self, hidden_states):
        norm_hidden_states = self.norm1(hidden_states)
        attn_output = self.attn(norm_hidden_states)

        hidden_states = attn_output + hidden_states
        norm_hidden_states = self.norm2(hidden_states)

        ff_output = self.ff(norm_hidden_states)
        hidden_states = ff_output + hidden_states

        return hidden_states

def unpatchify(x, out_channels, patch_size):
    """
    x: (N, T, patch_size**2 * C)
    imgs: (N, C, H, W)
    """
    c = out_channels
    p = patch_size
    h = w = int(x.shape[1] ** 0.5)
    assert h * w == x.shape[1]

    x = x.reshape(shape=(x.shape[0], h, w, p, p, c))
    x = torch.einsum('nhwpqc->nchpwq', x)
    imgs = x.reshape(shape=(x.shape[0], c, h * p, h * p))
    return imgs

class FinalLayer(nn.Module):
    """
    The final layer of DiT.
    """
    def __init__(self, hidden_size, patch_size, out_channels):
        super().__init__()
        self.norm_final = nn.LayerNorm(hidden_size, elementwise_affine=False, eps=1e-6)
        self.linear = nn.Linear(hidden_size, patch_size * patch_size * out_channels, bias=True)
        self.adaLN_modulation = nn.Sequential(
            nn.SiLU(),
            nn.Linear(hidden_size, 2 * hidden_size, bias=True)
        )

    def forward(self, x, c):
        shift, scale = self.adaLN_modulation(c).chunk(2, dim=1)
        x = self.linear(x)
        return x

