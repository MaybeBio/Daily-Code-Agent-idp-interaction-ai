"""
Building blocks of FlexBind (see Section "Method" of the paper).

    ProtT5 embedding H (frozen, precomputed)
        │
        ├─► MultiScaleComplexRepresentation          (shared across tasks)
        │       complex-Morlet CWT  →  Re / Im streams  →  learnable 1x1 conv (d1 → d2)
        │
        └─► TaskBranch  ×4 (IDR / protein / RNA / DNA, independently parameterised)
                ResidueWiseScaleAttention (RWSA)     — real & imaginary streams
                ComplexAxisAlignment (CAA)           — learnable complex scalar (a_k + i b_k)
                DualStreamIntegration                — W_F [axis1 || axis2], H + W_O C̄
                MambaBlock × M                       — padding-aware selective SSM
                MLP head                             — residue-level logit
"""
import math
from typing import Dict, Optional, Sequence, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.checkpoint import checkpoint

try:  # Optional fused CUDA kernel from the official Mamba implementation.
    from mamba_ssm.ops.selective_scan_interface import selective_scan_fn as _mamba_selective_scan_fn
    HAS_MAMBA_SSM = True
except Exception:  # pragma: no cover - depends on the environment
    _mamba_selective_scan_fn = None
    HAS_MAMBA_SSM = False


def _mask_like(mask: Optional[torch.Tensor], x: torch.Tensor) -> torch.Tensor:
    """Returns a (B, L) float mask; all-ones if mask is None. x is (B, L, ...)."""
    if mask is None:
        return torch.ones(x.shape[0], x.shape[1], device=x.device, dtype=x.dtype)
    return mask.to(dtype=x.dtype)


# --------------------------------------------------------------------------- #
# 1. Continuous wavelet transform with a complex Morlet mother wavelet
# --------------------------------------------------------------------------- #
class MorletCWT(nn.Module):
    """Fixed (non-learnable) continuous wavelet transform along the sequence axis.

    Mother wavelet (complex Morlet, with the admissibility correction term):
        psi(x) = pi^{-1/4} (exp(i w0 x) - exp(-w0^2 / 2)) exp(-x^2 / 2)

    For a scale s the daughter wavelet psi(t / s) / sqrt(s) is sampled on the integer grid
    t in [-floor(K s / 2), floor(K s / 2)], i.e. its support covers ~K*s residues
    (K = kernel_size = 5 → scales {1,...,32} span ~5–160 residues, as in the paper).

    W(s, b) = sum_t x(t) conj(psi((t - b) / s)) / sqrt(s)  — implemented as a depth-wise
    cross-correlation (F.conv1d) with the conjugated kernel.

    Padded positions are zeroed before the transform (so they contribute nothing) and the
    coefficients at padded positions are zeroed afterwards.
    """

    def __init__(self, scales: Sequence[int], kernel_size: int = 5, omega0: float = 6.0):
        super().__init__()
        self.scales = [int(s) for s in scales]
        self.kernel_size = int(kernel_size)
        self.omega0 = float(omega0)
        for i, s in enumerate(self.scales):
            k_re, k_im = self.morlet_kernel(s, self.kernel_size, self.omega0)
            # Non-persistent: deterministic, recomputed at construction, not stored in checkpoints.
            self.register_buffer(f"kernel_{i}", torch.stack([k_re, k_im], dim=0), persistent=False)

    @staticmethod
    def morlet_kernel(scale: int, kernel_size: int, omega0: float) -> Tuple[torch.Tensor, torch.Tensor]:
        half = (kernel_size * scale) // 2
        t = torch.arange(-half, half + 1, dtype=torch.float64)
        x = t / scale
        gauss = torch.exp(-0.5 * x ** 2)
        norm = math.pi ** -0.25 / math.sqrt(scale)
        correction = math.exp(-0.5 * omega0 ** 2)
        re = norm * (torch.cos(omega0 * x) - correction) * gauss
        im = norm * torch.sin(omega0 * x) * gauss
        # cross-correlation with conj(psi): real part unchanged, imaginary part negated
        return re.float(), (-im).float()

    def forward_scale(self, x: torch.Tensor, i: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """x: (B, C, L) already masked. Returns Re, Im coefficients at scale i, each (B, C, L)."""
        B, C, L = x.shape
        k = getattr(self, f"kernel_{i}").to(dtype=x.dtype)          # (2, K)
        K = k.shape[-1]
        weight = k.unsqueeze(1).repeat(C, 1, 1)                     # (2C, 1, K): [re_0, im_0, re_1, ...]
        out = F.conv1d(x, weight, padding=K // 2, groups=C)         # (B, 2C, L)
        out = out.view(B, C, 2, L)
        return out[:, :, 0], out[:, :, 1]

    def forward(self, x: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        """x: (B, C, L) → complex tensor (B, C, S, L)."""
        m = None if mask is None else mask.to(x.dtype).unsqueeze(1)
        if m is not None:
            x = x * m
        re, im = [], []
        for i in range(len(self.scales)):
            r, j = self.forward_scale(x, i)
            if m is not None:
                r, j = r * m, j * m
            re.append(r)
            im.append(j)
        return torch.complex(torch.stack(re, dim=2), torch.stack(im, dim=2))


class MultiScaleComplexRepresentation(nn.Module):
    """Shared stage: H (B, L, d1) → CWT → (Re, Im) → learnable 1x1 convs → (B, d2, S, L) each.

    The 1x1 convolution is applied scale by scale (identical to a Conv2d with a 1x1 kernel over
    the (S, L) plane) so that the full d1-channel complex tensor never has to be materialised.
    """

    def __init__(self, d1: int, d2: int, scales: Sequence[int], kernel_size: int = 5, omega0: float = 6.0):
        super().__init__()
        self.cwt = MorletCWT(scales, kernel_size, omega0)
        self.conv_real = nn.Conv1d(d1, d2, kernel_size=1)
        self.conv_imag = nn.Conv1d(d1, d2, kernel_size=1)

    @property
    def num_scales(self) -> int:
        return len(self.cwt.scales)

    def forward(self, H: torch.Tensor, mask: Optional[torch.Tensor] = None) -> Tuple[torch.Tensor, torch.Tensor]:
        x = H.transpose(1, 2)                                       # H' : (B, d1, L)
        m = None if mask is None else mask.to(x.dtype).unsqueeze(1)
        if m is not None:
            x = x * m
        real, imag = [], []
        for i in range(self.num_scales):
            re, im = self.cwt.forward_scale(x, i)
            r = self.conv_real(re)
            j = self.conv_imag(im)
            if m is not None:
                r, j = r * m, j * m
            real.append(r)
            imag.append(j)
        return torch.stack(real, dim=2), torch.stack(imag, dim=2)   # (B, d2, S, L) ×2


# --------------------------------------------------------------------------- #
# 2. Residue-wise scale attention (RWSA)
# --------------------------------------------------------------------------- #
class ResidueWiseScaleAttention(nn.Module):
    """Each residue attends over the S wavelet scales.

    Q_l = H_l W_Q (shared by both streams); K/V are stream-specific projections of the
    wavelet coefficients. Attention is a softmax over the scale axis; padded residues
    receive zero attention weight (and thus a zero context vector).
    """

    def __init__(self, d1: int, d2: int):
        super().__init__()
        self.d2 = d2
        self.W_Q = nn.Linear(d1, d2, bias=False)
        self.W_K_real = nn.Linear(d2, d2, bias=False)
        self.W_V_real = nn.Linear(d2, d2, bias=False)
        self.W_K_imag = nn.Linear(d2, d2, bias=False)
        self.W_V_imag = nn.Linear(d2, d2, bias=False)

    def _attend(self, Q, Ht, W_K, W_V, m):
        X = Ht.permute(0, 3, 2, 1)                                  # (B, L, S, d2)
        K = W_K(X)
        V = W_V(X)
        scores = torch.einsum("bld,blsd->bls", Q, K) / math.sqrt(self.d2)
        alpha = torch.softmax(scores, dim=-1)
        # alpha_{s,l} = Mask_l * softmax_s(scores): the softmax runs over the scale axis and
        # the residue mask is applied multiplicatively, so padded residues get zero weight
        # on every scale (an additive -inf mask would make the whole softmax row NaN).
        alpha = alpha * m.unsqueeze(-1)
        ctx = torch.einsum("bls,blsd->bld", alpha, V)               # (B, L, d2)
        return ctx, alpha

    def forward(self, H, H_real, H_imag, mask=None):
        m = _mask_like(mask, H)
        Q = self.W_Q(H)                                             # (B, L, d2)
        c_real, a_real = self._attend(Q, H_real, self.W_K_real, self.W_V_real, m)
        c_imag, a_imag = self._attend(Q, H_imag, self.W_K_imag, self.W_V_imag, m)
        return c_real, c_imag, (a_real, a_imag)


# --------------------------------------------------------------------------- #
# 3. Complex axis alignment (CAA)
# --------------------------------------------------------------------------- #
class ComplexAxisAlignment(nn.Module):
    """[axis1; axis2] = [[a_k, -b_k], [b_k, a_k]] [c^R; c^I]  (per channel k).

    Equivalent to multiplying (c^R + i c^I) by the learnable complex scalar (a_k + i b_k):
    a channel-wise phase rotation + amplitude scaling. Initialised to the identity.
    """

    def __init__(self, d2: int):
        super().__init__()
        self.a = nn.Parameter(torch.ones(d2))
        self.b = nn.Parameter(torch.zeros(d2))

    def forward(self, c_real, c_imag):
        axis1 = self.a * c_real - self.b * c_imag
        axis2 = self.b * c_real + self.a * c_imag
        return axis1, axis2


# --------------------------------------------------------------------------- #
# 4. Dual-stream integration
# --------------------------------------------------------------------------- #
class DualStreamIntegration(nn.Module):
    """C̄ = W_F [C^axis1 || C^axis2];  H_fused = H + W_O C̄."""

    def __init__(self, d1: int, d2: int):
        super().__init__()
        self.W_F = nn.Linear(2 * d2, d2)
        self.W_O = nn.Linear(d2, d1)

    def forward(self, H, axis1, axis2, mask=None):
        C_bar = self.W_F(torch.cat([axis1, axis2], dim=-1))
        H_fused = H + self.W_O(C_bar)
        return H_fused * _mask_like(mask, H).unsqueeze(-1)


# --------------------------------------------------------------------------- #
# 5. Mamba (selective state-space) block with explicit padding masking
# --------------------------------------------------------------------------- #
def _selective_scan_chunk(u, delta, A, Bm, Cm, h):
    """Sequential reference scan for one chunk.
    u, delta: (B, T, D); A: (D, N); Bm, Cm: (B, T, N); h: (B, D, N)."""
    ys = []
    for t in range(u.shape[1]):
        dt = delta[:, t].unsqueeze(-1)                              # (B, D, 1)
        h = torch.exp(dt * A) * h + dt * Bm[:, t].unsqueeze(1) * u[:, t].unsqueeze(-1)
        ys.append(torch.einsum("bdn,bn->bd", h, Cm[:, t]))
    return torch.stack(ys, dim=1), h


def selective_scan_reference(u, delta, A, Bm, Cm, D, chunk_size: int = 64):
    """Pure-PyTorch selective scan (Gu & Dao, 2024). Memory is bounded by gradient
    checkpointing every `chunk_size` steps. u, delta: (B, L, D) → y: (B, L, D)."""
    Bsz, L, Dm = u.shape
    h = u.new_zeros(Bsz, Dm, A.shape[1])
    ys = []
    for s in range(0, L, chunk_size):
        args = (u[:, s:s + chunk_size], delta[:, s:s + chunk_size], A,
                Bm[:, s:s + chunk_size], Cm[:, s:s + chunk_size], h)
        if torch.is_grad_enabled() and u.requires_grad:
            y, h = checkpoint(_selective_scan_chunk, *args, use_reentrant=False)
        else:
            y, h = _selective_scan_chunk(*args)
        ys.append(y)
    return torch.cat(ys, dim=1) + u * D


class MambaMixer(nn.Module):
    """Mamba selective-SSM mixer (Gu & Dao, 2024).

    Padding handling ("padded positions are explicitly masked out during Mamba state updates"):
      * inputs at padded positions are zeroed before the causal convolution,
      * the step size Δ is set to 0 at padded positions, so exp(ΔA) = 1 and ΔB·u = 0, i.e.
        the hidden state is left untouched by padded tokens,
      * outputs at padded positions are zeroed.
    Uses the fused `mamba_ssm` CUDA kernel when available, otherwise a pure-PyTorch scan.
    The parameterisation is identical in both cases, so checkpoints are interchangeable.
    """

    def __init__(self, d_model: int, d_state: int = 16, d_conv: int = 4, expand: int = 2,
                 dt_rank: Optional[int] = None, dt_min: float = 1e-3, dt_max: float = 1e-1,
                 dt_init_floor: float = 1e-4, use_fused_kernel: bool = True):
        super().__init__()
        self.d_model = d_model
        self.d_state = d_state
        self.d_inner = expand * d_model
        self.dt_rank = dt_rank or math.ceil(d_model / 16)
        self.use_fused_kernel = use_fused_kernel

        self.in_proj = nn.Linear(d_model, 2 * self.d_inner, bias=False)
        self.conv1d = nn.Conv1d(self.d_inner, self.d_inner, kernel_size=d_conv,
                                groups=self.d_inner, padding=d_conv - 1)
        self.x_proj = nn.Linear(self.d_inner, self.dt_rank + 2 * d_state, bias=False)
        self.dt_proj = nn.Linear(self.dt_rank, self.d_inner, bias=True)
        self.out_proj = nn.Linear(self.d_inner, d_model, bias=False)

        # Standard Mamba initialisation
        dt_init_std = self.dt_rank ** -0.5
        nn.init.uniform_(self.dt_proj.weight, -dt_init_std, dt_init_std)
        dt = torch.exp(torch.rand(self.d_inner) * (math.log(dt_max) - math.log(dt_min)) + math.log(dt_min))
        dt = dt.clamp(min=dt_init_floor)
        with torch.no_grad():
            self.dt_proj.bias.copy_(dt + torch.log(-torch.expm1(-dt)))   # inverse softplus
        A = torch.arange(1, d_state + 1, dtype=torch.float32).repeat(self.d_inner, 1)
        self.A_log = nn.Parameter(torch.log(A))
        self.D = nn.Parameter(torch.ones(self.d_inner))

    def forward(self, x: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        Bsz, L, _ = x.shape
        m = _mask_like(mask, x).unsqueeze(-1)                       # (B, L, 1)

        xz = self.in_proj(x)
        u, z = xz.chunk(2, dim=-1)                                  # (B, L, d_inner)
        u = u * m
        u = self.conv1d(u.transpose(1, 2))[..., :L].transpose(1, 2)  # causal depth-wise conv
        u = F.silu(u) * m

        dbc = self.x_proj(u)
        dt, Bm, Cm = torch.split(dbc, [self.dt_rank, self.d_state, self.d_state], dim=-1)
        delta = F.softplus(self.dt_proj(dt)) * m                    # Δ = 0 at padded positions
        A = -torch.exp(self.A_log.float())

        if self.use_fused_kernel and HAS_MAMBA_SSM and x.is_cuda:
            y = _mamba_selective_scan_fn(
                u.transpose(1, 2).contiguous(), delta.transpose(1, 2).contiguous(), A,
                Bm.transpose(1, 2).contiguous(), Cm.transpose(1, 2).contiguous(),
                self.D.float(), z=None, delta_bias=None, delta_softplus=False,
            ).transpose(1, 2)
        else:
            y = selective_scan_reference(u, delta, A, Bm, Cm, self.D)

        y = y * F.silu(z)
        return self.out_proj(y) * m


class MambaBlock(nn.Module):
    """Pre-norm residual Mamba block: E_m = E_{m-1} + Mixer(LayerNorm(E_{m-1}))."""

    def __init__(self, d_model: int, **mixer_kwargs):
        super().__init__()
        self.norm = nn.LayerNorm(d_model)
        self.mixer = MambaMixer(d_model, **mixer_kwargs)

    def forward(self, x, mask=None):
        out = x + self.mixer(self.norm(x), mask)
        return out * _mask_like(mask, x).unsqueeze(-1)


# --------------------------------------------------------------------------- #
# 6. Task-specific branch
# --------------------------------------------------------------------------- #
class TaskBranch(nn.Module):
    """RWSA → CAA → dual-stream integration → M Mamba blocks → MLP (logit)."""

    def __init__(self, d1: int, d2: int, n_mamba: int = 3, d_state: int = 16, d_conv: int = 4,
                 expand: int = 2, head_hidden: int = 256, dropout: float = 0.1,
                 use_fused_kernel: bool = True):
        super().__init__()
        self.rwsa = ResidueWiseScaleAttention(d1, d2)
        self.caa = ComplexAxisAlignment(d2)
        self.integration = DualStreamIntegration(d1, d2)
        self.mamba = nn.ModuleList([
            MambaBlock(d1, d_state=d_state, d_conv=d_conv, expand=expand,
                       use_fused_kernel=use_fused_kernel)
            for _ in range(n_mamba)
        ])
        self.norm = nn.LayerNorm(d1)
        self.head = nn.Sequential(
            nn.Linear(d1, head_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(head_hidden, 1),
        )

    def forward(self, H, H_real, H_imag, mask=None, return_attention: bool = False):
        c_real, c_imag, att = self.rwsa(H, H_real, H_imag, mask)
        axis1, axis2 = self.caa(c_real, c_imag)
        E = self.integration(H, axis1, axis2, mask)                 # E_0 = H_fused
        for blk in self.mamba:
            E = blk(E, mask)
        logits = self.head(self.norm(E)).squeeze(-1)                # (B, L); sigmoid applied outside
        if return_attention:
            return logits, att
        return logits
