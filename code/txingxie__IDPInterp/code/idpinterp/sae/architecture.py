"""
Sparse autoencoder architectures: TopK and L1.

Both architectures follow the standard SAE formulation:
    Encode:  z = activate(W_enc @ (x - b_dec) + b_enc)
    Decode:  x_hat = W_dec @ z + b_dec

The decoder columns W_dec[:, j] are constrained to unit norm after each
optimizer step to prevent magnitude cheating (the network learning to
shrink activations and scale up weights to artificially reduce loss).
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class BaseSAE(nn.Module):
    """
    Abstract base class for Sparse Autoencoders.
    
    Provides common interface: encode(), decode(), forward(), loss(),
    normalize_decoder(), and the d_hidden property.
    """

    d_in: int
    d_latent: int

    @property
    def d_hidden(self) -> int:
        """Alias for d_latent — used by trainer and normalizer."""
        return self.d_latent

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        z = self.encode(x)
        x_hat = self.decode(z)
        return x_hat, z

    def loss(
        self, x: torch.Tensor, l1_mult: float = 1.0,
        dead_mask: torch.Tensor | None = None,
    ) -> tuple[dict, torch.Tensor]:
        """
        Compute reconstruction + sparsity loss.

        `dead_mask` (bool, d_latent) marks currently-dead features; used only by
        TopKSAE for the AuxK auxiliary loss. Ignored here (and by L1SAE).

        Returns:
            loss_dict: Dict with keys 'loss', 'mse_loss', 'l1_loss', 'fve'.
            z: Latent activations (for dead-feature tracking).
        """
        x_hat, z = self.forward(x)

        # MSE reconstruction loss
        mse_loss = F.mse_loss(x_hat, x)

        # L1 sparsity penalty (overridden by subclasses if needed)
        l1_loss = self._l1_penalty(z, l1_mult)

        total_loss = mse_loss + l1_loss

        # Fraction of variance explained
        with torch.no_grad():
            total_var = ((x - x.mean(dim=0)) ** 2).sum()
            resid_var = ((x - x_hat) ** 2).sum()
            fve = 1.0 - (resid_var / (total_var + 1e-10))

        return {
            "loss": total_loss,
            "mse_loss": mse_loss.item(),
            "l1_loss": l1_loss.item() if isinstance(l1_loss, torch.Tensor) else l1_loss,
            "aux_loss": 0.0,
            "fve": fve.item(),
        }, z

    def _l1_penalty(self, z: torch.Tensor, l1_mult: float) -> torch.Tensor:
        """Default: no L1 penalty. Overridden by L1SAE."""
        return torch.tensor(0.0, device=z.device)

    @torch.no_grad()
    def normalize_decoder(self) -> None:
        """Constrain decoder columns to unit norm. Call after each optimizer step."""
        norms = self.W_dec.norm(dim=1, keepdim=True).clamp(min=1e-8)
        self.W_dec.data = self.W_dec.data / norms


class TopKSAE(BaseSAE):
    """
    Sparse Autoencoder with explicit Top-K sparsity constraint.
    
    Only the top-k activation values are retained; all others are zeroed.
    This provides exact control over the L0 sparsity level.
    
    Args:
        d_in: Input dimensionality (PLM hidden_dim).
        expansion_factor: Ratio of latent dims to input dims.
        k: Number of active features per input.
    """

    def __init__(
        self, d_in: int, expansion_factor: int = 8, k: int = 32,
        aux_k: int = 512, aux_coeff: float = 0.0,
    ):
        super().__init__()
        self.d_in = d_in
        self.d_latent = d_in * expansion_factor
        self.k = k
        # AuxK (Gao et al. 2024): reconstruct the residual from the top-`aux_k`
        # *dead* features, weighted by `aux_coeff`. Off when aux_coeff == 0.
        self.aux_k = aux_k
        self.aux_coeff = aux_coeff

        # Encoder: input → latent
        self.W_enc = nn.Parameter(torch.empty(d_in, self.d_latent))
        self.b_enc = nn.Parameter(torch.zeros(self.d_latent))

        # Decoder: latent → input
        self.W_dec = nn.Parameter(torch.empty(self.d_latent, d_in))
        self.b_dec = nn.Parameter(torch.zeros(d_in))

        # Initialize with small random values
        nn.init.kaiming_uniform_(self.W_enc)
        nn.init.kaiming_uniform_(self.W_dec)

    def _pre_acts(self, x: torch.Tensor) -> torch.Tensor:
        return (x - self.b_dec) @ self.W_enc + self.b_enc

    def _topk_z(self, pre_acts: torch.Tensor) -> torch.Tensor:
        """Keep only the top-k pre-activations (ReLU'd); zero the rest."""
        topk_values, topk_indices = pre_acts.topk(self.k, dim=-1)
        z = torch.zeros_like(pre_acts)
        z.scatter_(-1, topk_indices, F.relu(topk_values))
        return z

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode input to sparse latent representation."""
        return self._topk_z(self._pre_acts(x))

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """Decode from latent representation back to input space."""
        return z @ self.W_dec + self.b_dec

    def loss(
        self, x: torch.Tensor, l1_mult: float = 1.0,
        dead_mask: torch.Tensor | None = None,
    ) -> tuple[dict, torch.Tensor]:
        """MSE reconstruction + AuxK dead-feature revival loss."""
        pre = self._pre_acts(x)
        z = self._topk_z(pre)
        x_hat = self.decode(z)
        mse_loss = F.mse_loss(x_hat, x)

        # AuxK: let the dead features explain the residual the live ones missed.
        aux_loss = torch.zeros((), device=x.device, dtype=mse_loss.dtype)
        if (
            self.aux_coeff > 0
            and dead_mask is not None
            and bool(dead_mask.any())
        ):
            n_dead = int(dead_mask.sum())
            kaux = min(self.aux_k, n_dead)
            if kaux > 0:
                resid = (x - x_hat).detach()
                # Restrict the top-k_aux selection to dead features only.
                pre_dead = pre.masked_fill(~dead_mask.unsqueeze(0), float("-inf"))
                av, ai = pre_dead.topk(kaux, dim=-1)
                z_aux = torch.zeros_like(pre)
                z_aux.scatter_(-1, ai, F.relu(av))
                e_hat = z_aux @ self.W_dec  # no b_dec: predicting the residual
                aux_loss = self.aux_coeff * F.mse_loss(e_hat, resid)

        total_loss = mse_loss + aux_loss

        with torch.no_grad():
            total_var = ((x - x.mean(dim=0)) ** 2).sum()
            resid_var = ((x - x_hat) ** 2).sum()
            fve = 1.0 - (resid_var / (total_var + 1e-10))

        return {
            "loss": total_loss,
            "mse_loss": mse_loss.item(),
            "l1_loss": 0.0,
            "aux_loss": aux_loss.item() if isinstance(aux_loss, torch.Tensor) else aux_loss,
            "fve": fve.item(),
        }, z


class L1SAE(BaseSAE):
    """
    Sparse autoencoder with an L1 penalty on activations.
    
    Uses ReLU activation instead of Top-K. Sparsity is encouraged via
    an L1 penalty term in the loss function. The L1 coefficient is
    ramped linearly from 0 to its target value over the first 5% of
    training steps.
    
    Args:
        d_in: Input dimensionality (PLM hidden_dim).
        expansion_factor: Ratio of latent dims to input dims.
        l1_coeff: Target L1 penalty coefficient.
    """

    def __init__(self, d_in: int, expansion_factor: int = 8, l1_coeff: float = 0.1):
        super().__init__()
        self.d_in = d_in
        self.d_latent = d_in * expansion_factor
        self.l1_coeff_target = l1_coeff

        # Encoder
        self.W_enc = nn.Parameter(torch.empty(d_in, self.d_latent))
        self.b_enc = nn.Parameter(torch.zeros(self.d_latent))

        # Decoder
        self.W_dec = nn.Parameter(torch.empty(self.d_latent, d_in))
        self.b_dec = nn.Parameter(torch.zeros(d_in))

        nn.init.kaiming_uniform_(self.W_enc)
        nn.init.kaiming_uniform_(self.W_dec)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode with ReLU activation."""
        return F.relu((x - self.b_dec) @ self.W_enc + self.b_enc)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return z @ self.W_dec + self.b_dec

    def _l1_penalty(self, z: torch.Tensor, l1_mult: float) -> torch.Tensor:
        """L1 penalty on activations, scaled by the current warmup multiplier."""
        return self.l1_coeff_target * l1_mult * z.abs().mean()


def build_sae(sae_cfg, model_cfg) -> BaseSAE:
    """
    Factory function to build an SAE from Hydra config.
    
    Args:
        sae_cfg: SAE config (type, expansion_factor, k/l1_coeff).
        model_cfg: Model config (hidden_dim).
    
    Returns:
        Initialized SAE module.
    """
    d_in = model_cfg.hidden_dim
    expansion = sae_cfg.expansion_factor

    if sae_cfg.type == "topk":
        return TopKSAE(
            d_in=d_in, expansion_factor=expansion, k=sae_cfg.k,
            aux_k=sae_cfg.get("aux_k", 512),
            aux_coeff=sae_cfg.get("aux_coeff", 0.0),
        )
    elif sae_cfg.type == "l1":
        return L1SAE(d_in=d_in, expansion_factor=expansion, l1_coeff=sae_cfg.l1_coeff)
    else:
        raise ValueError(f"Unknown SAE type: {sae_cfg.type}")
