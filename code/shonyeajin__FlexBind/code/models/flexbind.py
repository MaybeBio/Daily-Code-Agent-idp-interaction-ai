from typing import Dict, Optional, Sequence

import torch
import torch.nn as nn

from .layers import MultiScaleComplexRepresentation, TaskBranch

TASKS = ("disorder", "protein", "rna", "dna")   # IDR identification + protein/RNA/DNA binding


class FlexBind(nn.Module):
    """FlexBind: multi-scale complex learning for joint IDR and binding-residue prediction.

    Input : frozen ProtT5 residue embeddings H (B, L, d1) and a binary mask (B, L).
    Shared: complex-Morlet CWT + learnable 1x1 convolutions (multi-scale complex representation).
    Per task (independently parameterised): RWSA → CAA → dual-stream integration
              → M Mamba blocks → MLP.
    Output: dict task → logits (B, L).  Probabilities are Sigmoid(logits) (see `predict_proba`).
    """

    def __init__(
        self,
        d_plm: int = 1024,
        d2: int = 128,
        scales: Sequence[int] = (1, 3, 4, 5, 8, 16, 32),
        kernel_size: int = 5,
        omega0: float = 6.0,
        n_mamba: int = 3,
        d_state: int = 16,
        d_conv: int = 4,
        expand: int = 2,
        head_hidden: int = 256,
        dropout: float = 0.1,
        tasks: Sequence[str] = TASKS,
        use_fused_kernel: bool = True,
    ):
        super().__init__()
        self.config = dict(
            d_plm=d_plm, d2=d2, scales=list(scales), kernel_size=kernel_size, omega0=omega0,
            n_mamba=n_mamba, d_state=d_state, d_conv=d_conv, expand=expand,
            head_hidden=head_hidden, dropout=dropout, tasks=list(tasks),
            use_fused_kernel=use_fused_kernel,
        )
        self.tasks = list(tasks)
        self.shared = MultiScaleComplexRepresentation(d_plm, d2, scales, kernel_size, omega0)
        self.branches = nn.ModuleDict({
            t: TaskBranch(d_plm, d2, n_mamba=n_mamba, d_state=d_state, d_conv=d_conv,
                          expand=expand, head_hidden=head_hidden, dropout=dropout,
                          use_fused_kernel=use_fused_kernel)
            for t in self.tasks
        })

    def forward(self, emb: torch.Tensor, mask: Optional[torch.Tensor] = None,
                return_attention: bool = False) -> Dict[str, torch.Tensor]:
        H_real, H_imag = self.shared(emb, mask)
        out, att = {}, {}
        for t in self.tasks:
            if return_attention:
                out[t], att[t] = self.branches[t](emb, H_real, H_imag, mask, return_attention=True)
            else:
                out[t] = self.branches[t](emb, H_real, H_imag, mask)
        return (out, att) if return_attention else out

    @torch.no_grad()
    def predict_proba(self, emb: torch.Tensor, mask: Optional[torch.Tensor] = None) -> Dict[str, torch.Tensor]:
        return {t: torch.sigmoid(v) for t, v in self.forward(emb, mask).items()}


# Backward-compatible name used by the previous version of the repository.
IDRModel = FlexBind
