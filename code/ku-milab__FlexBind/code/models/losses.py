"""Training objective (Section "Training Objective"):
task-specific asymmetric focal loss + uncertainty-based task weighting."""
from typing import Dict, Sequence

import torch
import torch.nn as nn
import torch.nn.functional as F


def asymmetric_focal_loss(logits: torch.Tensor, targets: torch.Tensor, valid: torch.Tensor,
                          gamma_pos: float = 0.0, gamma_neg: float = 2.0) -> torch.Tensor:
    """L_t = -mean_{valid residues} [ y (1-p)^{g+} log p + (1-y) p^{g-} log(1-p) ].

    Computed from logits with log-sigmoid for numerical stability. Padded / unlabelled
    residues are excluded through `valid`.
    """
    if valid.sum() == 0:
        return logits.sum() * 0.0
    x = logits[valid].float()
    y = targets[valid].float()
    p = torch.sigmoid(x)
    loss_pos = y * (1.0 - p).pow(gamma_pos) * F.logsigmoid(x)
    loss_neg = (1.0 - y) * p.pow(gamma_neg) * F.logsigmoid(-x)
    return -(loss_pos + loss_neg).mean()


class UncertaintyWeightedLoss(nn.Module):
    """L = sum_t ( L_t / (2 sigma_t^2) + log sigma_t )   (Kendall et al., 2018).

    sigma_t is parameterised as exp(log_sigma_t) for positivity; sigma_t is initialised to
    `init_sigma` (1.0 in the paper) and optimised jointly with the network.
    """

    def __init__(self, tasks: Sequence[str], init_sigma: float = 1.0):
        super().__init__()
        self.tasks = list(tasks)
        self.log_sigma = nn.Parameter(torch.full((len(self.tasks),), float(init_sigma)).log())

    def forward(self, task_losses: Dict[str, torch.Tensor]) -> torch.Tensor:
        total = 0.0
        for i, t in enumerate(self.tasks):
            s = self.log_sigma[i]
            total = total + 0.5 * torch.exp(-2.0 * s) * task_losses[t] + s
        return total

    def sigmas(self) -> Dict[str, float]:
        return {t: float(self.log_sigma[i].exp()) for i, t in enumerate(self.tasks)}
