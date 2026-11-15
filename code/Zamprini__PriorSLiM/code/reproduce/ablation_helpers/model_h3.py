from pathlib import Path
import sys
import torch
from torch import nn
import model_jwarm as baseline
FEATURE_NAMES = ['serine_fraction_29', 'tyrosine_fraction_29', 'phenylalanine_fraction_29']

def extend_density(intrinsic):
    if intrinsic.chemical_mean.numel() != 61:
        raise ValueError('H3 requires chemical58 + S/Y/F fractions')
    intrinsic.chemistry_indices = (*baseline.CHEMISTRY_INDICES, 58, 59, 60)
    chemistry = intrinsic.chemistry
    chemistry.group_indices[0] = torch.cat((chemistry.group_indices[0], torch.arange(55, 58)))
    with torch.random.fork_rng(devices=[]):
        for name in ('value', 'gate'):
            old = getattr(chemistry.encoders[0], name)
            new = nn.Linear(old.in_features + 3, old.out_features)
            with torch.no_grad():
                new.weight.zero_()
                new.weight[:, :old.in_features].copy_(old.weight)
                new.bias.copy_(old.bias)
            setattr(chemistry.encoders[0], name, new)

class H3Intrinsic(baseline.IntrinsicEncoder):

    def __init__(self, scalers, motif_dropout=0.4):
        super().__init__(scalers, motif_dropout)
        extend_density(self)

class H3Joint(baseline.JointModel):

    def __init__(self, scalers, feature_count, seed, motif_dropout=0.4):
        super().__init__(scalers, feature_count, seed, motif_dropout)
        extend_density(self.intrinsic)
