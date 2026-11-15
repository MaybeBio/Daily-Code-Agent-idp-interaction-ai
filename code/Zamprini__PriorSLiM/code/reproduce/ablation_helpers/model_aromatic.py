from pathlib import Path
import sys
import torch
from torch import nn
ROOT = Path(__import__('os').environ['PRIORSLIM_DATA_ROOT'])
import model_jwarm as baseline
EXTRA_NAMES = [f'{metric}_{aa}{suffix}' for aa in 'FYW' for metric in ('density', 'X', 'center_enrich') for suffix in ('', '_percentile_IDR')]
EXTRA_GROUPS = [0, 4, 1, 4, 2, 4] * 3

def configure(model, scalers):
    layout = list(scalers['aromatic_layout'])
    old_chem = list(baseline.CHEMISTRY_INDICES)
    chemistry_sources = [v for v in layout if v not in baseline.MOTIF_INDICES]
    model.chemistry_indices = tuple((layout.index(v) for v in chemistry_sources))
    model.motif_indices = tuple((layout.index(v) for v in baseline.MOTIF_INDICES))
    with torch.random.fork_rng(devices=[]):
        for g, old_indices in enumerate(baseline.CHEMISTRY_GROUPS):
            retained = [(col, old_chem[loc]) for col, loc in enumerate(old_indices) if old_chem[loc] in layout]
            added = [58 + j for j, group in enumerate(EXTRA_GROUPS) if group == g]
            source = [src for _, src in retained] + added
            model.chemistry.group_indices[g] = torch.tensor([chemistry_sources.index(v) for v in source])
            for name in ('value', 'gate'):
                old = getattr(model.chemistry.encoders[g], name)
                new = nn.Linear(len(source), old.out_features)
                with torch.no_grad():
                    new.weight.zero_()
                    new.bias.copy_(old.bias)
                    new.weight[:, :len(retained)].copy_(old.weight[:, [col for col, _ in retained]])
                setattr(model.chemistry.encoders[g], name, new)

class AromaticIntrinsic(baseline.IntrinsicEncoder):

    def __init__(self, scalers, motif_dropout=0.4):
        super().__init__({k: v for k, v in scalers.items() if k != 'aromatic_layout'}, motif_dropout)
        configure(self, scalers)

class AromaticJoint(baseline.JointModel):

    def __init__(self, scalers, feature_count, seed, motif_dropout=0.4):
        super().__init__({k: v for k, v in scalers.items() if k != 'aromatic_layout'}, feature_count, seed, motif_dropout)
        configure(self.intrinsic, scalers)
