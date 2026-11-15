from __future__ import annotations
import torch
from torch import nn
from torch.nn import functional as F

class BinaryGate(nn.Module):

    def __init__(self, features):
        super().__init__()
        self.logits = nn.Parameter(torch.zeros(features))

    def forward(self):
        probability = self.logits.sigmoid()
        hard = (probability > 0.5).to(probability.dtype)
        return (hard + (probability - probability.detach()), probability)

class SparseFrozenDecode(torch.autograd.Function):

    @staticmethod
    def forward(ctx, gates, vectors, indices, values, residues, n_residues, chunk):
        result = vectors.new_zeros((n_residues, vectors.shape[1]))
        for start in range(0, len(indices), chunk):
            j, r = (indices[start:start + chunk], residues[start:start + chunk])
            a = values[start:start + chunk] * gates[j]
            result.index_add_(0, r, vectors[j] * a[:, None])
        ctx.save_for_backward(vectors, indices, values, residues)
        ctx.n_features, ctx.chunk = (len(gates), chunk)
        return result

    @staticmethod
    def backward(ctx, grad):
        vectors, indices, values, residues = ctx.saved_tensors
        gate_grad = grad.new_zeros(ctx.n_features)
        for start in range(0, len(indices), ctx.chunk):
            j, r = (indices[start:start + ctx.chunk], residues[start:start + ctx.chunk])
            contribution = (grad[r] * vectors[j]).sum(1) * values[start:start + ctx.chunk]
            gate_grad.index_add_(0, j, contribution)
        return (gate_grad, None, None, None, None, None, None)

class DecoderReadout(nn.Module):

    def __init__(self, decoder_vectors, decoder_bias, chunk=2048):
        super().__init__()
        self.register_buffer('decoder_vectors', decoder_vectors.detach().clone(), persistent=False)
        self.register_buffer('decoder_bias', decoder_bias.detach().clone(), persistent=False)
        self.projection = nn.Linear(decoder_vectors.shape[1], 64)
        self.conv3 = nn.Conv1d(64, 32, 3, padding=1)
        self.conv5 = nn.Conv1d(64, 32, 5, padding=2)
        self.attention = nn.Sequential(nn.Linear(64, 32), nn.Tanh(), nn.Linear(32, 1))
        self.head = nn.Linear(64, 1)
        nn.init.zeros_(self.head.weight)
        nn.init.zeros_(self.head.bias)
        self.chunk = chunk
        self.register_buffer('projected_vectors', None, persistent=False)
        self.register_buffer('projected_bias', None, persistent=False)

    def freeze_for_gates(self):
        self.eval().requires_grad_(False)
        with torch.no_grad():
            self.projected_vectors = self.decoder_vectors @ self.projection.weight.T
            self.projected_bias = F.linear(self.decoder_bias, self.projection.weight, self.projection.bias)

    def reconstruct(self, batch, gates):
        indices, values, residues, valid = batch
        x = SparseFrozenDecode.apply(gates, self.decoder_vectors, indices, values, residues, valid.numel(), self.chunk)
        x = x.view(len(valid), 29, -1) + self.decoder_bias
        return x * valid[..., None].to(x.dtype)

    def forward(self, batch, gates=None, fused=False, return_residue=False):
        indices, values, residues, valid = batch
        if gates is None:
            gates = self.decoder_vectors.new_ones(len(self.decoder_vectors))
        if fused:
            if self.projected_vectors is None or any((p.requires_grad for p in self.parameters())):
                raise RuntimeError('Frozen composition requires freeze_for_gates()')
            x = SparseFrozenDecode.apply(gates, self.projected_vectors, indices, values, residues, valid.numel(), self.chunk)
            h = x.view(len(valid), 29, 64) + self.projected_bias
        else:
            h = self.projection(self.reconstruct(batch, gates))
        h = h * valid[..., None].to(h.dtype)
        channels = h.transpose(1, 2)
        local = torch.cat([F.gelu(self.conv3(channels)), F.gelu(self.conv5(channels))], 1).transpose(1, 2)
        local = local * valid[..., None].to(local.dtype)
        score = self.attention(local).squeeze(-1).masked_fill(~valid, -torch.inf)
        attention = score.softmax(dim=1)
        pooled = (local * attention[..., None]).sum(1)
        z = self.head(pooled).squeeze(-1)
        return (z, h, attention) if return_residue else z

def profile_sum(events, gates):
    indices, values, owners, count = events
    return gates.new_zeros(count).index_add(0, owners, values * gates[indices])
