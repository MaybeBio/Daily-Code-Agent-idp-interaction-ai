from __future__ import annotations

from .charge import charge_features
from .complexity import longest_single_run, seg_lowcomplexity_frac
from .composition import composition_features, mean_hydropathy


def segment_sequence_features(seg: str) -> dict[str, float]:
    f: dict[str, float] = {}
    f.update(composition_features(seg))
    f.update(charge_features(seg))
    f["hydropathy_mean"] = mean_hydropathy(seg)
    f["lowcomplexity_frac"] = seg_lowcomplexity_frac(seg)
    f["longest_single_run"] = float(longest_single_run(seg))
    return f


def apply_mutation(seq: str, pos: int, alt: str) -> str:
    if pos < 1 or pos > len(seq):
        raise ValueError("pos out of range")
    if len(alt) != 1:
        raise ValueError("alt must be a single residue")
    return seq[: pos - 1] + alt + seq[pos:]


def variant_delta(seq: str, pos: int, alt: str, start: int, end: int) -> dict[str, float]:
    if not (start <= pos <= end):
        raise ValueError("pos must fall inside the segment")
    wt_seg = seq[start - 1 : end]
    mut_seg = apply_mutation(seq, pos, alt)[start - 1 : end]
    wt = segment_sequence_features(wt_seg)
    mu = segment_sequence_features(mut_seg)
    return {f"d_{k}": mu[k] - wt[k] for k in wt}
