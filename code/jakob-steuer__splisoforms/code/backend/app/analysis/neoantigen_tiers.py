"""
neoantigen_tiers.py — the single source of truth for neoepitope ranking (v2).

Rationale (2026-07): the original scheme was built around a single MHCflurry run,
so the headline tier was gated on splice/PTM origin + pLDDT. We now have three
independent predictors — MHCflurry, NetMHCpan-4.2, BigMHC — giving multi-tool
PRESENTATION consensus plus a distinct IMMUNOGENICITY signal (BigMHC-IM). v2 ranks
on immunological evidence and demotes splice/PTM/expression/fold to orthogonal
DECORATING flags (they annotate, they do not gate).

Two axes drive the tier:
  1. presentation_tier  — multi-tool consensus that the peptide is MHC-presented.
  2. is_immunogenic     — BigMHC-IM predicts T-cell reactivity.

priority_tier (P1–P4) is the transparent conjunction of those two. Everything else
(is_expressed, splice-enabled, PTM-disrupting, high-fidelity fold) is a flag layered
on top so the paper can say "of P1 epitopes, X% are splice-enabled" rather than
gating on it.

All thresholds are field-standard EL %rank cutoffs; the only genuinely tunable
immunology knob is BIGMHC_IM_THRESH. Keep this module import-only (no DB, no I/O) so
both the loader and the re-tier backfill share identical logic.
"""
from __future__ import annotations

from typing import Optional, Tuple

# ── thresholds ───────────────────────────────────────────────────────────────
# Eluted-ligand %rank (lower = better) — the modern presentation standard.
STRONG_RANK = 0.5   # strong presented
WEAK_RANK = 2.0     # weak presented
# BigMHC EL is a 0–1 presentation score (higher = better); no natural %rank.
BIGMHC_EL_STRONG = 0.50
BIGMHC_EL_WEAK = 0.25
# BigMHC IM is a 0–1 immunogenicity score (higher = better). The one tunable knob.
BIGMHC_IM_THRESH = 0.50
# high-fidelity fold gate on the epitope window (pLDDT of the novel structure)
FOLD_MEAN_MIN = 70.0
FOLD_MIN_MIN = 50.0

# ── descriptive tier labels (kept here so API/frontend/docs stay in sync) ─────
PRESENTATION_LABELS = {
    "STRONG": "Strong presentation",
    "MODERATE": "Moderate presentation",
    "WEAK": "Weak presentation",
    "NONE": "Not presented",
}
PRIORITY_LABELS = {
    "P1": "High-confidence",
    "P2": "Strong presentation",
    "P3": "Supported",
    "P4": "Candidate",
}
PRIORITY_DESCRIPTIONS = {
    "P1": "Strong multi-tool presentation and predicted immunogenic.",
    "P2": "Strong multi-tool presentation; not predicted immunogenic.",
    "P3": "Moderate presentation, predicted immunogenic.",
    "P4": "Presentation evidence from a single tool or without immunogenicity support.",
}
PRIORITY_ORDER = {"P1": 1, "P2": 2, "P3": 3, "P4": 4}  # for sorting (lower = better)


def presentation_calls(
    mhcflurry_percentile: Optional[float],
    netmhcpan_rank_el: Optional[float],
    bigmhc_el_score: Optional[float],
) -> Tuple[int, int]:
    """Return (n_strong, n_weak): how many of the three tools call the peptide
    presented at the strong / weak stringency. Missing tools simply don't vote."""
    n_strong = n_weak = 0
    for strong, weak in (
        (mhcflurry_percentile is not None and mhcflurry_percentile <= STRONG_RANK,
         mhcflurry_percentile is not None and mhcflurry_percentile <= WEAK_RANK),
        (netmhcpan_rank_el is not None and netmhcpan_rank_el <= STRONG_RANK,
         netmhcpan_rank_el is not None and netmhcpan_rank_el <= WEAK_RANK),
        (bigmhc_el_score is not None and bigmhc_el_score >= BIGMHC_EL_STRONG,
         bigmhc_el_score is not None and bigmhc_el_score >= BIGMHC_EL_WEAK),
    ):
        n_strong += int(strong)
        n_weak += int(weak)
    return n_strong, n_weak


def presentation_tier(n_strong: int, n_weak: int) -> str:
    """STRONG (≥2 strong) · MODERATE (1 strong or ≥2 weak) · WEAK (1 weak) · NONE."""
    if n_strong >= 2:
        return "STRONG"
    if n_strong == 1 or n_weak >= 2:
        return "MODERATE"
    if n_weak == 1:
        return "WEAK"
    return "NONE"


def is_immunogenic(bigmhc_im_score: Optional[float]) -> bool:
    return bigmhc_im_score is not None and bigmhc_im_score >= BIGMHC_IM_THRESH


def priority_tier(pres_tier: str, immunogenic: bool) -> Optional[str]:
    """The headline P1–P4 conjunction (presentation × immunogenicity).
    Returns None when the peptide is not presented (not a candidate)."""
    if pres_tier == "STRONG":
        return "P1" if immunogenic else "P2"
    if pres_tier == "MODERATE":
        return "P3" if immunogenic else "P4"
    if pres_tier == "WEAK":
        return "P4"
    return None


def high_fidelity_fold(avg_plddt: Optional[float], min_plddt: Optional[float]) -> bool:
    return (avg_plddt is not None and min_plddt is not None
            and avg_plddt > FOLD_MEAN_MIN and min_plddt > FOLD_MIN_MIN)


def is_splice_enabled(ptm_impact_status: Optional[str]) -> bool:
    """The epitope exists only because of the splice event."""
    return bool(ptm_impact_status) and ptm_impact_status.startswith("Splice-Enabled")


def classify(
    *,
    mhcflurry_percentile: Optional[float],
    netmhcpan_rank_el: Optional[float],
    bigmhc_el_score: Optional[float],
    bigmhc_im_score: Optional[float],
) -> dict:
    """Convenience: full presentation/immunogenicity classification for one row.
    (Decorating flags that need DB context — is_expressed — are set by the caller.)"""
    n_strong, n_weak = presentation_calls(mhcflurry_percentile, netmhcpan_rank_el, bigmhc_el_score)
    pres = presentation_tier(n_strong, n_weak)
    immuno = is_immunogenic(bigmhc_im_score)
    return {
        "n_strong": n_strong,
        "n_weak": n_weak,
        "presentation_tier": pres,
        "is_immunogenic": immuno,
        "priority_tier": priority_tier(pres, immuno),
    }
