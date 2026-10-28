"""Assemble one feature row per IDR segment into a tidy table.

The column set is fixed and documented in ``docs/features.md``; ``TABLE_COLUMNS`` is the single
source of truth for both the output order and the docs. Sequence-based physicochemical features
come from localCIDER. Predictor, low-complexity, and annotation features are filled with missing
values when their model, binary, or database is not supplied, so the table records what was
actually computed rather than implying a zero where the source was simply absent.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .annotations import (
    PTM_BUCKETS,
    drllps_member,
    elm_segment_features,
    phasepro_overlap,
    ptm_segment_features,
)
from .charge import charge_features
from .complexity import longest_single_run
from .composition import composition_features, mean_hydropathy
from .constants import AA
from .disorder import (
    aiupred_available,
    aiupred_binding_scores,
    aiupred_disorder_scores,
    default_backend,
    metapredict_available,
    metapredict_scores,
    primary_disorder,
    segments_from_scores,
)

ID_COLUMNS = [
    "accession",
    "idr_index",
    "idr_start",
    "idr_end",
    "seg_length",
    "sequence",
    "disorder_backend",
]

_COMPOSITION_COLUMNS = (
    [f"aa_frac_{aa}" for aa in AA]
    + ["frac_charged", "frac_polar", "frac_aromatic", "frac_proline", "frac_glycine"]
    + ["frac_disorder_promoting"]
)
_CHARGE_COLUMNS = ["fcr", "ncpr", "net_charge", "kappa", "omega", "scd", "frac_positive", "frac_negative"]
_HYDROPATHY_COLUMNS = ["hydropathy_mean"]
_COMPLEXITY_COLUMNS = ["lowcomplexity_frac", "longest_single_run"]
_DISORDER_COLUMNS = [
    "disorder_mean",
    "disorder_max",
    "disorder_frac",
    "metapredict_mean",
    "metapredict_max",
    "metapredict_frac",
    "aiupred_disorder_mean",
    "plddt_mean",
]
_BINDING_COLUMNS = ["aiupred_binding_mean", "aiupred_binding_max", "aiupred_binding_frac"]
_MOTIF_COLUMNS = [
    "elm_motif_residues",
    "elm_motif_density",
    "elm_classes_distinct",
    "overlaps_binding_motif",
]
_PTM_COLUMNS = ["ptm_count", "ptm_density"] + [
    f"ptm_{b}_{stat}" for b in PTM_BUCKETS for stat in ("count", "density")
]
_PHASESEP_COLUMNS = ["phasepro_overlap", "drllps_member"]
_CONSERVATION_COLUMNS = ["phylop_mean", "phastcons_mean"]

FEATURE_COLUMNS = (
    _COMPOSITION_COLUMNS
    + _CHARGE_COLUMNS
    + _HYDROPATHY_COLUMNS
    + _COMPLEXITY_COLUMNS
    + _DISORDER_COLUMNS
    + _BINDING_COLUMNS
    + _MOTIF_COLUMNS
    + _PTM_COLUMNS
    + _PHASESEP_COLUMNS
    + _CONSERVATION_COLUMNS
)
TABLE_COLUMNS = ID_COLUMNS + FEATURE_COLUMNS


@dataclass
class Annotations:
    """Optional annotation indexes. A None field leaves its columns as missing values."""

    elm: dict | None = None
    ptm: dict | None = None
    phasepro: dict | None = None
    drllps: set | None = None
    plddt: dict | None = None  # accession -> per-residue pLDDT array
    phylop: dict | None = None  # accession -> per-residue phyloP array
    phastcons: dict | None = None  # accession -> per-residue phastCons array
    lowcomplexity: dict | None = None  # accession -> per-residue SEG low-complexity mask


def _segment_stats(scores: np.ndarray | None, start: int, end: int, threshold: float):
    if scores is None:
        return np.nan, np.nan, np.nan
    seg = np.asarray(scores)[start - 1 : end]
    if seg.size == 0:
        return np.nan, np.nan, np.nan
    return float(seg.mean()), float(seg.max()), float((seg >= threshold).mean())


def _segment_mean(values: np.ndarray | None, start: int, end: int) -> float:
    if values is None:
        return np.nan
    seg = np.asarray(values)[start - 1 : end]
    return float(seg.mean()) if seg.size else np.nan


def _resolve_primary(cfg: dict) -> str:
    backend = cfg["disorder"].get("primary", "metapredict")
    if backend == "metapredict" and not metapredict_available():
        raise RuntimeError(
            "metapredict is required to call IDR segments; install the 'predictors' extra"
        )
    if backend == "aiupred" and not aiupred_available():
        raise RuntimeError("aiupred backend requested but aiupred is not installed")
    return backend


def _segment_row(
    acc: str,
    idr_index: int,
    start: int,
    end: int,
    full_seq: str,
    backend: str,
    primary: np.ndarray,
    meta_scores: np.ndarray | None,
    aiu_disorder: np.ndarray | None,
    aiu_binding: np.ndarray | None,
    cfg: dict,
    anns: Annotations,
) -> dict:
    seg = full_seq[start - 1 : end]
    seg_len = len(seg)
    threshold = cfg["disorder"]["score_threshold"]
    binding_threshold = cfg["disorder"]["binding_threshold"]

    row: dict = {
        "accession": acc,
        "idr_index": idr_index,
        "idr_start": start,
        "idr_end": end,
        "seg_length": seg_len,
        "sequence": seg,
        "disorder_backend": backend,
    }
    row.update(composition_features(seg))
    row.update(charge_features(seg))
    row["hydropathy_mean"] = mean_hydropathy(seg)
    row["lowcomplexity_frac"] = _segment_mean(
        anns.lowcomplexity.get(acc) if anns.lowcomplexity else None, start, end
    )
    row["longest_single_run"] = longest_single_run(seg)

    d_mean, d_max, d_frac = _segment_stats(primary, start, end, threshold)
    row["disorder_mean"], row["disorder_max"], row["disorder_frac"] = d_mean, d_max, d_frac
    m_mean, m_max, m_frac = _segment_stats(meta_scores, start, end, threshold)
    row["metapredict_mean"], row["metapredict_max"], row["metapredict_frac"] = m_mean, m_max, m_frac
    row["aiupred_disorder_mean"] = _segment_mean(aiu_disorder, start, end)
    row["plddt_mean"] = _segment_mean(anns.plddt.get(acc) if anns.plddt else None, start, end)

    b_mean, b_max, b_frac = _segment_stats(aiu_binding, start, end, binding_threshold)
    row["aiupred_binding_mean"] = b_mean
    row["aiupred_binding_max"] = b_max
    row["aiupred_binding_frac"] = b_frac

    if anns.elm is not None:
        row.update(elm_segment_features(acc, start, end, anns.elm, seg_len))
    else:
        for col in _MOTIF_COLUMNS:
            row[col] = np.nan
    if anns.ptm is not None:
        row.update(ptm_segment_features(acc, start, end, anns.ptm, seg_len))
    else:
        for col in _PTM_COLUMNS:
            row[col] = np.nan
    row["phasepro_overlap"] = (
        phasepro_overlap(acc, start, end, anns.phasepro) if anns.phasepro is not None else np.nan
    )
    row["drllps_member"] = drllps_member(acc, anns.drllps) if anns.drllps is not None else np.nan

    row["phylop_mean"] = _segment_mean(anns.phylop.get(acc) if anns.phylop else None, start, end)
    row["phastcons_mean"] = _segment_mean(
        anns.phastcons.get(acc) if anns.phastcons else None, start, end
    )
    return row


def build_feature_table(
    seqs: dict[str, str],
    cfg: dict,
    annotations: Annotations | None = None,
) -> pd.DataFrame:
    """One row per IDR segment across all input sequences, in ``TABLE_COLUMNS`` order."""
    anns = annotations or Annotations()
    backend = _resolve_primary(cfg)
    definitions = cfg["disorder"].get("definitions", [])
    threshold = cfg["disorder"]["score_threshold"]
    min_len = cfg["disorder"]["min_segment_length"]
    use_meta = "metapredict" in definitions and metapredict_available()
    use_aiupred = "aiupred" in definitions and aiupred_available()

    rows: list[dict] = []
    for acc in sorted(seqs):
        seq = seqs[acc]
        primary, _ = primary_disorder(seq, backend=backend)
        meta_scores = metapredict_scores(seq) if use_meta else None
        if backend == "metapredict" and meta_scores is None and metapredict_available():
            meta_scores = primary
        aiu_disorder = aiupred_disorder_scores(seq) if use_aiupred else None
        aiu_binding = aiupred_binding_scores(seq) if use_aiupred else None

        for idr_index, (start, end) in enumerate(
            segments_from_scores(primary, threshold, min_len), start=1
        ):
            rows.append(
                _segment_row(
                    acc,
                    idr_index,
                    start,
                    end,
                    seq,
                    backend,
                    primary,
                    meta_scores,
                    aiu_disorder,
                    aiu_binding,
                    cfg,
                    anns,
                )
            )

    df = pd.DataFrame(rows, columns=TABLE_COLUMNS)
    return df


# Re-exported so callers do not need to know the backend resolution lives here.
__all__ = [
    "Annotations",
    "FEATURE_COLUMNS",
    "TABLE_COLUMNS",
    "build_feature_table",
    "default_backend",
]
