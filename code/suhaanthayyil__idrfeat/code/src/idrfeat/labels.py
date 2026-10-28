"""Loader and IDR aggregation for a per-residue protein fitness screen.

Expected per-residue columns (names configurable): Uniprot_ID, AA_loc, AA, IDR,
Fitness_depleted, LFC, P_value, Peptide_count. IDR is UniProt-derived. Residues with
Peptide_count == 0 have NaN LFC and P_value. The IDR source, the p-value correction, the
minimum peptide coverage, and the residue-to-segment aggregation rule are all configurable."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_COLS = {
    "accession": "Uniprot_ID",
    "aa_loc": "AA_loc",
    "aa": "AA",
    "idr": "IDR",
    "depleted": "Fitness_depleted",
    "lfc": "LFC",
    "pvalue": "P_value",
    "peptide_count": "Peptide_count",
}


def _read_any(path: Path):
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    sep = "\t" if path.suffix in {".txt", ".tsv"} else ","
    return pd.read_csv(path, sep=sep)


def load_residue_table(path: str | Path, cols: dict[str, str] | None = None) -> pd.DataFrame:
    mapping = {**DEFAULT_COLS, **(cols or {})}
    path = Path(path)
    raw = _read_any(path)
    out = pd.DataFrame(
        {
            "accession": raw[mapping["accession"]].astype(str),
            "aa_loc": raw[mapping["aa_loc"]].astype(int),
            "aa": raw[mapping["aa"]].astype(str),
            "idr": raw[mapping["idr"]].astype(int),
            "depleted": pd.to_numeric(raw[mapping["depleted"]], errors="coerce"),
            "lfc": pd.to_numeric(raw[mapping["lfc"]], errors="coerce"),
            "pvalue": pd.to_numeric(raw[mapping["pvalue"]], errors="coerce"),
            "peptide_count": pd.to_numeric(raw[mapping["peptide_count"]], errors="coerce").fillna(0),
        }
    )
    return out


def check_coordinate_concordance(df: pd.DataFrame, seqs: dict[str, str]) -> dict:
    checked = 0
    match = 0
    bad: set[str] = set()
    for row in df.itertuples(index=False):
        seq = seqs.get(row.accession)
        if seq is None or row.aa_loc < 1 or row.aa_loc > len(seq):
            bad.add(row.accession)
            continue
        checked += 1
        if seq[row.aa_loc - 1] == row.aa:
            match += 1
        else:
            bad.add(row.accession)
    return {
        "checked": checked,
        "match": match,
        "match_frac": match / checked if checked else float("nan"),
        "mismatched_accessions": sorted(bad),
    }


def idr_segments_from_residues(df: pd.DataFrame, min_len: int = 1) -> pd.DataFrame:
    out = []
    for acc, sub in df.groupby("accession"):
        sub = sub.sort_values("aa_loc")
        start = None
        prev = None
        for loc, flag in zip(sub["aa_loc"].tolist(), sub["idr"].tolist(), strict=True):
            if flag == 1 and start is None:
                start = prev = loc
            elif flag == 1 and loc == prev + 1:
                prev = loc
            elif flag == 1:
                if prev - start + 1 >= min_len:
                    out.append((acc, start, prev))
                start = prev = loc
            else:
                if start is not None and prev - start + 1 >= min_len:
                    out.append((acc, start, prev))
                start = prev = None
        if start is not None and prev - start + 1 >= min_len:
            out.append((acc, start, prev))
    seg = pd.DataFrame(out, columns=["accession", "idr_start", "idr_end"])
    seg["seg_length"] = seg["idr_end"] - seg["idr_start"] + 1
    return seg


def aggregate_segment_labels(
    df: pd.DataFrame,
    segments: pd.DataFrame,
    min_covered: int = 3,
    frac_cutoff: float = 0.5,
    lfc_agg: str = "mean",
) -> pd.DataFrame:
    reducer = {"mean": np.mean, "median": np.median, "min": np.min}[lfc_agg]
    by_acc = {acc: sub for acc, sub in df.groupby("accession")}
    recs = []
    for seg in segments.itertuples(index=False):
        sub = by_acc.get(seg.accession)
        if sub is None:
            recs.append((0, np.nan, np.nan, np.nan))
            continue
        within = sub[(sub["aa_loc"] >= seg.idr_start) & (sub["aa_loc"] <= seg.idr_end)]
        covered = within[within["peptide_count"] > 0]
        n = len(covered)
        if n < min_covered:
            recs.append((n, np.nan, np.nan, np.nan))
            continue
        frac = float(covered["depleted"].mean())
        mean_lfc = float(reducer(covered["lfc"].to_numpy(dtype=float)))
        label = 1.0 if frac >= frac_cutoff else 0.0
        recs.append((n, mean_lfc, frac, label))
    res = pd.DataFrame(recs, columns=["n_covered", "mean_lfc", "frac_depleted", "fitness_label"])
    return pd.concat([segments.reset_index(drop=True), res], axis=1)
