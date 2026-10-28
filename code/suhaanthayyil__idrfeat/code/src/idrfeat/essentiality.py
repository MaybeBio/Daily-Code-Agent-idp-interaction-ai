from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_depmap_gene_effect(
    path: str | Path,
    essential_cutoff: float = -0.5,
    cell_line: str | None = None,
) -> pd.DataFrame:
    path = Path(path)
    raw = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    first = raw.columns[0]
    genes = [c for c in raw.columns if c != first]
    if cell_line is not None:
        row = raw[raw[first] == cell_line]
        if row.empty:
            raise ValueError(f"cell line {cell_line} not found")
        effect = row[genes].iloc[0].to_numpy(dtype=float)
    else:
        effect = raw[genes].to_numpy(dtype=float).mean(axis=0)
    names = [g.split(" (")[0] for g in genes]
    out = pd.DataFrame({"gene": names, "gene_effect": effect})
    out["gene_essential"] = (out["gene_effect"] <= essential_cutoff).astype(int)
    return out


def attach_gene_effect(
    df: pd.DataFrame,
    depmap: pd.DataFrame,
    gene_col: str = "gene",
) -> pd.DataFrame:
    eff = dict(zip(depmap["gene"], depmap["gene_effect"], strict=True))
    ess = dict(zip(depmap["gene"], depmap["gene_essential"], strict=True))
    out = df.copy()
    out["gene_effect"] = out[gene_col].map(eff)
    out["gene_essential"] = out[gene_col].map(ess)
    return out


def depletion_by_region_stratified(
    df: pd.DataFrame,
    essential_col: str = "gene_essential",
) -> pd.DataFrame:
    covered = df[df["peptide_count"] > 0]
    rows = []
    for ess_val, ess_name in [(0, "non_essential"), (1, "essential")]:
        stratum = covered[covered[essential_col] == ess_val]
        for region, sub in [("IDR", stratum[stratum["idr"] == 1]), ("ordered", stratum[stratum["idr"] == 0])]:
            rows.append(
                {
                    "gene_group": ess_name,
                    "region": region,
                    "n": int(len(sub)),
                    "pct_depleted": float(sub["depleted"].mean() * 100) if len(sub) else float("nan"),
                    "mean_lfc": float(sub["lfc"].mean()) if len(sub) else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def idr_enrichment_within_protein(df: pd.DataFrame) -> pd.DataFrame:
    covered = df[df["peptide_count"] > 0]
    rows = []
    for acc, sub in covered.groupby("accession"):
        idr = sub[sub["idr"] == 1]
        ordered = sub[sub["idr"] == 0]
        if len(idr) == 0 or len(ordered) == 0:
            continue
        rows.append(
            {
                "accession": acc,
                "idr_depleted": float(idr["depleted"].mean()),
                "ordered_depleted": float(ordered["depleted"].mean()),
                "diff": float(idr["depleted"].mean() - ordered["depleted"].mean()),
                "n_idr": int(len(idr)),
                "n_ordered": int(len(ordered)),
            }
        )
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    return res.sort_values("diff", ascending=False).reset_index(drop=True)


def sign_test(diffs) -> dict:
    d = np.asarray([x for x in diffs if not np.isnan(x) and x != 0.0], dtype=float)
    n = len(d)
    pos = int((d > 0).sum())
    return {
        "n": n,
        "n_idr_higher": pos,
        "n_ordered_higher": n - pos,
        "frac_idr_higher": pos / n if n else float("nan"),
    }
