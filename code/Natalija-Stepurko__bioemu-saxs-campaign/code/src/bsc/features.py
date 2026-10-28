"""Pre-acquisition features for the validation experiments (DESIGN §7, D1 and D2): quantities known
before a profile is measured. Sequence-level descriptors, the disorder score and class, the radius
of gyration expected from chain length by the scaling laws for folded and for disordered chains, and
the model's own ensemble Rg and spread from its back-calculated curves. Nothing from the experimental
curve enters.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from bsc import config as C

# Rg (Å) against chain length N: folded proteins ~ 2.2 N^0.38 (Skolnick & Kolinski); chemically
# denatured / disordered chains 1.927 N^0.598 (Kohn et al. 2004, PNAS)
FOLDED_LAW = (2.2, 0.38)
DISORDERED_LAW = (1.927, 0.598)

GROUPS = {"hydrophobic": "AVILMFWY", "aromatic": "FWY", "positive": "KR", "negative": "DE", "proline": "P",
          "glycine": "G", "polar_st": "ST", "amide_qn": "QN", "cysteine": "C", "histidine": "H"}
TARGET_COLS = ["chi2_raw", "phi_at_chi2_1", "phi_at_chi2_2", "kind"]


def scaling_rg(length: np.ndarray, law: tuple[float, float]) -> np.ndarray:
    return law[0] * np.asarray(length, float) ** law[1]


def composition(seq: str) -> dict[str, float]:
    s = seq.upper()
    n = max(1, len(s))
    out = {f"frac_{k}": sum(s.count(a) for a in aa) / n for k, aa in GROUPS.items()}
    pos, neg = out["frac_positive"], out["frac_negative"]
    out["fcr"] = pos + neg                   # fraction of charged residues
    out["ncpr"] = pos - neg                  # net charge per residue
    return out


def build(d: pd.DataFrame, ent: pd.DataFrame) -> pd.DataFrame:
    """One row per clean entry of the model under study: features, targets and the sequence group."""
    m = d[(d["model"] == C.MODEL) & d["clean"]].copy()
    seq = ent.set_index("label")["sequence"]
    comp = pd.DataFrame([composition(seq[lab]) for lab in m["label"]], index=m.index)
    X = pd.DataFrame(index=m.index)
    X["length"] = m["length"].astype(float)
    X["log_length"] = np.log10(X["length"])
    X["disorder_mean"] = m["disorder_mean"]
    for cls, _, _ in C.DISORDER_CLASSES:
        X[f"class_{cls.split()[0]}"] = (m["disorder_class"] == cls).astype(float)
    X = pd.concat([X, comp], axis=1)
    X["rg_law_folded"] = scaling_rg(X["length"], FOLDED_LAW)
    X["rg_law_disordered"] = scaling_rg(X["length"], DISORDERED_LAW)
    f = X["disorder_mean"]
    X["rg_law_blend"] = (1 - f) * X["rg_law_folded"] + f * X["rg_law_disordered"]
    X["rg_model"] = m["rg_model"]
    X["rg_model_over_law"] = m["rg_model"] / X["rg_law_blend"]
    X["rg_model_over_law_disordered"] = m["rg_model"] / X["rg_law_disordered"]
    X["rg_conformer_spread"] = m["rg_conformer_iqr"] / m["rg_conformer_median"]
    out = X.copy()
    out.insert(0, "label", m["label"].to_numpy())
    for c in TARGET_COLS:
        out[c] = m[c].to_numpy()
    out["disorder_class"] = m["disorder_class"].to_numpy()
    return out.reset_index(drop=True)


def feature_columns(table: pd.DataFrame) -> list[str]:
    skip = {"label", "disorder_class", "group", "dup_cluster", *TARGET_COLS}
    return [c for c in table.columns if c not in skip]
