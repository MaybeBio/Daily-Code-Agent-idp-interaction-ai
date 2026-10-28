from __future__ import annotations

import numpy as np
import pandas as pd

from idrfeat.essentiality import (
    attach_gene_effect,
    depletion_by_region_stratified,
    idr_enrichment_within_protein,
    load_depmap_gene_effect,
    sign_test,
)


def test_load_depmap_mean_and_essential_flag(tmp_path) -> None:
    raw = pd.DataFrame(
        {
            "cell_line": ["A", "B"],
            "TP53 (7157)": [-1.0, -1.2],
            "XYZ (999)": [0.1, -0.1],
        }
    )
    path = tmp_path / "depmap.csv"
    raw.to_csv(path, index=False)
    out = load_depmap_gene_effect(path, essential_cutoff=-0.5).set_index("gene")
    assert out.loc["TP53", "gene_essential"] == 1
    assert out.loc["XYZ", "gene_essential"] == 0


def test_load_depmap_specific_cell_line(tmp_path) -> None:
    raw = pd.DataFrame({"cell_line": ["A", "B"], "TP53 (7157)": [-1.0, 0.4]})
    path = tmp_path / "depmap.csv"
    raw.to_csv(path, index=False)
    out = load_depmap_gene_effect(path, essential_cutoff=-0.5, cell_line="B").set_index("gene")
    assert out.loc["TP53", "gene_essential"] == 0


def _screen() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene": ["G1", "G1", "G2", "G2"],
            "accession": ["P1", "P1", "P2", "P2"],
            "idr": [1, 0, 1, 0],
            "depleted": [0.0, 1.0, 1.0, 0.0],
            "lfc": [0.2, -1.0, -1.0, 0.2],
            "peptide_count": [3, 3, 3, 3],
        }
    )


def test_attach_and_stratify() -> None:
    depmap = pd.DataFrame(
        {"gene": ["G1", "G2"], "gene_effect": [-1.0, 0.0], "gene_essential": [1, 0]}
    )
    joined = attach_gene_effect(_screen(), depmap)
    assert joined["gene_essential"].tolist() == [1, 1, 0, 0]
    strat = depletion_by_region_stratified(joined)
    assert set(strat["gene_group"]) == {"essential", "non_essential"}
    assert set(strat["region"]) == {"IDR", "ordered"}


def test_idr_enrichment_within_protein_and_sign_test() -> None:
    enr = idr_enrichment_within_protein(_screen())
    assert set(enr["accession"]) == {"P1", "P2"}
    st = sign_test(enr["diff"])
    assert st["n"] == 2
    assert st["n_idr_higher"] + st["n_ordered_higher"] == 2


def test_sign_test_ignores_zero_and_nan() -> None:
    st = sign_test([0.0, np.nan, 0.3, -0.1])
    assert st["n"] == 2
    assert st["n_idr_higher"] == 1
