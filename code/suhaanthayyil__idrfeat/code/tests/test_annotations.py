"""Tests for annotation loaders and per-segment annotation aggregation."""

from __future__ import annotations

import gzip
import json

from idrfeat.annotations import (
    drllps_member,
    elm_segment_features,
    load_dbptm,
    load_drllps,
    load_elm,
    load_phasepro,
    load_plddt_pdb,
    phasepro_overlap,
    ptm_segment_features,
)

ELM_INDEX = {
    "P1": [
        (5, 10, "LIG", "LIG_SH3_3"),
        (8, 20, "MOD", "MOD_CK1_1"),
        (100, 110, "DOC", "DOC_MAPK_1"),
    ]
}
PTM_INDEX = {
    "P1": [(7, "phospho"), (9, "phospho"), (12, "acetyl"), (200, "ubiq")],
}


def test_elm_segment_features_overlap() -> None:
    feats = elm_segment_features("P1", 1, 15, ELM_INDEX, seg_length=15)
    assert feats["elm_motif_residues"] == 11  # positions 5..15 covered
    assert feats["elm_motif_density"] == 11 / 15
    assert feats["elm_classes_distinct"] == 2  # LIG_SH3_3 and MOD_CK1_1
    assert feats["overlaps_binding_motif"] is True  # a LIG motif overlaps


def test_elm_segment_features_no_overlap() -> None:
    feats = elm_segment_features("P1", 50, 60, ELM_INDEX, seg_length=11)
    assert feats["elm_motif_residues"] == 0
    assert feats["elm_classes_distinct"] == 0
    assert feats["overlaps_binding_motif"] is False


def test_elm_segment_features_missing_accession() -> None:
    feats = elm_segment_features("ZZZ", 1, 15, ELM_INDEX, seg_length=15)
    assert feats["elm_motif_residues"] == 0
    assert feats["overlaps_binding_motif"] is False


def test_ptm_segment_features_by_type() -> None:
    feats = ptm_segment_features("P1", 1, 15, PTM_INDEX, seg_length=15)
    assert feats["ptm_count"] == 3
    assert feats["ptm_density"] == 3 / 15
    assert feats["ptm_phospho_count"] == 2
    assert feats["ptm_acetyl_count"] == 1
    assert feats["ptm_ubiq_count"] == 0
    assert feats["ptm_methyl_count"] == 0


def test_phasepro_overlap() -> None:
    index = {"P1": [(10, 50)]}
    assert phasepro_overlap("P1", 1, 15, index) is True
    assert phasepro_overlap("P1", 60, 70, index) is False
    assert phasepro_overlap("ZZZ", 1, 15, index) is False


def test_drllps_member() -> None:
    accs = {"P1", "P2"}
    assert drllps_member("P1", accs) is True
    assert drllps_member("ZZZ", accs) is False


def test_load_elm(tmp_path) -> None:
    path = tmp_path / "elm.tsv"
    path.write_text(
        "#comment line\n"
        "Primary_Acc\tELMType\tELMIdentifier\tStart\tEnd\tInstanceLogic\n"
        "P1\tLIG\tLIG_SH3_3\t5\t10\ttrue positive\n"
        "P1\tMOD\tMOD_CK1_1\t8\t20\tfalse positive\n"
    )
    index = load_elm(path)
    assert index["P1"] == [(5, 10, "LIG", "LIG_SH3_3")]  # false positive dropped


def test_load_dbptm(tmp_path) -> None:
    path = tmp_path / "dbptm_phos.gz"
    rows = "GENE1\tP1\t7\tPhosphorylation\t12345\tSEQ\nGENE1\tP1\t9\tPhosphorylation\t12345\tSEQ\n"
    with gzip.open(path, "wt") as fh:
        fh.write(rows)
    index = load_dbptm([path])
    assert sorted(index["P1"]) == [(7, "phospho"), (9, "phospho")]


def test_load_phasepro(tmp_path) -> None:
    path = tmp_path / "phasepro.json"
    path.write_text(
        json.dumps(
            {
                "e1": {"accession": "P1", "organism": "Homo sapiens", "boundaries": "10-50;60-70"},
                "e2": {"accession": "P2", "organism": "Mus musculus", "boundaries": "1-5"},
            }
        )
    )
    index = load_phasepro(path)
    assert index["P1"] == [(10, 50), (60, 70)]
    assert "P2" not in index  # non-human dropped


def test_load_drllps(tmp_path) -> None:
    path = tmp_path / "drllps.txt"
    path.write_text(
        "UniProt ID\tSpecies\tLLPS Type\n"
        "P1\tHomo sapiens\tScaffold\n"
        "P2\tHomo sapiens\tClient\n"
        "P3\tMus musculus\tScaffold\n"
    )
    accs = load_drllps(path)
    assert accs == {"P1"}  # only human scaffold/regulator


def test_load_plddt_pdb(tmp_path) -> None:
    path = tmp_path / "model.pdb"
    # Two CA atoms with B-factors (pLDDT) 80.0 and 40.0.
    path.write_text(
        "ATOM      1  N   MET A   1      0.000   0.000   0.000  1.00 80.00           N\n"
        "ATOM      2  CA  MET A   1      1.000   0.000   0.000  1.00 80.00           C\n"
        "ATOM      3  CA  GLU A   2      2.000   0.000   0.000  1.00 40.00           C\n"
    )
    plddt = load_plddt_pdb(path)
    assert list(plddt) == [80.0, 40.0]
