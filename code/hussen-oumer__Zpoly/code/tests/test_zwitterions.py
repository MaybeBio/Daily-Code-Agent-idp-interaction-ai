"""Tests for the zwitterion monomer library (rqm/zwitterion_library.md)."""

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
import yaml
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from zpoly.zwitterions import (
    AmbiguousPolymerizableGroupError,
    ChainError,
    ChargeError,
    DuplicateEntryError,
    EndGroupError,
    Ends,
    ForceField,
    IncompatibleBlocksError,
    InvalidStructureError,
    LibraryFileError,
    NoPolymerizableGroupError,
    UnknownBlockError,
    UnknownEntryError,
    UnknownInitiatorError,
    VariantError,
    build_chain,
    compose,
    custom_initiator,
    default_library,
    from_smiles,
    load_library,
    make_variant,
)
from zpoly.zwitterions._chem import assemble

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def lib():
    return default_library()


def write_yaml(tmp_path, body, name="user.yaml"):
    path = tmp_path / name
    path.write_text(textwrap.dedent(body), encoding="utf-8")
    return path


def charges(z):
    return [a.GetFormalCharge() for a in z.to_mol().GetAtoms()]


# rq-ae25abea
def test_builtin_library_entries_pass_conformance(lib):
    assert len(lib) > 0
    for z in lib:
        if not z.zwitterionic:
            continue  # comonomers: test_builtin_comonomers
        q = charges(z)
        assert sum(q) == 0, z.name
        assert any(c > 0 for c in q) and any(c < 0 for c in q), z.name
        if z.blocks is not None:
            built = Chem.MolToSmiles(assemble([lib.block(k) for k in z.blocks]))
            assert built == z.smiles, z.name
    for z in lib.entries(polymerizable=True):
        assert z.polymerizable_group is not None and z.vinyl_atoms is not None, z.name
    small = lib.entries(polymerizable=False)
    assert small and all(z.family.startswith("Zwitterionic liquid") for z in small)


# rq-83c0765f
def test_builtin_library_contains_catalog_and_paper(lib):
    for name in ["DMAPS", "MPC", "CBAA", "SPV", "ZB", "ZB2", "ZB3", "ZB4", "ZB4-C3", "ZB5",
                 "TMAO", "MIm-PS", *[f"M{i}" for i in range(1, 20) if i != 5]]:
        assert name in lib, name
    assert "M5" not in lib
    assert not lib.get("TMAO").polymerizable
    assert not lib.get("MIm-PS").polymerizable


# rq-471ad711
def test_lookup_by_alias_with_normalization(lib):
    assert lib.get("sb") is lib.get("SBMA") is lib.get("DMAPS")
    assert lib.get("cbma_1") is lib.get("CBMA-1")


# rq-881ae9e4
def test_unknown_name_gives_suggestions(lib):
    with pytest.raises(UnknownEntryError) as exc:
        lib.get("DMAP")
    assert exc.value.name == "DMAP"
    assert "DMAPS" in exc.value.suggestions
    assert "DMAPS" in str(exc.value)


# rq-ad1e99c1
def test_filter_by_family_and_status(lib):
    result = lib.entries(family="carboxybetaine", status="common")
    assert result
    assert all("Carboxybetaine" in z.family and z.status == "common" for z in result)
    assert "CBMA" in [z.name for z in result]


# rq-c9865bcb
def test_filter_by_anion_class(lib):
    result = lib.entries(anion="OSO3")
    names = [z.name for z in result]
    assert "ZB" in names and "ZB5" in names
    assert all("OSO3" in z.anion_classes for z in result)


# rq-6171c7c9
def test_filter_by_force_field(lib):
    result = lib.entries(force_field="oplsaa")
    assert "DMAPS" in [z.name for z in result]
    assert all(ForceField.OPLSAA in z.parameters for z in result)
    assert lib.entries(force_field=ForceField.GAFF2) == []


# rq-abf15540
def test_force_field_names_case_insensitive():
    assert ForceField("OPLS-AA") is ForceField.OPLSAA
    assert ForceField("GAFF2") is ForceField.GAFF2
    with pytest.raises(ValueError):
        ForceField("charmm")


# rq-82d1c963
def test_paper_properties_attached(lib):
    z = lib.get("ZB5")
    assert z.properties["fmin_kbt"] == pytest.approx(2.77)
    assert z.source == ""  # built-in references live once in REFERENCES.md
    assert "Schneider" in (Path(__file__).parents[1] / "REFERENCES.md").read_text()


# rq-e5af84ba
def test_diallyl_is_cyclopolymerizing(lib):
    assert lib.get("DAMAPS").polymerizable_group == "Al"


# rq-9769c65a
def test_compose_dmaps(lib):
    z = compose("MA", "C2", "NMe2", "C3", "SO3")
    assert z.smiles == lib.get("DMAPS").smiles
    assert z.status == "user"
    assert z.name == "MA-C2-NMe2-C3-SO3"


# rq-7141c236
def test_compose_spacer_shorthand():
    z = compose("MA", "C2", "NMe2", "C6", "SO3")
    patt = Chem.MolFromSmarts("[NX4+][CH2][CH2][CH2][CH2][CH2][CH2]S(=O)(=O)[O-]")
    assert z.to_mol(add_hs=False).HasSubstructMatch(patt)
    assert sum(charges(z)) == 0


# rq-f26c0a78
def test_compose_inverted_phosphorylcholine(lib):
    assert compose("MA", "C2", "PO4i", "C2", "NMe3").smiles == lib.get("MPC").smiles


# rq-e3d5290e
def test_compose_without_linker(lib):
    assert compose("V", "Py2", "C3", "SO3").smiles == lib.get("SPV").smiles


# rq-30c6d511
def test_compose_single_group_zwitterion():
    q = charges(compose("MA", "C2", "NO"))
    assert sum(q) == 0 and q.count(1) == 1 and q.count(-1) == 1


# rq-fb61858f
def test_compose_unknown_block():
    with pytest.raises(UnknownBlockError) as exc:
        compose("MA", "C2", "XYZ", "C3", "SO3")
    assert exc.value.key == "XYZ"


# rq-bf757833
def test_compose_non_neutral():
    with pytest.raises(ChargeError) as exc:
        compose("MA", "C2", "NMe3")
    assert exc.value.net_charge == 1


# rq-3c28b422
def test_compose_terminal_block_in_middle():
    with pytest.raises(IncompatibleBlocksError) as exc:
        compose("MA", "SO3", "C3", "NMe3")
    assert exc.value.position == 1


# rq-372fced8
def test_compose_must_start_polymerizable():
    with pytest.raises(IncompatibleBlocksError) as exc:
        compose("C2", "NMe2", "C3", "SO3")
    assert exc.value.position == 0


# rq-e70c4dde
def test_from_smiles_classifies():
    z = from_smiles("C=C(C)C(=O)OCC[N+]1(CCCS(=O)(=O)[O-])CCCCC1", "MyPip")
    assert z.polymerizable_group == "MA"
    assert "Pip" in z.cation_classes and "SO3" in z.anion_classes
    assert z.status == "user"


# rq-878e77a9
def test_from_smiles_unparsable():
    with pytest.raises(InvalidStructureError):
        from_smiles("C=C(C", "Bad")


# rq-3e28fc65
def test_from_smiles_charged():
    with pytest.raises(ChargeError) as exc:
        from_smiles("C=C(C)C(=O)OCC[N+](C)(C)C", "METAC")
    assert exc.value.net_charge == 1


# rq-bf11ad76
def test_from_smiles_small_molecule():
    with pytest.raises(NoPolymerizableGroupError):
        from_smiles("C[N+](C)(C)CC(=O)[O-]", "Betaine")
    z = from_smiles("C[N+](C)(C)CC(=O)[O-]", "Betaine", require_polymerizable=False)
    assert not z.polymerizable and z.vinyl_atoms is None


# rq-9ea65b16
def test_from_smiles_two_vinyl_groups():
    smi = "C=C(C)C(=O)OCC[N+](C)(Cc1ccc(C=C)cc1)CCCS(=O)(=O)[O-]"
    with pytest.raises(AmbiguousPolymerizableGroupError) as exc:
        from_smiles(smi, "Two")
    assert len(exc.value.candidates) == 2
    z = from_smiles(smi, "Two", vinyl_atoms=(0, 1))
    assert z.polymerizable_group == "MA"
    head, tail = z.vinyl_atoms
    mol = z.to_mol(add_hs=False)
    assert mol.GetBondBetweenAtoms(head, tail).GetBondType() == Chem.BondType.DOUBLE
    assert mol.GetAtomWithIdx(head).GetTotalNumHs() == 2


# rq-456dad55
def test_repeat_unit_dmaps(lib):
    expected = Chem.CanonSmiles("[*:1]CC([*:2])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]")
    assert lib.get("DMAPS").repeat_unit_smiles() == expected


# rq-052fdbed
def test_repeat_unit_small_molecule(lib):
    with pytest.raises(NoPolymerizableGroupError):
        lib.get("TMAO").repeat_unit_smiles()


# rq-3ae76166
def test_load_user_library(tmp_path, lib):
    path = write_yaml(tmp_path, """
        version: 1
        monomers:
          - name: MYSB
            full_name: my sulfobetaine
            family: Sulfobetaine
            status: specialty
            blocks: [AAm, C3, Pip, C4, SO3]
    """)
    user = load_library(path)
    assert "MYSB" in user
    assert all(z.name in user for z in lib)
    assert user.get("MYSB").polymerizable_group == "AAm"


# rq-bd5de7cc
def test_user_file_redefines_builtin_alias(tmp_path):
    path = write_yaml(tmp_path, """
        version: 1
        monomers:
          - name: NEWCB
            aliases: [cbma-1]
            full_name: x
            family: Carboxybetaine
            status: specialty
            blocks: [MA, C2, NMe2, C1, COO]
    """)
    with pytest.raises(DuplicateEntryError) as exc:
        load_library(path)
    assert exc.value.name == "cbma-1"


# rq-6698acbe
def test_user_file_invalid_yaml(tmp_path):
    path = write_yaml(tmp_path, "monomers: [unclosed\n")
    with pytest.raises(LibraryFileError) as exc:
        load_library(path)
    assert exc.value.path == str(path)


# rq-a50b7bb2
def test_user_file_missing():
    with pytest.raises(LibraryFileError) as exc:
        load_library("missing.yaml")
    assert exc.value.path == "missing.yaml"


# rq-63630b82
def test_user_file_unknown_field(tmp_path):
    path = write_yaml(tmp_path, """
        version: 1
        monomers:
          - name: TYPO
            full_name: x
            family: Sulfobetaine
            status: specialty
            smile: C=C
    """)
    with pytest.raises(LibraryFileError) as exc:
        load_library(path)
    assert "smile" in exc.value.reason


# rq-134bdc27
def test_user_file_unsupported_version(tmp_path):
    with pytest.raises(LibraryFileError):
        load_library(write_yaml(tmp_path, "version: 2\n"))


# rq-cf1496b1
def test_user_entry_blocks_and_smiles_disagree(tmp_path):
    path = write_yaml(tmp_path, """
        version: 1
        monomers:
          - name: MISMATCH
            full_name: x
            family: Sulfobetaine
            status: specialty
            blocks: [MA, C2, NMe2, C3, SO3]
            smiles: C=C(C)C(=O)OCC[N+](C)(C)CCC(=O)[O-]
    """)
    with pytest.raises(LibraryFileError) as exc:
        load_library(path)
    assert "MISMATCH" in exc.value.reason


# rq-3008e140
def test_user_entry_not_neutral(tmp_path):
    path = write_yaml(tmp_path, """
        version: 1
        monomers:
          - name: CATIONIC
            full_name: x
            family: Cationic
            status: specialty
            smiles: C=C(C)C(=O)OCC[N+](C)(C)C
    """)
    with pytest.raises(LibraryFileError) as exc:
        load_library(path)
    assert "CATIONIC" in exc.value.reason


# rq-33491788
def test_catalog_document_up_to_date():
    result = subprocess.run(
        [sys.executable, str(REPO / "devtools" / "build_zwitterion_catalog.py")],
        capture_output=True, text=True, check=True,
    )
    committed = (REPO / "docs" / "zwitterion_catalog.md").read_text(encoding="utf-8")
    assert result.stdout == committed, "regenerate docs/zwitterion_catalog.md"


# Initiators and end groups
SULFATE_ON_CH2 = Chem.MolFromSmarts("[CH2]OS(=O)(=O)[O-]")
SULFATE_ON_QUAT_C = Chem.MolFromSmarts("[CX4H0]OS(=O)(=O)[O-]")


def chain_mol(chain):
    return Chem.MolFromSmiles(chain.smiles)


# rq-30331843
def test_builtin_kps_initiator(lib):
    kps = lib.initiator("SO4")
    assert kps is lib.initiator("kso4")
    assert kps.name == "KPS"
    assert Chem.CanonSmiles(kps.head_smiles) == Chem.CanonSmiles("[*:1]OS(=O)(=O)[O-]")
    assert Chem.CanonSmiles(kps.tail_smiles) == Chem.CanonSmiles("[*:1]OS(=O)(=O)[O-]")
    assert kps.head_charge == kps.tail_charge == -1
    assert kps.counterion == "K+"


# rq-f5dfd922
def test_unknown_initiator(lib):
    with pytest.raises(UnknownInitiatorError) as exc:
        lib.initiator("KPSS")
    assert "KPS" in exc.value.suggestions


# rq-9218ad81
def test_kps_caps_both_ends_by_default():
    chain = build_chain(["DMAPS"] * 3, initiator="KPS")
    assert chain.ends is Ends.BOTH and chain.net_charge == -2
    mol = chain_mol(chain)
    assert len(mol.GetSubstructMatches(SULFATE_ON_CH2)) == 1
    assert len(mol.GetSubstructMatches(SULFATE_ON_QUAT_C)) == 1


# rq-94a7ac1f
def test_kps_caps_head_only():
    chain = build_chain(["DMAPS"] * 3, initiator="KPS", ends="head", cap="H")
    assert chain.net_charge == -1
    mol = chain_mol(chain)
    assert mol.HasSubstructMatch(SULFATE_ON_CH2)
    assert not mol.HasSubstructMatch(SULFATE_ON_QUAT_C)
    # Tail alpha-carbon gained a hydrogen: an alpha-methyl tertiary CH with the ester.
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("[CH3][CH1]([CH2])C(=O)O"))


# rq-4fa24446
def test_kps_caps_tail_only():
    chain = build_chain(["DMAPS"] * 3, initiator="KPS", ends="tail")
    assert chain.net_charge == -1
    mol = chain_mol(chain)
    assert mol.HasSubstructMatch(SULFATE_ON_QUAT_C)
    assert not mol.HasSubstructMatch(SULFATE_ON_CH2)
    # Head CH2 became CH3 bonded to the quaternary alpha carbon of the first unit.
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("[CH3][CX4H0]([CH3])([CH2])C(=O)O"))


# rq-f9862ff3
def test_chain_without_end_groups(lib):
    chain = build_chain(["DMAPS"] * 3, cap="H")
    assert chain.ends is Ends.NONE and chain.net_charge == 0 and chain.counterion is None
    expected = Chem.MolFromSmiles(".".join([lib.get("DMAPS").smiles] * 3 + ["[H][H]"]))
    assert chain.formula == rdMolDescriptors.CalcMolFormula(expected)
    methyl = Chem.MolFromSmiles(".".join([lib.get("DMAPS").smiles] * 3 + ["C"]))  # H + CH3 = CH4
    assert build_chain(["DMAPS"] * 3).formula == rdMolDescriptors.CalcMolFormula(methyl)


# rq-0c9b91ab
def test_end_cap_without_initiator():
    with pytest.raises(EndGroupError):
        build_chain(["DMAPS"] * 3, ends="both")


# rq-3ed1d11c
def test_raft_head_and_tail_differ():
    mol = chain_mol(build_chain(["DMAPS"] * 2, initiator="CPA"))
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("N#CC(C)(CCC(=O)O)[CH2]"))
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4]SC(=S)c1ccccc1"))


# rq-21f96ecb
def test_copolymer_sequence():
    chain = build_chain(["DMAPS", "CBMA", "DMAPS"], initiator="KPS")
    assert chain.units == ("DMAPS", "CBMA", "DMAPS")
    assert chain.net_charge == -2


# rq-8d1fb6d8
def test_custom_initiator_in_code():
    init = custom_initiator("MYINIT", head="[*:1]OC(=O)c1ccccc1")
    chain = build_chain(["DMAPS"] * 2, initiator=init)
    assert chain.ends is Ends.HEAD and chain.initiator == "MYINIT"
    with pytest.raises(EndGroupError):
        build_chain(["DMAPS"] * 2, initiator=init, ends="both")


# rq-b042d5d7
def test_custom_initiator_invalid_fragment():
    with pytest.raises(InvalidStructureError):
        custom_initiator("BAD", head="[*:1]OC(=O)[*:2]")


# rq-82181c68
def test_custom_initiator_without_fragments():
    with pytest.raises(InvalidStructureError):
        custom_initiator("EMPTY")


# rq-c62fdb4d
def test_initiator_in_user_file(tmp_path):
    path = write_yaml(tmp_path, """
        version: 1
        initiators:
          - name: VA044
            full_name: imidazolinium azo initiator
            head: '[*:1]C(C)(C)C1=[NH+]CCN1'
            tail: '[*:1]C(C)(C)C1=[NH+]CCN1'
            counterion: Cl-
    """)
    user = load_library(path)
    init = user.initiator("VA044")
    assert init.head_charge == 1 and init.tail_charge == 1
    assert build_chain(["MPC"] * 2, initiator="VA044", library=user).net_charge == 2


# rq-dec70936
def test_user_initiator_redefines_builtin(tmp_path):
    path = write_yaml(tmp_path, """
        version: 1
        initiators:
          - name: MYKPS
            aliases: [k2s2o8]
            full_name: x
            head: '[*:1]O'
    """)
    with pytest.raises(DuplicateEntryError) as exc:
        load_library(path)
    assert exc.value.name == "k2s2o8"


# rq-42202c0e
def test_chain_from_small_molecule():
    with pytest.raises(ChainError):
        build_chain(["TMAO"])


# rq-cffe1741
def test_empty_chain():
    with pytest.raises(ChainError):
        build_chain([])


# Variants: spacer length and N-alkyl length
def ch2_between(z, a_smarts, b_smarts, n):
    patt = Chem.MolFromSmarts(a_smarts + "[CH2]" * n + b_smarts)
    longer = Chem.MolFromSmarts(a_smarts + "[CH2]" * (n + 1) + b_smarts)
    mol = z.to_mol(add_hs=False)
    return mol.HasSubstructMatch(patt) and not mol.HasSubstructMatch(longer)


# rq-9fa59470
def test_variant_longer_spacer():
    v = make_variant("DMAPS", spacer=5)
    assert ch2_between(v, "[NX4+]", "S(=O)(=O)[O-]", 5)
    assert v.name == "DMAPS-sp5" and v.status == "user"
    assert sum(a.GetFormalCharge() for a in v.to_mol().GetAtoms()) == 0


# rq-1c79d863
def test_variant_inverted_spacer():
    v = make_variant("MPC", spacer=4)
    assert ch2_between(v, "P(=O)([O-])O", "[NX4+]", 4)


# rq-10b3fee1
def test_variant_alkyl(lib):
    v = make_variant("DMAPS", alkyl=2)
    assert v.to_mol(add_hs=False).HasSubstructMatch(Chem.MolFromSmarts("[CH3][CH2][NX4+][CH2][CH3]"))
    assert v.smiles == lib.get("SBMA-Et").smiles


# rq-f083b2fa
def test_variant_spacer_and_alkyl():
    v = make_variant("CBMA", spacer=3, alkyl=3)
    assert v.name == "CBMA-sp3-alk3"
    mol = v.to_mol(add_hs=False)
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("[CH3][CH2][CH2][NX4+][CH2][CH2][CH3]"))
    assert ch2_between(v, "[NX4+]", "C(=O)[O-]", 3)


# rq-046183e6
def test_variant_noop(lib):
    assert make_variant("DMAPS", spacer=3, alkyl=1) is lib.get("DMAPS")


# rq-942f1853
def test_variant_non_plain_spacer():
    with pytest.raises(VariantError) as exc:
        make_variant("HSBMA", spacer=4)
    assert exc.value.option == "spacer"


# rq-33cd927c
def test_variant_alkyl_non_ammonium():
    with pytest.raises(VariantError) as exc:
        make_variant("SPVI", alkyl=2)
    assert exc.value.option == "alkyl"


# rq-9ca58870
def test_variant_without_blocks():
    with pytest.raises(VariantError):
        make_variant("M1", spacer=2)


# rq-5ec68583
def test_variant_invalid_length():
    with pytest.raises(ValueError):
        make_variant("DMAPS", spacer=0)


# rq-5769e593
def test_builtin_comonomers(lib):
    for name, charge, ion in (("NIPAM", 0, None), ("IBOA", 0, None), ("VBD", 2, "Cl-"), ("BVD", 2, "Cl-")):
        z = lib.get(name)
        assert z.polymerizable and not z.zwitterionic and z.counterion == ion, name
        assert Chem.GetFormalCharge(z.to_mol(add_hs=False)) == charge, name
    with pytest.raises(VariantError):
        make_variant("VBD", alkyl=2)


# rq-914103fa
def test_comonomer_charge_rules(tmp_path):
    for smiles, ion in (("C=Cc1ccc(C[N+](C)(C)C)cc1", None), ("C=CC(=O)NC(C)C", "Cl-"),
                        ("C=Cc1ccc(C[N+](C)(C)C)cc1", "K+")):
        path = tmp_path / "c.yaml"
        rec = {"name": "X", "full_name": "X", "family": "test", "status": "specialty", "smiles": smiles}
        if ion:
            rec["counterion"] = ion
        path.write_text(yaml.safe_dump({"version": 1, "comonomers": [rec]}))
        with pytest.raises(LibraryFileError):
            load_library(path)


# rq-b7a9f7b4
def test_copolymer_with_cationic_comonomer():
    chain = build_chain(["DMAPS", "VBD", "DMAPS"], initiator="KPS")
    assert chain.net_charge == 0 and chain.counterion is None
    chain = build_chain(["VBD", "VBD"], initiator="KPS")
    assert chain.net_charge == 2 and chain.counterion == "Cl-"
    assert build_chain(["DMAPS"] * 2, initiator="KPS").counterion == "K+"


# rq-21cde3c2
def test_symmetric_divinyl():
    mol = Chem.MolFromSmiles(build_chain(["BVD"] * 3).smiles)
    double = [b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE
              and b.GetBeginAtom().GetSymbol() == "C" == b.GetEndAtom().GetSymbol() and not b.GetIsAromatic()]
    assert len(double) == 3
