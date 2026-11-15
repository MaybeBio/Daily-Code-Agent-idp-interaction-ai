from __future__ import annotations
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from ..utils.logging_utils import Logger
from ..utils.structure_utils import get_single_chain, get_residues_by_chain
from typing import Dict, List, Tuple, Any
from ..utils.clean_utils import (
    standardize_resname,
    choose_atom_altloc,
    normalize_atom_name,
    is_hydrogen_atom,
    normalize_aa_name_to_one_letter,
    get_single_chain_from_pdbfixer,
    count_single_chain_heterogens_from_pdbfixer_chain,
    build_nonstandard_residue_count_from_fixer,
    build_missing_atom_stats_from_fixer,
    get_single_chain_protein_residue_info_from_pdbfixer_chain,
    count_hydrogen_atoms_in_fixer,
    renumber_single_chain_fixer_residues,
    postprocess_clean_report_to_schema
)
from ..resources.aa_resources import (
    AA3_STANDARD,
    AA3_EXPECTED_HEAVY_ATOM_SET,
    AA3_ALLOWED_HEAVY_ATOM_SET_WITH_OXT,
    BACKBONE_REQUIRED_ATOMS,
)
from Bio.PDB.Atom import Atom
from openmm.app import Modeller, ForceField, PDBFile
import math
from pdbfixer import PDBFixer
import numpy as np
from openmm import unit

def clean_pdbfixer_to_single_chain_A(fixer: PDBFixer,add_H: bool,logger: Logger,pH: float = 7.0,force_field_file: str = "charmm36.xml") -> Tuple[PDBFixer, List[Dict[str, Dict[str, Any]]], Dict[str, int]] | None:
    topology = fixer.topology
    if topology is None:
        logger.print(f"[ERROR] No topology found in PDBFixer")
        return None

    chain_list = list(topology.chains())
    num_chains = len(chain_list)
    if num_chains == 0:
        logger.print(f"[ERROR] No chain found in PDBFixer topology")
        return None

    # keep only the first chain
    if num_chains > 1:
        try:
            fixer.removeChains(range(1, num_chains))
        except Exception as e:
            logger.print(f"[ERROR] Exception in removing extra chains: {e}")
            return None

    # obtain single old chain
    old_chain = get_single_chain_from_pdbfixer(fixer, logger)
    if old_chain is None:
        return None

    # mapping
    mapping_old_to_new: List[Dict[str, Dict[str, Any]]] = []

    # locate non-standard residues
    try:
        fixer.findNonstandardResidues()
    except Exception as e:
        logger.print(f"[ERROR] Exception in findNonstandardResidues: {e}")
        return None

    # old statistics before cleaning
    removed_heterogen = count_single_chain_heterogens_from_pdbfixer_chain(old_chain)
    changed_resname = build_nonstandard_residue_count_from_fixer(fixer)

    # build old mapping side
    old_residue_info_list = get_single_chain_protein_residue_info_from_pdbfixer_chain(old_chain)
    for old_residue_info in old_residue_info_list:
        mapping_old_to_new.append(
            {
                "old_residue": old_residue_info,
                "new_residue": {},
            }
        )

    # execute cleaning
    try:
        fixer.removeHeterogens(keepWater=False)
    except Exception as e:
        logger.print(f"[ERROR] Exception in removeHeterogens: {e}")
        return None

    try:
        fixer.replaceNonstandardResidues()
    except Exception as e:
        logger.print(f"[ERROR] Exception in replaceNonstandardResidues: {e}")
        return None

    fixer.missingResidues = {}

    try:
        fixer.findMissingAtoms()
    except Exception as e:
        logger.print(f"[ERROR] Exception in findMissingAtoms: {e}")
        return None

    # statistics after replacement and after missing atoms are located
    fixed_residues, added_heavy_atoms = build_missing_atom_stats_from_fixer(fixer)

    try:
        fixer.addMissingAtoms(seed=202602)
    except Exception as e:
        logger.print(f"[ERROR] Exception in addMissingAtoms: {e}")
        return None

    hydrogen_count_before = count_hydrogen_atoms_in_fixer(fixer)

    if add_H:
        try:
            ff = ForceField(force_field_file) if force_field_file else None
            fixer.addMissingHydrogens(pH=pH, forcefield=ff)
        except Exception as e:
            logger.print(f"[ERROR] Exception in addMissingHydrogens: {e}")
            return None

    hydrogen_count_after = count_hydrogen_atoms_in_fixer(fixer)
    added_hydrogen_atoms = hydrogen_count_after - hydrogen_count_before
    if added_hydrogen_atoms < 0:
        added_hydrogen_atoms = 0

    # check NaN coordinates before renumbering / writing
    try:
        positions_array = np.array(fixer.positions.value_in_unit(unit.nanometer))
        if np.isnan(positions_array).any():
            logger.print("[ERROR] NaN coordinates detected after PDBFixer cleaning")
            return None
    except Exception as e:
        logger.print(f"[ERROR] Exception in checking coordinates after cleaning: {e}")
        return None

    # renumber after fixing
    fixer = renumber_single_chain_fixer_residues(fixer, logger)
    if fixer is None:
        return None

    # update new mapping side
    new_chain = get_single_chain_from_pdbfixer(fixer, logger)
    if new_chain is None:
        return None

    new_residue_info_list = get_single_chain_protein_residue_info_from_pdbfixer_chain(new_chain)

    if len(old_residue_info_list) != len(new_residue_info_list):
        logger.print(
            f"[ERROR] Old and new protein residue counts are inconsistent: "
            f"{len(old_residue_info_list)} vs {len(new_residue_info_list)}"
        )
        return None

    for i, new_residue_info in enumerate(new_residue_info_list):
        mapping_old_to_new[i]["new_residue"] = new_residue_info

    kept_residues = len(new_residue_info_list)

    stats = {
        "removed_heterogen": removed_heterogen,
        "changed_resname": changed_resname,
        "fixed_residues": fixed_residues,
        "added_heavy_atoms": added_heavy_atoms,
        "added_hydrogen_atoms": added_hydrogen_atoms,
        "kept_residues": kept_residues,
    }

    return fixer, mapping_old_to_new, stats


def generate_clean_report(structure: Structure,cleaned_structure: Structure,mapping_old_to_new: List[Dict[str, Dict[str, Any]]],stats: Dict[str, int],logger: Logger,) -> dict | None:

    amino_acid_mapping_old_2_new: List[Dict[str, Dict[str, Any]]] = []

    old_chain = get_single_chain(structure, logger)
    if old_chain is None:
        return None

    new_chain = get_single_chain(cleaned_structure, logger)
    if new_chain is None:
        return None

    old_residue_list: List[Residue] = []
    for res in old_chain.get_residues():
        hetflag, _, _ = res.id
        if str(hetflag).strip():
            continue
        old_residue_list.append(res)

    new_residue_list: List[Residue] = []
    for res in new_chain.get_residues():
        hetflag, _, _ = res.id
        if str(hetflag).strip():
            continue
        new_residue_list.append(res)

    if len(mapping_old_to_new) != len(old_residue_list):
        logger.print(
            f"[ERROR] Mapping length and old residue count are inconsistent: "
            f"{len(mapping_old_to_new)} vs {len(old_residue_list)}"
        )
        return None

    if len(mapping_old_to_new) != len(new_residue_list):
        logger.print(
            f"[ERROR] Mapping length and cleaned residue count are inconsistent: "
            f"{len(mapping_old_to_new)} vs {len(new_residue_list)}"
        )
        return None

    for i, mapping_item in enumerate(mapping_old_to_new):
        old_res = old_residue_list[i]
        new_res = new_residue_list[i]

        old_resname = old_res.get_resname().strip()
        new_resname = new_res.get_resname().strip()

        old_h_count = sum(
            1
            for atom in old_res.get_atoms()
            if is_hydrogen_atom(atom)
        )

        new_h_count = sum(
            1
            for atom in new_res.get_atoms()
            if is_hydrogen_atom(atom)
        )

        old_residue_info = mapping_item.get("old_residue", {})
        new_residue_info = mapping_item.get("new_residue", {})

        amino_acid_mapping_old_2_new.append(
            {
                "old_residue": {
                    "aa_id": old_residue_info.get("aa_index"),
                    "aa_name": old_residue_info.get("aa_name", normalize_aa_name_to_one_letter(old_resname)),
                    "hydrogen_atom_count": old_h_count,
                },
                "new_residue": {
                    "aa_id": new_residue_info.get("aa_index"),
                    "aa_name": new_residue_info.get("aa_name", normalize_aa_name_to_one_letter(new_resname)),
                    "hydrogen_atom_count": new_h_count,
                },
            }
        )

    raw_report = {
        "output_type": "enzywizard_clean",
        "amino_acid_mapping_old_to_new": amino_acid_mapping_old_2_new,
        "clean_statistics": stats,
    }

    report = postprocess_clean_report_to_schema(raw_report)
    if report is None:
        return None

    return report

def check_cleaned_structure(struct: Structure, logger: Logger) -> bool:
    model_count = len(list(struct.get_models()))
    if model_count != 1:
        logger.print("[ERROR] Structure must contain exactly one model. Please run 'enzywizard clean' first.")
        return False

    model = next(struct.get_models())
    chains = list(model.get_chains())

    if len(chains) != 1:
        logger.print("[ERROR] Structure must contain exactly one chain. Please run 'enzywizard clean' first.")
        return False

    chain = chains[0]
    if chain.id != "A":
        logger.print("[ERROR] Cleaned structure must use chain ID 'A'. Please run 'enzywizard clean' first.")
        return False

    expected_resseq = 1

    for res in chain.get_residues():
        hetflag, resseq, icode = res.id

        if str(hetflag).strip():
            logger.print(f"[ERROR] Non-protein or hetero residue detected at residue {resseq}{icode}. Please run 'enzywizard clean' first.")
            return False

        if str(icode).strip():
            logger.print(f"[ERROR] Insertion code detected at residue {resseq}{icode}. Please run 'enzywizard clean' first.")
            return False

        resname = res.get_resname().strip()
        resname_std = standardize_resname(resname)
        if resname != resname_std:
            logger.print(f"[ERROR] Residue name '{resname}' at residue {resseq} is not standardized. Please run 'enzywizard clean' first.")
            return False

        if resname not in AA3_STANDARD:
            logger.print(f"[ERROR] Non-standard residue '{resname}' detected at residue {resseq}. Please run 'enzywizard clean' first.")
            return False

        if int(resseq) != expected_resseq:
            logger.print(
                f"[ERROR] Residue numbering is not continuous: expected {expected_resseq}, got {resseq}. Please run 'enzywizard clean' first."
            )
            return False

        atoms_by_name: Dict[str, List[Atom]] = {}
        for atom in res.get_atoms():
            if is_hydrogen_atom(atom):
                continue

            atom_name = normalize_atom_name(atom.get_name())
            atoms_by_name.setdefault(atom_name, []).append(atom)

        for atom_name in BACKBONE_REQUIRED_ATOMS:
            if atom_name not in atoms_by_name:
                logger.print(
                    f"[ERROR] Missing backbone atom '{atom_name}' at residue {resseq}. Please run 'enzywizard clean' first.")
                return False

        expected_heavy_atom_set = AA3_EXPECTED_HEAVY_ATOM_SET.get(resname)
        if expected_heavy_atom_set is None:
            logger.print(f"[ERROR] No expected heavy atom definition for residue '{resname}' at residue {resseq}.")
            return False

        actual_heavy_atom_set = set(atoms_by_name.keys())

        if not expected_heavy_atom_set.issubset(actual_heavy_atom_set):
            missing_atom_name_list = sorted(expected_heavy_atom_set - actual_heavy_atom_set)
            logger.print(
                f"[ERROR] Missing heavy atoms at residue {resseq} ({resname}): {missing_atom_name_list}. "
                f"Please run 'enzywizard clean' first."
            )
            return False

        allowed_heavy_atom_set = AA3_ALLOWED_HEAVY_ATOM_SET_WITH_OXT[resname]
        if not actual_heavy_atom_set.issubset(allowed_heavy_atom_set):
            unexpected_atom_name_list = sorted(actual_heavy_atom_set - allowed_heavy_atom_set)
            logger.print(
                f"[ERROR] Unexpected heavy atoms at residue {resseq} ({resname}): {unexpected_atom_name_list}. "
                f"Please run 'enzywizard clean' first."
            )
            return False

        expected_resseq += 1

    return True