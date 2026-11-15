from __future__ import annotations
from pathlib import Path
from ..utils.logging_utils import Logger
from typing import Optional, List, Dict, Set,Any
import requests
from ..utils.substrate_utils import clean_compound_name, build_retry_session, chebi_search_exact, choose_best_chebi_smiles, pubchem_name_to_cid,pubchem_cid_to_smiles, pubchem_cid_to_synonyms, expand_synonyms_with_brackets, is_valid_smiles, get_mol_from_smiles, get_mol_h_from_mol_2d, get_fingerprint_from_mol_2d, get_minimized_mol_3d_list_from_mol_3d_list, get_mol_3d_list_from_mol_h, get_uff_energy_from_mol_3d ,get_2d_descriptor_dict_from_mol_2d, get_3d_descriptor_dict_from_mol_3d, postprocess_substrate_report_to_schema
from ..utils.common_utils import get_clean_filename, get_optimized_filename


def get_smiles_from_substrate_name(substrate_name: str, logger: Logger, session: Optional[requests.Session] = None, timeout: int = 15, max_pubchem_synonyms_to_retry_chebi: int = 20) -> str | None:
    if is_valid_smiles(substrate_name):
        logger.print("[ERROR] SMILES was provided not substrate name.")
        return None

    try:
        cleaned_name = clean_compound_name(substrate_name)
        if not cleaned_name:
            logger.print("[ERROR] Empty/abnormal substrate name was provided.")
            return None

        used_session = build_retry_session(session)

        chebi_hits = chebi_search_exact(compound_name=cleaned_name, session=used_session, logger=logger,timeout=timeout)
        chebi_smiles = choose_best_chebi_smiles(query_name=cleaned_name,hit_list=chebi_hits,session=used_session,logger=logger,timeout=timeout)
        if chebi_smiles is not None and is_valid_smiles(chebi_smiles):
            return chebi_smiles

        pubchem_cid = pubchem_name_to_cid(compound_name=cleaned_name,session=used_session,logger=logger,timeout=timeout)
        if pubchem_cid is None:
            logger.print(f"[WARNING] Failed to obtain SMILES for substrate: {cleaned_name}")
            return ""

        pubchem_smiles = pubchem_cid_to_smiles(cid=pubchem_cid,session=used_session,logger=logger,timeout=timeout)
        if pubchem_smiles is not None and is_valid_smiles(pubchem_smiles):
            return pubchem_smiles

        synonym_list = pubchem_cid_to_synonyms(cid=pubchem_cid,session=used_session,logger=logger,timeout=timeout)

        expanded_synonym_list: List[str] = []
        seen_keys: Set[str] = set()

        for synonym in synonym_list:
            for expanded_synonym in expand_synonyms_with_brackets(synonym):
                synonym_key = expanded_synonym.casefold()
                if synonym_key in seen_keys:
                    continue
                seen_keys.add(synonym_key)
                expanded_synonym_list.append(expanded_synonym)

        expanded_synonym_list = expanded_synonym_list[:max_pubchem_synonyms_to_retry_chebi]

        for synonym in expanded_synonym_list:
            chebi_hits = chebi_search_exact(compound_name=synonym,session=used_session,logger=logger,timeout=timeout)
            chebi_smiles = choose_best_chebi_smiles(query_name=synonym,hit_list=chebi_hits,session=used_session,logger=logger,timeout=timeout)
            if chebi_smiles is not None and is_valid_smiles(chebi_smiles):
                return chebi_smiles

        logger.print(f"[WARNING] Failed to obtain SMILES for substrate: {cleaned_name}")
        return ""

    except Exception as e:
        logger.print(f"[ERROR] Unexpected error when resolving substrate: {substrate_name}. Reason: {e}")
        return None

def get_substrate_dict_list_from_input(substrate_names: str, logger: Logger) -> List[Dict[str, str]] | None:
    if not isinstance(substrate_names, str):
        logger.print("[ERROR] substrate_names must be a string.")
        return None

    try:
        names = [x.strip() for x in substrate_names.split(";") if x.strip()]
        if len(names) == 0:
            logger.print("[ERROR] substrate_names is empty after parsing.")
            return None

        result: List[Dict[str, str]] = []
        smiles_count = 0
        substrate_name_set: Set[str] = set()

        for name in names:
            if is_valid_smiles(name):
                smiles_count += 1
                substrate_name = f"smiles{smiles_count}"
                smiles = name
            else:
                substrate_name = name
                smiles = ""

            if substrate_name in substrate_name_set:
                logger.print(f"[ERROR] Duplicate substrate name detected: {substrate_name}")
                return None

            substrate_name_set.add(substrate_name)
            result.append({
                "substrate_name": substrate_name,
                "smiles": smiles
            })


        return result

    except Exception as e:
        logger.print(f"[ERROR] Failed to parse substrate_names. Reason: {e}")
        return None

def get_completed_smiles_list(substrate_list: List[Dict[str, str]],logger: Logger,session: Optional[requests.Session] = None, max_synonyms=20) -> List[Dict[str, str]] | None:

    if not isinstance(substrate_list, list):
        logger.print("[ERROR] substrate_list must be a list.")
        return None

    try:
        session = build_retry_session(session)

        for item in substrate_list:
            if item.get("smiles", "") != "":
                continue

            name = item.get("substrate_name", "")
            if not name:
                logger.print("[ERROR] Missing substrate_name.")
                return None

            smiles = get_smiles_from_substrate_name(name, logger, session=session, max_pubchem_synonyms_to_retry_chebi=max_synonyms)

            if smiles is None:
                logger.print(f"[ERROR] Failed to fetch SMILES for substrate: {name}")
                return None

            if not smiles:
                logger.print(f"[ERROR] Failed to obtain SMILES for substrate: {name}")
                return None

            item["smiles"] = smiles

        return substrate_list

    except Exception as e:
        logger.print(f"[ERROR] Failed to complete SMILES list. Reason: {e}")
        return None

def get_substrate_feature_list(substrate_list: List[Dict[str, str]],logger: Logger,fp_radius: int = 2, n_bits: int = 512,num_confs: int = 5, prune_rms: float = 0.5) -> List[Dict[str, Any]] | None:

    if not isinstance(substrate_list, list):
        logger.print("[ERROR] substrate_list must be a list.")
        return None

    result: List[Dict[str, Any]] = []

    try:
        for item in substrate_list:
            name = item.get("substrate_name", "")
            smiles = item.get("smiles", "")

            if not isinstance(name, str) or not name.strip():
                logger.print("[ERROR] Missing substrate_name.")
                return None

            if not isinstance(smiles, str) or not smiles.strip():
                logger.print(f"[ERROR] Missing SMILES for substrate: {name}")
                return None

            out: Dict[str, Any] = {
                "substrate_name": name,
                "smiles": smiles,
                "fingerprint": [],
                "num_atoms": 0,
                "mol_weight": 0.0,
                "logp": 0.0,
                "tpsa": 0.0,
                "heavy_atom_count": 0,
                "hbond_donor_count": 0,
                "hbond_acceptor_count": 0,
                "rotatable_bond_count": 0,
                "molar_refractivity": 0.0,
                "structures": []
            }

            mol_2d = get_mol_from_smiles(smiles, logger)
            if mol_2d is None:
                logger.print(f"[ERROR] Failed to generate Mol(2D) for substrate: {name}")
                return None

            fp = get_fingerprint_from_mol_2d(mol_2d, logger, radius=fp_radius,n_bits=n_bits)
            desc = get_2d_descriptor_dict_from_mol_2d(mol_2d, logger)

            if fp is None:
                logger.print(f"[ERROR] Failed to generate fingerprint for substrate: {name}")
                return None

            if desc is None:
                logger.print(f"[ERROR] Failed to generate 2D descriptors for substrate: {name}")
                return None

            out["fingerprint"] = fp
            out["num_atoms"] = desc.get("substrate_num_atoms")
            out["mol_weight"] = desc.get("substrate_molecular_weight")
            out["logp"] = desc.get("substrate_mol_logp")
            out["tpsa"] = desc.get("substrate_tpsa")
            out["heavy_atom_count"] = desc.get("substrate_num_heavy_atoms")
            out["hbond_donor_count"] = desc.get("substrate_hbond_donor_count")
            out["hbond_acceptor_count"] = desc.get("substrate_hbond_acceptor_count")
            out["rotatable_bond_count"] = desc.get("substrate_rotatable_bond_count")
            out["molar_refractivity"] = desc.get("substrate_molar_refractivity")

            mol_h = get_mol_h_from_mol_2d(mol_2d, logger)
            if mol_h is None:
                result.append(out)
                continue

            mol_3d_list = get_mol_3d_list_from_mol_h(mol_h, logger,num_confs=num_confs,prune_rms=prune_rms)
            if mol_3d_list is None:
                result.append(out)
                continue

            mol_3d_list = get_minimized_mol_3d_list_from_mol_3d_list(mol_3d_list, logger)
            if mol_3d_list is None:
                result.append(out)
                continue

            structures: List[Dict[str, Any]] = []

            for i, mol_3d in enumerate(mol_3d_list, start=1):
                energy = get_uff_energy_from_mol_3d(mol_3d, logger)
                descriptor_3d = get_3d_descriptor_dict_from_mol_3d(mol_3d, logger)

                structure_name = f"{name}_{i}"

                if energy is None or descriptor_3d is None:
                    logger.print(f"[WARNING] Skipping incomplete substrate structure: {structure_name}")
                    continue

                structure_dict = {
                    "structure_name": structure_name,
                    "structure_energy": energy,
                    "structure_mol": mol_3d
                }
                structure_dict.update(descriptor_3d)

                structures.append(structure_dict)

            structures.sort(key=lambda x: x["structure_energy"])

            for i, structure_dict in enumerate(structures, start=1):
                structure_dict["structure_name"] = f"{name}_{i}"

            out["structures"] = structures
            result.append(out)

        return result

    except Exception as e:
        logger.print(f"[ERROR] Failed to generate substrate feature list. Reason: {e}")
        return None

def generate_substrate_report(substrate_feature_list: List[Dict[str, Any]], logger: Logger) -> Dict[str, Any] | None:

    if not isinstance(substrate_feature_list, list):
        logger.print("[ERROR] substrate_feature_list must be a list.")
        return None

    try:
        cleaned_list = []

        for item in substrate_feature_list:
            if not isinstance(item, dict):
                logger.print("[ERROR] Invalid substrate feature entry.")
                return None

            required_substrate_fields = [
                "substrate_name",
                "smiles",
                "fingerprint",
                "num_atoms",
                "mol_weight",
                "logp",
                "tpsa",
                "heavy_atom_count",
                "hbond_donor_count",
                "hbond_acceptor_count",
                "rotatable_bond_count",
                "molar_refractivity",
                "structures",
            ]
            for field_name in required_substrate_fields:
                if field_name not in item:
                    logger.print(f"[ERROR] Missing substrate feature field: {field_name}")
                    return None

            new_item = dict(item)

            new_structures = []
            for s in item.get("structures", []):
                if not isinstance(s, dict):
                    logger.print("[ERROR] Invalid substrate structure entry.")
                    return None

                required_structure_fields = [
                    "structure_name",
                    "structure_energy",
                    "structure_max_3d_diameter",
                    "structure_mean_pairwise_atom_distance",
                    "structure_std_pairwise_atom_distance",
                    "structure_asphericity",
                    "structure_spherocity",
                    "structure_principal_moment_ratio",
                    "structure_radius_of_gyration",
                ]
                for field_name in required_structure_fields:
                    if field_name not in s:
                        logger.print(f"[ERROR] Missing substrate structure field: {field_name}")
                        return None

                new_structures.append({
                    "structure_name": Path(
                        get_optimized_filename(
                            f"{get_clean_filename(s['structure_name'])}.sdf"
                        )
                    ).stem,
                    "structure_energy": s["structure_energy"],
                    "structure_max_3d_diameter": s["structure_max_3d_diameter"],
                    "structure_mean_pairwise_atom_distance": s["structure_mean_pairwise_atom_distance"],
                    "structure_std_pairwise_atom_distance": s["structure_std_pairwise_atom_distance"],
                    "structure_asphericity": s["structure_asphericity"],
                    "structure_spherocity": s["structure_spherocity"],
                    "structure_principal_moment_ratio": s["structure_principal_moment_ratio"],
                    "structure_radius_of_gyration": s["structure_radius_of_gyration"],
                })

            new_item["structures"] = new_structures
            cleaned_list.append(new_item)

        raw_report = {
            "output_type": "enzywizard_substrate",
            "substrates": cleaned_list
        }

        return postprocess_substrate_report_to_schema(raw_report, logger)

    except Exception as e:
        logger.print(f"[ERROR] Failed to generate substrate report. Reason: {e}")
        return None
