from __future__ import annotations
from Bio.PDB.Structure import Structure
from ..utils.logging_utils import Logger
from typing import Dict, List, Tuple, Optional, Union, Any
from ..utils.structure_utils import get_single_chain, get_residues_by_chain
from ..utils.aaprops_utils import get_dssp_fields, postprocess_aaprops_report_to_schema
from Bio.PDB.DSSP import DSSP
from Bio.Data.IUPACData import protein_letters_3to1
from ..utils.common_utils import one_hot_vec_generator, multi_hot_vec_generator
from ..resources.aa_resources import AA_20NAME_INDEX, AA_8CLASSES, AA_8CLASSES_TO_INDICES, DSSP_8STATE_INDEX
from ..resources.aa_physicochemical_props import net_charge_dict, pka_dict, volume_dict, hydrophobicity_dict, molecular_weight_dict, pi_dict
from ..utils.sequence_utils import normalize_aa_name_to_one_letter

def calculate_aa_props(struct: Structure, dssp: DSSP, logger: Logger) -> List[Dict[str,Any]] | None:
    return_list: List[Dict[str,Any]] =[]

    chain=get_single_chain(struct,logger)
    if chain is None:
        return None

    residue_list=get_residues_by_chain(chain,logger)
    if residue_list is None:
        return None

    for residue in residue_list:
        return_dict: Dict[str,Any] = {}
        (hetflag, resseq, icode), resname, coord = residue
        chain_id = chain.get_id()
        dssp_fields = get_dssp_fields(dssp, (chain_id, (hetflag, resseq, icode)), logger)
        if dssp_fields is None:
            logger.print(f"[ERROR] DSSP fields missing for residue {resseq}")
            return None
        ss, rsa, phi, psi = dssp_fields

        return_dict['aa_id']=resseq
        aa=protein_letters_3to1.get(resname.capitalize())
        if aa is None:
            logger.print(f"[ERROR] Unrecognized residue name: {resname}")
            return None
        return_dict['aa_name']=normalize_aa_name_to_one_letter(aa)
        return_dict['aa_name_one_hot']=one_hot_vec_generator(len(AA_20NAME_INDEX), aa, AA_20NAME_INDEX)
        return_dict['aa_class'] = "/".join(AA_8CLASSES[i][0] for i in AA_8CLASSES_TO_INDICES.get(aa, []))
        return_dict['aa_class_one_hot'] = multi_hot_vec_generator(len(AA_8CLASSES),aa,AA_8CLASSES_TO_INDICES)
        return_dict['aa_ss']=ss
        return_dict['aa_ss_one_hot']=one_hot_vec_generator(len(DSSP_8STATE_INDEX), return_dict['aa_ss'],DSSP_8STATE_INDEX)
        return_dict['aa_rsa'] = rsa
        return_dict['aa_phi'] = phi
        return_dict['aa_psi'] = psi
        return_dict['aa_net_charge']=net_charge_dict[aa]
        return_dict['aa_pka'] = pka_dict[aa]
        return_dict['aa_volume'] = volume_dict[aa]
        return_dict['aa_hydrophobicity'] = hydrophobicity_dict[aa]
        return_dict['aa_molecular_weight'] = molecular_weight_dict[aa]
        return_dict['aa_pi'] = pi_dict[aa]
        return_dict['aa_coord'] = coord
        return_list.append(return_dict)
    return return_list

def calculate_aa_props_statistics(aa_props_list: List[Dict[str, Any]], logger: Logger) -> Dict[str, Dict[str, int]] | None:
    if not isinstance(aa_props_list, list):
        logger.print("[ERROR] aa_props_list must be a list.")
        return None

    aa_name_statistics: Dict[str, int] = {aa: 0 for aa in AA_20NAME_INDEX}
    aa_class_statistics: Dict[str, int] = {class_name: 0 for class_name, _ in AA_8CLASSES}
    aa_ss_statistics: Dict[str, int] = {ss: 0 for ss in DSSP_8STATE_INDEX}

    for i, aa_dict in enumerate(aa_props_list):
        if not isinstance(aa_dict, dict):
            logger.print(f"[ERROR] Invalid aa_props_list item at index {i}: not a dict.")
            return None

        aa_name = aa_dict.get("aa_name")
        aa_ss = aa_dict.get("aa_ss")
        aa_class = aa_dict.get("aa_class")

        if not isinstance(aa_name, str) or aa_name not in aa_name_statistics:
            logger.print(f"[ERROR] Invalid aa_name at index {i}: {aa_name}")
            return None
        aa_name_statistics[aa_name] += 1

        if not isinstance(aa_ss, str) or aa_ss not in aa_ss_statistics:
            logger.print(f"[ERROR] Invalid aa_ss at index {i}: {aa_ss}")
            return None
        aa_ss_statistics[aa_ss] += 1

        if not isinstance(aa_class, str):
            logger.print(f"[ERROR] Invalid aa_class at index {i}: {aa_class}")
            return None

        for cls in aa_class.split("/"):
            cls = cls.strip()
            if cls not in aa_class_statistics:
                logger.print(f"[ERROR] Unknown aa_class '{cls}' at index {i}")
                return None
            aa_class_statistics[cls] += 1

    return {
        "aa_name_statistics": aa_name_statistics,
        "aa_class_statistics": aa_class_statistics,
        "aa_ss_statistics": aa_ss_statistics,
    }

def generate_aaprops_report(aa_props: List[Dict[str,Any]], aa_props_statistics: Dict[str, Dict[str, int]]) -> dict:

    raw_report = {
        "output_type": "enzywizard_aaprops",
        "aa_props": aa_props,
        "aa_props_statistics": aa_props_statistics,
    }

    return postprocess_aaprops_report_to_schema(raw_report)


