from __future__ import annotations
from typing import Dict, List, Any
from Bio.PDB.Structure import Structure
from ..utils.logging_utils import Logger
from ..utils.structure_utils import get_single_chain,get_residues_by_chain
import numpy as np
import prody as pd
from ..utils.sequence_utils import normalize_aa_name_to_one_letter

def compute_protein_rmsf(struct: Structure, logger: Logger, cutoff: float = 15.0, n_modes: int = 20, method: str = "ANM") -> List[Dict[str, Any]] | None:

    results: List[Dict[str, Any]] = []

    chain = get_single_chain(struct, logger)
    if chain is None:
        return None

    residue_list = get_residues_by_chain(chain, logger)
    if residue_list is None:
        return None

    if len(residue_list) < 3:
        logger.print("[ERROR] At least 3 residues with CA atoms are required for RMSF calculation.")
        return None

    if cutoff <= 0:
        logger.print("[ERROR] cutoff must be greater than 0.")
        return None

    try:
        ca_coords = np.array([coord for _, _, coord in residue_list], dtype=float)
    except Exception as e:
        logger.print(f"[ERROR] Failed to collect CA coordinates: {e}")
        return None
    if not np.all(np.isfinite(ca_coords)):
        logger.print("[ERROR] Non-finite CA coordinates detected.")
        return None

    method = str(method).upper().strip()

    try:
        if method == "ANM":
            max_modes = max(1, 3 * len(ca_coords) - 6)
            use_modes = max(1, min(int(n_modes), max_modes))

            model = pd.ANM("Protein_ANM")
            model.buildHessian(ca_coords, cutoff=float(cutoff))
            model.calcModes(n_modes=use_modes)
            sqf = pd.calcSqFlucts(model[:use_modes])

        elif method == "GNM":
            max_modes = max(1, len(ca_coords) - 1)
            use_modes = max(1, min(int(n_modes), max_modes))

            model = pd.GNM("Protein_GNM")
            model.buildKirchhoff(ca_coords, cutoff=float(cutoff))
            model.calcModes(n_modes=use_modes)
            sqf = pd.calcSqFlucts(model[:use_modes])

        else:
            logger.print(f"[ERROR] Unsupported RMSF method: {method}. Use 'ANM' or 'GNM'.")
            return None

    except Exception as e:
        logger.print(f"[ERROR] Failed to calculate RMSF by ProDy: {e}")
        return None

    try:
        rmsf = np.sqrt(np.asarray(sqf, dtype=float))
    except Exception as e:
        logger.print(f"[ERROR] Failed to convert square fluctuations to RMSF: {e}")
        return None

    if len(rmsf) != len(residue_list):
        logger.print("[ERROR] RMSF result length does not match residue count.")
        return None
    if not np.all(np.isfinite(rmsf)):
        logger.print("[ERROR] Non-finite RMSF values detected.")
        return None

    try:
        for (res_id, resname, coord), rmsf_value in zip(residue_list, rmsf):
            _, resseq, _ = res_id

            result_dict: Dict[str, Any] = {
                "aa_id": resseq,
                "aa_name": normalize_aa_name_to_one_letter(resname),
                "rmsf": float(rmsf_value),
            }
            results.append(result_dict)
    except Exception as e:
        logger.print(f"[ERROR] Failed to assemble RMSF results: {e}")
        return None

    return results


def postprocess_flexibility_report_to_schema(raw_report: Dict[str, Any]) -> Dict[str, Any]:
    """
    Postprocess the raw EnzyWizard-Flexibility report into the new JSON Schema format.

    """

    protein_rmsf = raw_report.get("protein_rmsf", [])

    protein_flexibility: List[Dict[str, Any]] = []

    for residue_rmsf in protein_rmsf:
        protein_flexibility.append(
            {
                "residue_index": residue_rmsf.get("aa_id"),
                "residue_name": residue_rmsf.get("aa_name"),
                "residue_root_mean_square_fluctuation": residue_rmsf.get("rmsf"),
            }
        )

    schema_report: Dict[str, Any] = {
        "report_type": raw_report.get("output_type", "enzywizard_flexibility"),
        "protein_flexibility": protein_flexibility,
    }

    return schema_report

def generate_flexibility_report(protein_rmsf: List[Dict[str, Any]]) -> Dict[str, Any]:

    raw_report: Dict[str, Any] = {
        "output_type": "enzywizard_flexibility",
        "protein_rmsf": protein_rmsf,
    }

    return postprocess_flexibility_report_to_schema(raw_report)
