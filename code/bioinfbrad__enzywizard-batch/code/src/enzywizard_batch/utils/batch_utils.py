from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, List, Tuple

from ..utils.logging_utils import Logger
from ..utils.integrate_utils import save_integrate_json, split_integrated_graph_entries
from ..utils.common_utils import get_optimized_filename
from ..utils.structure_utils import get_single_chain
from ..utils.sequence_utils import normalize_aa_name_to_one_letter
from Bio.PDB.Structure import Structure


def validate_batch_parameter_ranges(
    logger: Logger,
    cutoff_area: float = 10.0,
    minimize_energy: bool = True,
    minimization_iteration: int = 100,
    flexibility_cutoff: float = 15.0,
    n_modes: int = 20,
    window_size: int = 11,
    min_region_length: int = 5,
    pocket_min_rad: float = 1.8,
    pocket_max_rad: float = 6.2,
    pocket_min_volume: int = 50,
    max_synonyms: int = 20,
    fp_radius: int = 2,
    n_bits: int = 512,
    num_confs: int = 5,
    prune_rms: float = 0.5,
    max_docking_attempt_num: int = 20,
    exhaustiveness: int = 8,
    dock_min_rad: float = 1.8,
    dock_max_rad: float = 6.2,
    dock_min_volume: int = 50,
    use_manual_docking_box: bool = False,
    bonded_h_min_distance_A: float = 0.8,
    bonded_h_max_distance_A: float = 1.3,
    da_max_distance_A: float = 3.9,
    ha_max_distance_A: float = 2.5,
    dha_min_angle_deg: float = 90.0,
    ionic_distance_cutoff_A: float = 4.0,
    mu: float = 0.01,
    ring_center_distance_cutoff_A: float = 6.5,
    ring_cation_distance_cutoff_A: float = 5.0,
    ring_cation_angle_cutoff_deg: float = 45.0,
    ss_max_distance_A: float = 2.5,
    docked_heavy_atom_distance_cutoff_A: float = 6.5,
    min_residue_index_gap: int = 3,
) -> bool:

    if cutoff_area <= 0:
        logger.print(f"[ERROR] Invalid cutoff_area: {cutoff_area}. Must be a positive number.")
        return False

    if minimization_iteration <= 0:
        logger.print(f"[ERROR] Invalid minimization_iteration: {minimization_iteration}. Must be a positive integer.")
        return False

    if flexibility_cutoff <= 0 or n_modes <= 0:
        logger.print(f"[ERROR] Invalid cutoff or n_modes. Must be a positive number.")
        return False

    if not isinstance(window_size, int) or window_size < 3 or window_size > 50 or window_size % 2 == 0:
        logger.print(f"[ERROR] Invalid window_size: {window_size}. Must be odd integer in [3, 50].")
        return False

    if min_region_length < 3 or min_region_length > 50 or min_region_length > window_size:
        logger.print(
            f"[ERROR] Invalid min_region_length: {min_region_length}. Must be integer in [3, 50] and ≤ window_size ({window_size})."
        )
        return False

    if pocket_min_rad < 1.2 or pocket_max_rad <= 0 or pocket_min_volume <= 20 or pocket_min_rad >= pocket_max_rad:
        logger.print(
            f"[ERROR] Invalid pocket parameters. Require: min_rad ≥ 1.2, max_rad > min_rad, min_volume > 20."
        )
        return False

    if (
        max_synonyms <= 0 or max_synonyms > 200
        or fp_radius <= 0 or fp_radius > 5
        or n_bits <= 0 or n_bits > 2048
        or num_confs <= 0 or num_confs > 20
        or prune_rms <= 0 or prune_rms > 5.0
    ):
        logger.print(
            f"[ERROR] Invalid substrate generation parameters. Require: max_synonyms (1–200), fp_radius (1–5), n_bits (1–2048), num_confs (1–20), prune_rms (0–5]."
        )
        return False

    if (
        max_docking_attempt_num <= 0 or max_docking_attempt_num > 100
        or exhaustiveness <= 0 or exhaustiveness > 64
    ):
        logger.print(
            f"[ERROR] Invalid docking parameters. Require: max_docking_attempt_num (1–100), exhaustiveness (1–64)."
        )
        return False

    if not use_manual_docking_box and (
        dock_min_rad < 1.2
        or dock_min_volume <= 20
        or dock_min_rad >= dock_max_rad
    ):
        logger.print(
            f"[ERROR] Invalid docking pocket parameters. Require: min_rad ≥ 1.2, max_rad > min_rad, min_volume > 20."
        )
        return False

    if bonded_h_min_distance_A < 0.5 or bonded_h_min_distance_A > 1.1:
        logger.print(f"[ERROR] bonded_h_min_distance_A out of range [0.5, 1.1]: {bonded_h_min_distance_A}")
        return False

    if bonded_h_max_distance_A < 1.0 or bonded_h_max_distance_A > 1.5:
        logger.print(f"[ERROR] bonded_h_max_distance_A out of range [1.0, 1.5]: {bonded_h_max_distance_A}")
        return False

    if bonded_h_min_distance_A >= bonded_h_max_distance_A:
        logger.print(f"[ERROR] bonded_h_min_distance_A must be < bonded_h_max_distance_A.")
        return False

    if da_max_distance_A < 2.5 or da_max_distance_A > 4.5:
        logger.print(f"[ERROR] da_max_distance_A out of range [2.5, 4.5]: {da_max_distance_A}")
        return False

    if ha_max_distance_A < 1.5 or ha_max_distance_A > 3.0:
        logger.print(f"[ERROR] ha_max_distance_A out of range [1.5, 3.0]: {ha_max_distance_A}")
        return False

    if dha_min_angle_deg < 60 or dha_min_angle_deg > 180:
        logger.print(f"[ERROR] dha_min_angle_deg out of range [60, 180]: {dha_min_angle_deg}")
        return False

    if ionic_distance_cutoff_A < 2.0 or ionic_distance_cutoff_A > 6.0:
        logger.print(f"[ERROR] ionic_distance_cutoff_A out of range [2.0, 6.0]: {ionic_distance_cutoff_A}")
        return False

    if mu < 0.001 or mu > 0.1:
        logger.print(f"[ERROR] mu out of range [0.001, 0.1]: {mu}")
        return False

    if ring_center_distance_cutoff_A < 4.0 or ring_center_distance_cutoff_A > 8.0:
        logger.print(
            f"[ERROR] ring_center_distance_cutoff_A out of range [4.0, 8.0]: {ring_center_distance_cutoff_A}"
        )
        return False

    if ring_cation_distance_cutoff_A < 3.0 or ring_cation_distance_cutoff_A > 7.0:
        logger.print(
            f"[ERROR] ring_cation_distance_cutoff_A out of range [3.0, 7.0]: {ring_cation_distance_cutoff_A}"
        )
        return False

    if ring_cation_angle_cutoff_deg < 0 or ring_cation_angle_cutoff_deg > 90:
        logger.print(
            f"[ERROR] ring_cation_angle_cutoff_deg out of range [0, 90]: {ring_cation_angle_cutoff_deg}"
        )
        return False

    if ss_max_distance_A < 1.8 or ss_max_distance_A > 3.0:
        logger.print(f"[ERROR] ss_max_distance_A out of range [1.8, 3.0]: {ss_max_distance_A}")
        return False

    if docked_heavy_atom_distance_cutoff_A < 4.0 or docked_heavy_atom_distance_cutoff_A > 10.0:
        logger.print(
            f"[ERROR] docked_heavy_atom_distance_cutoff_A out of range [4.0, 10.0]: {docked_heavy_atom_distance_cutoff_A}"
        )
        return False

    if min_residue_index_gap < 1 or min_residue_index_gap > 5:
        logger.print(f"[ERROR] min_residue_index_gap out of range [1, 5]: {min_residue_index_gap}")
        return False

    return True


def build_batch_output_paths(protein_name: str, msa_name: str, output_dir: str | Path) -> Dict[str, Path]:
    output_dir = Path(output_dir)

    return {
        "cleaned_sto": output_dir / get_optimized_filename(f"cleaned_{msa_name}.sto"),
        "hmm": output_dir / get_optimized_filename(f"hmm_profile_{msa_name}.hmm"),
    }


def save_batch_integrate_outputs(
    integrate_report: Dict[str, Any],
    output_dir: str | Path,
    protein_name: str,
    logger: Logger,
) -> bool:
    output_dir = Path(output_dir)

    report_path = output_dir / get_optimized_filename(f"integrate_report_{protein_name}.json")
    nodes_path = output_dir / get_optimized_filename(f"integrate_nodes_{protein_name}.json")
    edges_path = output_dir / get_optimized_filename(f"integrate_edges_{protein_name}.json")

    if not save_integrate_json(integrate_report, report_path, logger):
        return False
    logger.print(f"[INFO] Report JSON saved: {report_path}")

    integrated_graph = integrate_report.get("integrated_graph")
    if not isinstance(integrated_graph, list):
        logger.print("[ERROR] integrated_graph missing in batch integrate report.")
        return False

    split_result = split_integrated_graph_entries(integrated_graph, logger)
    if split_result is None:
        return False

    node_list, edge_list = split_result

    if not save_integrate_json(node_list, nodes_path, logger):
        return False
    logger.print(f"[INFO] Node list JSON saved: {nodes_path}")

    if not save_integrate_json(edge_list, edges_path, logger):
        return False
    logger.print(f"[INFO] Edge list JSON saved: {edges_path}")

    return True

def build_identity_clean_mapping_from_structure(structure: Structure,logger: Logger) -> Tuple[List[Dict[str, Dict[str, Any]]], Dict[str, int]] | None:

    chain = get_single_chain(structure, logger)
    if chain is None:
        return None

    mapping_old_to_new: List[Dict[str, Dict[str, Any]]] = []

    kept_residues = 0

    for res in chain.get_residues():
        hetflag, resseq, icode = res.id

        # 跳过非蛋白
        if str(hetflag).strip():
            continue

        resname = res.get_resname().strip()

        h_count = sum(
            1 for atom in res.get_atoms()
            if atom.element and atom.element.strip().upper() == "H"
        )

        residue_info = {
            "aa_index": int(resseq),
            "aa_name": normalize_aa_name_to_one_letter(resname),
            "hydrogen_atom_count": h_count,
        }

        mapping_old_to_new.append(
            {
                "old_residue": residue_info,
                "new_residue": residue_info.copy(),  # identity mapping
            }
        )

        kept_residues += 1

    if kept_residues == 0:
        logger.print("[ERROR] No valid protein residues found in structure.")
        return None

    stats = {
        "removed_heterogen": 0,
        "changed_resname": 0,
        "fixed_residues": 0,
        "added_heavy_atoms": 0,
        "added_hydrogen_atoms": 0,
        "kept_residues": kept_residues,
    }

    return mapping_old_to_new, stats
